"""
Fraggler Diagnostics — Job Execution

AsyncExecutor for background threads + pipeline/QC/DIT job wrappers.
Compatible with Python 3.10+.
"""
from __future__ import annotations

import os
import shutil
import tempfile
import threading
from pathlib import Path
from typing import Callable, List, Optional, Any

import param

from core.log import log
from config import APP_SETTINGS

# We lazy import the analysis engine inside the job functions to prevent 
# `fraggler` from immediately executing `pn.extension(template="fast")` 
# and polluting the global Panel state during application startup.

# ============================================================
# CONSTANTS
# ============================================================

CHUNK_SIZE = 25
SAFE_MAX_FILES_PER_PATIENT = 5000


# ============================================================
# ASYNC EXECUTOR
# ============================================================

class AsyncExecutor(param.Parameterized):
    """Run a callable in a background daemon thread."""
    running = param.Boolean(default=False)

    def __init__(self, **params):
        super().__init__(**params)
        self._thread: Optional[threading.Thread] = None

    def run_background(self, target: Callable, *args, **kwargs) -> bool:
        if self.running:
            log("[WARN] A job is already running.")
            return False
        self.running = True

        def _worker():
            try:
                target(*args, **kwargs)
            except Exception as e:
                log(f"[ERROR] Job exception: {e}")
            finally:
                self.running = False
                log("[INFO] Job finished.")

        self._thread = threading.Thread(target=_worker, daemon=True)
        self._thread.start()
        return True


executor = AsyncExecutor()


# ============================================================
# FILE STAGING HELPERS
# ============================================================

import hashlib


def _safe_link_name(p: Path, idx: int) -> str:
    h = hashlib.md5(str(p).encode("utf-8")).hexdigest()[:8]
    return f"{idx:05d}_{h}_{p.name}"


def stage_files(files: List[Path], use_symlink: bool = True) -> Path:
    """Create a temp directory with symlinks (or copies) of all files."""
    tdir = Path(tempfile.mkdtemp(prefix="fraggler_stage_"))
    for i, src in enumerate(sorted(files), start=1):
        dst = tdir / _safe_link_name(src, i)
        if use_symlink:
            try:
                os.symlink(src, dst)
            except Exception:
                shutil.copy2(src, dst)
        else:
            shutil.copy2(src, dst)
    return tdir


def cleanup_temp(p: Optional[Path]) -> None:
    if p and p.exists():
        shutil.rmtree(p, ignore_errors=True)


def build_filtered_input(src: Path, needle: str) -> Optional[Path]:
    """Create a temp directory with staged .fsa files matching needle."""
    matched = [p for p in sorted(src.glob("*.fsa")) if needle.lower() in p.name.lower()]
    if not matched:
        return None
    log(f"[INFO] Custom filter: {len(matched)} files matched '{needle}'.")
    return stage_files(matched)


def _should_stage_explicit_files_once(fsa_dir: Optional[Path], files: List[Path], *, chunk_files: bool) -> bool:
    """Stage a complete folder-sized explicit file set in one shot when safe."""
    if not files:
        return False
    if not chunk_files:
        return True

    if fsa_dir is None:
        return False

    folder = Path(fsa_dir)
    if not folder.is_dir():
        return False

    resolved_files = [Path(p) for p in files]
    if any(p.parent != folder for p in resolved_files):
        return False

    from core.utils import is_water_file

    folder_files = [
        p
        for p in sorted(folder.glob("*.fsa"))
        if p.is_file() and not is_water_file(p.name)
    ]
    if len(folder_files) != len(resolved_files):
        return False

    return {p.name for p in folder_files} == {p.name for p in resolved_files}


def _can_use_exact_source_dir(fsa_dir: Optional[Path], files: List[Path], *, chunk_files: bool) -> bool:
    """Reuse the original folder directly when the explicit file list already matches it exactly."""
    if not _should_stage_explicit_files_once(fsa_dir, files, chunk_files=chunk_files):
        return False
    return fsa_dir is not None and Path(fsa_dir).is_dir()


def _emit_progress(progress_callback, **event) -> None:
    if progress_callback is None:
        return
    payload = {
        "folder_name": "",
        "job_name": "",
        "phase": "",
        "file_name": "",
        "files_done": 0,
        "files_total": 0,
        "jobs_done": 0,
        "jobs_total": 0,
        "heartbeat_at": "",
        "note": "",
    }
    payload.update(event)
    progress_callback(payload)


# ============================================================
# PIPELINE JOB
# ============================================================

def run_pipeline_job(
    fsa_dir: Optional[Path],
    base_outdir: Path,
    out_folder_name: str,
    scope: str,
    needle: str,
    files: Optional[List[Path]] = None,
    *,
    update_tracking_workbook: bool = True,
) -> Optional[list]:
    """
    Run pipeline on a folder or an explicit file list.
    Returns entries if return_entries behaviour is needed.
    """
    effective_mode = "all"
    if scope == "controls":
        effective_mode = "controls"
    elif scope == "custom" and not needle.strip():
        raise ValueError("scope=custom requires an assay filter.")

    if not files:
        tmp_input = None
        try:
            if scope == "custom":
                tmp_input = build_filtered_input(fsa_dir, needle)
                if not tmp_input:
                    raise ValueError(f"No .fsa files matched '{needle}'.")
                fsa_to_run = tmp_input
            else:
                fsa_to_run = fsa_dir

            from core.pipeline import run_pipeline
            run_pipeline(
                fsa_dir=fsa_to_run,
                base_outdir=base_outdir,
                assay_folder_name=out_folder_name,
                mode=effective_mode,
                update_tracking_workbook=update_tracking_workbook,
            )
        finally:
            cleanup_temp(tmp_input)
        return None

    # Chunked mode for explicit file lists
    total = len(files)
    if total > SAFE_MAX_FILES_PER_PATIENT:
        raise ValueError(
            f"File count ({total}) exceeds SAFE_MAX={SAFE_MAX_FILES_PER_PATIENT}."
        )

    if scope == "custom":
        files = [p for p in files if needle.lower() in p.name.lower()]
        if not files:
            raise ValueError(f"No .fsa files matched '{needle}'.")

    if APP_SETTINGS.get("active_analysis") == "general":
        tmp_input = None
        try:
            tmp_input = stage_files(files)
            from core.pipeline import run_pipeline
            run_pipeline(
                fsa_dir=tmp_input,
                base_outdir=base_outdir,
                assay_folder_name=out_folder_name,
                mode=effective_mode,
                update_tracking_workbook=update_tracking_workbook,
            )
        finally:
            cleanup_temp(tmp_input)
        return None

    if _should_stage_explicit_files_once(fsa_dir, files, chunk_files=True):
        tmp_input = None
        try:
            run_input = fsa_dir if _can_use_exact_source_dir(fsa_dir, files, chunk_files=True) else stage_files(files)
            if run_input is not fsa_dir:
                tmp_input = run_input
            from core.pipeline import run_pipeline
            run_pipeline(
                fsa_dir=run_input,
                base_outdir=base_outdir,
                assay_folder_name=out_folder_name,
                mode=effective_mode,
                return_entries=True,
                make_dit_reports=False,
                update_tracking_workbook=update_tracking_workbook,
            )
        finally:
            cleanup_temp(tmp_input)
        return None

    ok_chunks = 0
    failed_chunks = 0
    collected_entries = []
    for offset in range(0, len(files), CHUNK_SIZE):
        chunk = files[offset: offset + CHUNK_SIZE]
        log(f"[INFO] Chunk {offset // CHUNK_SIZE + 1} ({len(chunk)} files)")
        tmp_input = None
        try:
            tmp_input = stage_files(chunk)
            from core.pipeline import run_pipeline
            chunk_entries = run_pipeline(
                fsa_dir=tmp_input,
                base_outdir=base_outdir,
                assay_folder_name=out_folder_name,
                mode=effective_mode,
                return_entries=True,
                make_dit_reports=False,
                update_tracking_workbook=update_tracking_workbook,
            )
            collected_entries.extend(chunk_entries or [])
            ok_chunks += 1
        except Exception as e:
            failed_chunks += 1
            log(f"[ERROR] Chunk failed: {e}")
        finally:
            cleanup_temp(tmp_input)

    if failed_chunks:
        raise RuntimeError(
            f"Pipeline completed with errors: {ok_chunks} ok, {failed_chunks} failed."
        )
    if collected_entries and effective_mode != "controls":
        from core.assay_config import OUTDIR_NAME
        from core.html_reports import build_dit_html_reports

        assay_outdir = base_outdir / (out_folder_name or OUTDIR_NAME)
        build_dit_html_reports(collected_entries, assay_outdir)
    return None


# ============================================================
# PIPELINE JOB (with return_entries for DIT aggregation)
# ============================================================

def run_pipeline_job_collect(
    fsa_dir: Optional[Path],
    base_outdir: Path,
    out_folder_name: str,
    scope: str,
    needle: str,
    files: Optional[List[Path]] = None,
    *,
    chunk_files: bool = True,
    tracking_excel_path: Path | None = None,
    progress_callback=None,
    update_tracking_workbook: bool = False,
) -> list:
    """
    Like run_pipeline_job but returns entries for cross-folder DIT aggregation.
    """
    effective_mode = "all"
    if scope == "controls":
        effective_mode = "controls"
    elif scope == "custom" and not needle.strip():
        raise ValueError("scope=custom requires an assay filter.")

    effective_in = fsa_dir
    tmp_input = None

    try:
        if files:
            if len(files) > SAFE_MAX_FILES_PER_PATIENT:
                raise ValueError(
                    f"File count ({len(files)}) exceeds SAFE_MAX={SAFE_MAX_FILES_PER_PATIENT}."
                )

            selected_files = files
            if scope == "custom" and needle.strip():
                selected_files = [p for p in files if needle.lower() in p.name.lower()]
                if not selected_files:
                    raise ValueError(f"No .fsa files matched '{needle}'.")

            from core.pipeline import run_pipeline
            if _should_stage_explicit_files_once(fsa_dir, selected_files, chunk_files=chunk_files):
                _emit_progress(
                    progress_callback,
                    phase="stage_files",
                    files_done=0,
                    files_total=len(selected_files),
                    note="using_exact_source_dir" if _can_use_exact_source_dir(fsa_dir, selected_files, chunk_files=chunk_files) else "staging_explicit_files_once",
                )
                run_input = fsa_dir if _can_use_exact_source_dir(fsa_dir, selected_files, chunk_files=chunk_files) else stage_files(selected_files)
                if run_input is not fsa_dir:
                    tmp_input = run_input
                try:
                    entries = run_pipeline(
                        fsa_dir=run_input,
                        base_outdir=base_outdir,
                        assay_folder_name=out_folder_name,
                        mode=effective_mode,
                        return_entries=True,
                        make_dit_reports=False,
                        tracking_excel_path=tracking_excel_path,
                        update_tracking_workbook=update_tracking_workbook,
                        progress_callback=progress_callback,
                    )
                    return entries or []
                finally:
                    cleanup_temp(tmp_input)
                    tmp_input = None

            if not chunk_files:
                _emit_progress(
                    progress_callback,
                    phase="stage_files",
                    files_done=0,
                    files_total=len(selected_files),
                    note="staging_explicit_files",
                )
                tmp_input = stage_files(selected_files)
                try:
                    entries = run_pipeline(
                        fsa_dir=tmp_input,
                        base_outdir=base_outdir,
                        assay_folder_name=out_folder_name,
                        mode=effective_mode,
                        return_entries=True,
                        make_dit_reports=False,
                        tracking_excel_path=tracking_excel_path,
                        update_tracking_workbook=update_tracking_workbook,
                        progress_callback=progress_callback,
                    )
                    return entries or []
                finally:
                    cleanup_temp(tmp_input)
                    tmp_input = None

            all_entries = []
            total_files = len(selected_files)
            for offset in range(0, len(selected_files), CHUNK_SIZE):
                chunk = selected_files[offset: offset + CHUNK_SIZE]
                chunk_index = (offset // CHUNK_SIZE) + 1
                _emit_progress(
                    progress_callback,
                    phase="stage_files",
                    files_done=offset,
                    files_total=total_files,
                    note=f"staging_chunk_{chunk_index}",
                )
                tmp_input = stage_files(chunk)
                try:
                    def _chunk_progress(event):
                        files_done = int(event.get("files_done", 0)) + offset
                        _emit_progress(
                            progress_callback,
                            **{
                                **event,
                                "files_done": min(total_files, files_done),
                                "files_total": total_files,
                            },
                        )

                    entries = run_pipeline(
                        fsa_dir=tmp_input,
                        base_outdir=base_outdir,
                        assay_folder_name=out_folder_name,
                        mode=effective_mode,
                        return_entries=True,
                        make_dit_reports=False,
                        tracking_excel_path=tracking_excel_path,
                        update_tracking_workbook=update_tracking_workbook,
                        progress_callback=_chunk_progress,
                    )
                    all_entries.extend(entries or [])
                finally:
                    cleanup_temp(tmp_input)
                    tmp_input = None
            return all_entries

        if scope == "custom" and needle.strip():
            filtered = build_filtered_input(effective_in, needle)
            if not filtered:
                raise ValueError(f"No .fsa files matched '{needle}'.")
            if tmp_input:
                cleanup_temp(tmp_input)
            tmp_input = filtered
            effective_in = filtered

        from core.pipeline import run_pipeline
        entries = run_pipeline(
            fsa_dir=effective_in,
            base_outdir=base_outdir,
            assay_folder_name=out_folder_name,
            mode=effective_mode,
            return_entries=True,
            make_dit_reports=False,  # We build DIT reports ourselves
            tracking_excel_path=tracking_excel_path,
            update_tracking_workbook=update_tracking_workbook,
            progress_callback=progress_callback,
        )
        return entries or []
    finally:
        cleanup_temp(tmp_input)


# ============================================================
# QC JOB
# ============================================================

def run_qc_job(
    fsa_dir: Optional[Path],
    base_outdir: Path,
    out_html_name: str,
    excel_name: str,
    rules: Any, # Use Any here to avoid top-level import of QCRules
    files: Optional[List[Path]] = None,
    *,
    tracking_excel_path: Path | None = None,
    update_tracking_workbook: bool = True,
    return_entries: bool = False,
    progress_callback=None,
) -> Optional[Path] | tuple[Optional[Path], list[dict]]:
    """Run QC analysis and return the path to the HTML report."""
    from datetime import datetime
    
    def _empty_or_unreadable_fsa_summary(input_dir: Optional[Path]) -> tuple[int, int]:
        if input_dir is None or not Path(input_dir).exists():
            return 0, 0
        total = 0
        empty_or_bad = 0
        for candidate in Path(input_dir).glob("*.fsa"):
            total += 1
            try:
                if candidate.stat().st_size <= 0:
                    empty_or_bad += 1
            except OSError:
                empty_or_bad += 1
        return total, empty_or_bad

    tmp_input = None
    try:
        effective_in = fsa_dir
        if files:
            tmp_input = stage_files(files)
            effective_in = tmp_input
            
        from core.pipeline import run_pipeline
        from core.assay_config import OUTDIR_NAME
        from core.qc.qc_html import build_qc_html
        from core.qc.qc_excel import update_excel_trends

        qc_outdir = base_outdir / OUTDIR_NAME
        qc_outdir.mkdir(parents=True, exist_ok=True)
        out_html = qc_outdir / out_html_name
        excel_path = base_outdir / excel_name

        entries = run_pipeline(
            fsa_dir=effective_in,
            base_outdir=base_outdir,
            return_entries=True,
            make_dit_reports=False,
            mode="controls",
            tracking_excel_path=tracking_excel_path,
            update_tracking_workbook=update_tracking_workbook,
            progress_callback=progress_callback,
        )
        if not entries:
            total_fsa, empty_or_bad = _empty_or_unreadable_fsa_summary(effective_in)
            if total_fsa > 0 and empty_or_bad == total_fsa:
                raise RuntimeError("All QC .fsa files were empty or unreadable.")
            raise RuntimeError("No QC entries found (check file names or skipped unreadable files).")

        build_qc_html(entries, out_html, rules, excel_path)
        run_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        update_excel_trends(excel_path, entries, rules, run_ts)

        if return_entries:
            return out_html, entries
        return out_html
    finally:
        cleanup_temp(tmp_input)


# ============================================================
# DIT JOB
# ============================================================

def run_dit_job(
    fsa_dir: Optional[Path],
    base_outdir: Path,
    out_folder_name: str,
    scope: str,
    needle: str,
    files: Optional[List[Path]] = None,
) -> None:
    """Run DIT report generation."""
    tmp_input = None
    try:
        effective_in = fsa_dir
        effective_mode = scope if scope != "controls" else "all"

        if files:
            tmp_input = stage_files(files)
            effective_in = tmp_input

        if scope == "custom":
            if not needle.strip():
                raise ValueError("scope=custom requires assay_filter_substring.")
            filtered = build_filtered_input(effective_in, needle)
            if not filtered:
                raise ValueError(f"No .fsa files matched '{needle}'.")
            if tmp_input:
                cleanup_temp(tmp_input)
            tmp_input = filtered
            effective_in = filtered

        from core.pipeline import run_pipeline
        run_pipeline(
            fsa_dir=effective_in,
            base_outdir=base_outdir,
            assay_folder_name=out_folder_name,
            mode=effective_mode if effective_mode != "controls" else "all",
        )
    finally:
        cleanup_temp(tmp_input)
