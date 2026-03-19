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
    """Create a temp directory with symlinks to .fsa files matching needle."""
    tmpdir = Path(tempfile.mkdtemp(prefix="fraggler_filter_"))
    count = 0
    for p in sorted(src.glob("*.fsa")):
        if needle.lower() in p.name.lower():
            os.symlink(p, tmpdir / p.name)
            count += 1
    if count == 0:
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None
    log(f"[INFO] Custom filter: {count} files matched '{needle}'.")
    return tmpdir


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
) -> list:
    """
    Like run_pipeline_job but returns entries for cross-folder DIT aggregation.
    """
    effective_mode = "all"
    if scope == "controls":
        effective_mode = "controls"

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
                    return []

            from core.pipeline import run_pipeline
            if not chunk_files:
                tmp_input = stage_files(selected_files)
                try:
                    entries = run_pipeline(
                        fsa_dir=tmp_input,
                        base_outdir=base_outdir,
                        assay_folder_name=out_folder_name,
                        mode=effective_mode,
                        return_entries=True,
                        make_dit_reports=False,
                    )
                    return entries or []
                finally:
                    cleanup_temp(tmp_input)
                    tmp_input = None

            all_entries = []
            for offset in range(0, len(selected_files), CHUNK_SIZE):
                chunk = selected_files[offset: offset + CHUNK_SIZE]
                tmp_input = stage_files(chunk)
                try:
                    entries = run_pipeline(
                        fsa_dir=tmp_input,
                        base_outdir=base_outdir,
                        assay_folder_name=out_folder_name,
                        mode=effective_mode,
                        return_entries=True,
                        make_dit_reports=False,
                    )
                    all_entries.extend(entries or [])
                finally:
                    cleanup_temp(tmp_input)
                    tmp_input = None
            return all_entries

        if scope == "custom" and needle.strip():
            filtered = build_filtered_input(effective_in, needle)
            if not filtered:
                return []
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
) -> Optional[Path]:
    """Run QC analysis and return the path to the HTML report."""
    from datetime import datetime

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
        )
        if not entries:
            raise RuntimeError("No QC entries found (check file names).")

        build_qc_html(entries, out_html, rules, excel_path)
        run_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        update_excel_trends(excel_path, entries, rules, run_ts)

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
