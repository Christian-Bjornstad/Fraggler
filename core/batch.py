"""
Fraggler Diagnostics — Batch Processing

Logic for automated scanning of subfolders / YAML files and executing
jobs (pipeline/qc/dit). Also implements cross-folder DIT aggregation.
Compatible with Python 3.10+.
"""
from __future__ import annotations

import os
import re
import time
import yaml
import threading
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

from config import resolve_analysis_excel_output_path
from core.log import log
from core.runner import run_pipeline_job, run_pipeline_job_collect, run_qc_job, run_dit_job

# Lazy load fraggler modules to prevent global Panel state pollution on import


# ============================================================
# SCANNING UTILITIES
# ============================================================

def scan_jobs_from_subfolders(base_dir: Path, target_depth: int = 1) -> List[Path]:
    """Scan base_dir for subdirectories at target_depth."""
    if not base_dir.is_dir():
        log(f"[WARN] Base directory not found: {base_dir}")
        return []
    
    if target_depth == 1:
        folders = [d for d in base_dir.iterdir() if d.is_dir()]
        return sorted(folders)
        
    log(f"[WARN] target_depth != 1 is not fully implemented for basic scan.")
    return []


def scan_jobs_from_yaml(yaml_path: Path) -> List[Path]:
    """Read a list of directories from a YAML file."""
    if not yaml_path.is_file():
        log(f"[WARN] YAML file not found: {yaml_path}")
        return []
        
    try:
        with open(yaml_path, "r", encoding="utf-8", errors="replace") as f:
            data = yaml.safe_load(f)
    except Exception as e:
        log(f"[ERROR] Could not parse YAML: {e}")
        return []

    if isinstance(data, dict):
        if "directories" in data:
            data = data.get("directories", [])
        elif "folders" in data:
            data = data.get("folders", [])
        else:
            log("[ERROR] YAML file must contain 'directories' or 'folders' list.")
            return []

    if not isinstance(data, list):
        log("[ERROR] YAML root must be a list or a dict with a list of directories.")
        return []

    valid_dirs = []
    base_dir = yaml_path.parent
    for item in data:
        p = Path(item)
        if not p.is_absolute():
            p = (base_dir / p).resolve()
        if p.is_dir():
            valid_dirs.append(p)
        else:
            log(f"[WARN] Skipping invalid path from YAML: {p}")
            
    return sorted(valid_dirs)


# ============================================================
# PATIENT AGGREGATION UTILITIES
# ============================================================

def _scan_folder_fsa_files(path: Path, folder_files: Dict[Path, List[Path]]) -> List[Path]:
    """Return cached .fsa files for a folder/file, scanning only once per call site."""
    if path not in folder_files:
        from core.utils import is_water_file

        def _is_usable_fsa(candidate: Path) -> bool:
            if candidate.suffix.lower() != ".fsa" or is_water_file(candidate.name):
                return False
            try:
                if candidate.stat().st_size <= 0:
                    log(f"[WARN] Skipping empty .fsa file during job scan: {candidate.name}")
                    return False
            except OSError as exc:
                log(f"[WARN] Skipping unreadable .fsa file during job scan: {candidate.name} ({exc})")
                return False
            return True

        if path.is_file():
            if _is_usable_fsa(path):
                folder_files[path] = [path]
            else:
                folder_files[path] = []
        elif not path.is_dir():
            folder_files[path] = []
        else:
            raw_candidates = [f for f in sorted(path.glob("*.fsa")) if not is_water_file(f.name)]
            folder_files[path] = [
                f for f in raw_candidates
                if _is_usable_fsa(f)
            ]
            if raw_candidates and not folder_files[path]:
                log(f"[WARN] Folder has only empty/unreadable .fsa files: {path}")
    return folder_files[path]


def find_all_fsa_files(
    paths: List[Path],
    folder_files: Dict[Path, List[Path]] | None = None,
) -> List[Path]:
    """Find all .fsa files in the given folders or individual paths."""
    folder_files = folder_files or {}

    fsa_files = []
    for p in paths:
        fsa_files.extend(_scan_folder_fsa_files(p, folder_files))
    return fsa_files


def group_files_by_patient(fsa_files: List[Path], regex_pattern: str) -> Dict[str, List[Path]]:
    """Group a list of .fsa files by patient ID extracted via regex.
    QC files starting with PK/NK/RK are separated into their own folder-based groups.
    """
    grouped = {}
    pattern = None
    if regex_pattern:
        try:
            pattern = re.compile(regex_pattern)
        except Exception as e:
            log(f"[ERROR] Invalid regex '{regex_pattern}': {e}")
            
    from core.utils import CONTROL_PREFIX_RE, is_water_file, strip_stage_prefix
    qc_pattern = CONTROL_PREFIX_RE

    for f in sorted(fsa_files):
        if is_water_file(f.name):
            continue

        # Use the stripped name for identification logic
        clean_name = strip_stage_prefix(f.name)
        
        if qc_pattern.match(clean_name):
            grouped.setdefault("QC", []).append(f)
            continue
            
        pid = None
        if pattern:
            match = pattern.search(clean_name)
            if match:
                pid = match.group()
        
        if not pid:
            # Fallback 1: Split by underscore and take first part
            stem = Path(clean_name).stem
            parts = stem.split("_")
            if len(parts) > 1:
                pid = parts[0]
            else:
                # Fallback 2: Just use the stem
                pid = stem
            
        grouped.setdefault(pid, []).append(f)
    
    return grouped


def generate_jobs(
    input_paths: List[Path], 
    aggregate_patients: bool = True,
    patient_regex: str = r"\d{2}OUM\d{5}"
) -> List[Dict[str, Any]]:
    """
    Generate a list of standard job dicts from a list of folders.
    Each job has:
      - 'name': str (display name)
      - 'type': str ("pipeline" or "qc")
      - 'path': Path (for subfolder mode) OR None
      - 'files': list[Path] (if aggregated or filtered)
    """
    from config import APP_SETTINGS

    folders_to_scan = []
    folder_files: Dict[Path, List[Path]] = {}
    active_analysis = APP_SETTINGS.get("active_analysis", "clonality")
    for p in input_paths:
        if p.is_file():
            if p.suffix.lower() == ".yaml":
                folders_to_scan.extend(scan_jobs_from_yaml(p))
                continue
            if p.suffix.lower() == ".fsa":
                folders_to_scan.append(p)
                continue

        if _scan_folder_fsa_files(p, folder_files):
            folders_to_scan.append(p)

        if p.is_dir():
            subfolders = [d for d in p.iterdir() if d.is_dir()]
            for sub in subfolders:
                if _scan_folder_fsa_files(sub, folder_files):
                    folders_to_scan.append(sub)

    folders_to_scan = list(dict.fromkeys(folders_to_scan))
        
    if not folders_to_scan:
        log(f"[WARN] No folders with .fsa data found.")
        return []

        
    jobs = []
    
    if aggregate_patients:
        all_fsa = find_all_fsa_files(folders_to_scan, folder_files)
        all_fsa = sorted(list(set(all_fsa)))
        
        grouped = group_files_by_patient(all_fsa, patient_regex)
        
        for pid, files in sorted(grouped.items()):
            jtype = "qc" if pid == "QC" else "pipeline"
            jobs.append({
                "name": pid,
                "type": jtype,
                "path": None,
                "files": files
            })
        if jobs:
            log(f"[INFO] Aggregated {len(all_fsa)} files into {len(jobs)} jobs.")
    else:
        from core.utils import CONTROL_PREFIX_RE, strip_stage_prefix
        qc_pattern = CONTROL_PREFIX_RE
        all_qc_files = []
        for folder in folders_to_scan:
            fsa_files = _scan_folder_fsa_files(folder, folder_files)
            qc_files = [f for f in fsa_files if qc_pattern.match(strip_stage_prefix(f.name))]
            pat_files = [f for f in fsa_files if not qc_pattern.match(strip_stage_prefix(f.name))]
            
            all_qc_files.extend(qc_files)
            
            if pat_files:
                jobs.append({
                    "name": folder.name,
                    "type": "pipeline",
                    "path": folder,
                    "files": pat_files
                })
        
        if all_qc_files:
            all_qc_files = sorted(list(set(all_qc_files)))
            jobs.append({
                "name": "QC",
                "type": "qc",
                "path": None,
                "files": all_qc_files
            })
            
    return jobs


def _extract_patient_ids_from_files(files: List[Path], regex_pattern: str) -> set[str]:
    from core.utils import CONTROL_PREFIX_RE, strip_stage_prefix

    try:
        pattern = re.compile(regex_pattern) if regex_pattern else None
    except re.error:
        return set()

    patient_ids: set[str] = set()
    for file_path in files or []:
        clean_name = strip_stage_prefix(Path(file_path).name)
        if CONTROL_PREFIX_RE.match(clean_name):
            continue
        if not pattern:
            return set()
        match = pattern.search(clean_name)
        if not match:
            return set()
        patient_ids.add(match.group())
    return patient_ids


def _can_stream_aggregated_dit_reports(jobs: List[Dict[str, Any]], patient_regex: str) -> bool:
    """
    Stream DIT report generation job-by-job when each pipeline/DIT job already maps
    to exactly one unique patient ID and no patient appears in more than one job.
    """
    seen_ids: set[str] = set()
    relevant_jobs = [job for job in jobs if job.get("type") in {"pipeline", "dit"}]
    if not relevant_jobs:
        return False

    for job in relevant_jobs:
        patient_ids = _extract_patient_ids_from_files(job.get("files", []), patient_regex)
        if len(patient_ids) != 1:
            return False
        patient_id = next(iter(patient_ids))
        if patient_id in seen_ids:
            return False
        seen_ids.add(patient_id)
    return True

# ============================================================
# BATCH EXECUTION RUNNER
# ============================================================

def run_batch_jobs(
    jobs: List[Dict[str, Any]],
    output_base: Path,
    out_folder_tmpl: str,
    outfile_html_tmpl: str,
    excel_name_tmpl: str,
    pipeline_scope: str,
    assay_filter: str,
    aggregate_dit_reports: bool,
    continue_on_error: bool,
    update_callback: Any = None,
    progress_callback: Any = None,
    max_workers: int | None = None,
    tracking_excel_path: Path | None = None,
    aggregate_outdir_name: str | None = None,
    defer_tracking_workbook_refresh: bool = False,
    defer_dit_html_reports: bool | None = None,
) -> dict[str, Any]:
    """
    Run all generated jobs.
    Reads QC parameters from APP_SETTINGS and uses explicit per-analysis batch flags.
    """
    from config import APP_SETTINGS
    from core.qc.qc_rules import QCRules
    from core.assay_config import OUTDIR_NAME
    s_qc = APP_SETTINGS.get("qc", {})
    active_analysis = APP_SETTINGS.get("active_analysis", "clonality")
    aggregate_dit_reports = bool(aggregate_dit_reports) and active_analysis != "general"
    analysis_batch = APP_SETTINGS.get("analyses", {}).get(active_analysis, {}).get("batch", {})
    patient_regex = analysis_batch.get("patient_id_regex", r"\d{2}OUM\d{5}")
    stream_aggregated_dit = aggregate_dit_reports and _can_stream_aggregated_dit_reports(jobs, patient_regex)
    defer_dit_html_reports = defer_tracking_workbook_refresh if defer_dit_html_reports is None else bool(defer_dit_html_reports)
    sample_window = s_qc.get("sample_peak_window_bp", s_qc.get("w_sample", 3.0))
    ladder_window = s_qc.get("ladder_peak_window_bp", s_qc.get("w_ladder", 3.0))
    qc_rules = QCRules(
        min_r2_ok=s_qc.get("min_r2_ok", 0.999),
        min_r2_warn=s_qc.get("min_r2_warn", 0.995),
        sample_peak_window_bp=sample_window,
        sample_peak_window_bp_fallback=s_qc.get("sample_peak_window_bp_fallback", max(float(sample_window) + 4.0, 8.0)),
        ladder_peak_window_bp=ladder_window,
    )
    
    from concurrent.futures import ThreadPoolExecutor

    total = len(jobs)
    log(f"[BATCH] Starting batch run of {total} jobs.")

    agg_outdir = output_base / (aggregate_outdir_name or OUTDIR_NAME) if aggregate_dit_reports else None
    if agg_outdir is not None:
        agg_outdir.mkdir(exist_ok=True, parents=True)
    
    # Storage for cross-folder aggregation
    all_collected_entries_by_job: dict[int, list[Any]] = {}
    deferred_tracking_entries_by_job: dict[int, list[Any]] | None = {} if defer_tracking_workbook_refresh else None
    failed_jobs = []
    completed_jobs = []
    aggregation_failed = False
    
    # Thread safety locks
    data_lock = threading.Lock()
    callback_lock = threading.Lock()

    def _extend_entries(
        bucket: dict[int, list[Any]] | None,
        job_index: int,
        entries: list[Any],
    ) -> None:
        if bucket is None or not entries:
            return
        bucket.setdefault(job_index, []).extend(entries)

    def _materialize_entries(bucket: dict[int, list[Any]] | None) -> list[Any]:
        if not bucket:
            return []
        ordered: list[Any] = []
        for job_index in sorted(bucket):
            ordered.extend(bucket[job_index])
        return sorted(ordered, key=_entry_sort_key)

    def _entry_sort_key(entry: Any) -> tuple[str, str, str, str]:
        if not isinstance(entry, dict):
            return ("", "", "", "")
        fsa = entry.get("fsa")
        file_name = str(getattr(fsa, "file_name", "") or entry.get("File") or "")
        return (
            str(entry.get("dit") or entry.get("IdentityKey") or ""),
            file_name,
            str(entry.get("MarkerName") or ""),
            str(entry.get("Assay") or entry.get("assay") or ""),
        )

    def _emit_progress(
        job_name: str,
        phase: str,
        *,
        file_name: str = "",
        files_done: int | None = None,
        files_total: int | None = None,
        jobs_done: int | None = None,
        note: str = "",
        folder_name: str = "",
    ) -> None:
        if progress_callback is None:
            return
        if jobs_done is None:
            with data_lock:
                jobs_done = len(completed_jobs) + len(failed_jobs)
        payload = {
            "folder_name": folder_name,
            "job_name": job_name,
            "phase": phase,
            "file_name": file_name,
            "files_done": None if files_done is None else int(files_done),
            "files_total": None if files_total is None else int(files_total),
            "jobs_done": int(jobs_done or 0),
            "jobs_total": int(total),
            "heartbeat_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "note": note,
        }
        with callback_lock:
            progress_callback(payload)

    def _with_heartbeat(job_name: str, phase: str, label: str, func, *args, **kwargs):
        stop_event = threading.Event()
        started = time.monotonic()
        slow_job_logged = False

        def _heartbeat() -> None:
            nonlocal slow_job_logged
            while not stop_event.wait(30.0):
                elapsed = int(time.monotonic() - started)
                log(f"[BATCH] Still running: {label}")
                note = "heartbeat"
                if elapsed >= 120:
                    note = "slow_job"
                    if not slow_job_logged:
                        log(f"[BATCH] Slow job warning: {job_name} has remained in {phase} for {elapsed}s.")
                        slow_job_logged = True
                _emit_progress(job_name, phase, note=note)

        thread = threading.Thread(target=_heartbeat, daemon=True)
        thread.start()
        try:
            return func(*args, **kwargs)
        finally:
            stop_event.set()
            thread.join(timeout=0.1)
    
    def process_job(i, job):
        job_name = job["name"]
        job_path = job["path"]
        job_files = job["files"]
        job_file_total = len(job_files or [])
        started = time.monotonic()
        job_status = "success"
        
        log(f"[BATCH] Processing job {i+1}/{total}: {job_name}")
        with callback_lock:
            if update_callback:
                update_callback(i, total, job_name, "running")
        _emit_progress(
            job_name,
            "job_start",
            files_done=0,
            files_total=job_file_total,
            note="job_started",
        )
            
        try:
            # Format output params
            resolved_out_folder = out_folder_tmpl.replace("{name}", job_name)
            resolved_html = outfile_html_tmpl.replace("{name}", job_name)
            resolved_excel = excel_name_tmpl.replace("{name}", job_name)
            
            job_type = job.get("type", "pipeline")

            def _job_progress(event: dict) -> None:
                _emit_progress(
                    job_name,
                    event.get("phase", "analyze"),
                    file_name=str(event.get("file_name", "") or ""),
                    files_done=int(event.get("files_done", 0) or 0),
                    files_total=int(event.get("files_total", job_file_total) or job_file_total),
                    note=str(event.get("note", "") or ""),
                )
            
            # Execute based on job type
            if job_type == "pipeline":
                if aggregate_dit_reports:
                    # Collect mode
                    entries = _with_heartbeat(
                        job_name,
                        "collect_entries",
                        f"collecting entries for {job_name}",
                        run_pipeline_job_collect,
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files,
                        chunk_files=(active_analysis != "flt3"),
                        tracking_excel_path=tracking_excel_path if active_analysis == "clonality" else None,
                        progress_callback=_job_progress,
                    )
                    if stream_aggregated_dit and agg_outdir is not None:
                        from core.html_reports import build_dit_html_reports
                        if defer_dit_html_reports:
                            with data_lock:
                                _extend_entries(all_collected_entries_by_job, i, entries)
                            _emit_progress(
                                job_name,
                                "build_report",
                                files_done=job_file_total,
                                files_total=job_file_total,
                                note="build_report_deferred",
                            )
                        else:
                            _emit_progress(
                                job_name,
                                "build_report",
                                files_done=job_file_total,
                                files_total=job_file_total,
                                note="build_report_started",
                            )
                            _with_heartbeat(
                                job_name,
                                "build_report",
                                f"building DIT report for {job_name}",
                                build_dit_html_reports,
                                entries,
                                agg_outdir,
                            )
                        if active_analysis == "clonality":
                            from core.analyses.clonality.tracking_excel import (
                                CLONALITY_TRACKING_FILENAME,
                                update_clonality_tracking_workbook,
                            )

                            if defer_tracking_workbook_refresh and deferred_tracking_entries_by_job is not None:
                                with data_lock:
                                    _extend_entries(deferred_tracking_entries_by_job, i, entries)
                            elif not defer_dit_html_reports:
                                _emit_progress(
                                    job_name,
                                    "write_tracking_excel",
                                    files_done=job_file_total,
                                    files_total=job_file_total,
                                    note="write_tracking_excel_started",
                                )
                                with data_lock:
                                    _with_heartbeat(
                                        job_name,
                                        "write_tracking_excel",
                                        f"updating clonality tracking workbook for {job_name}",
                                        update_clonality_tracking_workbook,
                                        tracking_excel_path
                                        or resolve_analysis_excel_output_path(
                                            "clonality",
                                            agg_outdir,
                                            CLONALITY_TRACKING_FILENAME,
                                            ),
                                        entries,
                                    )
                        if defer_dit_html_reports:
                            log(f"[BATCH] Deferred aggregated DIT report for {job_name} with {len(entries)} entries.")
                        else:
                            log(f"[BATCH] Built aggregated DIT report for {job_name} with {len(entries)} entries.")
                    else:
                        with data_lock:
                            _extend_entries(all_collected_entries_by_job, i, entries)
                            if defer_tracking_workbook_refresh and deferred_tracking_entries_by_job is not None:
                                _extend_entries(deferred_tracking_entries_by_job, i, entries)
                        log(f"[BATCH] Collected {len(entries)} entries from {job_name}.")
                else:
                    # Normal mode
                    _emit_progress(
                        job_name,
                        "collect_entries",
                        files_done=0,
                        files_total=job_file_total,
                        note="job_started",
                    )
                    run_pipeline_job(
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files
                    )
            
            elif job_type == "qc":
                qc_entries: list[dict[str, Any]] = []
                if active_analysis == "clonality" and aggregate_dit_reports:
                    _, qc_entries = _with_heartbeat(
                        job_name,
                        "collect_entries",
                        f"collecting QC entries for {job_name}",
                        run_qc_job,
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_html_name=resolved_html,
                        excel_name=resolved_excel,
                        rules=qc_rules,
                        files=job_files,
                        tracking_excel_path=tracking_excel_path,
                        update_tracking_workbook=False,
                        return_entries=True,
                        progress_callback=_job_progress,
                    )
                    if qc_entries:
                        if defer_tracking_workbook_refresh and deferred_tracking_entries_by_job is not None:
                            with data_lock:
                                _extend_entries(deferred_tracking_entries_by_job, i, qc_entries)
                            log(f"[BATCH] Deferred {len(qc_entries)} tracking entries from QC job {job_name}.")
                        else:
                            from core.analyses.clonality.tracking_excel import (
                                CLONALITY_TRACKING_FILENAME,
                                update_clonality_tracking_workbook,
                            )

                            _emit_progress(
                                job_name,
                                "write_tracking_excel",
                                files_done=job_file_total,
                                files_total=job_file_total,
                                note="write_tracking_excel_started",
                            )
                            with data_lock:
                                _with_heartbeat(
                                    job_name,
                                    "write_tracking_excel",
                                    f"updating clonality tracking workbook for {job_name}",
                                    update_clonality_tracking_workbook,
                                    tracking_excel_path
                                    or resolve_analysis_excel_output_path(
                                        "clonality",
                                        output_base,
                                        CLONALITY_TRACKING_FILENAME,
                                    ),
                                    qc_entries,
                                )
                        log(f"[BATCH] Collected {len(qc_entries)} tracking entries from QC job {job_name}.")
                else:
                    run_qc_job(
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_html_name=resolved_html,
                        excel_name=resolved_excel,
                        rules=qc_rules,
                        files=job_files,
                    )
                
            elif job_type == "dit":
                if aggregate_dit_reports:
                    entries = _with_heartbeat(
                        job_name,
                        "collect_entries",
                        f"collecting DIT entries for {job_name}",
                        run_pipeline_job_collect,
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files,
                        chunk_files=(active_analysis != "flt3"),
                        tracking_excel_path=tracking_excel_path if active_analysis == "clonality" else None,
                        progress_callback=_job_progress,
                    )
                    if stream_aggregated_dit and agg_outdir is not None:
                        from core.html_reports import build_dit_html_reports
                        if defer_dit_html_reports:
                            with data_lock:
                                _extend_entries(all_collected_entries_by_job, i, entries)
                            _emit_progress(
                                job_name,
                                "build_report",
                                files_done=job_file_total,
                                files_total=job_file_total,
                                note="build_report_deferred",
                            )
                        else:
                            _emit_progress(
                                job_name,
                                "build_report",
                                files_done=job_file_total,
                                files_total=job_file_total,
                                note="build_report_started",
                            )
                            _with_heartbeat(
                                job_name,
                                "build_report",
                                f"building aggregated DIT report for {job_name}",
                                build_dit_html_reports,
                                entries,
                                agg_outdir,
                            )
                        if active_analysis == "clonality":
                            from core.analyses.clonality.tracking_excel import (
                                CLONALITY_TRACKING_FILENAME,
                                update_clonality_tracking_workbook,
                            )

                            if defer_tracking_workbook_refresh and deferred_tracking_entries_by_job is not None:
                                with data_lock:
                                    _extend_entries(deferred_tracking_entries_by_job, i, entries)
                            elif not defer_dit_html_reports:
                                _emit_progress(
                                    job_name,
                                    "write_tracking_excel",
                                    files_done=job_file_total,
                                    files_total=job_file_total,
                                    note="write_tracking_excel_started",
                                )
                                with data_lock:
                                    _with_heartbeat(
                                        job_name,
                                        "write_tracking_excel",
                                        f"updating clonality tracking workbook for {job_name}",
                                        update_clonality_tracking_workbook,
                                        tracking_excel_path
                                        or resolve_analysis_excel_output_path(
                                            "clonality",
                                            agg_outdir,
                                            CLONALITY_TRACKING_FILENAME,
                                        ),
                                        entries,
                                    )
                        if defer_dit_html_reports:
                            log(f"[BATCH] Deferred aggregated DIT report for {job_name} with {len(entries)} entries.")
                        else:
                            log(f"[BATCH] Built aggregated DIT report for {job_name} with {len(entries)} entries.")
                    else:
                        with data_lock:
                            _extend_entries(all_collected_entries_by_job, i, entries)
                            if defer_tracking_workbook_refresh and deferred_tracking_entries_by_job is not None:
                                _extend_entries(deferred_tracking_entries_by_job, i, entries)
                        log(f"[BATCH] Collected {len(entries)} entries from {job_name} for DIT aggregation.")
                else:
                    run_dit_job(
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files
                    )
            
            else:
                raise ValueError(f"Unknown job_type: {job_type}")
                    
            with data_lock:
                completed_jobs.append(job_name)
            with callback_lock:
                if update_callback:
                    update_callback(i + 1, total, job_name, "success")
            _emit_progress(
                job_name,
                "done",
                files_done=job_file_total,
                files_total=job_file_total,
                note="job_complete",
            )
            job_status = "success"
                
        except Exception as e:
            log(f"[ERROR] Job '{job_name}' failed: {e}")
            job_status = "failed"
            with data_lock:
                failed_jobs.append(job_name)
            with callback_lock:
                if update_callback:
                    update_callback(i + 1, total, job_name, f"error: {e}")
            _emit_progress(
                job_name,
                "failed",
                files_done=job_file_total,
                files_total=job_file_total,
                note=str(e),
            )
            if not continue_on_error:
                log("[BATCH] Stopping batch due to error.")
                raise  # Stop early by propagating error to executor
        finally:
            log(f"[BATCH] Job {job_name} finished in {time.monotonic() - started:.1f}s ({job_status}).")

    # Perform multi-threaded patient processing
    # Max workers = 3 (modest to prevent over-subscription of child Pool processes)
    max_patient_workers = int(max_workers) if max_workers is not None else min(3, max(1, os.cpu_count() // 2 or 1))
    max_patient_workers = max(1, max_patient_workers)
    
    try:
        with ThreadPoolExecutor(max_workers=max_patient_workers) as executor:
            futures = [
                executor.submit(process_job, i, job)
                for i, job in enumerate(jobs)
            ]
            from concurrent.futures import as_completed
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as ex:
                    if not continue_on_error:
                        raise
                    log(f"[BATCH] Worker future failed: {ex}")
    except Exception as ex:
        if not continue_on_error:
            log(f"[BATCH] Batch execution halted after error: {ex}")
        else:
            log(f"[BATCH] Batch execution completed with some errors.")
                
    # --- CROSS-FOLDER DIT AGGREGATION ---
    all_collected_entries = _materialize_entries(all_collected_entries_by_job)
    deferred_tracking_entries = _materialize_entries(deferred_tracking_entries_by_job)

    if aggregate_dit_reports and all_collected_entries and not defer_dit_html_reports:
        log("\n[BATCH] Final step: Building aggregated DIT HTML reports...")
        
        from core.html_reports import build_dit_html_reports
        
        try:
            if (not stream_aggregated_dit) or defer_dit_html_reports:
                _with_heartbeat(
                    "DIT aggregation",
                    "build_report",
                    "building final aggregated DIT reports",
                    build_dit_html_reports,
                    all_collected_entries,
                    agg_outdir,
                )
            if active_analysis == "clonality" and not defer_tracking_workbook_refresh:
                from core.analyses.clonality.tracking_excel import (
                    CLONALITY_TRACKING_FILENAME,
                    update_clonality_tracking_workbook,
                )

                _with_heartbeat(
                    "DIT aggregation",
                    "write_tracking_excel",
                    "updating final clonality tracking workbook",
                    update_clonality_tracking_workbook,
                    tracking_excel_path
                    or resolve_analysis_excel_output_path(
                        "clonality",
                        agg_outdir,
                        CLONALITY_TRACKING_FILENAME,
                    ),
                    all_collected_entries,
                )
            log(f"[BATCH] Successfully built aggregated DIT reports in {agg_outdir}")
        except Exception as e:
            aggregation_failed = True
            failed_jobs.append("DIT aggregation")
            log(f"[ERROR] Failed to build aggregated DIT reports: {e}")
    elif aggregate_dit_reports and all_collected_entries and defer_dit_html_reports:
        log("[BATCH] Skipped aggregated DIT HTML report generation; entries were retained for deferred handling.")
    elif aggregate_dit_reports and stream_aggregated_dit:
        log("[BATCH] Aggregated DIT reports were streamed job-by-job to reduce memory pressure.")

    if aggregation_failed:
        log("[BATCH] Batch run complete with aggregation errors.")
    else:
        log("[BATCH] Batch run complete.")
    if update_callback:
        update_callback(total, total, "Done", "done")
    return {
        "total_jobs": total,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
        "collected_entries": deferred_tracking_entries if defer_tracking_workbook_refresh else all_collected_entries,
    }
