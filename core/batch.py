"""
Fraggler Diagnostics — Batch Processing

Logic for automated scanning of subfolders / YAML files and executing
jobs (pipeline/qc/dit). Also implements cross-folder DIT aggregation.
Compatible with Python 3.10+.
"""
from __future__ import annotations

import os
import re
import yaml
from pathlib import Path
from typing import Dict, List, Any

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
        with open(yaml_path, "r", encoding="utf-8") as f:
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

def find_all_fsa_files(folders: List[Path]) -> List[Path]:
    """Find all .fsa files in the given folders."""
    from core.utils import is_water_file

    fsa_files = []
    for f in folders:
        if f.is_dir():
            fsa_files.extend(
                path for path in f.glob("*.fsa")
                if not is_water_file(path.name)
            )
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
    active_analysis = APP_SETTINGS.get("active_analysis", "clonality")
    for p in input_paths:
        if p.is_file() and p.name.endswith(".yaml"):
            folders_to_scan.extend(scan_jobs_from_yaml(p))
            continue
            
        if any(p.glob("*.fsa")):
            folders_to_scan.append(p)
            
        if p.is_dir():
            subfolders = [d for d in p.iterdir() if d.is_dir()]
            for sub in subfolders:
                if any(sub.glob("*.fsa")):
                    folders_to_scan.append(sub)

    folders_to_scan = list(dict.fromkeys(folders_to_scan))
        
    if not folders_to_scan:
        log(f"[WARN] No folders with .fsa data found.")
        return []

    if active_analysis == "general":
        jobs = []
        for folder in folders_to_scan:
            jobs.append({
                "name": folder.name,
                "type": "pipeline",
                "path": folder,
                "files": [],
            })
        if jobs:
            log(f"[INFO] General analysis: prepared {len(jobs)} folder job(s) without patient aggregation.")
        return jobs
        
    jobs = []
    
    if aggregate_patients:
        all_fsa = find_all_fsa_files(folders_to_scan)
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
        from core.utils import CONTROL_PREFIX_RE, is_water_file, strip_stage_prefix
        qc_pattern = CONTROL_PREFIX_RE
        all_qc_files = []
        for folder in folders_to_scan:
            fsa_files = [
                f for f in folder.glob("*.fsa")
                if not is_water_file(f.name)
            ]
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
    update_callback: Any = None
) -> None:
    """
    Run all generated jobs.
    Reads QC parameters from APP_SETTINGS and uses explicit per-analysis batch flags.
    """
    from config import APP_SETTINGS
    from core.qc.qc_rules import QCRules
    s_qc = APP_SETTINGS.get("qc", {})
    active_analysis = APP_SETTINGS.get("active_analysis", "clonality")
    aggregate_dit_reports = bool(aggregate_dit_reports) and active_analysis != "general"
    qc_rules = QCRules(
        min_r2_ok=s_qc.get("min_r2_ok", 0.999),
        min_r2_warn=s_qc.get("min_r2_warn", 0.995),
        sample_peak_window_bp=s_qc.get("sample_peak_window_bp", 2.0)
    )
    
    total = len(jobs)
    log(f"[BATCH] Starting batch run of {total} jobs.")
    
    # Storage for cross-folder aggregation
    all_collected_entries = []
    failed_jobs = []
    completed_jobs = []
    
    for i, job in enumerate(jobs):
        job_name = job["name"]
        job_path = job["path"]
        job_files = job["files"]
        
        log(f"[BATCH] Processing job {i+1}/{total}: {job_name}")
        if update_callback:
            update_callback(i, total, job_name, "running")
            
        try:
            # Format output params
            resolved_out_folder = out_folder_tmpl.replace("{name}", job_name)
            resolved_html = outfile_html_tmpl.replace("{name}", job_name)
            resolved_excel = excel_name_tmpl.replace("{name}", job_name)
            
            job_type = job.get("type", "pipeline")
            
            # Execute based on job type
            if job_type == "pipeline":
                if aggregate_dit_reports:
                    # Collect mode: run pipeline but don't build DIT reports yet
                    entries = run_pipeline_job_collect(
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files,
                        chunk_files=(active_analysis != "flt3"),
                    )
                    all_collected_entries.extend(entries)
                    log(f"[BATCH] Collected {len(entries)} entries from {job_name}.")
                else:
                    # Normal mode: run pipeline and build DIT reports per-job
                    run_pipeline_job(
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files
                    )
            
            elif job_type == "qc":
                run_qc_job(
                    fsa_dir=job_path,
                    base_outdir=output_base,
                    out_html_name=resolved_html,
                    excel_name=resolved_excel,
                    rules=qc_rules,
                    files=job_files
                )
                
            elif job_type == "dit":
                if aggregate_dit_reports:
                    entries = run_pipeline_job_collect(
                        fsa_dir=job_path,
                        base_outdir=output_base,
                        out_folder_name=resolved_out_folder,
                        scope=pipeline_scope,
                        needle=assay_filter,
                        files=job_files,
                        chunk_files=(active_analysis != "flt3"),
                    )
                    all_collected_entries.extend(entries)
                    log(f"[BATCH] Collected {len(entries)} entries from {job_name} for DIT aggregation.")
                else:
                    # Explicit standalone DIT job (just rebuilds reports per folder)
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
                    
            if update_callback:
                update_callback(i + 1, total, job_name, "success")
            completed_jobs.append(job_name)
                
        except Exception as e:
            log(f"[ERROR] Job '{job_name}' failed: {e}")
            failed_jobs.append(job_name)
            if update_callback:
                update_callback(i + 1, total, job_name, f"error: {e}")
            if not continue_on_error:
                log("[BATCH] Stopping batch due to error.")
                break
                
    # --- CROSS-FOLDER DIT AGGREGATION ---
    if aggregate_dit_reports and all_collected_entries:
        log("\n[BATCH] Final step: Building aggregated DIT HTML reports...")
        
        # The reports will be placed in base_outdir / OUTDIR_NAME (default)
        # We'll just define a generic container for them
        from core.html_reports import build_dit_html_reports
        from core.assay_config import OUTDIR_NAME
        agg_outdir = output_base / OUTDIR_NAME
        agg_outdir.mkdir(exist_ok=True, parents=True)
        
        try:
            build_dit_html_reports(all_collected_entries, agg_outdir)
            log(f"[BATCH] Successfully built aggregated DIT reports in {agg_outdir}")
        except Exception as e:
            log(f"[ERROR] Failed to build aggregated DIT reports: {e}")

    log("[BATCH] Batch run complete.")
    if update_callback:
        update_callback(total, total, "Done", "done")
    return {
        "total_jobs": total,
        "completed_jobs": completed_jobs,
        "failed_jobs": failed_jobs,
    }
