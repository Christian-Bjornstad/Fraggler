"""
Fraggler Diagnostics — Rust Engine Bridge.

Provides a hybrid mode where the fast Rust engine is used to detect
and fit the size standard peaks, while Python maintains the rest of
the pipeline for full compatibility with existing Plotly HTML reports
and QC log tracking.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from core.log import log
from fraggler.fraggler import FsaFile, fit_size_standard_to_ladder


ROX_PREFERRED_TIME_MIN = 1500.0
ROX_PREFERRED_TIME_MAX = 4000.0
ROX_HARD_TIME_MIN = 1300.0
ROX_HARD_TIME_MAX = 4300.0
ROX_MAX_FIRST_ANCHOR = 1900.0
ROX_MIN_SPAN = 1100.0
ROX_MIN_MEDIAN_GAP = 26.0
ROX_MIN_HARD_WINDOW_FRACTION = 0.75

LIZ_HARD_TIME_MIN = 1150.0
LIZ_HARD_TIME_MAX = 4300.0
LIZ_MAX_FIRST_ANCHOR = 1700.0
LIZ_MIN_SPAN = 900.0
LIZ_MIN_MEDIAN_GAP = 22.0
LIZ_MIN_HARD_WINDOW_FRACTION = 0.80


def _is_rox_ladder(fsa: FsaFile, expected_bps: list[float]) -> bool:
    ladder_name = str(getattr(fsa, "ladder", "") or "").upper()
    if "ROX" in ladder_name:
        return True
    return len(expected_bps) >= 20


def _anchor_intensity(trace: np.ndarray, scan_idx: int) -> float:
    if trace.size == 0:
        return float("nan")
    idx = int(np.clip(scan_idx, 0, trace.size - 1))
    return float(trace[idx])


def _validate_rust_anchor_selection(
    fsa: FsaFile,
    scan_indices: list[int],
    expected_bps: list[float],
) -> tuple[bool, str]:
    if not scan_indices or len(scan_indices) != len(expected_bps):
        return False, "scan/step length mismatch"

    scans = np.asarray(scan_indices, dtype=float)
    if scans.size < 3:
        return False, "too few anchor points"
    if np.any(np.diff(scans) <= 0):
        return False, "anchors are not strictly increasing"

    span = float(scans[-1] - scans[0])
    gaps = np.diff(scans)
    median_gap = float(np.median(gaps)) if gaps.size else 0.0

    size_standard = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    anchor_signal = np.asarray([_anchor_intensity(size_standard, int(v)) for v in scan_indices], dtype=float)
    median_signal = float(np.nanmedian(anchor_signal)) if anchor_signal.size else float("nan")

    is_rox = _is_rox_ladder(fsa, expected_bps)
    if is_rox:
        in_hard = np.logical_and(scans >= ROX_HARD_TIME_MIN, scans <= ROX_HARD_TIME_MAX)
        hard_fraction = float(np.mean(in_hard)) if scans.size else 0.0
        if hard_fraction < ROX_MIN_HARD_WINDOW_FRACTION:
            return False, f"ROX anchors mostly outside expected time window ({hard_fraction:.2f})"
        if scans[0] > ROX_MAX_FIRST_ANCHOR:
            return False, f"ROX first anchor too late ({scans[0]:.0f})"
        if not (ROX_PREFERRED_TIME_MIN <= float(np.median(scans)) <= ROX_PREFERRED_TIME_MAX):
            return False, f"ROX median anchor outside preferred window ({np.median(scans):.0f})"
        if span < ROX_MIN_SPAN:
            return False, f"ROX anchor span too small ({span:.0f})"
        if median_gap < ROX_MIN_MEDIAN_GAP:
            return False, f"ROX anchors too tightly clustered (median gap {median_gap:.1f})"
        if np.isfinite(median_signal) and median_signal < 45.0:
            return False, f"ROX anchor signal too weak (median {median_signal:.1f})"
        return True, "ROX anchor checks passed"

    in_hard = np.logical_and(scans >= LIZ_HARD_TIME_MIN, scans <= LIZ_HARD_TIME_MAX)
    hard_fraction = float(np.mean(in_hard)) if scans.size else 0.0
    if hard_fraction < LIZ_MIN_HARD_WINDOW_FRACTION:
        return False, f"LIZ anchors mostly outside expected time window ({hard_fraction:.2f})"
    if scans[0] > LIZ_MAX_FIRST_ANCHOR:
        return False, f"LIZ first anchor too late ({scans[0]:.0f})"
    if span < LIZ_MIN_SPAN:
        return False, f"LIZ anchor span too small ({span:.0f})"
    if median_gap < LIZ_MIN_MEDIAN_GAP:
        return False, f"LIZ anchors too tightly clustered (median gap {median_gap:.1f})"
    if np.isfinite(median_signal) and median_signal < 35.0:
        return False, f"LIZ anchor signal too weak (median {median_signal:.1f})"
    return True, "LIZ anchor checks passed"

def run_ladder_fit_hybrid(fsa: FsaFile, analysis_kind: str) -> FsaFile | None:
    """
    Passes the FSA file to the fraggler-cli to perform baseline correction,
    peak detection, and ladder fitting. Retrieves the mapped ladder steps 
    and applies them directly to the Python FsaFile.
    """
    root = Path(__file__).resolve().parent.parent
    
    # Resolve the CLI binary
    if getattr(sys, 'frozen', False):
        cli_bin = Path(sys._MEIPASS) / "fraggler-cli"
        if not cli_bin.exists():
            cli_bin = Path(sys.executable).parent / "fraggler-cli"
    else:
        search_paths = [
            root / "fraggler-v2" / "target" / "release" / "fraggler-cli",
            root / "fraggler-v2" / "target" / "debug" / "fraggler-cli",
            root / "bin" / "fraggler-cli",
        ]
        cli_bin = next((p for p in search_paths if p.exists()), None)

    if not cli_bin or not cli_bin.exists():
        log("[RUST ERROR] Could not find fraggler-cli. Falling back to Python engine.")
        return None

    fsa_path = Path(fsa.file)
    with tempfile.TemporaryDirectory(prefix="fraggler_hybrid_") as tdir:
        tdir_path = Path(tdir)
        from config import APP_SETTINGS
        skip_html_reports = APP_SETTINGS.get("engine", {}).get("skip_html_reports", False)

        req = {
            "contract_version": {"major": 1, "minor": 0},
            "run_kind": "analyze",
            "analysis_kind": analysis_kind,
            "correlation_id": "00000000-0000-0000-0000-000000000000",
            "inputs": {
                "paths": [str(fsa_path)],
                "manifest_path": None,
                "report_source_path": None
            },
            "output": {
                "root_dir": str(tdir_path),
                "report_dir": None,
                "artifacts_dir": None
            },
            "options": {
                "max_workers": 1,
                "deterministic": True,
                "emit_compact_json": True,
                "open_reports_in_browser": False,
                "shadow_reference_python": False,
                "skip_html_reports": skip_html_reports,
                "extra": {}
            }
        }
        
        req_path = tdir_path / "req.json"
        with open(req_path, "w") as f:
            json.dump(req, f)

        cmd = [str(cli_bin), "analyze", "--json-request", str(req_path)]
        import time
        start_time = time.monotonic()
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60, stdin=subprocess.DEVNULL)
            elapsed = time.monotonic() - start_time
            log(f"[RUST] Engine finished in {elapsed:.1f}s for {fsa.file_name}")
        except subprocess.TimeoutExpired as e:
            log(f"[RUST ERROR] CLI timed out after 60s for {fsa.file_name}. Stderr: {e.stderr}")
            return None
        except subprocess.CalledProcessError as e:
            log(f"[RUST ERROR] CLI failed with code {e.returncode} for {fsa.file_name}. Stderr: {e.stderr}")
            return None

        summary_path = tdir_path / "analyze_summary.json"
        if not summary_path.exists():
            log("[RUST ERROR] Missing analyze_summary.json from Rust engine.")
            return None

        try:
            with open(summary_path, "r", encoding="utf-8") as f:
                results = json.load(f)
        except Exception as e:
            log(f"[RUST ERROR] Failed to parse analyze_summary.json: {e}")
            return None
            
        if not results:
            log("[RUST ERROR] results is empty")
            return None
            
        res = results[0]
        fit_preview = res.get("ladder_fit_preview")
        if not fit_preview:
            log("[RUST ERROR] ladder_fit_preview missing")
            return None

        # Extract the scan indices chosen by Rust
        refinement = fit_preview.get("refinement")
        if refinement and refinement.get("refined_scan_indices"):
            scan_indices = refinement["refined_scan_indices"]
        else:
            scan_indices = fit_preview.get("best_scan_indices", [])
            
        if not scan_indices:
            log(f"[RUST ERROR] scan_indices empty for {fsa.file_name}")
            if res and isinstance(res, dict) and res.get("stderr"):
                log(f"[RUST DIAG] Stderr: {res.get('stderr').strip()}")
            return None

        # Extract the expected ladder basepairs 
        model_preview = fit_preview.get("sizing_model")
        if not model_preview or not model_preview.get("predicted_ladder_basepairs"):
            log("[RUST ERROR] model_preview missing predicted_ladder_basepairs")
            return None
            
        expected_bps = model_preview["predicted_ladder_basepairs"]
        
        if len(scan_indices) != len(expected_bps):
            log(f"[RUST ERROR] Mismatch between selected scan indices ({len(scan_indices)}) and predicted basepairs ({len(expected_bps)}).")
            return None

        ok, reason = _validate_rust_anchor_selection(fsa, scan_indices, expected_bps)
        if not ok:
            log(
                f"[RUST GUARDRAIL] Rejected anchor set for {fsa.file_name}: {reason}. "
                "Falling back to Python ladder search."
            )
            return None

        # Transfer to the Python FsaFile
        fsa.best_size_standard = np.array(scan_indices, dtype=float)
        fsa.ladder_steps = np.array(expected_bps, dtype=float)
        fsa.expected_ladder_steps = np.array(expected_bps, dtype=float)

        # Let the standard python function map traces using these chosen anchor points
        try:
            fsa = fit_size_standard_to_ladder(fsa)
        except Exception as e:
            log(f"[HYBRID ERROR] fit_size_standard_to_ladder failed: {e}")
            return None
            
        return fsa
