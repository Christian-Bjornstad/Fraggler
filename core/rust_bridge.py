"""
HemaFrag Diagnostics — Rust Engine Bridge.

Provides a hybrid mode where the fast Rust engine is used to detect
and fit the size standard peaks, while Python maintains the rest of
the pipeline for full compatibility with existing Plotly HTML reports
and QC log tracking.
"""
from __future__ import annotations

import json
import os
import select
import subprocess
import sys
import tempfile
import threading
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.log import log
from fraggler.fraggler import FsaFile, baseline_arPLS, fit_size_standard_to_ladder


ROX_PREFERRED_TIME_MIN = 1500.0
ROX_PREFERRED_TIME_MAX = 4000.0
ROX_HARD_TIME_MIN = 1300.0
ROX_HARD_TIME_MAX = 4300.0
ROX_MAX_FIRST_ANCHOR = 1900.0
ROX_MIN_SPAN = 1100.0
ROX_MIN_MEDIAN_GAP = 26.0
ROX_MIN_HARD_WINDOW_FRACTION = 0.75

GS500ROX_PREFERRED_TIME_MIN = 1400.0
GS500ROX_PREFERRED_TIME_MAX = 4200.0
GS500ROX_HARD_TIME_MIN = 1180.0
GS500ROX_HARD_TIME_MAX = 4550.0
GS500ROX_MAX_FIRST_ANCHOR = 1700.0
GS500ROX_MIN_SPAN = 2500.0
GS500ROX_MIN_MEDIAN_GAP = 36.0
GS500ROX_MIN_HARD_WINDOW_FRACTION = 0.60

LIZ_HARD_TIME_MIN = 1150.0
LIZ_HARD_TIME_MAX = 4300.0
LIZ_MAX_FIRST_ANCHOR = 1700.0
LIZ_MIN_SPAN = 900.0
LIZ_MIN_MEDIAN_GAP = 22.0
LIZ_MIN_HARD_WINDOW_FRACTION = 0.80

_CLI_BIN_CACHE: Path | None = None
_RUST_WORKER: "_RustPrimitiveWorker | None" = None
_RUST_WORKER_LOCK = threading.Lock()


class _RustSizingModel:
    def __init__(
        self,
        *,
        strategy: str,
        coefficients: list[float],
        scan_indices: np.ndarray,
        ladder_steps: np.ndarray,
    ) -> None:
        self.strategy = str(strategy or "")
        self.coefficients = np.asarray(coefficients, dtype=float)
        self.scan_indices = np.asarray(scan_indices, dtype=float)
        self.ladder_steps = np.asarray(ladder_steps, dtype=float)

    def predict(self, x_values: np.ndarray) -> np.ndarray:
        x_array = np.asarray(x_values, dtype=float).reshape(-1)
        if self.strategy == "willros_monotone_spline":
            return np.asarray(
                [
                    _eval_monotone_cubic_spline(
                        self.scan_indices,
                        self.ladder_steps,
                        self.coefficients,
                        float(xq),
                    )
                    for xq in x_array
                ],
                dtype=float,
            )
        if self.strategy == "polynomial_fallback":
            return np.asarray(
                [_eval_polynomial(self.coefficients, float(xq)) for xq in x_array],
                dtype=float,
            )
        raise ValueError(f"Unsupported Rust sizing strategy: {self.strategy}")


def _eval_polynomial(coefficients: np.ndarray, x_value: float) -> float:
    return float(
        sum(float(coefficient) * (x_value ** power) for power, coefficient in enumerate(coefficients))
    )


def _eval_monotone_cubic_spline(
    x: np.ndarray,
    y: np.ndarray,
    tangents: np.ndarray,
    x_query: float,
) -> float:
    if x.size == 1:
        return float(y[0])
    if x_query <= float(x[0]):
        return float(y[0] + tangents[0] * (x_query - x[0]))
    if x_query >= float(x[-1]):
        return float(y[-1] + tangents[-1] * (x_query - x[-1]))

    upper = int(np.searchsorted(x, x_query, side="right"))
    index = min(max(upper - 1, 0), x.size - 2)
    step = float(x[index + 1] - x[index])
    if step <= 0.0:
        return float(y[index])
    t = (x_query - float(x[index])) / step
    t2 = t * t
    t3 = t2 * t

    h00 = 2.0 * t3 - 3.0 * t2 + 1.0
    h10 = t3 - 2.0 * t2 + t
    h01 = -2.0 * t3 + 3.0 * t2
    h11 = t3 - t2

    return float(
        h00 * y[index]
        + h10 * step * tangents[index]
        + h01 * y[index + 1]
        + h11 * step * tangents[index + 1]
    )


def _apply_rust_sizing_model_to_fsa(
    fsa: FsaFile,
    scan_indices: list[int],
    expected_bps: list[float],
    model_preview: dict[str, Any],
) -> FsaFile | None:
    strategy = str(model_preview.get("strategy") or "")
    coefficients = model_preview.get("coefficients")
    if strategy not in {"willros_monotone_spline", "polynomial_fallback"} or not isinstance(coefficients, list):
        return None

    scan_array = np.asarray(scan_indices, dtype=float)
    ladder_array = np.asarray(expected_bps, dtype=float)
    model = _RustSizingModel(
        strategy=strategy,
        coefficients=[float(value) for value in coefficients],
        scan_indices=scan_array,
        ladder_steps=ladder_array,
    )

    sample_trace = np.asarray(getattr(fsa, "sample_data", []), dtype=float)
    if sample_trace.size == 0:
        return None

    time_values = np.arange(sample_trace.size, dtype=float)
    basepairs = model.predict(time_values).round(2)
    sample_df = (
        pd.DataFrame({"time": time_values.astype(int), "peaks": sample_trace, "basepairs": basepairs})
        .loc[lambda df: df["basepairs"] >= 0]
        .reset_index(drop=True)
    )
    if sample_df.empty:
        return None

    fsa.ladder_model = model
    fsa.sample_data_with_basepairs = sample_df
    fsa.fitted_to_model = True
    rust_qc = model_preview.get("qc_metrics")
    if isinstance(rust_qc, dict):
        fsa.rust_ladder_qc_metrics = rust_qc
    return fsa


def _is_rox_ladder(fsa: FsaFile, expected_bps: list[float]) -> bool:
    ladder_name = str(getattr(fsa, "ladder", "") or "").upper()
    if "ROX" in ladder_name:
        return True
    return len(expected_bps) >= 20


def _is_gs500rox_ladder(fsa: FsaFile, expected_bps: list[float]) -> bool:
    ladder_name = str(getattr(fsa, "ladder", "") or "").upper()
    analysis_id = str(getattr(fsa, "analysis_id", "") or "").lower()
    if "GS500ROX" in ladder_name:
        return True
    if analysis_id == "flt3":
        return True
    if len(expected_bps) == 16:
        rounded = {int(round(float(value))) for value in expected_bps}
        return 500 in rounded and 490 in rounded and 35 in rounded
    return False


def _resolve_cli_bin() -> Path | None:
    global _CLI_BIN_CACHE
    if _CLI_BIN_CACHE is not None and _CLI_BIN_CACHE.exists():
        return _CLI_BIN_CACHE

    root = Path(__file__).resolve().parent.parent
    if getattr(sys, 'frozen', False):
        cli_bin = Path(sys._MEIPASS) / "fraggler-cli"
        if not cli_bin.exists():
            cli_bin = Path(sys.executable).parent / "fraggler-cli"
        if cli_bin.exists():
            _CLI_BIN_CACHE = cli_bin
            return cli_bin
        return None

    preferred_paths = [
        root / "fraggler-v2" / "target" / "release" / "fraggler-cli",
        root / "fraggler-v2" / "target" / "debug" / "fraggler-cli",
        root / "bin" / "fraggler-cli",
    ]
    cli_bin = next((p for p in preferred_paths if p.exists()), None)
    if cli_bin is not None:
        _CLI_BIN_CACHE = cli_bin
    return cli_bin


class _RustPrimitiveWorker:
    def __init__(self, cli_bin: Path) -> None:
        self.cli_bin = cli_bin
        self._proc = subprocess.Popen(
            [str(cli_bin), "serve-primitives", "--log-filter", "error"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        self._lock = threading.Lock()

    def close(self) -> None:
        proc = self._proc
        if proc.poll() is None:
            try:
                proc.terminate()
                proc.wait(timeout=1)
            except Exception:
                try:
                    proc.kill()
                except Exception:
                    pass

    def request(
        self,
        fsa_path: Path,
        analysis_kind: str,
        timeout_seconds: int,
    ) -> dict[str, Any] | None:
        if self._proc.poll() is not None:
            return None
        if self._proc.stdin is None or self._proc.stdout is None:
            return None

        payload = {
            "input": str(fsa_path),
            "analysis": str(analysis_kind or "").lower() or None,
        }

        with self._lock:
            if self._proc.poll() is not None:
                return None
            try:
                self._proc.stdin.write(json.dumps(payload) + "\n")
                self._proc.stdin.flush()
            except Exception:
                return None

            fd = self._proc.stdout.fileno()
            ready, _, _ = select.select([fd], [], [], max(timeout_seconds, 1))
            if not ready:
                return {"error": f"worker timeout after {timeout_seconds}s"}

            line = self._proc.stdout.readline()
            if not line:
                stderr = ""
                if self._proc.stderr is not None:
                    try:
                        stderr = self._proc.stderr.read()[-1000:]
                    except Exception:
                        stderr = ""
                return {"error": f"worker closed unexpectedly: {stderr.strip()}"}
            try:
                response = json.loads(line)
            except Exception as exc:
                return {"error": f"invalid worker response: {exc}"}
            return response


def _get_rust_worker() -> _RustPrimitiveWorker | None:
    global _RUST_WORKER
    with _RUST_WORKER_LOCK:
        if _RUST_WORKER is not None and _RUST_WORKER._proc.poll() is None:
            return _RUST_WORKER
        cli_bin = _resolve_cli_bin()
        if cli_bin is None or not cli_bin.exists():
            return None
        _RUST_WORKER = _RustPrimitiveWorker(cli_bin)
        return _RUST_WORKER


def _invalidate_rust_worker() -> None:
    global _RUST_WORKER
    with _RUST_WORKER_LOCK:
        if _RUST_WORKER is not None:
            _RUST_WORKER.close()
            _RUST_WORKER = None


def _rust_timeout_seconds(analysis_kind: str) -> int:
    from config import APP_SETTINGS

    engine_settings = APP_SETTINGS.get("engine", {})
    kind = str(analysis_kind or "").lower()
    if kind == "rox":
        return int(engine_settings.get("rust_timeout_seconds_rox", engine_settings.get("rust_timeout_seconds", 60)))
    if kind == "liz":
        return int(engine_settings.get("rust_timeout_seconds_liz", engine_settings.get("rust_timeout_seconds", 60)))
    return int(engine_settings.get("rust_timeout_seconds", 60))


def _anchor_intensity(trace: np.ndarray, scan_idx: int) -> float:
    if trace.size == 0:
        return float("nan")
    idx = int(np.clip(scan_idx, 0, trace.size - 1))
    return float(trace[idx])


def _baseline_correct_for_validation(trace: np.ndarray) -> np.ndarray:
    if trace.size == 0:
        return trace
    try:
        baseline = np.asarray(baseline_arPLS(trace), dtype=float)
        if baseline.shape != trace.shape:
            return np.maximum(trace, 0.0)
        return np.maximum(trace - baseline, 0.0)
    except Exception:
        return np.maximum(trace, 0.0)


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
    validation_trace = _baseline_correct_for_validation(size_standard)
    anchor_signal = np.asarray(
        [_anchor_intensity(validation_trace, int(v)) for v in scan_indices],
        dtype=float,
    )
    median_signal = float(np.nanmedian(anchor_signal)) if anchor_signal.size else float("nan")

    is_gs500rox = _is_gs500rox_ladder(fsa, expected_bps)
    is_rox = _is_rox_ladder(fsa, expected_bps)
    if is_gs500rox:
        in_hard = np.logical_and(scans >= GS500ROX_HARD_TIME_MIN, scans <= GS500ROX_HARD_TIME_MAX)
        hard_fraction = float(np.mean(in_hard)) if scans.size else 0.0
        if hard_fraction < GS500ROX_MIN_HARD_WINDOW_FRACTION:
            return False, f"GS500ROX anchors mostly outside expected time window ({hard_fraction:.2f})"
        if scans[0] > GS500ROX_MAX_FIRST_ANCHOR:
            return False, f"GS500ROX first anchor too late ({scans[0]:.0f})"
        if not (GS500ROX_PREFERRED_TIME_MIN <= float(np.median(scans)) <= GS500ROX_PREFERRED_TIME_MAX):
            return False, f"GS500ROX median anchor outside preferred window ({np.median(scans):.0f})"
        if span < GS500ROX_MIN_SPAN:
            return False, f"GS500ROX anchor span too small ({span:.0f})"
        if median_gap < GS500ROX_MIN_MEDIAN_GAP:
            return False, f"GS500ROX anchors too tightly clustered (median gap {median_gap:.1f})"
        if np.isfinite(median_signal) and median_signal < 20.0:
            return False, f"GS500ROX anchor signal too weak (median {median_signal:.1f})"
        return True, "GS500ROX anchor checks passed"

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
    cli_bin = _resolve_cli_bin()
    if not cli_bin or not cli_bin.exists():
        log("[RUST ERROR] Could not find fraggler-cli. Rust runtime analysis cannot continue.")
        return None

    fsa_path = Path(fsa.file)
    timeout_seconds = max(_rust_timeout_seconds(analysis_kind), 1)

    worker = _get_rust_worker()
    if worker is not None:
        worker_response = worker.request(fsa_path, analysis_kind, timeout_seconds)
        if worker_response and worker_response.get("ok") and worker_response.get("result"):
            res = worker_response["result"]
        else:
            if worker_response and worker_response.get("error"):
                log(f"[RUST ERROR] Worker failed for {fsa.file_name}: {worker_response['error']}")
            _invalidate_rust_worker()
            res = None
    else:
        res = None

    if res is None:
        res = _run_cli_once(cli_bin, fsa_path, analysis_kind, fsa.file_name)
        if res is None:
            return None

    return _apply_rust_result_to_fsa(fsa, res)


def _run_cli_once(cli_bin: Path, fsa_path: Path, analysis_kind: str, file_name: str) -> dict[str, Any] | None:
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
        timeout_seconds = max(_rust_timeout_seconds(analysis_kind), 1)
        try:
            subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=timeout_seconds,
                stdin=subprocess.DEVNULL,
            )
            elapsed = time.monotonic() - start_time
            log(f"[RUST] Engine finished in {elapsed:.1f}s for {file_name}")
        except subprocess.TimeoutExpired as e:
            log(
                f"[RUST ERROR] CLI timed out after {timeout_seconds}s for {file_name}. "
                f"Stderr: {e.stderr}"
            )
            return None
        except subprocess.CalledProcessError as e:
            log(f"[RUST ERROR] CLI failed with code {e.returncode} for {file_name}. Stderr: {e.stderr}")
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
        return results[0]


def _apply_rust_result_to_fsa(fsa: FsaFile, res: dict[str, Any]) -> FsaFile | None:
    fit_preview = res.get("ladder_fit_preview")
    if not fit_preview:
        log("[RUST ERROR] ladder_fit_preview missing")
        return None

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

    model_preview = fit_preview.get("sizing_model")
    if not model_preview or not model_preview.get("predicted_ladder_basepairs"):
        log("[RUST ERROR] model_preview missing predicted_ladder_basepairs")
        return None

    expected_bps = model_preview["predicted_ladder_basepairs"]
    if len(scan_indices) != len(expected_bps):
        log(
            f"[RUST ERROR] Mismatch between selected scan indices ({len(scan_indices)}) "
            f"and predicted basepairs ({len(expected_bps)})."
        )
        return None

    ok, reason = _validate_rust_anchor_selection(fsa, scan_indices, expected_bps)
    if not ok:
        log(
            f"[RUST GUARDRAIL] Rejected anchor set for {fsa.file_name}: {reason}. "
            "Returning control to the runtime without Python rescue."
        )
        return None

    fsa.best_size_standard = np.array(scan_indices, dtype=float)
    fsa.ladder_steps = np.array(expected_bps, dtype=float)
    fsa.expected_ladder_steps = np.array(expected_bps, dtype=float)

    hydrated = _apply_rust_sizing_model_to_fsa(fsa, scan_indices, expected_bps, model_preview)
    if hydrated is not None:
        return hydrated

    try:
        fsa = fit_size_standard_to_ladder(fsa)
    except Exception as e:
        log(f"[HYBRID ERROR] fit_size_standard_to_ladder failed: {e}")
        return None

    return fsa
