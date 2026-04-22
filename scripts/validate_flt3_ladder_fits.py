from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import math
import os
import re
import signal
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


from core.analysis import compute_ladder_qc_metrics  # noqa: E402
from core.analyses.flt3.classification import classify_fsa  # noqa: E402
from core.analyses.flt3.config import PREFERRED_INJECTION_TIME  # noqa: E402
from core.analyses.flt3.pipeline import (  # noqa: E402
    FLT3_LADDER_QC_THRESHOLD,
    FLT3_REVIEW_MAX_RESIDUAL_BP,
    _analyse_fsa_candidate,
)
from core.analyses.shared_pipeline import scan_fsa_files  # noqa: E402


DEFAULT_DATA_ROOT = Path("/Volumes/T7 Shield/DATA/flt3")
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "validation_outputs"
FLT3_ASSAYS = {"FLT3-ITD", "FLT3-D835"}
DIT_SPECIMEN_RE = re.compile(r"^\d{2}OUM\d{5}$", re.IGNORECASE)
FLT3_TEMPLATE_RESCUE_AUTO_ACCEPT_MAX_RESIDUAL_BP = 4.0
FLT3_TEMPLATE_RESCUE_AUTO_ACCEPT_MIN_R2 = 0.9997


class Flt3ValidationTimeout(TimeoutError):
    pass

METRIC_FIELDS = [
    "path",
    "relative_path",
    "year",
    "run_dir",
    "file",
    "assay",
    "group",
    "analysis_type",
    "specimen_id",
    "well_id",
    "selection_key",
    "injection_time",
    "preferred_injection_time",
    "is_preferred_injection",
    "run_date",
    "run_time",
    "run_name",
    "source_run_dir",
    "status",
    "needs_manual_review",
    "review_reason",
    "ladder_fit_strategy",
    "ladder_fit_note",
    "ladder_review_required",
    "ladder_r2",
    "mean_abs_error_bp",
    "max_abs_error_bp",
    "n_ladder_steps",
    "n_size_standard_peaks",
    "ladder_expected_step_count",
    "ladder_fitted_step_count",
    "ladder_missing_expected_steps",
    "sizing_method",
    "error",
]

REVIEW_FIELDS = [
    "path",
    "queue_order",
    "queue_group",
    "status",
    "review_reason",
    "assay",
    "year",
    "run_dir",
    "file",
    "specimen_id",
    "well_id",
    "injection_time",
    "analysis_type",
    "r2",
    "max_bp_err",
    "mean_bp_err",
    "ladder_fit_strategy",
    "ladder_fit_note",
    "missing_expected_steps",
    "expected_step_count",
    "fitted_step_count",
    "relative_path",
]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate FLT3 ladder fits across a large .fsa tree and export a manual-review manifest."
        )
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help=f"Root containing FLT3 year folders. Default: {DEFAULT_DATA_ROOT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to validation_outputs/flt3_ladder_validation_<timestamp>.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, min(4, (os.cpu_count() or 2) // 2)),
        help="Parallel workers. Use 1 for fully sequential debugging.",
    )
    parser.add_argument(
        "--year",
        dest="years",
        action="append",
        default=[],
        help="Optional year folder to include, e.g. 2024. Repeat to include multiple years.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional max number of FSA files to validate after sorting.",
    )
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=None,
        help="Optional CSV containing a path/file_path/fsa_path column. When provided, validate only those files.",
    )
    parser.add_argument(
        "--include-npm1",
        action="store_true",
        help="Also validate NPM1 files found in the FLT3/NPM1 folders.",
    )
    parser.add_argument(
        "--dit-only",
        action="store_true",
        help="Keep only ordinary DIT patient specimens matching NN OUM NNNNN; skips NTC, IVS, MP/MQ, and other labels.",
    )
    parser.add_argument(
        "--require-run-name-contains",
        default="",
        help="Optional substring required in ABI run_name, e.g. 3730DNA.",
    )
    parser.add_argument(
        "--no-suppress-worker-output",
        action="store_true",
        help="Do not suppress verbose per-file Fraggler output from workers.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=45,
        help="Per-file validation timeout. Timed out files are added to the manual-review queue.",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=250,
        help="Write partial CSV/JSON/Markdown outputs after this many completed files. Use 0 to disable.",
    )
    parser.add_argument(
        "--use-rust",
        action="store_true",
        help="Use the Rust Engine for ladder fitting.",
    )
    parser.add_argument(
        "--skip-html-reports",
        action="store_true",
        help="Skip generation of individual HTML reports.",
    )
    return parser


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _safe_float(value: Any, default: float = float("nan")) -> float:
    try:
        result = float(value)
    except Exception:
        return default
    return result


def _json_safe(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def _csv_value(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    if isinstance(value, (list, tuple)):
        return ";".join(str(_csv_value(v)) for v in value)
    return value


def _preferred_injection(meta: dict[str, Any]) -> int:
    analysis_type = str(meta.get("analysis_type") or "")
    assay = str(meta.get("assay") or "")
    if analysis_type == "ratio_quant":
        return 1
    if assay == "FLT3-D835":
        return 3
    if analysis_type in PREFERRED_INJECTION_TIME:
        return int(PREFERRED_INJECTION_TIME[analysis_type])
    return int(PREFERRED_INJECTION_TIME.get(assay, meta.get("protocol_injection_time") or meta.get("injection_time") or 0))


def _should_auto_accept_flt3_review_case(
    *,
    ladder_fit_strategy: str,
    missing_steps: list[float],
    max_abs_error_bp: float,
    r2: float,
    fitted_step_count: int,
    expected_step_count: int,
) -> bool:
    if ladder_fit_strategy != "flt3_template_rescue":
        return False
    if missing_steps:
        return False
    if expected_step_count <= 0 or fitted_step_count != expected_step_count:
        return False
    if not math.isfinite(max_abs_error_bp) or not math.isfinite(r2):
        return False
    return bool(
        r2 >= FLT3_TEMPLATE_RESCUE_AUTO_ACCEPT_MIN_R2
        and max_abs_error_bp <= FLT3_TEMPLATE_RESCUE_AUTO_ACCEPT_MAX_RESIDUAL_BP
    )


def _year_for(path: Path, data_root: Path) -> str:
    try:
        rel = path.relative_to(data_root)
    except ValueError:
        return ""
    return rel.parts[0] if rel.parts else ""


def _review_reason(
    *,
    status: str,
    ladder_fit_strategy: str,
    missing_steps: list[float],
    max_abs_error_bp: float,
    r2: float,
    error: str = "",
) -> str:
    if error:
        lowered_error = error.lower()
        if "timed out" in lowered_error or "timeout" in lowered_error:
            return "Per-file validation timed out; treat as short/incomplete-trace review, not a ladder-fit failure."
        return error
    if status == "ladder_fit_failed":
        return "No ladder fit returned by the FLT3 analysis path."
    if status == "manual_adjustment":
        return "Saved manual ladder adjustment already applied."
    if ladder_fit_strategy == "short_trace":
        return "Short ROX trace: " + (
            "expected ladder steps beyond trace: " + ", ".join(f"{step:.0f}" for step in missing_steps)
            if missing_steps
            else "full GS500ROX ladder assignment is not reliable."
        )
    if missing_steps:
        return "Missing expected ladder steps: " + ", ".join(f"{step:.0f}" for step in missing_steps)
    if ladder_fit_strategy == "auto_partial":
        return "Auto partial fit selected."
    if math.isfinite(max_abs_error_bp) and max_abs_error_bp > FLT3_REVIEW_MAX_RESIDUAL_BP:
        return f"Max residual {max_abs_error_bp:.2f} bp exceeds {FLT3_REVIEW_MAX_RESIDUAL_BP:.2f} bp review gate."
    if math.isfinite(r2) and r2 <= FLT3_LADDER_QC_THRESHOLD:
        return f"R2 {r2:.6f} is at/below {FLT3_LADDER_QC_THRESHOLD:.2f} QC gate."
    if status != "ok":
        return status
    return "ok"


def _queue_group(row: dict[str, Any]) -> str:
    status = str(row.get("status") or "")
    strategy = str(row.get("ladder_fit_strategy") or "")
    missing = str(row.get("ladder_missing_expected_steps") or "")
    reason = str(row.get("review_reason") or "")
    error = str(row.get("error") or "")
    lowered_reason = reason.lower()
    lowered_error = error.lower()
    if (
        strategy == "timeout_review"
        or "timed out" in lowered_reason
        or "timed out" in lowered_error
        or "timeout" in lowered_reason
        or "timeout" in lowered_error
    ):
        return "short_trace_timeout"
    if status in {"analysis_error", "ladder_fit_failed"}:
        return "fit_failed"
    if strategy == "short_trace":
        return "short_trace"
    if missing:
        return "missing_steps"
    if "Max residual" in reason:
        return "high_residual"
    if status == "ladder_qc_failed":
        return "low_r2"
    if strategy == "auto_partial":
        return "partial_fit"
    return status or "review"


def _needs_manual_review(status: str, ladder_fit_strategy: str) -> bool:
    if status in {"review_required", "ladder_qc_failed", "ladder_fit_failed", "analysis_error"}:
        return True
    return ladder_fit_strategy == "auto_partial"


def _row_for_failed_path(path: Path, data_root: Path, status: str, error: str) -> dict[str, Any]:
    try:
        relative = str(path.relative_to(data_root))
    except ValueError:
        relative = path.name
    error_lower = error.lower()
    ladder_fit_strategy = "timeout_review" if ("timed out" in error_lower or "timeout" in error_lower) else ""
    return {
        "path": str(path),
        "relative_path": relative,
        "year": _year_for(path, data_root),
        "run_dir": path.parent.name,
        "file": path.name,
        "status": status,
        "needs_manual_review": True,
        "ladder_fit_strategy": ladder_fit_strategy,
        "ladder_review_required": True,
        "review_reason": _review_reason(
            status=status,
            ladder_fit_strategy=ladder_fit_strategy,
            missing_steps=[],
            max_abs_error_bp=float("inf"),
            r2=float("nan"),
            error=error,
        ),
        "error": error,
    }


def _raise_timeout(_signum: int, _frame: Any) -> None:
    raise Flt3ValidationTimeout("per-file validation timeout")


def validate_one_path(payload: tuple[str, str, bool, bool, int, bool, str, bool, bool]) -> dict[str, Any]:
    path_str, data_root_str, include_npm1, suppress_output, timeout_seconds, dit_only, required_run_name, use_rust, skip_html_reports = payload
    path = Path(path_str)
    data_root = Path(data_root_str)

    if use_rust:
        from config import APP_SETTINGS
        APP_SETTINGS.setdefault("engine", {})["use_rust"] = True
        APP_SETTINGS.setdefault("engine", {})["skip_html_reports"] = skip_html_reports

    stream_cm = contextlib.nullcontext()
    if suppress_output:
        stream_cm = contextlib.ExitStack()
        stream_cm.enter_context(contextlib.redirect_stdout(io.StringIO()))
        stream_cm.enter_context(contextlib.redirect_stderr(io.StringIO()))

    old_handler = None
    timeout_enabled = timeout_seconds > 0 and hasattr(signal, "SIGALRM")
    if timeout_enabled:
        old_handler = signal.signal(signal.SIGALRM, _raise_timeout)
        signal.setitimer(signal.ITIMER_REAL, float(timeout_seconds))

    try:
        with stream_cm:
            return _validate_one_path_inner(
                path,
                data_root,
                include_npm1,
                dit_only=dit_only,
                required_run_name=required_run_name,
                skip_html_reports=skip_html_reports,
            )
    except Flt3ValidationTimeout:
        return _row_for_failed_path(path, data_root, "review_required", f"timed out after {timeout_seconds}s")
    except Exception as exc:
        return _row_for_failed_path(path, data_root, "analysis_error", str(exc))
    finally:
        if timeout_enabled:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, old_handler or signal.SIG_DFL)


def _validate_one_path_inner(
    path: Path,
    data_root: Path,
    include_npm1: bool,
    *,
    dit_only: bool,
    required_run_name: str,
    skip_html_reports: bool = False,
) -> dict[str, Any]:
    try:
        meta = classify_fsa(path)
    except Exception as exc:
        return _row_for_failed_path(path, data_root, "analysis_error", f"classification failed: {exc}")

    if meta is None:
        row = _row_for_failed_path(path, data_root, "skipped_unclassified", "")
        row["needs_manual_review"] = False
        row["review_reason"] = "unclassified"
        return row

    assay = str(meta.get("assay") or "")
    if assay not in FLT3_ASSAYS and not include_npm1:
        row = _row_for_failed_path(path, data_root, "skipped_non_flt3", "")
        row.update({key: meta.get(key, "") for key in ("assay", "group", "analysis_type", "specimen_id", "well_id", "selection_key")})
        row["needs_manual_review"] = False
        row["review_reason"] = f"skipped assay {assay}"
        return row

    specimen_id = str(meta.get("specimen_id") or "")
    if dit_only and not DIT_SPECIMEN_RE.fullmatch(specimen_id):
        row = _row_for_failed_path(path, data_root, "skipped_non_dit", "")
        row.update({key: meta.get(key, "") for key in ("assay", "group", "analysis_type", "specimen_id", "well_id", "selection_key", "run_name")})
        row["needs_manual_review"] = False
        row["review_reason"] = f"skipped non-DIT specimen {specimen_id}"
        return row

    required_run_name = required_run_name.strip()
    run_name = str(meta.get("run_name") or "")
    if required_run_name and required_run_name.lower() not in run_name.lower():
        row = _row_for_failed_path(path, data_root, "skipped_run_name", "")
        row.update({key: meta.get(key, "") for key in ("assay", "group", "analysis_type", "specimen_id", "well_id", "selection_key", "run_name")})
        row["needs_manual_review"] = False
        row["review_reason"] = f"skipped run_name {run_name}"
        return row

    preferred_injection = _preferred_injection(meta)
    try:
        fsa = _analyse_fsa_candidate(
            path,
            str(meta.get("primary_peak_channel") or "DATA1"),
            assay,
            str(meta.get("analysis_type") or ""),
        )
    except Flt3ValidationTimeout:
        raise
    except Exception as exc:
        row = _row_for_failed_path(path, data_root, "analysis_error", str(exc))
        row.update(meta)
        row["preferred_injection_time"] = preferred_injection
        row["is_preferred_injection"] = int(meta.get("injection_time", 0) or 0) == preferred_injection
        return row

    if fsa is None:
        row = _row_for_failed_path(path, data_root, "ladder_fit_failed", "")
        row.update(meta)
        row["preferred_injection_time"] = preferred_injection
        row["is_preferred_injection"] = int(meta.get("injection_time", 0) or 0) == preferred_injection
        return row

    try:
        metrics = compute_ladder_qc_metrics(fsa)
    except Exception as exc:
        row = _row_for_failed_path(path, data_root, "analysis_error", f"QC metrics failed: {exc}")
        row.update(meta)
        row["preferred_injection_time"] = preferred_injection
        row["is_preferred_injection"] = int(meta.get("injection_time", 0) or 0) == preferred_injection
        return row

    r2 = _safe_float(metrics.get("r2"))
    mean_abs_error_bp = _safe_float(metrics.get("mean_abs_error_bp"), float("inf"))
    max_abs_error_bp = _safe_float(metrics.get("max_abs_error_bp"), float("inf"))
    ladder_fit_strategy = str(getattr(fsa, "ladder_fit_strategy", "auto_full"))
    missing_steps = [float(v) for v in getattr(fsa, "ladder_missing_expected_steps", []) or []]
    expected_steps_arr = np.asarray(
        getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])),
        dtype=float,
    )
    fitted_steps_arr = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    expected_step_count = int(expected_steps_arr.size)
    fitted_step_count = int(fitted_steps_arr.size)
    ladder_review_required = bool(
        getattr(fsa, "ladder_review_required", bool(missing_steps))
        or (math.isfinite(max_abs_error_bp) and max_abs_error_bp > FLT3_REVIEW_MAX_RESIDUAL_BP)
    )
    if _should_auto_accept_flt3_review_case(
        ladder_fit_strategy=ladder_fit_strategy,
        missing_steps=missing_steps,
        max_abs_error_bp=max_abs_error_bp,
        r2=r2,
        fitted_step_count=fitted_step_count,
        expected_step_count=expected_step_count,
    ):
        ladder_review_required = False

    if ladder_fit_strategy == "manual_adjustment":
        status = "manual_adjustment"
    elif ladder_review_required:
        status = "review_required"
    elif math.isfinite(r2) and r2 > FLT3_LADDER_QC_THRESHOLD:
        status = "ok"
    else:
        status = "ladder_qc_failed"

    review_reason = _review_reason(
        status=status,
        ladder_fit_strategy=ladder_fit_strategy,
        missing_steps=missing_steps,
        max_abs_error_bp=max_abs_error_bp,
        r2=r2,
    )
    try:
        relative = str(path.relative_to(data_root))
    except ValueError:
        relative = path.name

    expected_steps = getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", []))
    fitted_steps = getattr(fsa, "ladder_steps", [])
    row = {
        "path": str(path),
        "relative_path": relative,
        "year": _year_for(path, data_root),
        "run_dir": path.parent.name,
        "file": path.name,
        "assay": assay,
        "group": meta.get("group") or "",
        "analysis_type": meta.get("analysis_type") or "",
        "specimen_id": meta.get("specimen_id") or "",
        "well_id": meta.get("well_id") or "",
        "selection_key": meta.get("selection_key") or "",
        "injection_time": int(meta.get("injection_time", 0) or 0),
        "preferred_injection_time": preferred_injection,
        "is_preferred_injection": int(meta.get("injection_time", 0) or 0) == preferred_injection,
        "run_date": meta.get("run_date") or "",
        "run_time": meta.get("run_time") or "",
        "run_name": meta.get("run_name") or "",
        "source_run_dir": meta.get("source_run_dir") or "",
        "status": status,
        "needs_manual_review": _needs_manual_review(status, ladder_fit_strategy),
        "review_reason": review_reason,
        "ladder_fit_strategy": ladder_fit_strategy,
        "ladder_fit_note": str(getattr(fsa, "ladder_fit_note", "")),
        "ladder_review_required": ladder_review_required,
        "ladder_r2": r2,
        "mean_abs_error_bp": mean_abs_error_bp,
        "max_abs_error_bp": max_abs_error_bp,
        "n_ladder_steps": int(metrics.get("n_ladder_steps") or 0),
        "n_size_standard_peaks": int(metrics.get("n_size_standard_peaks") or 0),
        "ladder_expected_step_count": int(len(expected_steps)),
        "ladder_fitted_step_count": int(len(fitted_steps)),
        "ladder_missing_expected_steps": missing_steps,
        "sizing_method": str(getattr(fsa, "_flt3_sizing_method", "")),
        "error": "",
    }
    return row


def _load_manifest_paths(manifest_path: Path) -> list[Path]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        path_col = next(
            (col for col in (reader.fieldnames or []) if col.lower() in {"path", "file_path", "fsa_path"}),
            None,
        )
        if path_col is None:
            raise ValueError(f"Manifest {manifest_path} must include a path, file_path, or fsa_path column.")
        files: list[Path] = []
        for row in reader:
            raw_path = str(row.get(path_col) or "").strip()
            if not raw_path:
                continue
            path = Path(raw_path).expanduser()
            if not path.is_absolute():
                path = (manifest_path.parent / path).resolve()
            files.append(path)
    return sorted(dict.fromkeys(files))


def _discover_files(data_root: Path, selected_years: list[str], limit: int, input_manifest: Path | None = None) -> list[Path]:
    files = _load_manifest_paths(input_manifest.expanduser().resolve()) if input_manifest else scan_fsa_files(data_root, recursive=True)
    if selected_years:
        allowed = set(str(year) for year in selected_years)
        files = [path for path in files if _year_for(path, data_root) in allowed]
    files = sorted(files)
    if limit > 0:
        files = files[:limit]
    return files


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field, "")) for field in fields})


def _review_sort_key(row: dict[str, Any]) -> tuple[int, float, float, str]:
    status = str(row.get("status") or "")
    status_rank = {
        "analysis_error": 0,
        "ladder_fit_failed": 1,
        "review_required": 2,
        "ladder_qc_failed": 3,
    }.get(status, 9)
    max_err = _safe_float(row.get("max_abs_error_bp"), -1.0)
    r2 = _safe_float(row.get("ladder_r2"), float("inf"))
    return (status_rank, -max_err if math.isfinite(max_err) else -999999.0, r2, str(row.get("path") or ""))


def _build_review_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    review_rows = [row for row in rows if bool(row.get("needs_manual_review"))]
    review_rows = sorted(review_rows, key=_review_sort_key)
    exported: list[dict[str, Any]] = []
    for idx, row in enumerate(review_rows, start=1):
        exported.append(
            {
                "path": row.get("path", ""),
                "queue_order": idx,
                "queue_group": _queue_group(row),
                "status": row.get("status", ""),
                "review_reason": row.get("review_reason", ""),
                "assay": row.get("assay", ""),
                "year": row.get("year", ""),
                "run_dir": row.get("run_dir", ""),
                "file": row.get("file", ""),
                "specimen_id": row.get("specimen_id", ""),
                "well_id": row.get("well_id", ""),
                "injection_time": row.get("injection_time", ""),
                "analysis_type": row.get("analysis_type", ""),
                "r2": row.get("ladder_r2", ""),
                "max_bp_err": row.get("max_abs_error_bp", ""),
                "mean_bp_err": row.get("mean_abs_error_bp", ""),
                "ladder_fit_strategy": row.get("ladder_fit_strategy", ""),
                "ladder_fit_note": row.get("ladder_fit_note", ""),
                "missing_expected_steps": row.get("ladder_missing_expected_steps", ""),
                "expected_step_count": row.get("ladder_expected_step_count", ""),
                "fitted_step_count": row.get("ladder_fitted_step_count", ""),
                "relative_path": row.get("relative_path", ""),
            }
        )
    return exported


def _actionable_review_sort_key(row: dict[str, Any]) -> tuple[int, float, float, str]:
    group = str(row.get("queue_group") or "")
    group_rank = {
        "missing_steps": 0,
        "high_residual": 1,
        "low_r2": 2,
        "partial_fit": 3,
    }.get(group, 9)
    max_err = _safe_float(row.get("max_bp_err"), -1.0)
    r2 = _safe_float(row.get("r2"), float("inf"))
    return (
        group_rank,
        -max_err if math.isfinite(max_err) else 0.0,
        r2 if math.isfinite(r2) else float("inf"),
        str(row.get("relative_path") or row.get("path") or ""),
    )


def _renumber_review_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    exported = [dict(row) for row in rows]
    for index, row in enumerate(exported, start=1):
        row["queue_order"] = index
    return exported


def _split_review_rows(review_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    fit_failed = [
        row
        for row in review_rows
        if row.get("queue_group") == "fit_failed" or row.get("status") == "ladder_fit_failed"
    ]
    actionable = [
        row
        for row in review_rows
        if row.get("queue_group") != "fit_failed" and row.get("status") != "ladder_fit_failed"
    ]
    return (
        _renumber_review_rows(sorted(actionable, key=_actionable_review_sort_key)),
        _renumber_review_rows(fit_failed),
    )


def _manual_adjustment_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    adjusted = [
        row
        for row in rows
        if row.get("status") == "manual_adjustment"
        or row.get("ladder_fit_strategy") == "manual_adjustment"
    ]
    return sorted(adjusted, key=lambda row: str(row.get("relative_path") or row.get("path") or ""))


def _counter_by(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    return dict(Counter(str(row.get(field) or "") for row in rows))


def _summarize(rows: list[dict[str, Any]], review_rows: list[dict[str, Any]], data_root: Path, output_dir: Path) -> dict[str, Any]:
    processed_rows = [row for row in rows if not str(row.get("status", "")).startswith("skipped")]
    finite_residuals = [
        _safe_float(row.get("max_abs_error_bp"))
        for row in processed_rows
        if math.isfinite(_safe_float(row.get("max_abs_error_bp")))
    ]
    finite_r2 = [
        _safe_float(row.get("ladder_r2"))
        for row in processed_rows
        if math.isfinite(_safe_float(row.get("ladder_r2")))
    ]
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "output_dir": str(output_dir),
        "total_rows": len(rows),
        "validated_count": len(processed_rows),
        "manual_review_count": len(review_rows),
        "status_counts": _counter_by(rows, "status"),
        "year_counts": _counter_by(processed_rows, "year"),
        "assay_counts": _counter_by(processed_rows, "assay"),
        "strategy_counts": _counter_by(processed_rows, "ladder_fit_strategy"),
        "queue_group_counts": dict(Counter(str(row.get("queue_group") or "") for row in review_rows)),
        "max_residual_bp": max(finite_residuals) if finite_residuals else None,
        "median_max_residual_bp": float(np.median(finite_residuals)) if finite_residuals else None,
        "median_r2": float(np.median(finite_r2)) if finite_r2 else None,
        "worst_cases": review_rows[:50],
    }


def _write_markdown_summary(summary: dict[str, Any], output_dir: Path) -> None:
    status_lines = "\n".join(
        f"- {key or '(blank)'}: {value}" for key, value in sorted(summary["status_counts"].items())
    )
    queue_lines = "\n".join(
        f"- {key or '(blank)'}: {value}" for key, value in sorted(summary["queue_group_counts"].items())
    )
    strategy_lines = "\n".join(
        f"- {key or '(blank)'}: {value}" for key, value in sorted(summary["strategy_counts"].items())
    )
    worst_lines = []
    for row in summary["worst_cases"][:20]:
        worst_lines.append(
            "- "
            f"{row.get('queue_order')}. {row.get('status')} / {row.get('queue_group')} | "
            f"{row.get('assay')} | R2={_csv_value(row.get('r2'))} | "
            f"max={_csv_value(row.get('max_bp_err'))} | {row.get('relative_path')}"
        )
    markdown = f"""# FLT3 Ladder Validation Summary

Generated UTC: {summary["generated_at_utc"]}

Data root: `{summary["data_root"]}`
Output dir: `{summary["output_dir"]}`

Validated files: **{summary["validated_count"]}**
Manual review queue: **{summary["manual_review_count"]}**
Median R2: **{summary["median_r2"]}**
Median max residual bp: **{summary["median_max_residual_bp"]}**
Worst max residual bp: **{summary["max_residual_bp"]}**

## Status Counts
{status_lines or "- none"}

## Queue Groups
{queue_lines or "- none"}

## Ladder Fit Strategies
{strategy_lines or "- none"}

## First 20 Review Cases
{chr(10).join(worst_lines) if worst_lines else "- none"}
"""
    (output_dir / "summary.md").write_text(markdown, encoding="utf-8")


def _write_output_bundle(
    output_dir: Path,
    rows: list[dict[str, Any]],
    data_root: Path,
    *,
    suffix: str = "",
) -> dict[str, Any]:
    ordered_rows = sorted(rows, key=lambda row: str(row.get("path") or ""))
    review_rows = _build_review_rows(ordered_rows)
    summary = _summarize(ordered_rows, review_rows, data_root, output_dir)

    metrics_csv = output_dir / f"flt3_ladder_metrics{suffix}.csv"
    review_csv = output_dir / f"flt3_ladder_review_manifest{suffix}.csv"
    summary_json = output_dir / f"summary{suffix}.json"
    summary_md = output_dir / f"summary{suffix}.md"

    _write_csv(metrics_csv, ordered_rows, METRIC_FIELDS)
    _write_csv(review_csv, review_rows, REVIEW_FIELDS)
    actionable_rows, fit_failed_rows = _split_review_rows(review_rows)
    _write_csv(
        output_dir / f"flt3_ladder_actionable_review_manifest{suffix}.csv",
        actionable_rows,
        REVIEW_FIELDS,
    )
    _write_csv(
        output_dir / f"flt3_ladder_fit_failed_manifest{suffix}.csv",
        fit_failed_rows,
        REVIEW_FIELDS,
    )
    _write_csv(
        output_dir / f"flt3_ladder_manual_adjustments_applied{suffix}.csv",
        _manual_adjustment_rows(ordered_rows),
        METRIC_FIELDS,
    )
    summary_json.write_text(json.dumps(_json_safe(summary), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_markdown_summary(summary, output_dir)
    if suffix:
        (output_dir / "summary.md").replace(summary_md)

    return summary


def run_validation(
    data_root: Path,
    output_dir: Path,
    *,
    workers: int,
    selected_years: list[str],
    limit: int,
    input_manifest: Path | None,
    include_npm1: bool,
    suppress_worker_output: bool,
    timeout_seconds: int,
    checkpoint_every: int,
    dit_only: bool,
    required_run_name: str,
    excluded_basenames: list[str] | None = None,
    use_rust: bool = False,
    skip_html_reports: bool = False,
    progress_callback=None,
    progress_max_callback=None,
    status_callback=None,
) -> dict[str, Any]:
    data_root = data_root.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    files = _discover_files(data_root, selected_years, limit, input_manifest=input_manifest)
    excluded_name_set = {
        str(name).strip()
        for name in (excluded_basenames or [])
        if str(name).strip()
    }
    if excluded_name_set:
        files = [path for path in files if path.name not in excluded_name_set]
    payloads = [
        (
            str(path),
            str(data_root),
            include_npm1,
            suppress_worker_output,
            int(timeout_seconds),
            bool(dit_only),
            str(required_run_name or ""),
            bool(use_rust),
            bool(skip_html_reports),
        )
        for path in files
    ]

    rows: list[dict[str, Any]] = []
    print(f"Discovered {len(files)} FSA files under {data_root}", flush=True)
    if progress_max_callback is not None:
        try:
            progress_max_callback.emit(len(payloads))
        except AttributeError:
            progress_max_callback(len(payloads))
    if status_callback is not None:
        message = f"Discovered {len(payloads)} FLT3 files to validate."
        try:
            status_callback.emit(message)
        except AttributeError:
            status_callback(message)

    def maybe_checkpoint(idx: int) -> None:
        if checkpoint_every <= 0 or idx <= 0 or idx % checkpoint_every != 0:
            return
        _write_output_bundle(output_dir, rows, data_root, suffix=".partial")
        print(f"Checkpoint saved after {idx}/{len(payloads)} files", flush=True)
        if status_callback is not None:
            message = f"Checkpoint saved after {idx}/{len(payloads)} files ({max(len(payloads) - idx, 0)} remaining)."
            try:
                status_callback.emit(message)
            except AttributeError:
                status_callback(message)

    if workers <= 1:
        iterator = map(validate_one_path, payloads)
        for idx, row in enumerate(iterator, start=1):
            rows.append(row)
            if idx == 1 or idx % 100 == 0 or idx == len(payloads):
                print(f"Validated {idx}/{len(payloads)} files", flush=True)
                if status_callback is not None:
                    message = f"Validated {idx}/{len(payloads)} files ({max(len(payloads) - idx, 0)} remaining)."
                    try:
                        status_callback.emit(message)
                    except AttributeError:
                        status_callback(message)
            if progress_callback is not None:
                try:
                    progress_callback.emit(idx)
                except AttributeError:
                    progress_callback(idx)
            maybe_checkpoint(idx)
    else:
        from multiprocessing import Pool

        with Pool(workers) as pool:
            for idx, row in enumerate(pool.imap_unordered(validate_one_path, payloads, chunksize=1), start=1):
                rows.append(row)
                if idx == 1 or idx % 100 == 0 or idx == len(payloads):
                    print(f"Validated {idx}/{len(payloads)} files", flush=True)
                    if status_callback is not None:
                        message = f"Validated {idx}/{len(payloads)} files ({max(len(payloads) - idx, 0)} remaining)."
                        try:
                            status_callback.emit(message)
                        except AttributeError:
                            status_callback(message)
                if progress_callback is not None:
                    try:
                        progress_callback.emit(idx)
                    except AttributeError:
                        progress_callback(idx)
                maybe_checkpoint(idx)

    summary = _write_output_bundle(output_dir, rows, data_root)
    summary["excluded_basenames"] = sorted(excluded_name_set)
    summary["excluded_count"] = len(excluded_name_set)

    print(json.dumps(_json_safe(summary), indent=2, sort_keys=True))
    print(f"Metrics CSV: {output_dir / 'flt3_ladder_metrics.csv'}")
    print(f"Review manifest: {output_dir / 'flt3_ladder_review_manifest.csv'}")
    print(f"Markdown summary: {output_dir / 'summary.md'}")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    output_dir = args.output_dir or (DEFAULT_OUTPUT_ROOT / f"flt3_ladder_validation_{_timestamp()}")
    run_validation(
        args.data_root,
        output_dir,
        workers=max(1, int(args.workers)),
        selected_years=[str(year) for year in args.years],
        limit=max(0, int(args.limit)),
        input_manifest=args.input_manifest,
        include_npm1=bool(args.include_npm1),
        suppress_worker_output=not bool(args.no_suppress_worker_output),
        timeout_seconds=max(0, int(args.timeout_seconds)),
        checkpoint_every=max(0, int(args.checkpoint_every)),
        dit_only=bool(args.dit_only),
        required_run_name=str(args.require_run_name_contains or ""),
        use_rust=getattr(args, "use_rust", False),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
