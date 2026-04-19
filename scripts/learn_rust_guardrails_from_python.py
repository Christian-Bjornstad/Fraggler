from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config import APP_SETTINGS
from core.analysis import (
    LADDER_FIT_PROFILE_CLONALITY_LIZ500,
    LADDER_FIT_PROFILE_CLONALITY_ROX400HD,
    analyse_fsa_liz,
    analyse_fsa_rox,
    compute_ladder_qc_metrics,
)
from core.analyses.clonality.config import LIZ_LADDER, ROX_LADDER
from core.analyses.clonality.pipeline import classify_fsa
from core.rust_bridge import run_ladder_fit_hybrid
from fraggler.fraggler import FsaFile, print_green, print_warning
from scripts.reanalyze_clonality_review_cases import resolve_fsa_path


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _fit_status(fsa: Any) -> str:
    if fsa is None:
        return "analysis_error"
    strategy = str(getattr(fsa, "ladder_fit_strategy", "") or "")
    review_required = bool(getattr(fsa, "ladder_review_required", False))
    if strategy == "manual_adjustment":
        return "manual_adjustment"
    if review_required:
        return "review_required"
    return "ok"


def _anchor_stats(fsa: Any) -> dict[str, float]:
    anchors = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    if anchors.size == 0:
        return {
            "anchor_count": 0.0,
            "anchor_first": float("nan"),
            "anchor_last": float("nan"),
            "anchor_span": float("nan"),
            "anchor_median_gap": float("nan"),
        }
    gaps = np.diff(anchors)
    return {
        "anchor_count": float(anchors.size),
        "anchor_first": float(anchors[0]),
        "anchor_last": float(anchors[-1]),
        "anchor_span": float(anchors[-1] - anchors[0]) if anchors.size > 1 else 0.0,
        "anchor_median_gap": float(np.median(gaps)) if gaps.size else float("nan"),
    }


def _build_base_fsa(fsa_path: Path, sample_channel: str, ladder_kind: str) -> FsaFile:
    is_liz = "LIZ" in ladder_kind.upper()
    if is_liz:
        return FsaFile(
            file=str(fsa_path),
            ladder=LIZ_LADDER,
            sample_channel=sample_channel,
            min_distance_between_peaks=30,
            min_size_standard_height=300,
            size_standard_channel="DATA105",
        )
    return FsaFile(
        file=str(fsa_path),
        ladder=ROX_LADDER,
        sample_channel=sample_channel,
        min_distance_between_peaks=30,
        min_size_standard_height=300,
        size_standard_channel="DATA4",
    )


def _run_rust_only(fsa_path: Path, sample_channel: str, ladder_kind: str) -> Any:
    base_fsa = _build_base_fsa(fsa_path, sample_channel, ladder_kind)
    base_fsa.analysis_id = "clonality"
    rust_fsa = run_ladder_fit_hybrid(base_fsa, "clonality")
    if rust_fsa is None:
        return None
    qc = compute_ladder_qc_metrics(rust_fsa)
    from core.analysis import _finalize_auto_fit_metadata, _annotate_fit_qc_review

    rust_fsa = _finalize_auto_fit_metadata(rust_fsa)
    rust_fsa = _annotate_fit_qc_review(rust_fsa, qc)
    return rust_fsa


def _run_python_only(fsa_path: Path, sample_channel: str, ladder_kind: str) -> Any:
    APP_SETTINGS.setdefault("engine", {})["use_rust"] = False
    if "LIZ" in ladder_kind.upper():
        return analyse_fsa_liz(
            fsa_path,
            sample_channel,
            ladder_name=LIZ_LADDER,
            ladder_fit_profile=LADDER_FIT_PROFILE_CLONALITY_LIZ500,
        )
    return analyse_fsa_rox(
        fsa_path,
        sample_channel,
        ladder_name=ROX_LADDER,
        ladder_fit_profile=LADDER_FIT_PROFILE_CLONALITY_ROX400HD,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Learn Rust guardrail priors from Python fits by comparing Rust-only and Python-only ladder fits."
        )
    )
    parser.add_argument("--month-run-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--sheet-name", type=str, default="Patient_Runs")
    parser.add_argument("--filter-status", type=str, default="review_required")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    return parser


def run(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    workbook = args.month_run_dir.expanduser().resolve() / "track-clonality.xlsx"
    data_root = args.data_root.expanduser().resolve()

    df = pd.read_excel(workbook, sheet_name=args.sheet_name)
    status = str(args.filter_status or "").strip().lower()
    if status:
        df = df[df.get("LadderQC", pd.Series(dtype=object)).astype(str).str.strip().str.lower().eq(status)]
    if args.limit and int(args.limit) > 0:
        df = df.head(int(args.limit))
    df = df.reset_index(drop=True)
    print_green(f"Training set rows: {len(df)}")

    rows: list[dict[str, Any]] = []
    for idx, row in df.iterrows():
        source_run_dir = str(row.get("SourceRunDir", "") or "")
        assay = str(row.get("Assay", "") or "")
        well = str(row.get("Well", "") or "")
        run_code = str(row.get("RunCode", "") or "")

        try:
            fsa_path = resolve_fsa_path(
                data_dir=data_root,
                source_run_dir=source_run_dir,
                assay=assay,
                well=well,
                run_code=run_code,
            )
            classified = classify_fsa(fsa_path)
            if not classified:
                raise RuntimeError("classify_fsa returned None")
            _assay, _group, ladder_kind, trace_channels, *_rest = classified
            sample_channel = trace_channels[0]

            rust_fsa = _run_rust_only(fsa_path, sample_channel, ladder_kind)
            python_fsa = _run_python_only(fsa_path, sample_channel, ladder_kind)

            rust_metrics = compute_ladder_qc_metrics(rust_fsa) if rust_fsa is not None else {}
            py_metrics = compute_ladder_qc_metrics(python_fsa) if python_fsa is not None else {}

            row_out = {
                "SourceRunDir": source_run_dir,
                "IdentityKey": str(row.get("IdentityKey", "") or ""),
                "Assay": assay,
                "Well": well,
                "Ladder": str(row.get("Ladder", "") or ladder_kind),
                "FsaPath": str(fsa_path),
                "RustStatus": _fit_status(rust_fsa),
                "PythonStatus": _fit_status(python_fsa),
                "RustR2": _safe_float(rust_metrics.get("r2", np.nan)),
                "PythonR2": _safe_float(py_metrics.get("r2", np.nan)),
                "RustMeanResidualBp": _safe_float(rust_metrics.get("mean_abs_error_bp", np.nan)),
                "PythonMeanResidualBp": _safe_float(py_metrics.get("mean_abs_error_bp", np.nan)),
                "RustMaxResidualBp": _safe_float(rust_metrics.get("max_abs_error_bp", np.nan)),
                "PythonMaxResidualBp": _safe_float(py_metrics.get("max_abs_error_bp", np.nan)),
                "RustMaxCurvature": _safe_float(rust_metrics.get("max_curvature", np.nan)),
                "PythonMaxCurvature": _safe_float(py_metrics.get("max_curvature", np.nan)),
                "RustStrategy": str(getattr(rust_fsa, "ladder_fit_strategy", "") or ""),
                "PythonStrategy": str(getattr(python_fsa, "ladder_fit_strategy", "") or ""),
            }
            row_out.update({f"Rust_{k}": v for k, v in _anchor_stats(rust_fsa).items()})
            row_out.update({f"Python_{k}": v for k, v in _anchor_stats(python_fsa).items()})

            row_out["PythonPreferredForRust"] = bool(
                row_out["PythonStatus"] == "ok"
                and row_out["RustStatus"] != "ok"
            )
            rows.append(row_out)
        except Exception as exc:
            rows.append(
                {
                    "SourceRunDir": source_run_dir,
                    "IdentityKey": str(row.get("IdentityKey", "") or ""),
                    "Assay": assay,
                    "Well": well,
                    "Ladder": str(row.get("Ladder", "") or ""),
                    "FsaPath": "",
                    "RustStatus": "analysis_error",
                    "PythonStatus": "analysis_error",
                    "Error": str(exc),
                    "PythonPreferredForRust": False,
                }
            )
            print_warning(f"[LEARN] Failed on {assay} {well} ({source_run_dir}): {exc}")

        if (idx + 1) % 10 == 0 or (idx + 1) == len(df):
            print_green(f"Progress: {idx + 1}/{len(df)}")

    out_df = pd.DataFrame(rows)
    args.output_csv.expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output_csv.expanduser().resolve(), index=False)

    preferred = out_df[out_df.get("PythonPreferredForRust", False).astype(bool)].copy()
    summary = {
        "total_cases": int(len(out_df)),
        "python_preferred_cases": int(len(preferred)),
        "rust_status_counts": out_df.get("RustStatus", pd.Series(dtype=object)).value_counts(dropna=False).to_dict(),
        "python_status_counts": out_df.get("PythonStatus", pd.Series(dtype=object)).value_counts(dropna=False).to_dict(),
        "recommended_rox_time_window": {},
        "recommended_liz_time_window": {},
    }

    if not preferred.empty:
        for ladder_name, key in [("ROX", "recommended_rox_time_window"), ("LIZ", "recommended_liz_time_window")]:
            subset = preferred[preferred.get("Ladder", pd.Series(dtype=object)).astype(str).str.upper().str.contains(ladder_name)]
            if subset.empty:
                continue
            first = pd.to_numeric(subset.get("Python_anchor_first"), errors="coerce").dropna()
            last = pd.to_numeric(subset.get("Python_anchor_last"), errors="coerce").dropna()
            span = pd.to_numeric(subset.get("Python_anchor_span"), errors="coerce").dropna()
            if first.empty or last.empty:
                continue
            summary[key] = {
                "first_q05": float(first.quantile(0.05)),
                "first_q50": float(first.quantile(0.50)),
                "last_q50": float(last.quantile(0.50)),
                "last_q95": float(last.quantile(0.95)),
                "span_q05": float(span.quantile(0.05)) if not span.empty else float("nan"),
                "span_q50": float(span.quantile(0.50)) if not span.empty else float("nan"),
            }

    args.summary_json.expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.expanduser().resolve().write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print_green(f"Wrote learning table: {args.output_csv}")
    print_green(f"Wrote learning summary: {args.summary_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
