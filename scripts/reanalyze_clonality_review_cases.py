from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


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
    get_ladder_candidates,
)
from core.analyses.clonality.candidate_artifacts import write_clonality_candidate_artifacts
from core.analyses.clonality.config import LIZ_LADDER, ROX_LADDER
from core.analyses.clonality.pipeline import classify_fsa
from fraggler.fraggler import print_green, print_warning


def _normalize_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


ASSAY_FILENAME_PATTERNS = {
    "IKZF1": ["ikzf1"],
    "Ktr-albumin": ["ktralbumin", "ktralbuminr", "ktralbuminr13"],
    "TCRbA": ["trbmixa", "tcrba", "tcrb_a"],
    "TCRbB": ["trbmixb", "tcrbb", "tcrb_b"],
    "TCRbC": ["trbmixc", "tcrbc", "tcrb_c"],
    "TCRgA": ["tcrga", "tcrg_a"],
    "TCRgB": ["tcrgb", "tcrg_b"],
    "FR1": ["fr1"],
    "FR2": ["fr2"],
    "FR3": ["fr3"],
    "DHJH_D": ["dhjhmixd", "dhjhd"],
    "DHJH_E": ["dhjhmixe", "dhjhe"],
    "IGK": ["igk"],
    "KDE": ["kde"],
}


def resolve_fsa_path(data_dir: Path, source_run_dir: str, assay: str, well: str, run_code: str) -> Path:
    run_dir = data_dir / source_run_dir
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Missing run directory: {run_dir}")

    assay_tokens = ASSAY_FILENAME_PATTERNS.get(assay, [_normalize_token(assay)])
    well_token = _normalize_token(well)
    run_code_token = _normalize_token(run_code)

    matches: list[Path] = []
    for candidate in sorted(run_dir.glob("*.fsa")):
        token = _normalize_token(candidate.name)
        if run_code_token and run_code_token not in token:
            continue
        if well_token and well_token not in token:
            continue
        if not any(pattern in token for pattern in assay_tokens):
            continue
        matches.append(candidate)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        return matches[0]

    fallback: list[Path] = []
    for candidate in sorted(run_dir.glob("*.fsa")):
        token = _normalize_token(candidate.name)
        if run_code_token and run_code_token not in token:
            continue
        if well_token and well_token not in token:
            continue
        fallback.append(candidate)
    if len(fallback) == 1:
        return fallback[0]
    if fallback:
        return fallback[0]

    raise FileNotFoundError(
        f"Could not resolve FSA for source_run_dir={source_run_dir}, assay={assay}, well={well}, run_code={run_code}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Re-analyze clonality ladder review-required cases from an existing monthly workbook, "
            "compare old/new ladder QC, and optionally write a fresh review bundle."
        )
    )
    parser.add_argument(
        "--month-run-dir",
        type=Path,
        required=True,
        help="Monthly run directory containing track-clonality.xlsx (for example .../month_runs/2025_01).",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        required=True,
        help="Root data directory that contains source run folders with .fsa files.",
    )
    parser.add_argument(
        "--sheet-name",
        type=str,
        default="Patient_Runs",
        help="Workbook sheet to read. Default: Patient_Runs.",
    )
    parser.add_argument(
        "--filter-status",
        type=str,
        default="review_required",
        help="Only re-analyze rows where LadderQC equals this value (case-insensitive). Default: review_required.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional max number of rows to re-analyze (0 means all matching rows).",
    )
    parser.add_argument(
        "--use-rust",
        action="store_true",
        help="Enable Rust hybrid ladder fit while re-analyzing.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Output CSV path. Default: <month-run-dir>/reanalysis/reanalysis_results.csv",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="Output summary JSON path. Default: <month-run-dir>/reanalysis/reanalysis_summary.json",
    )
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        default=None,
        help="Optional output directory for a fresh ladder review bundle built from remaining review_required cases.",
    )
    parser.add_argument(
        "--include-sl",
        action="store_true",
        help="Include SL rows when writing candidate artifacts (only used with --bundle-dir).",
    )
    return parser


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _fit_status_from_fsa(fsa: Any) -> str:
    strategy = str(getattr(fsa, "ladder_fit_strategy", "") or "")
    review_required = bool(getattr(fsa, "ladder_review_required", False))
    if strategy == "manual_adjustment":
        return "manual_adjustment"
    if review_required:
        return "review_required"
    return "ok"


def _row_filter(df: pd.DataFrame, filter_status: str) -> pd.DataFrame:
    status = str(filter_status or "").strip().lower()
    if not status:
        return df
    qc = df.get("LadderQC", pd.Series(dtype=object)).astype(str).str.strip().str.lower()
    return df[qc.eq(status)]


def _analyze_one_case(row: pd.Series, data_root: Path) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
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
    except Exception as exc:
        return (
            {
                "SourceRunDir": source_run_dir,
                "IdentityKey": str(row.get("IdentityKey", "") or ""),
                "Assay": assay,
                "Well": well,
                "Ladder": str(row.get("Ladder", "") or ""),
                "OldLadderQC": str(row.get("LadderQC", "") or ""),
                "OldLadderR2": _safe_float(row.get("LadderR2", np.nan)),
                "NewLadderQC": "analysis_error",
                "NewLadderR2": float("nan"),
                "ExpectedStepCount": int(_safe_float(row.get("LadderExpectedStepCount", 0)) or 0),
                "FittedStepCount": 0,
                "MissingStepCount": int(
                    max(
                        (_safe_float(row.get("LadderExpectedStepCount", 0)) or 0)
                        - (_safe_float(row.get("LadderFittedStepCount", 0)) or 0),
                        0,
                    )
                ),
                "LadderFitStrategy": "",
                "LadderFitNote": f"resolve_fsa_path failed: {exc}",
                "FsaPath": "",
            },
            None,
        )

    classified = classify_fsa(fsa_path)
    if not classified:
        return (
            {
                "SourceRunDir": source_run_dir,
                "IdentityKey": str(row.get("IdentityKey", "") or ""),
                "Assay": assay,
                "Well": well,
                "Ladder": str(row.get("Ladder", "") or ""),
                "OldLadderQC": str(row.get("LadderQC", "") or ""),
                "OldLadderR2": _safe_float(row.get("LadderR2", np.nan)),
                "NewLadderQC": "analysis_error",
                "NewLadderR2": float("nan"),
                "ExpectedStepCount": int(_safe_float(row.get("LadderExpectedStepCount", 0)) or 0),
                "FittedStepCount": 0,
                "MissingStepCount": int(
                    max(
                        (_safe_float(row.get("LadderExpectedStepCount", 0)) or 0)
                        - (_safe_float(row.get("LadderFittedStepCount", 0)) or 0),
                        0,
                    )
                ),
                "LadderFitStrategy": "",
                "LadderFitNote": "classify_fsa returned None",
                "FsaPath": str(fsa_path),
            },
            None,
        )

    (
        _assay,
        _group,
        ladder_name,
        trace_channels,
        _peak_channels,
        _primary_peak_channel,
        _bp_min,
        _bp_max,
    ) = classified

    sample_channel = trace_channels[0]
    try:
        if "LIZ" in str(ladder_name).upper():
            fsa = analyse_fsa_liz(
                fsa_path,
                sample_channel,
                ladder_name=LIZ_LADDER,
                ladder_fit_profile=LADDER_FIT_PROFILE_CLONALITY_LIZ500,
            )
        else:
            fsa = analyse_fsa_rox(
                fsa_path,
                sample_channel,
                ladder_name=ROX_LADDER,
                ladder_fit_profile=LADDER_FIT_PROFILE_CLONALITY_ROX400HD,
            )
    except Exception as exc:
        return (
            {
                "SourceRunDir": source_run_dir,
                "IdentityKey": str(row.get("IdentityKey", "") or ""),
                "Assay": assay,
                "Well": well,
                "Ladder": str(row.get("Ladder", "") or ""),
                "OldLadderQC": str(row.get("LadderQC", "") or ""),
                "OldLadderR2": _safe_float(row.get("LadderR2", np.nan)),
                "NewLadderQC": "analysis_error",
                "NewLadderR2": float("nan"),
                "ExpectedStepCount": int(_safe_float(row.get("LadderExpectedStepCount", 0)) or 0),
                "FittedStepCount": 0,
                "MissingStepCount": int(
                    max(
                        (_safe_float(row.get("LadderExpectedStepCount", 0)) or 0)
                        - (_safe_float(row.get("LadderFittedStepCount", 0)) or 0),
                        0,
                    )
                ),
                "LadderFitStrategy": "",
                "LadderFitNote": f"analysis failed: {exc}",
                "FsaPath": str(fsa_path),
            },
            None,
        )

    if fsa is None:
        return (
            {
                "SourceRunDir": source_run_dir,
                "IdentityKey": str(row.get("IdentityKey", "") or ""),
                "Assay": assay,
                "Well": well,
                "Ladder": str(row.get("Ladder", "") or ""),
                "OldLadderQC": str(row.get("LadderQC", "") or ""),
                "OldLadderR2": _safe_float(row.get("LadderR2", np.nan)),
                "NewLadderQC": "analysis_error",
                "NewLadderR2": float("nan"),
                "ExpectedStepCount": int(_safe_float(row.get("LadderExpectedStepCount", 0)) or 0),
                "FittedStepCount": 0,
                "MissingStepCount": int(
                    max(
                        (_safe_float(row.get("LadderExpectedStepCount", 0)) or 0)
                        - (_safe_float(row.get("LadderFittedStepCount", 0)) or 0),
                        0,
                    )
                ),
                "LadderFitStrategy": "",
                "LadderFitNote": "analysis returned None",
                "FsaPath": str(fsa_path),
            },
            None,
        )

    metrics = compute_ladder_qc_metrics(fsa)
    expected_steps = list(
        map(float, getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])))
    )
    fitted_steps = list(map(float, getattr(fsa, "ladder_steps", [])))
    missing_steps = list(map(float, getattr(fsa, "ladder_missing_expected_steps", [])))
    new_status = _fit_status_from_fsa(fsa)

    result = {
        "SourceRunDir": source_run_dir,
        "IdentityKey": str(row.get("IdentityKey", "") or ""),
        "Assay": assay,
        "Well": well,
        "Ladder": str(row.get("Ladder", "") or ""),
        "OldLadderQC": str(row.get("LadderQC", "") or ""),
        "OldLadderR2": _safe_float(row.get("LadderR2", np.nan)),
        "NewLadderQC": new_status,
        "NewLadderR2": _safe_float(metrics.get("r2", np.nan)),
        "NewMeanResidualBp": _safe_float(metrics.get("mean_abs_error_bp", np.nan)),
        "NewMaxResidualBp": _safe_float(metrics.get("max_abs_error_bp", np.nan)),
        "NewMaxCurvature": _safe_float(metrics.get("max_curvature", np.nan)),
        "ExpectedStepCount": int(len(expected_steps)),
        "FittedStepCount": int(len(fitted_steps)),
        "MissingStepCount": int(len(missing_steps)),
        "MissingStepsBp": ",".join(f"{bp:.1f}" for bp in missing_steps),
        "LadderFitStrategy": str(getattr(fsa, "ladder_fit_strategy", "") or ""),
        "LadderFitNote": str(getattr(fsa, "ladder_fit_note", "") or ""),
        "FsaPath": str(fsa_path),
    }

    candidate_entry = {
        "fsa": fsa,
        "file_name": fsa.file_name,
        "source_run_dir": source_run_dir,
        "assay": assay,
        "well": well,
        "run_code": run_code,
        "run_date": row.get("RunDate", ""),
        "identity_key": str(row.get("IdentityKey", "") or ""),
        "ladder": str(row.get("Ladder", "") or ""),
        "group": str(row.get("Group", "") or ""),
        "ladder_qc_status": new_status,
        "ladder_fit_strategy": str(getattr(fsa, "ladder_fit_strategy", "") or ""),
        "ladder_r2": _safe_float(metrics.get("r2", np.nan)),
    }
    return result, candidate_entry


def _summarize(results_df: pd.DataFrame) -> dict[str, Any]:
    old_status = results_df.get("OldLadderQC", pd.Series(dtype=object)).astype(str)
    new_status = results_df.get("NewLadderQC", pd.Series(dtype=object)).astype(str)

    changed_to_ok = int(((old_status.str.lower() != "ok") & (new_status.str.lower() == "ok")).sum())
    still_review = int((new_status.str.lower() == "review_required").sum())
    errors = int((new_status.str.lower() == "analysis_error").sum())
    total = int(len(results_df))

    by_ladder = Counter(results_df.get("Ladder", pd.Series(dtype=object)).astype(str).tolist())
    by_new_status = Counter(new_status.tolist())
    by_strategy = Counter(results_df.get("LadderFitStrategy", pd.Series(dtype=object)).astype(str).tolist())

    resolved = results_df[new_status.str.lower() == "ok"].copy()
    unresolved = results_df[new_status.str.lower() == "review_required"].copy()

    summary = {
        "total_reanalyzed": total,
        "resolved_to_ok": changed_to_ok,
        "still_review_required": still_review,
        "analysis_errors": errors,
        "resolution_rate": (changed_to_ok / total) if total > 0 else 0.0,
        "counts_by_ladder": dict(sorted(by_ladder.items())),
        "counts_by_new_status": dict(sorted(by_new_status.items())),
        "counts_by_strategy": dict(sorted(by_strategy.items())),
        "resolved_r2_quantiles": {},
        "unresolved_r2_quantiles": {},
    }
    if not resolved.empty:
        quantiles = resolved["NewLadderR2"].dropna().quantile([0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]).to_dict()
        summary["resolved_r2_quantiles"] = {str(k): float(v) for k, v in quantiles.items()}
    if not unresolved.empty:
        quantiles = unresolved["NewLadderR2"].dropna().quantile([0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]).to_dict()
        summary["unresolved_r2_quantiles"] = {str(k): float(v) for k, v in quantiles.items()}
    return summary


def _write_annotator_bundle(
    bundle_dir: Path,
    unresolved_entries: list[dict[str, Any]],
    results_df: pd.DataFrame,
    *,
    include_sl: bool = False,
) -> None:
    outputs = write_clonality_candidate_artifacts(
        bundle_dir,
        unresolved_entries,
        include_sl=include_sl,
        write_gold_label_template=True,
    )
    ladder_candidates_path = outputs["ladder_candidates"]
    try:
        cand_df = pd.read_csv(ladder_candidates_path) if ladder_candidates_path.exists() else pd.DataFrame()
    except EmptyDataError:
        cand_df = pd.DataFrame()

    if cand_df.empty:
        fallback_rows: list[dict[str, Any]] = []
        for entry in unresolved_entries:
            fsa = entry.get("fsa")
            if fsa is None:
                continue
            candidate_table = get_ladder_candidates(fsa)
            if candidate_table.empty:
                continue
            selected = list(
                zip(
                    np.asarray(getattr(fsa, "best_size_standard", []), dtype=float).tolist(),
                    np.asarray(getattr(fsa, "ladder_steps", []), dtype=float).tolist(),
                )
            )
            for _, raw in candidate_table.iterrows():
                time_value = float(raw.get("time", np.nan))
                step_bp = None
                for sel_t, sel_bp in selected:
                    if np.isclose(time_value, float(sel_t), atol=1e-6):
                        step_bp = float(sel_bp)
                        break
                fallback_rows.append(
                    {
                        "month": "",
                        "source_run_dir": str(entry.get("source_run_dir", "") or ""),
                        "assay": str(entry.get("assay", "") or ""),
                        "identity_key": str(entry.get("identity_key", "") or ""),
                        "run_date": str(entry.get("run_date", "") or ""),
                        "run_code": str(entry.get("run_code", "") or ""),
                        "well": str(entry.get("well", "") or ""),
                        "ladder": str(entry.get("ladder", "") or ""),
                        "artifact_row_key": "",
                        "join_key": "",
                        "ladder_join_key": "",
                        "candidate_index": int(raw.get("index", 0)),
                        "candidate_time": time_value,
                        "candidate_intensity": _safe_float(raw.get("intensity", np.nan)),
                        "candidate_source": str(raw.get("source", "") or "auto"),
                        "selected_for_fit": step_bp is not None,
                        "selected_step_bp": step_bp if step_bp is not None else np.nan,
                        "ladder_fit_strategy": str(entry.get("ladder_fit_strategy", "") or ""),
                        "ladder_r2": _safe_float(entry.get("ladder_r2", np.nan)),
                        "ladder_review_required": True,
                        "human_label": "",
                        "human_note": "",
                        "control": "",
                        "sample_kind": "",
                    }
                )
        cand_df = pd.DataFrame(fallback_rows)
        if cand_df.empty:
            print_warning(f"No ladder candidate rows found for annotator bundle: {bundle_dir}")
            return

    unresolved = results_df[results_df["NewLadderQC"].astype(str).str.lower() == "review_required"].copy()
    unresolved_key = set(unresolved["FsaPath"].astype(str).tolist())

    entry_rows: list[dict[str, Any]] = []
    for entry in unresolved_entries:
        fsa = entry.get("fsa")
        fsa_file = str(getattr(fsa, "file", "") or "")
        if fsa_file not in unresolved_key:
            continue
        entry_rows.append(
            {
                "fsa_path": fsa_file,
                "month": "",
                "scope": "Patient",
                "identity_key": str(entry.get("identity_key", "") or ""),
                "source_run_dir": str(entry.get("source_run_dir", "") or ""),
                "assay": str(entry.get("assay", "") or ""),
                "run_date": str(entry.get("run_date", "") or ""),
                "run_code": str(entry.get("run_code", "") or ""),
                "well": str(entry.get("well", "") or ""),
                "ladder": str(entry.get("ladder", "") or ""),
            }
        )
    case_base = pd.DataFrame(entry_rows)
    if case_base.empty:
        print_warning(f"No unresolved entry rows found for annotator bundle: {bundle_dir}")
        return

    merged = unresolved.merge(case_base, left_on="FsaPath", right_on="fsa_path", how="inner")

    case_rows = merged.loc[
        :,
        [
            "month",
            "scope",
            "identity_key",
            "source_run_dir",
            "assay",
            "run_date",
            "run_code",
            "well",
            "ladder",
            "LadderFitStrategy",
            "ExpectedStepCount",
            "FittedStepCount",
            "NewLadderR2",
            "NewMeanResidualBp",
            "NewMaxResidualBp",
            "LadderFitNote",
        ],
    ].copy()
    case_rows = case_rows.rename(
        columns={
            "LadderFitStrategy": "ladder_fit_strategy",
            "ExpectedStepCount": "ladder_expected_step_count",
            "FittedStepCount": "ladder_fitted_step_count",
            "NewLadderR2": "ladder_r2",
            "NewMeanResidualBp": "mean_residual",
            "NewMaxResidualBp": "max_residual",
            "LadderFitNote": "ladder_fit_note",
        }
    )
    case_rows["ladder_qc"] = "review_required"
    case_rows["artifact_row_key"] = ""
    case_rows["join_key"] = ""
    case_rows["ladder_join_key"] = ""
    case_rows["label"] = ""
    case_rows["label_note"] = ""
    case_rows["reviewed_at_utc"] = ""
    case_rows = case_rows.drop_duplicates(subset=["source_run_dir", "assay", "well"], keep="first")

    for col in [
        "month",
        "scope",
        "identity_key",
        "source_run_dir",
        "assay",
        "run_date",
        "run_code",
        "well",
        "ladder",
        "ladder_qc",
        "ladder_fit_strategy",
        "ladder_expected_step_count",
        "ladder_fitted_step_count",
        "ladder_r2",
        "artifact_row_key",
        "join_key",
        "ladder_join_key",
        "label",
        "label_note",
        "reviewed_at_utc",
        "mean_residual",
        "max_residual",
        "ladder_fit_note",
    ]:
        if col not in case_rows.columns:
            case_rows[col] = ""

    case_rows = case_rows[
        [
            "month",
            "scope",
            "identity_key",
            "source_run_dir",
            "assay",
            "run_date",
            "run_code",
            "well",
            "ladder",
            "ladder_qc",
            "ladder_fit_strategy",
            "ladder_expected_step_count",
            "ladder_fitted_step_count",
            "ladder_r2",
            "artifact_row_key",
            "join_key",
            "ladder_join_key",
            "label",
            "label_note",
            "reviewed_at_utc",
            "mean_residual",
            "max_residual",
            "ladder_fit_note",
        ]
    ]

    cand_rows = cand_df.copy()
    for col in [
        "month",
        "source_run_dir",
        "assay",
        "identity_key",
        "run_date",
        "run_code",
        "well",
        "ladder",
        "artifact_row_key",
        "join_key",
        "ladder_join_key",
        "candidate_index",
        "candidate_time",
        "candidate_intensity",
        "candidate_source",
        "selected_for_fit",
        "selected_step_bp",
        "ladder_fit_strategy",
        "ladder_r2",
        "ladder_review_required",
        "human_label",
        "human_note",
        "control",
        "sample_kind",
    ]:
        if col not in cand_rows.columns:
            cand_rows[col] = ""

    cand_rows["month"] = cand_rows["month"].astype(str)
    cand_rows["ladder_review_required"] = True
    cand_rows["human_label"] = cand_rows["human_label"].astype(str)
    cand_rows["human_note"] = cand_rows["human_note"].astype(str)
    cand_rows = cand_rows[
        [
            "month",
            "source_run_dir",
            "assay",
            "identity_key",
            "run_date",
            "run_code",
            "well",
            "ladder",
            "artifact_row_key",
            "join_key",
            "ladder_join_key",
            "candidate_index",
            "candidate_time",
            "candidate_intensity",
            "candidate_source",
            "selected_for_fit",
            "selected_step_bp",
            "ladder_fit_strategy",
            "ladder_r2",
            "ladder_review_required",
            "human_label",
            "human_note",
            "control",
            "sample_kind",
        ]
    ]

    case_path = bundle_dir / "ladder_review_cases.csv"
    cand_path = bundle_dir / "ladder_review_candidates.csv"
    case_rows.to_csv(case_path, index=False)
    cand_rows.to_csv(cand_path, index=False)
    print_green(f"Wrote annotator cases: {case_path}")
    print_green(f"Wrote annotator candidates: {cand_path}")


def run(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    month_run_dir = args.month_run_dir.expanduser().resolve()
    data_root = args.data_root.expanduser().resolve()
    workbook_path = month_run_dir / "track-clonality.xlsx"
    if not workbook_path.exists():
        print_warning(f"Missing workbook: {workbook_path}")
        return 2
    if not data_root.exists():
        print_warning(f"Missing data root: {data_root}")
        return 2

    out_root = month_run_dir / "reanalysis"
    out_root.mkdir(parents=True, exist_ok=True)
    output_csv = args.output_csv.expanduser().resolve() if args.output_csv else out_root / "reanalysis_results.csv"
    summary_json = args.summary_json.expanduser().resolve() if args.summary_json else out_root / "reanalysis_summary.json"

    print_green(f"Loading workbook: {workbook_path}")
    df = pd.read_excel(workbook_path, sheet_name=args.sheet_name)
    filtered = _row_filter(df, args.filter_status).copy()
    if args.limit and int(args.limit) > 0:
        filtered = filtered.head(int(args.limit)).copy()
    filtered = filtered.reset_index(drop=True)
    print_green(f"Rows selected for re-analysis: {len(filtered)}")

    APP_SETTINGS.setdefault("engine", {})["use_rust"] = bool(args.use_rust)
    if args.use_rust:
        print_green("Rust hybrid mode is ENABLED for this re-analysis.")
    else:
        print_warning("Rust hybrid mode is DISABLED for this re-analysis.")

    result_rows: list[dict[str, Any]] = []
    candidate_entries: list[dict[str, Any]] = []

    for idx, row in filtered.iterrows():
        result, candidate_entry = _analyze_one_case(row, data_root=data_root)
        if result is not None:
            result_rows.append(result)
        if candidate_entry is not None:
            candidate_entries.append(candidate_entry)
        if (idx + 1) % 10 == 0 or (idx + 1) == len(filtered):
            print_green(f"Progress: {idx + 1}/{len(filtered)}")

    results_df = pd.DataFrame(result_rows)
    results_df.to_csv(output_csv, index=False)
    print_green(f"Wrote re-analysis results: {output_csv}")

    summary = _summarize(results_df)
    summary["month_run_dir"] = str(month_run_dir)
    summary["workbook_path"] = str(workbook_path)
    summary["data_root"] = str(data_root)
    summary["filter_status"] = str(args.filter_status)
    summary["use_rust"] = bool(args.use_rust)
    summary["limit"] = int(args.limit or 0)
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print_green(f"Wrote summary: {summary_json}")

    if args.bundle_dir:
        bundle_dir = args.bundle_dir.expanduser().resolve()
        unresolved_paths = set(
            results_df.loc[
                results_df["NewLadderQC"].astype(str).str.lower() == "review_required",
                "FsaPath",
            ].astype(str).tolist()
        )
        unresolved_entries = [
            entry for entry in candidate_entries
            if str(getattr(entry.get("fsa"), "file", "") or "") in unresolved_paths
        ]
        bundle_dir.mkdir(parents=True, exist_ok=True)
        if unresolved_entries:
            _write_annotator_bundle(
                bundle_dir,
                unresolved_entries,
                results_df,
                include_sl=bool(args.include_sl),
            )
            print_green(f"Wrote review bundle for {len(unresolved_entries)} unresolved cases: {bundle_dir}")
        else:
            print_warning("No unresolved review_required cases remained; bundle was not written.")

    return 0


if __name__ == "__main__":
    raise SystemExit(run())
