from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from core.analysis import get_ladder_candidates
from core.analyses.clonality.tracking_excel import (
    build_clonality_qc_rules,
    build_tracking_join_fields,
    build_tracking_join_key,
    build_tracking_ladder_join_key,
    build_tracking_pk_join_key,
    build_tracking_row_key,
)
from core.qc.qc_markers import (
    evaluate_peak_near_bp_with_fallback,
    find_peak_near_bp,
    markers_for_entry,
)
from core.qc.qc_rules import QCRules


CANDIDATE_ARTIFACT_VERSION = "v1"
CANDIDATE_ARTIFACT_FILENAMES = {
    "combined": "clonality_candidate_artifacts.csv",
    "ladder_candidates": "clonality_ladder_candidates.csv",
    "pk_candidates": "clonality_pk_candidates.csv",
    "manifest": "clonality_candidate_manifest.json",
    "gold_labels": "clonality_gold_labels.csv",
}
GOLD_LABEL_COLUMNS = [
    "artifact_table",
    "artifact_row_key",
    "label",
    "label_source",
    "reviewer",
    "reviewed_at_utc",
    "notes",
]
_KEEP_ASSAYS = {"SL"}


def build_clonality_candidate_tables(
    entries: Sequence[dict[str, Any]],
    *,
    rules: QCRules | None = None,
    include_sl: bool = False,
) -> dict[str, pd.DataFrame]:
    rules = rules or build_clonality_qc_rules()
    ladder_rows: list[dict[str, Any]] = []
    pk_rows: list[dict[str, Any]] = []

    for entry in entries:
        join_fields = build_tracking_join_fields(entry)
        if not join_fields:
            continue

        assay = str(join_fields.get("assay") or "").strip()
        if not include_sl and assay.upper() in _KEEP_ASSAYS:
            continue

        ladder_rows.extend(_build_ladder_candidate_rows(entry, join_fields))
        if str(join_fields.get("control") or "").strip() in {"PK", "PK1", "PK2"}:
            pk_rows.extend(_build_pk_candidate_rows(entry, rules, join_fields))

    ladder_df = _finalize_frame(ladder_rows, table="ladder_candidates")
    pk_df = _finalize_frame(pk_rows, table="pk_candidates")
    combined_df = _finalize_frame(ladder_rows + pk_rows, table="combined")
    return {
        "combined": combined_df,
        "ladder_candidates": ladder_df,
        "pk_candidates": pk_df,
    }


def write_clonality_candidate_artifacts(
    output_dir: Path,
    entries: Sequence[dict[str, Any]],
    *,
    rules: QCRules | None = None,
    include_sl: bool = False,
    write_gold_label_template: bool = False,
) -> dict[str, Path]:
    tables = build_clonality_candidate_tables(
        entries,
        rules=rules,
        include_sl=include_sl,
    )
    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    combined_path = output_dir / CANDIDATE_ARTIFACT_FILENAMES["combined"]
    ladder_path = output_dir / CANDIDATE_ARTIFACT_FILENAMES["ladder_candidates"]
    pk_path = output_dir / CANDIDATE_ARTIFACT_FILENAMES["pk_candidates"]
    manifest_path = output_dir / CANDIDATE_ARTIFACT_FILENAMES["manifest"]

    tables["combined"].to_csv(combined_path, index=False)
    tables["ladder_candidates"].to_csv(ladder_path, index=False)
    tables["pk_candidates"].to_csv(pk_path, index=False)

    outputs = {
        "combined": combined_path,
        "ladder_candidates": ladder_path,
        "pk_candidates": pk_path,
        "manifest": manifest_path,
    }
    if write_gold_label_template:
        outputs["gold_labels"] = write_clonality_gold_label_template(output_dir)

    manifest = {
        "version": CANDIDATE_ARTIFACT_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "include_sl": bool(include_sl),
        "entry_count": int(len(entries)),
        "row_counts": {name: int(len(df)) for name, df in tables.items()},
        "output_files": {name: str(path) for name, path in outputs.items() if name != "manifest"},
        "assays": sorted(
            {
                str(value).strip()
                for value in pd.concat(
                    [
                        tables["ladder_candidates"].get("assay", pd.Series(dtype=str)),
                        tables["pk_candidates"].get("assay", pd.Series(dtype=str)),
                    ],
                    ignore_index=True,
                ).tolist()
                if str(value).strip()
            }
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return outputs


def write_clonality_gold_label_template(output_dir: Path) -> Path:
    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / CANDIDATE_ARTIFACT_FILENAMES["gold_labels"]
    pd.DataFrame(columns=GOLD_LABEL_COLUMNS).to_csv(path, index=False)
    return path


def _build_ladder_candidate_rows(entry: dict[str, Any], join_fields: dict[str, str]) -> list[dict[str, Any]]:
    fsa = entry.get("fsa")
    if fsa is None:
        return []

    candidates = get_ladder_candidates(fsa)
    if candidates.empty:
        return []

    identity_key = str(join_fields.get("identity_key") or "")
    assay = str(join_fields.get("assay") or "")
    run_code = str(join_fields.get("run_code") or "")
    well = str(join_fields.get("well") or "")
    ladder = str(join_fields.get("ladder") or "")
    join_key = build_tracking_join_key(
        identity_key=identity_key,
        assay=assay,
        run_code=run_code,
        well=well,
    )
    ladder_join_key = build_tracking_ladder_join_key(
        identity_key=identity_key,
        assay=assay,
        run_code=run_code,
        well=well,
        ladder=ladder,
    )
    selected_pairs = list(
        zip(
            np.asarray(getattr(fsa, "best_size_standard", []), dtype=float).tolist(),
            np.asarray(getattr(fsa, "ladder_steps", []), dtype=float).tolist(),
        )
    )

    rows: list[dict[str, Any]] = []
    for _, raw in candidates.iterrows():
        candidate_time = float(raw.get("time", np.nan))
        selected_step = _selected_step_bp(candidate_time, selected_pairs)
        candidate_index = int(raw.get("index", len(rows)))
        rows.append(
            {
                "artifact_kind": "ladder_candidate",
                "artifact_row_key": (
                    f"{build_tracking_row_key(artifact_kind='ladder', identity_key=identity_key)}"
                    f":candidate:{candidate_index}"
                ),
                "feature_source": "entry_candidates",
                "identity_key": identity_key,
                "join_key": join_key,
                "ladder_join_key": ladder_join_key,
                "source_run_dir": str(join_fields.get("source_run_dir") or ""),
                "assay": assay,
                "control": str(join_fields.get("control") or ""),
                "sample_kind": str(join_fields.get("sample_kind") or ""),
                "run_date": str(join_fields.get("run_date") or ""),
                "run_code": run_code,
                "well": well,
                "ladder": ladder,
                "candidate_index": candidate_index,
                "candidate_time": candidate_time,
                "candidate_intensity": _coerce_float(raw.get("intensity")),
                "candidate_source": str(raw.get("source") or ""),
                "selected_for_fit": selected_step is not None,
                "selected_step_bp": selected_step,
                "ladder_fit_strategy": str(entry.get("ladder_fit_strategy") or ""),
                "ladder_r2": _coerce_float(entry.get("ladder_r2")),
                "ladder_review_required": _ladder_review_required(entry, fsa),
            }
        )
    return rows


def _build_pk_candidate_rows(
    entry: dict[str, Any],
    rules: QCRules,
    join_fields: dict[str, str],
) -> list[dict[str, Any]]:
    fsa = entry.get("fsa")
    if fsa is None:
        return []

    identity_key = str(join_fields.get("identity_key") or "")
    assay = str(join_fields.get("assay") or "")
    run_code = str(join_fields.get("run_code") or "")
    well = str(join_fields.get("well") or "")
    join_key = build_tracking_join_key(
        identity_key=identity_key,
        assay=assay,
        run_code=run_code,
        well=well,
    )
    primary_channel = str(entry.get("primary_peak_channel") or "")
    rows: list[dict[str, Any]] = []

    for marker in markers_for_entry(entry, rules):
        channel = primary_channel if marker["channel"] == "primary" else str(marker["channel"])
        if marker["kind"] == "sample":
            evaluation = evaluate_peak_near_bp_with_fallback(
                fsa=fsa,
                channel=channel,
                target_bp=float(marker["expected_bp"]),
                window_bp=float(marker["window_bp"]),
                fallback_window_bp=float(getattr(rules, "sample_peak_window_bp_fallback", marker["window_bp"])),
                baseline_correct=True,
                name=str(marker.get("name") or ""),
            )
            candidates = list(evaluation.get("candidates") or [])
            selected_index = int(evaluation.get("selected_index", 0) or 0)
        else:
            candidates = [
                find_peak_near_bp(
                    fsa=fsa,
                    channel=channel,
                    target_bp=float(marker["expected_bp"]),
                    window_bp=float(marker["window_bp"]),
                    baseline_correct=True,
                )
            ]
            selected_index = 0

        pk_join_key = build_tracking_pk_join_key(
            identity_key=identity_key,
            assay=assay,
            run_code=run_code,
            well=well,
            marker_name=str(marker.get("name") or ""),
        )
        base_row_key = build_tracking_row_key(
            artifact_kind="pk",
            identity_key=identity_key,
            marker_name=str(marker.get("name") or ""),
        )
        for candidate_index, candidate in enumerate(candidates):
            candidate = dict(candidate or {})
            found_bp = _coerce_float(candidate.get("found_bp"))
            expected_bp = float(marker.get("expected_bp", np.nan))
            delta_bp = None if found_bp is None else float(found_bp - expected_bp)
            rows.append(
                {
                    "artifact_kind": "pk_candidate",
                    "artifact_row_key": f"{base_row_key}:candidate:{candidate_index}",
                    "feature_source": "entry_candidates",
                    "identity_key": identity_key,
                    "join_key": join_key,
                    "pk_join_key": pk_join_key,
                    "source_run_dir": str(join_fields.get("source_run_dir") or ""),
                    "assay": assay,
                    "control": str(join_fields.get("control") or ""),
                    "sample_kind": str(join_fields.get("sample_kind") or ""),
                    "run_date": str(join_fields.get("run_date") or ""),
                    "run_code": run_code,
                    "well": well,
                    "marker_name": str(marker.get("name") or ""),
                    "kind": str(marker.get("kind") or ""),
                    "channel": channel,
                    "expected_bp": expected_bp,
                    "window_bp": float(marker.get("window_bp", np.nan)),
                    "search_mode": str(candidate.get("search_mode") or ""),
                    "search_window_bp": _coerce_float(candidate.get("search_window_bp")),
                    "ok": bool(candidate.get("ok", False)),
                    "reason": str(candidate.get("reason") or ""),
                    "found_bp": found_bp,
                    "delta_bp": delta_bp,
                    "height": _coerce_float(candidate.get("height")),
                    "area": _coerce_float(candidate.get("area")),
                    "selection_score": _coerce_float(candidate.get("selection_score")),
                    "fallback_from_window_bp": _coerce_float(candidate.get("fallback_from_window_bp")),
                    "selected": candidate_index == selected_index,
                }
            )
    return rows


def _selected_step_bp(candidate_time: float, selected_pairs: list[tuple[float, float]]) -> float | None:
    for selected_time, step_bp in selected_pairs:
        if np.isclose(candidate_time, selected_time, atol=1e-6):
            return float(step_bp)
    return None


def _ladder_review_required(entry: dict[str, Any], fsa: Any) -> bool:
    if "ladder_review_required" in entry:
        return bool(entry.get("ladder_review_required"))
    if hasattr(fsa, "ladder_review_required"):
        return bool(getattr(fsa, "ladder_review_required"))
    status = str(entry.get("ladder_qc_status") or "").strip().lower()
    return bool(status and status != "ok")


def _finalize_frame(rows: list[dict[str, Any]], *, table: str) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df.insert(0, "artifact_version", CANDIDATE_ARTIFACT_VERSION)
    df.insert(1, "artifact_table", table)
    sort_cols = [col for col in ("assay", "identity_key", "marker_name", "artifact_row_key") if col in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols, kind="stable").reset_index(drop=True)
    return df


def _coerce_float(value: Any) -> float | None:
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


__all__ = [
    "CANDIDATE_ARTIFACT_FILENAMES",
    "CANDIDATE_ARTIFACT_VERSION",
    "GOLD_LABEL_COLUMNS",
    "build_clonality_candidate_tables",
    "write_clonality_candidate_artifacts",
    "write_clonality_gold_label_template",
]
