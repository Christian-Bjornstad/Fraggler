from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from core.analyses.clonality.tracking_excel import (
    build_tracking_join_key,
    build_tracking_ladder_join_key,
    build_tracking_pk_join_key,
    build_tracking_row_key,
)


FEATURE_ARTIFACT_VERSION = "v1"
FEATURE_ARTIFACT_FILENAMES = {
    "combined": "clonality_feature_artifacts.csv",
    "ladder": "clonality_ladder_features.csv",
    "pk": "clonality_pk_features.csv",
    "manifest": "clonality_feature_manifest.json",
}

_LADDER_KEEP_ASSAYS = {"SL"}


def write_clonality_feature_artifacts(
    excel_path: Path,
    output_dir: Path,
    *,
    entry_metadata: Mapping[str, Mapping[str, Any]] | Sequence[Mapping[str, Any]] | None = None,
    include_sl: bool = False,
) -> dict[str, Path]:
    tables = build_clonality_feature_tables(
        excel_path,
        entry_metadata=entry_metadata,
        include_sl=include_sl,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    combined_path = output_dir / FEATURE_ARTIFACT_FILENAMES["combined"]
    ladder_path = output_dir / FEATURE_ARTIFACT_FILENAMES["ladder"]
    pk_path = output_dir / FEATURE_ARTIFACT_FILENAMES["pk"]
    manifest_path = output_dir / FEATURE_ARTIFACT_FILENAMES["manifest"]

    tables["combined"].to_csv(combined_path, index=False)
    tables["ladder"].to_csv(ladder_path, index=False)
    tables["pk"].to_csv(pk_path, index=False)

    manifest = {
        "version": FEATURE_ARTIFACT_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "workbook": str(Path(excel_path).expanduser()),
        "include_sl": bool(include_sl),
        "row_counts": {name: int(len(df)) for name, df in tables.items()},
        "output_files": {
            "combined": str(combined_path),
            "ladder": str(ladder_path),
            "pk": str(pk_path),
        },
        "assays": sorted(
            {
                str(v).strip()
                for v in pd.concat(
                    [
                        tables["ladder"].get("assay", pd.Series(dtype=str)),
                        tables["pk"].get("assay", pd.Series(dtype=str)),
                    ],
                    ignore_index=True,
                ).tolist()
                if str(v).strip()
            }
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "combined": combined_path,
        "ladder": ladder_path,
        "pk": pk_path,
        "manifest": manifest_path,
    }


def build_clonality_feature_tables(
    excel_path: Path,
    *,
    entry_metadata: Mapping[str, Mapping[str, Any]] | Sequence[Mapping[str, Any]] | None = None,
    include_sl: bool = False,
) -> dict[str, pd.DataFrame]:
    workbook = Path(excel_path).expanduser()
    xls = pd.ExcelFile(workbook, engine="openpyxl")

    patient = _read_sheet(xls, "Patient_Runs")
    control = _read_sheet(xls, "Control_Runs")
    peaks = _read_sheet(xls, "PK_Peaks")
    metadata_index = _normalize_entry_metadata(entry_metadata)

    ladder_rows = _build_ladder_rows(patient, scope="Patient", metadata_index=metadata_index, include_sl=include_sl)
    ladder_rows.extend(_build_ladder_rows(control, scope="Control", metadata_index=metadata_index, include_sl=include_sl))
    pk_rows = _build_pk_rows(peaks, metadata_index=metadata_index, include_sl=include_sl)

    ladder_df = _finalize_frame(ladder_rows, kind="ladder")
    pk_df = _finalize_frame(pk_rows, kind="pk")
    combined_df = _finalize_frame(ladder_rows + pk_rows, kind="combined")

    return {
        "combined": combined_df,
        "ladder": ladder_df,
        "pk": pk_df,
    }


def _read_sheet(xls: pd.ExcelFile, sheet_name: str) -> pd.DataFrame:
    if sheet_name not in xls.sheet_names:
        return pd.DataFrame()
    return pd.read_excel(xls, sheet_name=sheet_name, engine="openpyxl").fillna("")


def _normalize_entry_metadata(
    entry_metadata: Mapping[str, Mapping[str, Any]] | Sequence[Mapping[str, Any]] | None,
) -> dict[str, dict[str, Any]]:
    if entry_metadata is None:
        return {}
    if isinstance(entry_metadata, Mapping):
        return {str(key): dict(value) for key, value in entry_metadata.items()}

    index: dict[str, dict[str, Any]] = {}
    for item in entry_metadata:
        if not isinstance(item, Mapping):
            continue
        payload = dict(item)
        for key in _metadata_keys_from_payload(payload):
            index.setdefault(key, payload)
    return index


def _metadata_keys_from_payload(payload: Mapping[str, Any]) -> list[str]:
    keys: list[str] = []
    for candidate in (
        payload.get("identity_key"),
        payload.get("IdentityKey"),
        payload.get("source_run_dir"),
        payload.get("SourceRunDir"),
        payload.get("file_name"),
        payload.get("File"),
    ):
        value = str(candidate or "").strip()
        if value:
            keys.append(value)

    source_run_dir = str(payload.get("source_run_dir") or payload.get("SourceRunDir") or "").strip()
    file_name = str(payload.get("file_name") or payload.get("File") or "").strip()
    if source_run_dir and file_name:
        keys.append(f"{source_run_dir}::{file_name}")
    return keys


def _build_ladder_rows(
    frame: pd.DataFrame,
    *,
    scope: str,
    metadata_index: Mapping[str, Mapping[str, Any]],
    include_sl: bool,
) -> list[dict[str, Any]]:
    if frame.empty:
        return []

    rows: list[dict[str, Any]] = []
    for _, raw in frame.iterrows():
        assay = _clean_text(raw, "Assay")
        if not include_sl and assay.upper() in _LADDER_KEEP_ASSAYS:
            continue

        identity_key = _clean_text(raw, "IdentityKey")
        run_code = _clean_text(raw, "RunCode")
        well = _clean_text(raw, "Well")
        ladder = _clean_text(raw, "Ladder")
        join_key = build_tracking_join_key(
            identity_key=identity_key,
            assay=assay,
            run_code=run_code,
            well=well,
        )
        row = {
            "artifact_kind": "ladder",
            "artifact_row_key": build_tracking_row_key(artifact_kind="ladder", identity_key=identity_key),
            "scope": scope,
            "identity_key": identity_key,
            "join_key": join_key,
            "ladder_join_key": build_tracking_ladder_join_key(
                identity_key=identity_key,
                assay=assay,
                run_code=run_code,
                well=well,
                ladder=ladder,
            ),
            "source_run_dir": _clean_text(raw, "SourceRunDir"),
            "assay": assay,
            "control": _clean_text(raw, "Control"),
            "sample_kind": _clean_text(raw, "SampleKind"),
            "run_date": _clean_text(raw, "RunDate"),
            "run_code": run_code,
            "well": well,
            "ladder": ladder,
            "ladder_qc": _clean_text(raw, "LadderQC"),
            "ladder_fit_strategy": _clean_text(raw, "LadderFitStrategy"),
            "ladder_expected_step_count": _coerce_int(raw.get("LadderExpectedStepCount")),
            "ladder_fitted_step_count": _coerce_int(raw.get("LadderFittedStepCount")),
            "ladder_step_gap": _step_gap(raw),
            "ladder_step_gap_abs": abs(_step_gap(raw)),
            "ladder_step_gap_ratio": _step_gap_ratio(raw),
            "ladder_r2": _coerce_float(raw.get("LadderR2")),
            "ladder_review_required": _ladder_review_required(raw),
            "is_sl": assay.upper() == "SL",
            "is_training_candidate": assay.upper() != "SL",
            "fit_is_partial": "partial" in _clean_text(raw, "LadderFitStrategy").lower(),
            "fit_is_rescue": "rescue" in _clean_text(raw, "LadderFitStrategy").lower(),
            "fit_is_manual": "manual" in _clean_text(raw, "LadderFitStrategy").lower(),
            "r2_band": _r2_band(_coerce_float(raw.get("LadderR2"))),
            "feature_source": "workbook",
        }
        row.update(_metadata_features(identity_key, metadata_index))
        rows.append(row)
    return rows


def _build_pk_rows(
    frame: pd.DataFrame,
    *,
    metadata_index: Mapping[str, Mapping[str, Any]],
    include_sl: bool,
) -> list[dict[str, Any]]:
    if frame.empty:
        return []

    rows: list[dict[str, Any]] = []
    for _, raw in frame.iterrows():
        assay = _clean_text(raw, "Assay")
        if not include_sl and assay.upper() in _LADDER_KEEP_ASSAYS:
            continue

        kind = _clean_text(raw, "Kind").lower()
        if kind and kind not in {"sample", "ladder"}:
            continue

        identity_key = _clean_text(raw, "IdentityKey")
        marker_name = _clean_text(raw, "MarkerName")
        run_code = _clean_text(raw, "RunCode")
        well = _clean_text(raw, "Well")
        join_key = build_tracking_join_key(
            identity_key=identity_key,
            assay=assay,
            run_code=run_code,
            well=well,
        )
        row = {
            "artifact_kind": "pk",
            "artifact_row_key": build_tracking_row_key(
                artifact_kind="pk",
                identity_key=identity_key,
                marker_name=marker_name,
            ),
            "scope": _artifact_scope(raw),
            "identity_key": identity_key,
            "join_key": join_key,
            "pk_join_key": build_tracking_pk_join_key(
                identity_key=identity_key,
                assay=assay,
                run_code=run_code,
                well=well,
                marker_name=marker_name,
            ),
            "source_run_dir": _clean_text(raw, "SourceRunDir"),
            "assay": assay,
            "control": _clean_text(raw, "Control"),
            "run_date": _clean_text(raw, "RunDate"),
            "run_code": run_code,
            "well": well,
            "marker_name": marker_name,
            "kind": kind,
            "channel": _clean_text(raw, "Channel"),
            "expected_bp": _coerce_float(raw.get("ExpectedBP")),
            "window_bp": _coerce_float(raw.get("WindowBP")),
            "search_mode": _clean_text(raw, "SearchMode"),
            "search_window_bp": _coerce_float(raw.get("SearchWindowBP")),
            "found_bp": _coerce_float(raw.get("FoundBP")),
            "delta_bp": _coerce_float(raw.get("DeltaBP")),
            "height": _coerce_float(raw.get("Height")),
            "area": _coerce_float(raw.get("Area")),
            "ok": _normalize_bool(raw.get("OK")),
            "reason": _clean_text(raw, "Reason"),
            "abs_delta_bp": _abs_delta_bp(raw),
            "delta_bucket": _delta_bucket(raw),
            "pk_review_required": _pk_review_required(raw),
            "is_sl": assay.upper() == "SL",
            "is_training_candidate": assay.upper() != "SL",
            "is_fallback": _clean_text(raw, "SearchMode").lower() == "fallback",
            "feature_source": "workbook",
        }
        row.update(_metadata_features(identity_key, metadata_index))
        rows.append(row)
    return rows


def _metadata_features(identity_key: str, metadata_index: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    if not identity_key:
        return {"metadata_source": ""}

    meta = metadata_index.get(identity_key)
    if meta is None:
        return {"metadata_source": ""}

    return {
        "metadata_source": "entry_metadata",
        "metadata_primary_peak_channel": _clean_value(meta.get("primary_peak_channel")),
        "metadata_ladder": _clean_value(meta.get("ladder")),
        "metadata_sample_kind": _clean_value(meta.get("sample_kind")),
    }


def _finalize_frame(rows: list[dict[str, Any]], *, kind: str) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df.insert(0, "artifact_version", FEATURE_ARTIFACT_VERSION)
    df.insert(1, "artifact_table", kind)
    sort_cols = [col for col in ("assay", "identity_key", "marker_name", "artifact_row_key") if col in df.columns]
    if sort_cols:
        df = df.sort_values(sort_cols, kind="stable").reset_index(drop=True)
    return df


def _artifact_scope(raw: Mapping[str, Any]) -> str:
    value = _clean_text(raw, "SampleKind")
    if value:
        return value
    control = _clean_text(raw, "Control")
    return "control" if control else "patient"


def _ladder_review_required(raw: Mapping[str, Any]) -> bool:
    status = _clean_text(raw, "LadderQC").strip().lower()
    return bool(status and status != "ok")


def _pk_review_required(raw: Mapping[str, Any]) -> bool:
    if not _normalize_bool(raw.get("OK")):
        return True
    abs_delta = _abs_delta_bp(raw)
    return bool(abs_delta and abs_delta > 2.0)


def _delta_bucket(raw: Mapping[str, Any]) -> str:
    abs_delta = _abs_delta_bp(raw)
    if abs_delta <= 2.0:
        return "within_2bp"
    if abs_delta <= 5.0:
        return "within_5bp"
    return "over_5bp"


def _step_gap(raw: Mapping[str, Any]) -> float:
    expected = _coerce_float(raw.get("LadderExpectedStepCount"))
    fitted = _coerce_float(raw.get("LadderFittedStepCount"))
    if expected is None or fitted is None:
        return 0.0
    return float(fitted - expected)


def _step_gap_ratio(raw: Mapping[str, Any]) -> float:
    expected = _coerce_float(raw.get("LadderExpectedStepCount"))
    if not expected:
        return 0.0
    return float(_step_gap(raw) / expected)


def _r2_band(value: float | None) -> str:
    if value is None:
        return "missing"
    if value >= 0.999:
        return "ok"
    if value >= 0.995:
        return "warn"
    return "fail"


def _abs_delta_bp(raw: Mapping[str, Any]) -> float:
    value = _coerce_float(raw.get("AbsDeltaBP"))
    if value is not None:
        return float(value)
    delta = _coerce_float(raw.get("DeltaBP"))
    if delta is not None:
        return abs(float(delta))
    found = _coerce_float(raw.get("FoundBP"))
    expected = _coerce_float(raw.get("ExpectedBP"))
    if found is not None and expected is not None:
        return abs(float(found) - float(expected))
    return 0.0


def _coerce_float(value: Any) -> float | None:
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_int(value: Any) -> int | None:
    if value in ("", None):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _normalize_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value in ("", None):
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y", "ok"}


def _clean_text(raw: Mapping[str, Any], key: str) -> str:
    return _clean_value(raw.get(key))


def _clean_value(value: Any) -> str:
    if value in ("", None):
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _make_row_key(prefix: str, value: str) -> str:
    safe = _clean_value(value).replace(" ", "_")
    return f"{prefix}:{safe}" if safe else prefix


__all__ = [
    "FEATURE_ARTIFACT_FILENAMES",
    "FEATURE_ARTIFACT_VERSION",
    "build_clonality_feature_tables",
    "write_clonality_feature_artifacts",
]
