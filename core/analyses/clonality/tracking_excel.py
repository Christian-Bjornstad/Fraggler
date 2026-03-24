from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from config import APP_SETTINGS
from fraggler.fraggler import print_green
from core.analyses.clonality.tracking_dashboard import refresh_clonality_tracking_dashboard
from core.qc.qc_markers import (
    control_id_from_filename,
    find_peak_near_bp,
    find_peak_near_bp_with_fallback,
    make_run_key,
    markers_for_entry,
    parse_batch_from_filename,
    parse_pcr_date_from_filename,
    parse_run_code_from_filename,
    parse_well_from_filename,
)
from core.qc.qc_rules import QCRules
from core.utils import strip_stage_prefix

import threading

CLONALITY_TRACKING_FILENAME = "Clonality_Tracking.xlsx"
_clonality_excel_lock = threading.Lock()
CONTROL_IDS = {"PK", "PK1", "PK2", "NK", "RK"}
RUN_SHEET_COLUMNS = [
    "IdentityKey",
    "File",
    "SourceRunDir",
    "DIT",
    "Assay",
    "SampleKind",
    "Group",
    "Control",
    "RunDate",
    "RunCode",
    "Well",
    "Batch",
    "Ladder",
    "LadderQC",
    "LadderFitStrategy",
    "LadderExpectedStepCount",
    "LadderFittedStepCount",
    "LadderR2",
]
PEAK_SHEET_COLUMNS = [
    "IdentityKey",
    "File",
    "SourceRunDir",
    "DIT",
    "Assay",
    "Control",
    "RunDate",
    "RunCode",
    "Well",
    "Batch",
    "MarkerName",
    "Kind",
    "Channel",
    "ExpectedBP",
    "WindowBP",
    "SearchMode",
    "SearchWindowBP",
    "FoundBP",
    "DeltaBP",
    "Height",
    "Area",
    "OK",
    "Reason",
    "AbsDeltaBP",
]


def build_clonality_qc_rules() -> QCRules:
    qc_settings = APP_SETTINGS.get("qc", {})
    sample_window = qc_settings.get("sample_peak_window_bp", qc_settings.get("w_sample", 3.0))
    ladder_window = qc_settings.get("ladder_peak_window_bp", qc_settings.get("w_ladder", 3.0))
    return QCRules(
        min_r2_ok=qc_settings.get("min_r2_ok", 0.999),
        min_r2_warn=qc_settings.get("min_r2_warn", 0.995),
        nk_ymax_floor=qc_settings.get("nk_ymax_floor", 250.0),
        sample_peak_window_bp=sample_window,
        sample_peak_window_bp_fallback=qc_settings.get("sample_peak_window_bp_fallback", max(float(sample_window) + 4.0, 8.0)),
        ladder_peak_window_bp=ladder_window,
        min_sl_total_area=qc_settings.get("min_sl_total_area", 1e4),
    )


def resolve_original_input_path(path_like: Path | str | None) -> Path | None:
    if not path_like:
        return None

    path = Path(path_like).expanduser()
    try:
        if path.is_symlink():
            return path.resolve(strict=False)
    except OSError:
        return None

    try:
        resolved = path.resolve(strict=False)
    except OSError:
        resolved = path

    if path.parent.name.startswith("fraggler_stage_") and resolved.parent == path.parent:
        return None
    return resolved


def resolve_source_run_dir(entry: dict) -> str:
    existing = str(entry.get("source_run_dir") or "").strip()
    if existing:
        return existing

    fsa = entry.get("fsa")
    original_path = resolve_original_input_path(getattr(fsa, "file", None))
    if original_path is not None and original_path.parent.name:
        return original_path.parent.name

    file_name = str(getattr(fsa, "file_name", "") or entry.get("file_name") or "")
    run_key = make_run_key(file_name)
    if run_key and run_key != "UNKNOWN":
        return run_key

    path = Path(getattr(fsa, "file", "")) if getattr(fsa, "file", None) else None
    if path is not None and path.parent.name and not path.parent.name.startswith("fraggler_stage_"):
        return path.parent.name

    return ""


def update_clonality_tracking_workbook(
    excel_path: Path,
    entries: list[dict],
    rules: QCRules | None = None,
) -> None:
    excel_path.parent.mkdir(parents=True, exist_ok=True)

    rules = rules or build_clonality_qc_rules()
    df_patient, df_control, df_peaks, pk_identity_keys = _build_tracking_frames(entries, rules)
    if df_patient.empty and df_control.empty and df_peaks.empty:
        return

    with _clonality_excel_lock:
        if excel_path.exists():
            try:
                with pd.ExcelFile(excel_path, engine="openpyxl") as xls:
                    has_patient = "Patient_Runs" in xls.sheet_names
                    has_control = "Control_Runs" in xls.sheet_names
                    has_peaks = "PK_Peaks" in xls.sheet_names
            except Exception:
                has_patient = False
                has_control = False
                has_peaks = False
                from fraggler.fraggler import print_warning # Import if not already globally available
                print_warning(f"Kunne ikke lese eksisterende {excel_path.name}, kanskje korrupt. Lager ny...")

            try:
                old_patient = pd.read_excel(excel_path, sheet_name="Patient_Runs", engine="openpyxl") if has_patient else pd.DataFrame(columns=RUN_SHEET_COLUMNS)
                old_control = pd.read_excel(excel_path, sheet_name="Control_Runs", engine="openpyxl") if has_control else pd.DataFrame(columns=RUN_SHEET_COLUMNS)
                old_peaks = pd.read_excel(excel_path, sheet_name="PK_Peaks", engine="openpyxl") if has_peaks else pd.DataFrame(columns=PEAK_SHEET_COLUMNS)
            except Exception:
                old_patient = pd.DataFrame(columns=RUN_SHEET_COLUMNS)
                old_control = pd.DataFrame(columns=RUN_SHEET_COLUMNS)
                old_peaks = pd.DataFrame(columns=PEAK_SHEET_COLUMNS)
        else:
            old_patient = pd.DataFrame(columns=RUN_SHEET_COLUMNS)
            old_control = pd.DataFrame(columns=RUN_SHEET_COLUMNS)
            old_peaks = pd.DataFrame(columns=PEAK_SHEET_COLUMNS)

        if not df_patient.empty and "IdentityKey" in old_patient.columns:
            old_patient = old_patient[~old_patient["IdentityKey"].isin(df_patient["IdentityKey"])]
        if not df_control.empty and "IdentityKey" in old_control.columns:
            old_control = old_control[~old_control["IdentityKey"].isin(df_control["IdentityKey"])]
        if pk_identity_keys and "IdentityKey" in old_peaks.columns:
            old_peaks = old_peaks[~old_peaks["IdentityKey"].isin(sorted(pk_identity_keys))]

        all_patient = _concat_frames(old_patient, df_patient)
        all_control = _concat_frames(old_control, df_control)
        all_peaks = _concat_frames(old_peaks, df_peaks)

        if not all_patient.empty and "IdentityKey" in all_patient.columns:
            all_patient = all_patient.drop_duplicates(subset=["IdentityKey"], keep="last")
        if not all_control.empty and "IdentityKey" in all_control.columns:
            all_control = all_control.drop_duplicates(subset=["IdentityKey"], keep="last")
        if not all_peaks.empty and {"IdentityKey", "MarkerName"}.issubset(all_peaks.columns):
            all_peaks = all_peaks.drop_duplicates(subset=["IdentityKey", "MarkerName"], keep="last")

        all_patient = _reindex_columns(all_patient, RUN_SHEET_COLUMNS)
        all_control = _reindex_columns(all_control, RUN_SHEET_COLUMNS)
        all_peaks = _reindex_columns(all_peaks, PEAK_SHEET_COLUMNS)

        writer_kwargs = {"engine": "openpyxl"}
        if excel_path.exists():
            writer_kwargs.update({"mode": "a", "if_sheet_exists": "replace"})

        with pd.ExcelWriter(excel_path, **writer_kwargs) as writer:
            all_patient.to_excel(writer, sheet_name="Patient_Runs", index=False)
            all_control.to_excel(writer, sheet_name="Control_Runs", index=False)
            all_peaks.to_excel(writer, sheet_name="PK_Peaks", index=False)
        refresh_clonality_tracking_dashboard(excel_path)
        print_green(f"Clonality tracking workbook updated in {excel_path}")


def _build_tracking_frames(entries: list[dict], rules: QCRules) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, set[str]]:
    patient_rows: list[dict] = []
    control_rows: list[dict] = []
    peak_rows: list[dict] = []
    pk_identity_keys: set[str] = set()

    for entry in entries:
        base_row = _build_run_row(entry)
        if not base_row:
            continue

        if base_row["SampleKind"] == "patient":
            patient_rows.append(base_row)
            continue

        control_rows.append(base_row)
        if base_row["Control"] not in {"PK", "PK1", "PK2"}:
            continue

        pk_identity_keys.add(base_row["IdentityKey"])
        peak_rows.extend(_build_pk_peak_rows(entry, rules, base_row))

    return (
        _reindex_columns(pd.DataFrame(patient_rows), RUN_SHEET_COLUMNS),
        _reindex_columns(pd.DataFrame(control_rows), RUN_SHEET_COLUMNS),
        _reindex_columns(pd.DataFrame(peak_rows), PEAK_SHEET_COLUMNS),
        pk_identity_keys,
    )


def _build_run_row(entry: dict) -> dict:
    fsa = entry.get("fsa")
    file_name = _resolve_entry_file_name(entry)
    if not file_name:
        return {}

    source_run_dir = resolve_source_run_dir(entry)
    control = control_id_from_filename(file_name)
    is_control = control in CONTROL_IDS
    identity_key = f"{source_run_dir}::{file_name}" if source_run_dir else file_name
    ladder_r2 = entry.get("ladder_r2")
    if ladder_r2 is None or not np.isfinite(ladder_r2):
        ladder_r2 = ""

    return {
        "IdentityKey": identity_key,
        "File": file_name,
        "SourceRunDir": source_run_dir,
        "DIT": entry.get("dit") or "",
        "Assay": entry.get("assay") or "",
        "SampleKind": "control" if is_control else "patient",
        "Group": entry.get("group") or "",
        "Control": control if is_control else "",
        "RunDate": parse_pcr_date_from_filename(file_name) or "",
        "RunCode": parse_run_code_from_filename(file_name) or "",
        "Well": parse_well_from_filename(file_name) or "",
        "Batch": parse_batch_from_filename(file_name) or "",
        "Ladder": entry.get("ladder") or "",
        "LadderQC": entry.get("ladder_qc_status") or "",
        "LadderFitStrategy": entry.get("ladder_fit_strategy") or "",
        "LadderExpectedStepCount": int(entry.get("ladder_expected_step_count", 0) or 0),
        "LadderFittedStepCount": int(entry.get("ladder_fitted_step_count", 0) or 0),
        "LadderR2": ladder_r2,
    }


def _build_pk_peak_rows(entry: dict, rules: QCRules, base_row: dict) -> list[dict]:
    marker_rows: list[dict] = []
    fsa = entry["fsa"]
    primary_channel = str(entry.get("primary_peak_channel") or "")

    for marker in markers_for_entry(entry, rules):
        channel = primary_channel if marker["channel"] == "primary" else str(marker["channel"])
        if marker["kind"] == "sample":
            result = find_peak_near_bp_with_fallback(
                fsa=fsa,
                channel=channel,
                target_bp=float(marker["expected_bp"]),
                window_bp=float(marker["window_bp"]),
                fallback_window_bp=float(getattr(rules, "sample_peak_window_bp_fallback", marker["window_bp"])),
                baseline_correct=True,
            )
        else:
            result = find_peak_near_bp(
                fsa=fsa,
                channel=channel,
                target_bp=float(marker["expected_bp"]),
                window_bp=float(marker["window_bp"]),
                baseline_correct=True,
            )

        row = {
            "IdentityKey": base_row["IdentityKey"],
            "File": base_row["File"],
            "SourceRunDir": base_row["SourceRunDir"],
            "DIT": base_row["DIT"],
            "Assay": base_row["Assay"],
            "Control": base_row["Control"],
            "RunDate": base_row["RunDate"],
            "RunCode": base_row["RunCode"],
            "Well": base_row["Well"],
            "Batch": base_row["Batch"],
            "MarkerName": marker.get("name") or "",
            "Kind": marker.get("kind") or "",
            "Channel": channel,
            "ExpectedBP": float(marker.get("expected_bp", np.nan)),
            "WindowBP": float(marker.get("window_bp", np.nan)),
            "SearchMode": result.get("search_mode", ""),
            "SearchWindowBP": float(result.get("search_window_bp", np.nan)) if result.get("search_window_bp") is not None else "",
            "FoundBP": "",
            "DeltaBP": "",
            "Height": "",
            "Area": "",
            "OK": bool(result.get("ok", False)),
            "Reason": result.get("reason") or "",
            "AbsDeltaBP": "",
        }
        if result.get("ok", False):
            found_bp = float(result["found_bp"])
            delta_bp = found_bp - float(marker.get("expected_bp", np.nan))
            row["FoundBP"] = found_bp
            row["DeltaBP"] = delta_bp
            row["Height"] = float(result.get("height", np.nan))
            row["Area"] = float(result.get("area", np.nan))
            row["AbsDeltaBP"] = abs(delta_bp)
        marker_rows.append(row)

    return marker_rows


def _reindex_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=columns)
    return df.reindex(columns=columns)


def _concat_frames(old_df: pd.DataFrame, new_df: pd.DataFrame) -> pd.DataFrame:
    frames = [df for df in (old_df, new_df) if not df.empty]
    if not frames:
        return pd.DataFrame(columns=(list(old_df.columns) if len(old_df.columns) else list(new_df.columns)))
    if len(frames) == 1:
        return frames[0].copy()
    return pd.concat(frames, ignore_index=True)


def _resolve_entry_file_name(entry: dict) -> str:
    fsa = entry.get("fsa")
    original_path = resolve_original_input_path(getattr(fsa, "file", None))
    if original_path is not None and original_path.name:
        return original_path.name
    return strip_stage_prefix(str(getattr(fsa, "file_name", "") or entry.get("file_name") or ""))
