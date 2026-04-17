from __future__ import annotations

import threading
from pathlib import Path

import pandas as pd

from fraggler.fraggler import print_green
from core.analyses.clonality.tracking_dashboard import refresh_clonality_tracking_dashboard
from core.qc.qc_markers import (
    make_run_key,
    parse_pcr_date_from_filename,
    parse_run_code_from_filename,
    parse_well_from_filename,
)


FLT3_NPM1_QC_TRACKER_FILENAME = "FLT3_NPM1_QC_TRACKER.xlsx"

PATIENT_RUN_SHEET_COLUMNS = [
    "Month",
    "IdentityKey",
    "File",
    "SourceRunDir",
    "DIT",
    "Assay",
    "SampleKind",
    "Group",
    "RunDate",
    "RunCode",
    "Well",
    "Ladder",
    "LadderQC",
    "LadderFitStrategy",
    "LadderExpectedStepCount",
    "LadderFittedStepCount",
    "LadderR2",
    "PeakQC",
]

CONTROL_RUN_SHEET_COLUMNS = [
    "Month",
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
    "Ladder",
    "LadderQC",
    "LadderFitStrategy",
    "LadderExpectedStepCount",
    "LadderFittedStepCount",
    "LadderR2",
    "PeakQC",
]

PEAK_SHEET_COLUMNS = [
    "Month",
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

_excel_lock = threading.Lock()

POSITIVE_CONTROL_IDS = {"IVS-P001", "IVS-P0001"}

MARKER_DEFINITIONS = {
    ("IVS-0000", "FLT3-ITD", None): [
        {
            "name": "IVS0000_ITD_329",
            "kind": "sample",
            "channel": "primary",
            "peak_label": "WT",
            "expected_bp": 329.0,
            "window_bp": 3.0,
            "delta_threshold_bp": 1.0,
            "analysis_label": "ITD",
        },
        {
            "name": "ITD_Ladder_350",
            "kind": "ladder",
            "channel": "size_standard",
            "expected_bp": 350.0,
            "window_bp": 2.5,
            "delta_threshold_bp": 1.5,
            "analysis_label": "ITD",
        },
    ],
    ("IVS-0000", "FLT3-ITD", "ratio_quant"): [
        {
            "name": "IVS0000_RATIO_329",
            "kind": "sample",
            "channel": "primary",
            "peak_label": "WT",
            "expected_bp": 329.0,
            "window_bp": 3.0,
            "delta_threshold_bp": 1.0,
            "analysis_label": "ITD-ratio",
        },
        {
            "name": "RATIO_Ladder_350",
            "kind": "ladder",
            "channel": "size_standard",
            "expected_bp": 350.0,
            "window_bp": 2.5,
            "delta_threshold_bp": 1.5,
            "analysis_label": "ITD-ratio",
        },
    ],
    ("IVS-0000", "FLT3-D835", None): [
        {
            "name": "IVS0000_D835_80",
            "kind": "sample",
            "channel": "primary",
            "peak_label": "WT",
            "expected_bp": 80.0,
            "window_bp": 2.5,
            "delta_threshold_bp": 1.0,
            "analysis_label": "D835",
        },
        {
            "name": "D835_Ladder_139",
            "kind": "ladder",
            "channel": "size_standard",
            "expected_bp": 139.0,
            "window_bp": 2.0,
            "delta_threshold_bp": 1.5,
            "analysis_label": "D835",
        },
    ],
    ("IVS-P001", "FLT3-D835", None): [
        {
            "name": "IVSP001_D835_128_129",
            "kind": "sample",
            "channel": "primary",
            "peak_label": "MUT",
            "expected_bp": 128.5,
            "window_bp": 2.5,
            "delta_threshold_bp": 1.0,
            "analysis_label": "D835",
            "expected_range": "128-129 bp",
        },
        {
            "name": "D835_Ladder_139",
            "kind": "ladder",
            "channel": "size_standard",
            "expected_bp": 139.0,
            "window_bp": 2.0,
            "delta_threshold_bp": 1.5,
            "analysis_label": "D835",
        },
    ],
    ("IVS-0000", "NPM1", None): [
        {
            "name": "IVS0000_NPM1_299",
            "kind": "sample",
            "channel": "primary",
            "peak_label": "WT",
            "expected_bp": 299.0,
            "window_bp": 2.5,
            "delta_threshold_bp": 1.0,
            "analysis_label": "NPM1",
        },
        {
            "name": "NPM1_Ladder_350",
            "kind": "ladder",
            "channel": "size_standard",
            "expected_bp": 350.0,
            "window_bp": 2.5,
            "delta_threshold_bp": 1.5,
            "analysis_label": "NPM1",
        },
    ],
}

def normalize_specimen_id(value: str | None) -> str:
    specimen = str(value or "").strip().upper()
    if specimen in POSITIVE_CONTROL_IDS:
        return "IVS-P001"
    return specimen


def control_code_for_entry(entry: dict) -> str:
    specimen = normalize_specimen_id(entry.get("specimen_id"))
    group = str(entry.get("group") or "")
    if specimen == "IVS-0000":
        return "RK"
    if specimen == "IVS-P001":
        return "PK"
    if group == "negative_control" or specimen == "NTC":
        return "NK"
    return ""


def is_tracking_control_entry(entry: dict) -> bool:
    return control_code_for_entry(entry) in {"RK", "PK", "NK"}


def marker_specs_for_entry(entry: dict) -> list[dict]:
    specimen = normalize_specimen_id(entry.get("specimen_id"))
    assay = str(entry.get("assay") or "")
    analysis_type = str(entry.get("analysis_type") or "")
    specs = MARKER_DEFINITIONS.get((specimen, assay, analysis_type))
    if specs is None:
        specs = MARKER_DEFINITIONS.get((specimen, assay, None), [])
    return [dict(spec) for spec in specs]


def build_tracking_base_row(entry: dict) -> dict:
    file_name = resolve_entry_file_name(entry)
    if not file_name:
        return {}

    control = control_code_for_entry(entry)
    sample_kind = "control" if control else "patient"
    source_run_dir = str(entry.get("source_run_dir") or "").strip()
    run_date = str(entry.get("run_date") or "").strip() or parse_pcr_date_from_filename(file_name) or ""
    run_code = parse_run_code_from_filename(file_name) or ""
    well = str(entry.get("well_id") or "").strip() or parse_well_from_filename(file_name) or ""
    identity_key = f"{source_run_dir}::{file_name}" if control else f"{source_run_dir}::{make_run_key(file_name)}::{well or file_name}"
    ladder_r2 = entry.get("ladder_r2")

    return {
        "Month": _month_bucket(run_date),
        "IdentityKey": identity_key,
        "File": file_name,
        "SourceRunDir": source_run_dir,
        "DIT": str(entry.get("dit") or ""),
        "Assay": str(entry.get("assay") or ""),
        "SampleKind": sample_kind,
        "Group": str(entry.get("group") or ""),
        "Control": control,
        "RunDate": run_date,
        "RunCode": run_code,
        "Well": well,
        "Batch": "",
        "Ladder": str(entry.get("ladder") or ""),
        "LadderQC": str(entry.get("ladder_qc_status") or ""),
        "LadderFitStrategy": str(entry.get("ladder_fit_strategy") or ""),
        "LadderExpectedStepCount": int(entry.get("ladder_expected_step_count", 0) or 0),
        "LadderFittedStepCount": int(entry.get("ladder_fitted_step_count", 0) or 0),
        "LadderR2": ladder_r2 if ladder_r2 is not None else "",
    }


def resolve_entry_file_name(entry: dict) -> str:
    fsa = entry.get("fsa")
    return str(getattr(fsa, "file_name", "") or entry.get("file_name") or "").strip()


def update_flt3_npm1_qc_tracker(
    excel_path: Path,
    patient_df: pd.DataFrame,
    control_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
) -> None:
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    patient_df = _reindex_columns(patient_df, PATIENT_RUN_SHEET_COLUMNS)
    control_df = _reindex_columns(control_df, CONTROL_RUN_SHEET_COLUMNS)
    peaks_df = _reindex_columns(peaks_df, PEAK_SHEET_COLUMNS)

    with _excel_lock:
        if excel_path.exists():
            try:
                with pd.ExcelFile(excel_path, engine="openpyxl") as xls:
                    old_patient = pd.read_excel(excel_path, sheet_name="Patient_Runs", engine="openpyxl") if "Patient_Runs" in xls.sheet_names else pd.DataFrame(columns=PATIENT_RUN_SHEET_COLUMNS)
                    old_control = pd.read_excel(excel_path, sheet_name="Control_Runs", engine="openpyxl") if "Control_Runs" in xls.sheet_names else pd.DataFrame(columns=CONTROL_RUN_SHEET_COLUMNS)
                    old_peaks = pd.read_excel(excel_path, sheet_name="PK_Peaks", engine="openpyxl") if "PK_Peaks" in xls.sheet_names else pd.DataFrame(columns=PEAK_SHEET_COLUMNS)
            except Exception:
                old_patient = pd.DataFrame(columns=PATIENT_RUN_SHEET_COLUMNS)
                old_control = pd.DataFrame(columns=CONTROL_RUN_SHEET_COLUMNS)
                old_peaks = pd.DataFrame(columns=PEAK_SHEET_COLUMNS)
        else:
            old_patient = pd.DataFrame(columns=PATIENT_RUN_SHEET_COLUMNS)
            old_control = pd.DataFrame(columns=CONTROL_RUN_SHEET_COLUMNS)
            old_peaks = pd.DataFrame(columns=PEAK_SHEET_COLUMNS)

        if not patient_df.empty:
            old_patient = old_patient[~old_patient["IdentityKey"].isin(patient_df["IdentityKey"])]
        if not control_df.empty:
            old_control = old_control[~old_control["IdentityKey"].isin(control_df["IdentityKey"])]
        if not peaks_df.empty and not old_peaks.empty:
            new_keys = peaks_df[["IdentityKey", "MarkerName"]].astype(str).agg("::".join, axis=1)
            old_keys = old_peaks[["IdentityKey", "MarkerName"]].astype(str).agg("::".join, axis=1)
            old_peaks = old_peaks[~old_keys.isin(set(new_keys.tolist()))]

        all_patient = _concat_frames(old_patient, patient_df).drop_duplicates(subset=["IdentityKey"], keep="last")
        all_control = _concat_frames(old_control, control_df).drop_duplicates(subset=["IdentityKey"], keep="last")
        all_peaks = _concat_frames(old_peaks, peaks_df).drop_duplicates(subset=["IdentityKey", "MarkerName"], keep="last")

        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            all_patient.to_excel(writer, sheet_name="Patient_Runs", index=False)
            all_control.to_excel(writer, sheet_name="Control_Runs", index=False)
            all_peaks.to_excel(writer, sheet_name="PK_Peaks", index=False)

        refresh_clonality_tracking_dashboard(
            excel_path,
            dashboard_title="Fraggler FLT3/NPM1 Tracking Dashboard",
        )

    print_green(f"FLT3/NPM1 QC tracker updated in {excel_path}")


def _reindex_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame(columns=columns)
    result = df.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = ""
    return result[columns]


def _concat_frames(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    old = old if old is not None else pd.DataFrame()
    new = new if new is not None else pd.DataFrame()
    if old.empty:
        return new.copy()
    if new.empty:
        return old.copy()
    all_columns = list(dict.fromkeys([*old.columns.tolist(), *new.columns.tolist()]))
    return pd.concat([old.reindex(columns=all_columns), new.reindex(columns=all_columns)], ignore_index=True)


def _month_bucket(run_date: str) -> str:
    value = str(run_date or "").strip()
    if len(value) >= 7 and value[4:5] == "-" and value[7:8] == "-":
        return value[:7].replace("-", "_")
    return ""
