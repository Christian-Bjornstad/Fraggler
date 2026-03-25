from __future__ import annotations

from collections import defaultdict
import os
import __main__
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from config import resolve_analysis_excel_output_path
from fraggler.fraggler import (
    FsaFile,
    calculate_best_combination_of_size_standard_peaks,
    find_size_standard_peaks,
    fit_size_standard_to_ladder,
    generate_combinations,
    print_green,
    print_warning,
    return_maxium_allowed_distance_between_size_standard_peaks,
)

from core.analysis import (
    MIN_R2_QUALITY,
    _select_best_ladder_candidate,
    analyse_fsa_rox,
    compute_ladder_qc_metrics,
    estimate_running_baseline,
)
from core.analyses.flt3.classification import classify_fsa
from core.analyses.flt3.config import ASSAY_CONFIG, BP_CORRECTION_OFFSETS, PREFERRED_INJECTION_TIME
from core.analyses.shared_pipeline import finalize_pipeline_run, normalize_pipeline_paths, scan_fsa_files
from core.html_reports import extract_dit_from_name


FLT3_LADDER_QC_THRESHOLD = 0.99
RELEVANT_PEAK_LABELS = {"WT", "MUT", "ITD"}
FLT3_QC_TRENDS_FILENAME = "FLT3_QC_TRENDS.xlsx"
FLT3_MANUAL_RATIO_VERSION = 2
MANUAL_RATIO_ASSAYS = {"FLT3-ITD", "FLT3-D835"}


def _scan_files(fsa_dir: Path, mode: str = "all") -> list[Path]:
    """Scan recursively for FLT3 .fsa files, excluding water/Vann files."""
    return scan_fsa_files(fsa_dir, mode=mode, recursive=True)


def _should_use_multiprocessing() -> bool:
    disabled = os.environ.get("FRAGGLER_DISABLE_MULTIPROCESSING", "").strip().lower()
    if disabled in {"1", "true", "yes", "on"}:
        return False
    if getattr(sys, "frozen", False):
        return False
    main_file = getattr(__main__, "__file__", "")
    if not main_file or str(main_file).startswith("<"):
        return False
    if not Path(main_file).exists():
        return False
    return True


def _preferred_injection_time(meta: dict) -> int:
    assay = meta.get("assay")
    analysis_type = meta.get("analysis_type")
    if analysis_type == "ratio_quant":
        return 1
    if assay == "FLT3-D835":
        return 3
    if analysis_type in PREFERRED_INJECTION_TIME:
        return int(PREFERRED_INJECTION_TIME[analysis_type])
    return int(PREFERRED_INJECTION_TIME.get(assay, meta.get("protocol_injection_time", meta.get("injection_time", 0))))


def _candidate_sort_key(item: tuple[Path, dict], preferred_injection: int) -> tuple[int, int, str, str]:
    path, meta = item
    injection = int(meta.get("injection_time", 0) or 0)
    return (
        0 if injection == preferred_injection else 1,
        abs(injection - preferred_injection),
        meta.get("source_run_dir", ""),
        path.name,
    )


def _calculate_auc(trace: np.ndarray, time_idx: np.ndarray) -> float:
    if time_idx.size == 0:
        return 0.0
    clipped = time_idx[(time_idx >= 0) & (time_idx < trace.size)]
    if clipped.size == 0:
        return 0.0
    return float(trace[clipped].sum())


def _peak_area_half_width_bp(assay: str, label: str, center_bp: float) -> float:
    if assay == "FLT3-D835":
        if label == "WT":
            return 1.2
        if label == "MUT":
            return 0.5
        if abs(center_bp - 150.0) <= 6.0:
            return 0.8
        return 0.8
    if assay == "FLT3-ITD":
        if label == "WT" or abs(center_bp - 330.0) <= 8.0:
            return 2.0
        if label in {"ITD", "MUT"} or center_bp >= 335.0:
            return 1.0
        return 2.0
    return 5.0


def _resolve_peak_area(assay: str, combined_area: float, channel_areas: dict[str, float]) -> float:
    if assay != "FLT3-ITD":
        return combined_area
    finite_channel_areas = [float(v) for v in channel_areas.values() if np.isfinite(v)]
    if not finite_channel_areas:
        return combined_area
    return max(finite_channel_areas)


def _peak_id_for_row(row: pd.Series | dict, ordinal: int | None = None) -> str:
    bp = float(row.get("basepairs", 0.0))
    height = float(row.get("peaks", 0.0))
    parts = [f"{int(round(bp * 10)):05d}", f"{int(round(height)):06d}"]
    if ordinal is not None:
        parts.append(f"{int(ordinal):03d}")
    return "pk_" + "_".join(parts)


def _ensure_peak_ids(peaks: pd.DataFrame) -> pd.DataFrame:
    if peaks.empty:
        if "peak_id" not in peaks.columns:
            peaks = peaks.copy()
            peaks["peak_id"] = pd.Series(dtype=str)
        return peaks
    ensured = peaks.copy()
    if "peak_id" not in ensured.columns:
        ensured["peak_id"] = [
            _peak_id_for_row(row, ordinal=index)
            for index, (_, row) in enumerate(ensured.iterrows(), start=1)
        ]
    return ensured


def _default_manual_ratio_selection() -> dict:
    return {
        "enabled": False,
        "version": FLT3_MANUAL_RATIO_VERSION,
        "wt": {"peak_id": None, "channel": None},
        "mutants": [],
    }


def _normalize_manual_peak_spec(spec: dict | None) -> dict:
    spec = spec if isinstance(spec, dict) else {}
    peak_id = spec.get("peak_id", spec.get("id", spec.get("peakId")))
    channel = spec.get("channel")
    if channel is not None:
        channel = str(channel).upper()
    return {
        "peak_id": peak_id,
        "channel": channel,
    }


def _normalize_manual_ratio_selection(raw: dict | None) -> dict:
    normalized = _default_manual_ratio_selection()
    if not isinstance(raw, dict):
        return normalized

    normalized["enabled"] = bool(raw.get("enabled", False))
    try:
        normalized["version"] = int(raw.get("version", FLT3_MANUAL_RATIO_VERSION))
    except (TypeError, ValueError):
        normalized["version"] = FLT3_MANUAL_RATIO_VERSION

    if isinstance(raw.get("wt"), dict):
        normalized["wt"] = _normalize_manual_peak_spec(raw.get("wt"))
    else:
        normalized["wt"] = _normalize_manual_peak_spec(
            {
                "peak_id": raw.get("wt_peak_id", raw.get("selected_wt_peak_id")),
                "channel": raw.get("wt_channel", raw.get("selected_wt_channel")),
            }
        )

    mutants = raw.get("mutants")
    normalized_mutants: list[dict] = []
    if isinstance(mutants, list):
        for item in mutants:
            if isinstance(item, dict):
                normalized_mutants.append(_normalize_manual_peak_spec(item))
    else:
        mutant_ids = raw.get("mutant_peak_ids", raw.get("selected_mutant_peak_ids", []))
        if mutant_ids is None:
            mutant_ids = []
        if not isinstance(mutant_ids, list):
            mutant_ids = [mutant_ids]
        mutant_channels = raw.get("mutant_channels", {})
        if not isinstance(mutant_channels, dict):
            mutant_channels = {}
        for peak_id in mutant_ids:
            normalized_mutants.append(
                _normalize_manual_peak_spec(
                    {
                        "peak_id": peak_id,
                        "channel": mutant_channels.get(peak_id),
                    }
                )
            )
    normalized["mutants"] = normalized_mutants
    return normalized


def _peak_area_for_channel(row: pd.Series, channel: str | None) -> float:
    if channel:
        value = row.get(f"area_{channel}", np.nan)
        if np.isfinite(value):
            return float(value)
    value = row.get("area", np.nan)
    return float(value) if np.isfinite(value) else float("nan")


def _peak_source_channel(row: pd.Series, fallback: str | None = None) -> str | None:
    channel = row.get("source_channel", fallback)
    if channel is None:
        return None
    channel = str(channel).upper()
    return channel if channel.startswith("DATA") else None


def _lookup_peak_row(peaks: pd.DataFrame, peak_id: str | None) -> pd.Series | None:
    if peaks.empty or not peak_id or "peak_id" not in peaks.columns:
        return None
    match = peaks[peaks["peak_id"].astype(str) == str(peak_id)]
    if match.empty:
        return None
    return match.iloc[0]


def _empty_manual_ratio_resolution(entry: dict, reason: str) -> dict:
    return {
        "ratio_mode": "manual_required",
        "manual_ratio_selection": _normalize_manual_ratio_selection(entry.get("manual_ratio_selection")),
        "manual_ratio_selection_valid": False,
        "manual_ratio_selection_reason": reason,
        "selected_wt_row": None,
        "selected_wt_rows": pd.DataFrame(),
        "selected_mut_rows": pd.DataFrame(),
        "selected_wt_peak_id": None,
        "selected_wt_peak_ids": [],
        "selected_mutant_peak_ids": [],
        "selected_wt_bp": np.nan,
        "selected_wt_bps": [],
        "selected_mutant_bps": [],
        "selected_wt_area": 0.0,
        "selected_wt_areas": [],
        "selected_mutant_area": 0.0,
        "selected_mutant_areas": [],
        "selected_wt_channel": None,
        "selected_wt_channels": [],
        "selected_mutant_channels": [],
        "ratio_numerator_area": 0.0,
        "ratio_denominator_area": 0.0,
        "ratio": 0.0,
        "mutant_fraction": 0.0,
    }


def _wt_candidates_for_assay(
    peaks: pd.DataFrame,
    assay: str,
    expected_wt_bp: float,
    *,
    channel: str | None = None,
) -> pd.DataFrame:
    if peaks.empty:
        return pd.DataFrame()

    channel_rows = peaks.copy()
    if channel and "source_channel" in channel_rows.columns:
        channel_rows = channel_rows[
            channel_rows["source_channel"].astype(str).str.upper() == str(channel).upper()
        ]
    if channel_rows.empty:
        return channel_rows

    if assay == "FLT3-D835":
        wt_min, wt_max = ASSAY_CONFIG.get("FLT3-D835", {}).get("wt_range", (expected_wt_bp - 4.0, expected_wt_bp + 4.0))
        wt_candidates = channel_rows[
            (channel_rows["basepairs"].astype(float) >= float(wt_min))
            & (channel_rows["basepairs"].astype(float) <= float(wt_max))
        ].copy()
    else:
        wt_candidates = channel_rows.assign(
            _wt_distance=(channel_rows["basepairs"].astype(float) - expected_wt_bp).abs()
        )
        wt_candidates = wt_candidates[wt_candidates["_wt_distance"] <= 8.0].copy()

    if wt_candidates.empty:
        return wt_candidates
    if "_wt_distance" not in wt_candidates.columns:
        wt_candidates["_wt_distance"] = (wt_candidates["basepairs"].astype(float) - expected_wt_bp).abs()
    return wt_candidates


def _peak_row_payload(row: pd.Series | None, *, channel: str | None = None) -> dict:
    if row is None:
        return {
            "peak_id": None,
            "basepairs": np.nan,
            "peaks": np.nan,
            "label": "",
            "area": 0.0,
            "channel": channel,
        }
    payload = {
        "peak_id": row.get("peak_id"),
        "basepairs": float(row.get("basepairs", np.nan)),
        "peaks": float(row.get("peaks", np.nan)),
        "label": row.get("label", ""),
        "area": float(row.get("area", 0.0)),
        "channel": channel,
    }
    if channel:
        payload["area"] = _peak_area_for_channel(row, channel)
    return payload


def _resolve_auto_ratio_selection(entry: dict, peaks: pd.DataFrame) -> dict:
    assay = entry.get("assay")
    wt_rows = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False) if not peaks.empty else pd.DataFrame()
    mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].copy() if not peaks.empty else pd.DataFrame()

    selected_mut_rows = mut_rows
    if assay == "FLT3-ITD":
        selected_mut_rows = _reportable_itd_mut_rows(entry, peaks, wt_rows=wt_rows, mut_rows=mut_rows)

    wt_main = wt_rows.iloc[0] if not wt_rows.empty else None
    if wt_main is None:
        return {
            "ratio_mode": "auto",
            "manual_ratio_selection": _normalize_manual_ratio_selection(entry.get("manual_ratio_selection")),
            "manual_ratio_selection_valid": False,
            "manual_ratio_selection_reason": "",
            "selected_wt_row": None,
            "selected_mut_rows": selected_mut_rows.iloc[0:0].copy(),
            "selected_wt_peak_id": None,
            "selected_mutant_peak_ids": [],
            "selected_wt_bp": np.nan,
            "selected_mutant_bps": [],
            "selected_wt_area": 0.0,
            "selected_mutant_area": 0.0,
            "selected_wt_channel": None,
            "selected_mutant_channels": [],
            "ratio_numerator_area": 0.0,
            "ratio_denominator_area": 0.0,
            "ratio": 0.0,
            "mutant_fraction": 0.0,
        }
    if assay == "FLT3-ITD" and not selected_mut_rows.empty:
        mut_area = float(selected_mut_rows.area.sum())
    elif assay in {"FLT3-D835", "NPM1"} and not mut_rows.empty:
        selected_mut_rows = mut_rows.sort_values("area", ascending=False).iloc[[0]]
        mut_area = float(selected_mut_rows.iloc[0].area)
    else:
        selected_mut_rows = selected_mut_rows.iloc[0:0].copy()
        mut_area = 0.0

    wt_area = float(wt_main.area) if wt_main is not None else 0.0
    ratio = (mut_area / wt_area) if wt_area > 0 else 0.0

    return {
        "ratio_mode": "auto",
        "manual_ratio_selection": _normalize_manual_ratio_selection(entry.get("manual_ratio_selection")),
        "manual_ratio_selection_valid": False,
        "manual_ratio_selection_reason": "",
        "selected_wt_row": wt_main,
        "selected_mut_rows": selected_mut_rows,
        "selected_wt_peak_id": wt_main.get("peak_id") if wt_main is not None else None,
        "selected_mutant_peak_ids": [row.peak_id for row in selected_mut_rows.itertuples(index=False)] if not selected_mut_rows.empty and "peak_id" in selected_mut_rows.columns else [],
        "selected_wt_bp": float(wt_main.basepairs) if wt_main is not None else np.nan,
        "selected_mutant_bps": [round(float(v), 2) for v in selected_mut_rows.basepairs.tolist()] if not selected_mut_rows.empty else [],
        "selected_wt_area": wt_area,
        "selected_mutant_area": mut_area,
        "selected_wt_channel": None,
        "selected_mutant_channels": [],
        "ratio_numerator_area": mut_area,
        "ratio_denominator_area": wt_area,
        "ratio": ratio,
        "mutant_fraction": (mut_area / (mut_area + wt_area)) if (mut_area + wt_area) > 0 else 0.0,
    }


def _resolve_manual_ratio_selection(entry: dict, peaks: pd.DataFrame) -> dict | None:
    assay = entry.get("assay")
    if assay not in MANUAL_RATIO_ASSAYS:
        return None
    if peaks.empty:
        return _empty_manual_ratio_resolution(entry, "Ingen manuelle peaks registrert")

    manual = _normalize_manual_ratio_selection(entry.get("manual_ratio_selection"))
    if not manual["enabled"]:
        return _empty_manual_ratio_resolution(entry, f"Manuelt peakvalg kreves for {assay}-ratio")

    wt_spec = manual.get("wt") or {}
    mut_specs = manual.get("mutants") or []
    if not mut_specs:
        return _empty_manual_ratio_resolution(entry, "Velg minst en mutantpeak manuelt")

    selected_mut_rows: list[pd.Series] = []
    selected_mut_ids: list[str] = []
    selected_mut_bps: list[float] = []
    selected_mut_areas: list[float] = []
    selected_mut_channels: list[str | None] = []
    mut_area = 0.0
    seen_pairs: set[tuple[str, str | None]] = set()

    for spec in mut_specs:
        peak_id = spec.get("peak_id")
        if not peak_id:
            continue
        mut_row = _lookup_peak_row(peaks, peak_id)
        if mut_row is None:
            continue
        mut_channel = spec.get("channel")
        if mut_channel is not None:
            mut_channel = str(mut_channel).upper()
        if mut_channel is None:
            mut_channel = entry.get("primary_peak_channel")
        mut_channel = _peak_source_channel(mut_row, fallback=mut_channel)
        selection_key = (str(peak_id), mut_channel)
        if selection_key in seen_pairs:
            continue
        mut_area_value = _peak_area_for_channel(mut_row, mut_channel)
        if not np.isfinite(mut_area_value) or mut_area_value <= 0:
            return _empty_manual_ratio_resolution(entry, "Valgt mutantpeak mangler brukbar kanal/area")
        seen_pairs.add(selection_key)
        selected_mut_rows.append(mut_row)
        selected_mut_ids.append(str(mut_row.get("peak_id")))
        selected_mut_bps.append(round(float(mut_row.get("basepairs", np.nan)), 2))
        selected_mut_channels.append(mut_channel)
        selected_mut_areas.append(float(mut_area_value))
        mut_area += float(mut_area_value)

    if not selected_mut_rows:
        return _empty_manual_ratio_resolution(entry, "Ingen gyldige manuelle mutantpeaks valgt")

    expected_wt_bp = float(ASSAY_CONFIG.get(assay, {}).get("wt_bp", entry.get("wt_bp", 330.0) or 330.0))
    selected_wt_rows: list[pd.Series] = []
    selected_wt_ids: list[str] = []
    selected_wt_bps: list[float] = []
    selected_wt_areas: list[float] = []
    selected_wt_channels: list[str] = []
    denominator_area = 0.0
    wt_row = _lookup_peak_row(peaks, wt_spec.get("peak_id"))
    wt_channel = wt_spec.get("channel")
    if wt_channel is not None:
        wt_channel = str(wt_channel).upper()
    if wt_channel is None:
        wt_channel = entry.get("primary_peak_channel")
    if wt_row is not None:
        wt_channel = _peak_source_channel(wt_row, fallback=wt_channel)
    if wt_row is not None or wt_channel:
        if wt_row is None:
            wt_candidates = _wt_candidates_for_assay(peaks, assay, expected_wt_bp, channel=wt_channel)
            if wt_candidates.empty:
                return _empty_manual_ratio_resolution(entry, "Mangler manuell WT-peak for valgt WT-kanal")
            wt_row = wt_candidates.sort_values(["_wt_distance", "peaks"], ascending=[True, False]).iloc[0]
        if wt_channel is None:
            wt_channel = _peak_source_channel(wt_row, fallback=entry.get("primary_peak_channel"))
        wt_area = _peak_area_for_channel(wt_row, wt_channel)
        if not np.isfinite(wt_area) or wt_area <= 0:
            return _empty_manual_ratio_resolution(entry, "Valgt WT-peak mangler brukbar area")
        if str(wt_row.get("peak_id")) in selected_mut_ids:
            return _empty_manual_ratio_resolution(entry, "WT-peaken kan ikke ogsa brukes som mutant")
        selected_wt_rows.append(wt_row)
        selected_wt_ids.append(str(wt_row.get("peak_id")))
        selected_wt_bps.append(round(float(wt_row.get("basepairs", np.nan)), 2))
        selected_wt_areas.append(float(wt_area))
        selected_wt_channels.append(wt_channel)
        denominator_area += float(wt_area)
    else:
        active_channels = [channel for channel in dict.fromkeys(selected_mut_channels) if channel]
        if assay == "FLT3-D835":
            active_channels = [entry.get("primary_peak_channel") or "DATA3"]
        for channel in active_channels:
            wt_candidates = _wt_candidates_for_assay(peaks, assay, expected_wt_bp, channel=channel)
            if wt_candidates.empty:
                return _empty_manual_ratio_resolution(entry, f"Mangler manuell WT-peak i {channel}")
            inferred_wt_row = wt_candidates.sort_values(["_wt_distance", "peaks"], ascending=[True, False]).iloc[0]
            wt_area = _peak_area_for_channel(inferred_wt_row, channel)
            if not np.isfinite(wt_area) or wt_area <= 0:
                return _empty_manual_ratio_resolution(entry, f"WT-peak i {channel} mangler brukbar area")
            if str(inferred_wt_row.get("peak_id")) in selected_mut_ids:
                return _empty_manual_ratio_resolution(entry, f"WT-peaken i {channel} er valgt som mutant")
            selected_wt_rows.append(inferred_wt_row)
            selected_wt_ids.append(str(inferred_wt_row.get("peak_id")))
            selected_wt_bps.append(round(float(inferred_wt_row.get("basepairs", np.nan)), 2))
            selected_wt_areas.append(float(wt_area))
            selected_wt_channels.append(channel)
            denominator_area += float(wt_area)

    if denominator_area <= 0:
        return _empty_manual_ratio_resolution(entry, "Ingen gyldig WT-area funnet for valgt mutantkanal")

    wt_row = selected_wt_rows[0] if selected_wt_rows else None
    wt_area = denominator_area

    return {
        "ratio_mode": "manual",
        "manual_ratio_selection": manual,
        "manual_ratio_selection_valid": True,
        "manual_ratio_selection_reason": "",
        "selected_wt_row": wt_row,
        "selected_wt_rows": pd.DataFrame(selected_wt_rows) if selected_wt_rows else pd.DataFrame(),
        "selected_mut_rows": pd.DataFrame(selected_mut_rows),
        "selected_wt_peak_id": wt_row.get("peak_id") if wt_row is not None else None,
        "selected_wt_peak_ids": selected_wt_ids,
        "selected_mutant_peak_ids": selected_mut_ids,
        "selected_wt_bp": float(wt_row.get("basepairs", np.nan)) if wt_row is not None else np.nan,
        "selected_wt_bps": selected_wt_bps,
        "selected_mutant_bps": selected_mut_bps,
        "selected_wt_area": wt_area,
        "selected_wt_areas": selected_wt_areas,
        "selected_mutant_area": mut_area,
        "selected_mutant_areas": selected_mut_areas,
        "selected_wt_channel": selected_wt_channels[0] if selected_wt_channels else None,
        "selected_wt_channels": selected_wt_channels,
        "selected_mutant_channels": selected_mut_channels,
        "ratio_numerator_area": mut_area,
        "ratio_denominator_area": wt_area,
        "ratio": (mut_area / wt_area) if wt_area > 0 else 0.0,
        "mutant_fraction": (mut_area / (mut_area + wt_area)) if (mut_area + wt_area) > 0 else 0.0,
    }


def _resolve_flt3_ratio_selection(entry: dict) -> dict:
    peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
    if entry.get("assay") in MANUAL_RATIO_ASSAYS:
        return _resolve_manual_ratio_selection(entry, peaks)
    return _resolve_auto_ratio_selection(entry, peaks)


def _calculate_peak_area(
    trace: np.ndarray,
    time_all: np.ndarray,
    bp_all: np.ndarray,
    center_bp: float,
    assay: str,
    label: str,
) -> float:
    half_width_bp = _peak_area_half_width_bp(assay, label, center_bp)
    auc_mask = (bp_all >= center_bp - half_width_bp) & (bp_all <= center_bp + half_width_bp)
    return _calculate_auc(trace, time_all[auc_mask])


def _correct_peak_channel_traces(
    fsa: FsaFile,
    channels: list[str],
    *,
    bin_size: int = 5000,
    quantile: float = 0.01,
) -> dict[str, np.ndarray]:
    """Baseline-correct each peak channel once and reuse the result."""
    corrected: dict[str, np.ndarray] = {}
    for ch in channels:
        if ch not in fsa.fsa:
            continue
        raw = np.asarray(fsa.fsa[ch]).astype(float)
        baseline = estimate_running_baseline(raw, bin_size=bin_size, quantile=quantile)
        corrected[ch] = np.maximum(raw - baseline, 0.0)
    return corrected


def _assay_positive_ratio(assay: str) -> float:
    return float(ASSAY_CONFIG.get(assay, {}).get("positive_ratio", 0.01))


def _bp_in_ranges(bp: float, ranges: list[tuple[float, float]] | None) -> bool:
    if not ranges:
        return False
    return any(start <= bp <= end for start, end in ranges)


def _apply_bp_offset(fsa: FsaFile, assay: str) -> None:
    offset = float(BP_CORRECTION_OFFSETS.get(assay, 0.0))
    sample_data = getattr(fsa, "sample_data_with_basepairs", None)
    if not offset or sample_data is None or sample_data.empty:
        return
    adjusted = sample_data.copy()
    adjusted["basepairs"] = adjusted["basepairs"].astype(float) + offset
    fsa.sample_data_with_basepairs = adjusted


def _infer_sizing_method(fsa: FsaFile) -> str:
    if hasattr(fsa, "_flt3_sizing_method"):
        return str(getattr(fsa, "_flt3_sizing_method"))
    model = getattr(fsa, "ladder_model", None)
    model_name = type(model).__name__
    if model_name == "Pipeline":
        return "spline"
    if model_name == "LinearRegression":
        return "polynomial_refinement"
    return "unknown"


def _detect_peaks(
    fsa: FsaFile,
    assay: str,
    wt_bp: float,
    trace: np.ndarray,
    mut_bp: float | None = None,
    analysis_type: str | None = None,
    corrected_channel_traces: dict[str, np.ndarray] | None = None,
) -> pd.DataFrame:
    """Detect WT and mutant peaks and estimate their corrected AUC."""
    sample_data = getattr(fsa, "sample_data_with_basepairs", None)
    if sample_data is None or sample_data.empty:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    time_all = sample_data["time"].astype(int).to_numpy()
    bp_all = sample_data["basepairs"].to_numpy()

    bp_min, bp_max = 50.0, 1000.0
    peaks: list[dict] = []
    assay_cfg = ASSAY_CONFIG.get(assay, {})

    from scipy.signal import find_peaks

    mask = (bp_all >= bp_min) & (bp_all <= bp_max)
    if not mask.any():
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    valid_time = time_all[mask]
    valid_time = valid_time[(valid_time >= 0) & (valid_time < trace.size)]
    if valid_time.size < 3:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    y_win = trace[time_all[mask]]
    bp_win = bp_all[mask]

    peak_height_min = float(assay_cfg.get("peak_height_min", 200))
    peak_distance = int(assay_cfg.get("peak_distance", 20))
    peak_idx, _ = find_peaks(y_win, height=peak_height_min, distance=peak_distance)

    wt_tol = 4.0 if assay == "FLT3-ITD" else 2.0 if assay in {"FLT3-D835", "NPM1"} else 5.0
    itd_min_bp = float(ASSAY_CONFIG.get("FLT3-ITD", {}).get("itd_min_bp", wt_bp + 4.9)) if wt_bp else None
    wt_range = assay_cfg.get("wt_range")
    mut_ranges = assay_cfg.get("mut_ranges")

    channel_traces = corrected_channel_traces or _correct_peak_channel_traces(
        fsa,
        assay_cfg.get("peak_channels", ["DATA1", "DATA2", "DATA3"]),
    )

    for idx in peak_idx:
        p_bp = float(bp_win[idx])
        p_h = float(y_win[idx])

        label = "unspecific"
        if _bp_in_ranges(p_bp, [wt_range] if wt_range else None) or (wt_bp and abs(p_bp - wt_bp) < wt_tol):
            label = "WT"
        elif assay == "FLT3-ITD" and itd_min_bp is not None and p_bp >= itd_min_bp:
            label = "ITD"
        elif _bp_in_ranges(p_bp, mut_ranges) or (mut_bp and abs(p_bp - mut_bp) < (wt_tol + 2)):
            label = "MUT"

        channel_areas = {
            ch: _calculate_peak_area(
                trace=corr_trace,
                time_all=time_all,
                bp_all=bp_all,
                center_bp=p_bp,
                assay=assay,
                label=label,
            )
            for ch, corr_trace in channel_traces.items()
        }
        combined_area = _calculate_peak_area(
            trace=trace,
            time_all=time_all,
            bp_all=bp_all,
            center_bp=p_bp,
            assay=assay,
            label=label,
        )
        p_area = _resolve_peak_area(assay=assay, combined_area=combined_area, channel_areas=channel_areas)

        peak_info = {
            "basepairs": p_bp,
            "peaks": p_h,
            "area": p_area,
            "label": label,
            "keep": True,
        }
        for ch, channel_area in channel_areas.items():
            peak_info[f"area_{ch}"] = channel_area
        peaks.append(peak_info)

    if not peaks:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label", "peak_id"])
    return _ensure_peak_ids(pd.DataFrame(peaks))


def _combine_peak_traces(
    fsa: FsaFile,
    peak_channels: list[str],
    primary_channel: str,
    corrected_channel_traces: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    combined_trace = None
    channel_traces = corrected_channel_traces or _correct_peak_channel_traces(fsa, peak_channels)
    for ch in peak_channels:
        corrected = channel_traces.get(ch)
        if corrected is None:
            continue
        combined_trace = corrected if combined_trace is None else combined_trace + corrected

    if combined_trace is not None:
        return combined_trace

    if primary_channel in fsa.fsa:
        return np.asarray(fsa.fsa[primary_channel]).astype(float)

    first_trace = next(iter(fsa.fsa.values()))
    return np.zeros(len(first_trace), dtype=float)


def _peak_qc_status(peaks: pd.DataFrame, group: str) -> tuple[bool, str]:
    if group == "negative_control":
        return True, "negative_control"
    if peaks.empty:
        return False, "no_peaks"
    relevant = peaks[peaks.label.isin(RELEVANT_PEAK_LABELS)]
    if relevant.empty:
        return False, "no_relevant_peaks"
    return True, "ok"


def _attempt_lenient_rox_fit(fsa_path: Path, sample_channel: str) -> FsaFile | None:
    configs = [{"min_h": 5, "min_d": 3}]
    for cfg in configs:
        try:
            fsa = FsaFile(
                file=str(fsa_path),
                ladder="GS500ROX",
                sample_channel=sample_channel,
                min_distance_between_peaks=cfg["min_d"],
                min_size_standard_height=cfg["min_h"],
                size_standard_channel="DATA4",
            )
            fsa = find_size_standard_peaks(fsa)
            ss_peaks = getattr(fsa, "size_standard_peaks", None)
            if ss_peaks is None or getattr(ss_peaks, "shape", [0])[0] < 3:
                continue

            fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=1.5)
            for _ in range(20):
                fsa = generate_combinations(fsa)
                best = getattr(fsa, "best_size_standard_combinations", None)
                if best is not None and best.shape[0] > 0:
                    break
                fsa.maxium_allowed_distance_between_size_standard_peaks += 10

            best = getattr(fsa, "best_size_standard_combinations", None)
            if best is None or best.shape[0] == 0:
                continue

            selected_fit = _select_best_ladder_candidate(fsa)
            if selected_fit is not None:
                fsa = selected_fit
            else:
                fsa = calculate_best_combination_of_size_standard_peaks(fsa)
                if not getattr(fsa, "fitted_to_model", False):
                    fsa = fit_size_standard_to_ladder(fsa)

            if not getattr(fsa, "fitted_to_model", False):
                continue

            qc = compute_ladder_qc_metrics(fsa)
            if qc["r2"] >= MIN_R2_QUALITY:
                setattr(fsa, "_flt3_sizing_method", "spline_lenient")
                return fsa

            setattr(fsa, "_flt3_sizing_method", "spline_lenient")
            return fsa
        except Exception:
            continue

    return None


def _analyse_fsa_candidate(fsa_path: Path, sample_channel: str, assay: str) -> FsaFile | None:
    fsa = analyse_fsa_rox(fsa_path, sample_channel)
    if fsa is not None:
        setattr(fsa, "_flt3_sizing_method", _infer_sizing_method(fsa))
        return fsa
    if assay == "FLT3-D835":
        return _attempt_lenient_rox_fit(fsa_path, sample_channel)
    return None


def _build_entry_from_candidate(fsa_path: Path, meta: dict) -> dict | None:
    fsa = _analyse_fsa_candidate(fsa_path, meta["primary_peak_channel"], meta["assay"])
    if fsa is None:
        return None

    _apply_bp_offset(fsa, meta["assay"])
    corrected_channel_traces = _correct_peak_channel_traces(
        fsa,
        meta.get("peak_channels", [meta["primary_peak_channel"]]),
    )
    combined_trace = _combine_peak_traces(
        fsa=fsa,
        peak_channels=meta.get("peak_channels", [meta["primary_peak_channel"]]),
        primary_channel=meta["primary_peak_channel"],
        corrected_channel_traces=corrected_channel_traces,
    )
    peaks = _detect_peaks(
        fsa=fsa,
        assay=meta["assay"],
        wt_bp=meta["wt_bp"],
        trace=combined_trace,
        mut_bp=meta.get("mut_bp"),
        analysis_type=meta.get("analysis_type"),
        corrected_channel_traces=corrected_channel_traces,
    )

    metrics = compute_ladder_qc_metrics(fsa)
    ladder_qc_status = "ok" if float(metrics.get("r2", float("nan"))) > FLT3_LADDER_QC_THRESHOLD else "ladder_qc_failed"
    peak_qc_pass, peak_qc_reason = _peak_qc_status(peaks, meta.get("group", "sample"))

    return {
        "fsa": fsa,
        "peaks_by_channel": {meta["primary_peak_channel"]: peaks},
        "trace_channels": meta["trace_channels"],
        "primary_peak_channel": meta["primary_peak_channel"],
        "ymax": float(np.max(combined_trace)) * 1.1 if combined_trace.size and np.any(combined_trace) else 1000.0,
        "assay": meta["assay"],
        "analysis_type": meta["analysis_type"],
        "parallel": meta.get("parallel"),
        "well_id": meta.get("well_id"),
        "specimen_id": meta.get("specimen_id"),
        "selection_key": meta.get("selection_key"),
        "group": meta["group"],
        "ladder": fsa.ladder,
        "bp_min": meta["bp_min"],
        "bp_max": meta["bp_max"],
        "dit": extract_dit_from_name(fsa.file_name),
        "ladder_qc_status": ladder_qc_status,
        "ladder_r2": float(metrics.get("r2", np.nan)),
        "n_ladder_steps": metrics.get("n_ladder_steps"),
        "n_size_standard_peaks": metrics.get("n_size_standard_peaks"),
        "injection_time": meta["injection_time"],
        "selected_injection": f"{int(meta['injection_time'])}s",
        "selected_injection_time": int(meta["injection_time"]),
        "preferred_injection_time": _preferred_injection_time(meta),
        "protocol_injection_time": meta.get("protocol_injection_time", meta["injection_time"]),
        "source_run_dir": meta.get("source_run_dir", ""),
        "run_name": meta.get("run_name", ""),
        "run_date": meta.get("run_date", ""),
        "run_time": meta.get("run_time", ""),
        "injection_protocol": meta.get("injection_protocol", ""),
        "selection_reason": "",
        "alternate_injections": [],
        "alternate_injections_summary": "",
        "sizing_method": _infer_sizing_method(fsa),
        "manual_ratio_selection": _default_manual_ratio_selection(),
        "ratio_mode": "auto",
        "manual_ratio_selection_valid": False,
        "manual_ratio_selection_reason": "",
        "selected_wt_peak_id": None,
        "selected_wt_peak_ids": [],
        "selected_mutant_peak_ids": [],
        "selected_wt_bp": np.nan,
        "selected_wt_bps": [],
        "selected_mutant_bps": [],
        "selected_wt_area": 0.0,
        "selected_wt_areas": [],
        "selected_mutant_area": 0.0,
        "selected_mutant_areas": [],
        "selected_wt_channel": None,
        "selected_wt_channels": [],
        "selected_mutant_channels": [],
        "peak_qc_pass": peak_qc_pass,
        "peak_qc_status": peak_qc_reason,
    }


def _candidate_audit_record(path: Path, meta: dict, status: str, reason: str) -> dict:
    return {
        "file": path.name,
        "injection_time": int(meta.get("injection_time", 0) or 0),
        "selected_injection": f"{int(meta.get('injection_time', 0) or 0)}s",
        "source_run_dir": meta.get("source_run_dir", ""),
        "status": status,
        "reason": reason,
    }


def _select_best_entry(candidates: list[tuple[Path, dict]]) -> dict | None:
    if not candidates:
        return None

    preferred_injection = _preferred_injection_time(candidates[0][1])
    ordered = sorted(candidates, key=lambda item: _candidate_sort_key(item, preferred_injection))

    audit_records: list[dict] = []
    best_available: tuple[dict, str] | None = None
    preferred_reason = "preferred injection unavailable"

    for index, (path, meta) in enumerate(ordered):
        same_as_preferred = int(meta.get("injection_time", 0) or 0) == preferred_injection
        entry = _build_entry_from_candidate(path, meta)
        if entry is None:
            reason = "ladder_fit_failed"
            if same_as_preferred:
                preferred_reason = reason
            audit_records.append(_candidate_audit_record(path, meta, "rejected", reason))
            continue

        acceptable = entry["ladder_qc_status"] == "ok" and entry["peak_qc_pass"]
        candidate_reason = "qc_pass"
        if not acceptable:
            if entry["ladder_qc_status"] != "ok":
                candidate_reason = entry["ladder_qc_status"]
            else:
                candidate_reason = entry["peak_qc_status"]
            if same_as_preferred:
                preferred_reason = candidate_reason
            if best_available is None:
                best_available = (entry, candidate_reason)
            audit_records.append(_candidate_audit_record(path, meta, "rejected", candidate_reason))
            continue

        if same_as_preferred:
            selection_reason = f"Preferred {preferred_injection}s injection selected"
        elif any(int(m.get("injection_time", 0) or 0) == preferred_injection for _, m in candidates):
            selection_reason = f"Preferred {preferred_injection}s failed ({preferred_reason}); selected {entry['selected_injection']} fallback"
        else:
            selection_reason = f"Preferred {preferred_injection}s unavailable; selected {entry['selected_injection']}"

        entry["selection_reason"] = selection_reason
        entry["preferred_injection_time"] = preferred_injection
        entry["alternate_injections"] = [
            record for record in audit_records
            if record["file"] != path.name
        ] + [
            _candidate_audit_record(other_path, other_meta, "not_selected", "preferred candidate passed")
            for other_path, other_meta in ordered[index + 1:]
        ]
        entry["alternate_injections_summary"] = "; ".join(
            f"{alt['selected_injection']} {alt['file']} ({alt['status']}: {alt['reason']})"
            for alt in entry["alternate_injections"]
        )
        return entry

    if best_available is None:
        return None

    entry, best_reason = best_available
    entry["selection_reason"] = f"No candidate passed QC; kept {entry['selected_injection']} ({best_reason})"
    entry["preferred_injection_time"] = preferred_injection
    entry["alternate_injections"] = audit_records
    entry["alternate_injections_summary"] = "; ".join(
        f"{alt['selected_injection']} {alt['file']} ({alt['status']}: {alt['reason']})"
        for alt in audit_records
        if alt["file"] != entry["fsa"].file_name
    )
    return entry


def _calculate_ratios(entries: list[dict]) -> None:
    """Calculate FLT3 mutant ratios and store explicit numerator/denominator fields."""
    for entry in entries:
        resolved = _resolve_flt3_ratio_selection(entry)
        entry["manual_ratio_selection"] = resolved.get("manual_ratio_selection", _default_manual_ratio_selection())
        entry["ratio_mode"] = resolved.get("ratio_mode", "auto")
        entry["manual_ratio_selection_valid"] = bool(resolved.get("manual_ratio_selection_valid", False))
        entry["manual_ratio_selection_reason"] = resolved.get("manual_ratio_selection_reason", "")
        entry["selected_wt_peak_id"] = resolved.get("selected_wt_peak_id")
        entry["selected_wt_peak_ids"] = resolved.get("selected_wt_peak_ids", [])
        entry["selected_mutant_peak_ids"] = resolved.get("selected_mutant_peak_ids", [])
        entry["selected_wt_bp"] = resolved.get("selected_wt_bp", np.nan)
        entry["selected_wt_bps"] = resolved.get("selected_wt_bps", [])
        entry["selected_mutant_bps"] = resolved.get("selected_mutant_bps", [])
        entry["selected_wt_area"] = float(resolved.get("selected_wt_area", 0.0))
        entry["selected_wt_areas"] = [float(v) for v in resolved.get("selected_wt_areas", [])]
        entry["selected_mutant_area"] = float(resolved.get("selected_mutant_area", 0.0))
        entry["selected_mutant_areas"] = [float(v) for v in resolved.get("selected_mutant_areas", [])]
        entry["selected_wt_channel"] = resolved.get("selected_wt_channel")
        entry["selected_wt_channels"] = resolved.get("selected_wt_channels", [])
        entry["selected_mutant_channels"] = resolved.get("selected_mutant_channels", [])
        entry["ratio_numerator_area"] = float(resolved.get("ratio_numerator_area", 0.0))
        entry["ratio_denominator_area"] = float(resolved.get("ratio_denominator_area", 0.0))
        entry["ratio"] = float(resolved.get("ratio", 0.0))
        entry["mutant_fraction"] = float(resolved.get("mutant_fraction", 0.0))


def _summarize_peak_areas(entry: dict) -> tuple[float, float]:
    return float(entry.get("ratio_denominator_area", 0.0)), float(entry.get("ratio_numerator_area", 0.0))


def _reportable_itd_mut_rows(
    entry: dict,
    peaks: pd.DataFrame,
    wt_rows: pd.DataFrame | None = None,
    mut_rows: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if entry.get("assay") != "FLT3-ITD" or peaks.empty:
        return mut_rows if mut_rows is not None else pd.DataFrame()

    wt_rows = wt_rows if wt_rows is not None else peaks[peaks.label == "WT"].sort_values("peaks", ascending=False)
    mut_rows = mut_rows if mut_rows is not None else peaks[peaks.label.isin(["MUT", "ITD"])].copy()
    if mut_rows.empty or wt_rows.empty:
        return mut_rows
    if entry.get("analysis_type") == "ratio_quant":
        return mut_rows

    wt_main = wt_rows.iloc[0]
    wt_bp = float(wt_main.basepairs)
    wt_area = float(wt_main.area)
    shoulder_bp_limit = wt_bp + 12.0
    shoulder_area_limit = max(4000.0, wt_area * 0.02)

    keep_mask = ~(
        (mut_rows.basepairs <= shoulder_bp_limit)
        & (mut_rows.area <= shoulder_area_limit)
    )
    return mut_rows[keep_mask].copy()


def _summarize_detected_peaks(entry: dict) -> dict:
    resolved = _resolve_flt3_ratio_selection(entry)
    wt_row = resolved.get("selected_wt_row")
    mut_rows = resolved.get("selected_mut_rows", pd.DataFrame())
    if not isinstance(mut_rows, pd.DataFrame):
        mut_rows = pd.DataFrame(mut_rows)
    mut_channels = resolved.get("selected_mutant_channels", [])
    wt_channel = resolved.get("selected_wt_channel")

    wt_bp = float(wt_row.basepairs) if wt_row is not None else np.nan
    wt_area = float(resolved.get("selected_wt_area", 0.0))
    mut_bps: list[float] = []
    mut_areas: list[float] = []
    for idx, (_, row) in enumerate(mut_rows.iterrows()):
        channel = mut_channels[idx] if idx < len(mut_channels) else None
        mut_bps.append(round(float(row.get("basepairs", np.nan)), 2))
        mut_areas.append(round(float(_peak_area_for_channel(row, channel)), 2))

    mut_main_bp = np.nan
    mut_main_area = 0.0
    if mut_areas:
        mut_main_idx = int(np.argmax(mut_areas))
        mut_main_bp = mut_bps[mut_main_idx]
        mut_main_area = float(mut_areas[mut_main_idx])

    return {
        "ratio_mode": resolved.get("ratio_mode", "auto"),
        "manual_ratio_selection_valid": bool(resolved.get("manual_ratio_selection_valid", False)),
        "manual_ratio_selection_reason": resolved.get("manual_ratio_selection_reason", ""),
        "selected_wt_peak_id": resolved.get("selected_wt_peak_id"),
        "selected_wt_peak_ids": resolved.get("selected_wt_peak_ids", []),
        "selected_mutant_peak_ids": resolved.get("selected_mutant_peak_ids", []),
        "selected_wt_channel": wt_channel,
        "selected_wt_channels": resolved.get("selected_wt_channels", []),
        "selected_mutant_channels": mut_channels,
        "wt_bp": wt_bp,
        "wt_area": wt_area,
        "wt_bps": resolved.get("selected_wt_bps", []),
        "wt_areas": [round(float(v), 2) for v in resolved.get("selected_wt_areas", [])],
        "mut_bps": mut_bps,
        "mut_areas": mut_areas,
        "mut_area_total": float(sum(mut_areas)),
        "mut_main_bp": mut_main_bp,
        "mut_main_area": mut_main_area,
    }


def _interpret_entry(entry: dict) -> str:
    assay = entry["assay"]
    ratio = float(entry.get("ratio", 0.0))
    peak_summary = _summarize_detected_peaks(entry)
    positive_ratio = _assay_positive_ratio(assay)

    if assay == "FLT3-ITD":
        if ratio >= positive_ratio:
            return "Positiv FLT3-ITD"
        if peak_summary["mut_bps"]:
            return "Negativ FLT3-ITD - lavniva dokumentert"
        return "Ingen FLT3-ITD pavist"
    if assay == "FLT3-D835":
        if ratio >= positive_ratio:
            return "Positiv FLT3-D835"
        if peak_summary["mut_bps"]:
            return "FLT3-D835 under positiv grense - dokumentert"
        return "Ingen FLT3-D835 pavist"
    if assay == "NPM1":
        if ratio >= positive_ratio:
            return "Positiv NPM1"
        if peak_summary["mut_bps"]:
            return "Mulig NPM1 - vurder manuelt"
        return "Ingen NPM1-mutasjon pavist"
    return "Ingen tolkning"


def generate_flt3_peak_report(entries: list[dict], outdir: Path) -> None:
    rows = []
    for entry in entries:
        peak_summary = _summarize_detected_peaks(entry)
        rows.append(
            {
                "DIT": entry.get("dit") or "",
                "File": entry["fsa"].file_name,
                "Assay": entry["assay"],
                "Group": entry.get("group") or "",
                "Parallel": entry.get("parallel") or "",
                "Well": entry.get("well_id") or "",
                "Treatment": entry.get("analysis_type") or "",
                "SelectedInjection": entry.get("selected_injection") or "",
                "PreferredInjection": f"{int(entry.get('preferred_injection_time', 0) or 0)}s" if entry.get("preferred_injection_time") else "",
                "SelectionReason": entry.get("selection_reason") or "",
                "SourceRunDir": entry.get("source_run_dir") or "",
                "AlternateInjections": entry.get("alternate_injections_summary") or "",
                "SizingMethod": entry.get("sizing_method") or "",
                "RatioMode": peak_summary.get("ratio_mode", "auto"),
                "ManualSelectionEnabled": bool(entry.get("manual_ratio_selection_valid", False)),
                "ManualSelectionReason": entry.get("manual_ratio_selection_reason") or "",
                "SelectedWT_PeakID": peak_summary.get("selected_wt_peak_id") or "",
                "SelectedWT_Channel": peak_summary.get("selected_wt_channel") or "",
                "SelectedMutant_PeakIDs": ", ".join(str(v) for v in peak_summary.get("selected_mutant_peak_ids", [])),
                "SelectedMutant_Channels": ", ".join(
                    str(v) for v in peak_summary.get("selected_mutant_channels", []) if v is not None
                ),
                "InjectionTime": entry.get("injection_time"),
                "WT_bp": round(float(peak_summary["wt_bp"]), 2) if not np.isnan(peak_summary["wt_bp"]) else "",
                "WT_Area": round(peak_summary["wt_area"], 2),
                "Mutant_bp": ", ".join(f"{bp:.2f}" for bp in peak_summary["mut_bps"]),
                "Mutant_Area": ", ".join(f"{area:.2f}" for area in peak_summary["mut_areas"]),
                "Mutant_Area_Total": round(peak_summary["mut_area_total"], 2),
                "RatioNumeratorArea": round(float(entry.get("ratio_numerator_area", 0.0)), 2),
                "RatioDenominatorArea": round(float(entry.get("ratio_denominator_area", 0.0)), 2),
                "Ratio": round(float(entry.get("ratio", 0.0)), 4),
                "MutantFractionMutPlusWT": round(float(entry.get("mutant_fraction", 0.0)), 4),
                "Interpretation": _interpret_entry(entry),
                "LadderQC": entry.get("ladder_qc_status", ""),
                "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
            }
        )

    if not rows:
        return

    df = pd.DataFrame(rows)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "Final_Detailed_Peak_Report.csv"
    df.to_csv(csv_path, index=False)
    print_green(f"FLT3 detailed peak report saved to {csv_path}")


def _build_control_qc_row(entry: dict) -> dict | None:
    group = entry.get("group")
    if group not in ["negative_control", "positive_control", "reactive_control"]:
        return None

    peak_summary = _summarize_detected_peaks(entry)
    if peak_summary.get("ratio_mode") == "manual_required":
        peaks = entry["peaks_by_channel"].get(entry["primary_peak_channel"], pd.DataFrame())
        if not peaks.empty:
            raw_wt_rows = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False)
            raw_mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].copy()
            if entry.get("assay") == "FLT3-ITD":
                raw_mut_rows = _reportable_itd_mut_rows(entry, peaks, wt_rows=raw_wt_rows, mut_rows=raw_mut_rows)

            if raw_wt_rows is not None and not raw_wt_rows.empty and (peak_summary.get("wt_bp") != peak_summary.get("wt_bp")):
                wt_main = raw_wt_rows.iloc[0]
                peak_summary["wt_bp"] = float(wt_main.basepairs)
                peak_summary["wt_area"] = float(wt_main.area)

            if raw_mut_rows is not None and not raw_mut_rows.empty and not peak_summary.get("mut_bps"):
                peak_summary["mut_bps"] = [round(float(v), 2) for v in raw_mut_rows.basepairs.tolist()]
                peak_summary["mut_areas"] = [round(float(v), 2) for v in raw_mut_rows.area.tolist()]
                peak_summary["mut_area_total"] = float(raw_mut_rows.area.sum())

    wt_area = float(peak_summary.get("wt_area", 0.0))
    mut_area = float(peak_summary.get("mut_area_total", 0.0))
    ratio = float(entry.get("ratio", 0.0))
    assay_cfg = ASSAY_CONFIG.get(entry["assay"], {})
    min_wt_area = float(assay_cfg.get("control_wt_min_area", 0.0))
    min_ratio = _assay_positive_ratio(entry["assay"])

    status = "FAIL"
    details = ""
    expectation = ""

    if group == "negative_control":
        expectation = "Ingen mutant/ITD-topper forventet"
        if not peak_summary["mut_bps"]:
            status = "PASS"
        else:
            details = f"Unexpected mutant peaks found: {peak_summary['mut_bps']}"
    elif group == "reactive_control":
        expectation = f"WT-topp forventet (min area {min_wt_area:.0f})"
        if peak_summary["wt_bp"] == peak_summary["wt_bp"] and wt_area >= min_wt_area:
            status = "PASS"
        elif peak_summary["wt_bp"] != peak_summary["wt_bp"]:
            details = "No WT peak detected"
        else:
            details = f"WT area below threshold ({wt_area:.0f} < {min_wt_area:.0f})"
    elif group == "positive_control":
        expectation = f"Mutantsignal forventet (ratio >= {min_ratio:.4f})"
        if peak_summary["mut_bps"] and (ratio >= min_ratio or peak_summary["wt_bp"] != peak_summary["wt_bp"]):
            status = "PASS"
        elif not peak_summary["mut_bps"]:
            details = "No mutant/ITD peak detected"
        else:
            details = f"Mutant ratio below threshold ({ratio:.4f} < {min_ratio:.4f})"

    return {
        "File": entry["fsa"].file_name,
        "ControlGroup": group,
        "Assay": entry["assay"],
        "Well": entry.get("well_id") or "",
        "SelectedInjection": entry.get("selected_injection") or "",
        "SelectionReason": entry.get("selection_reason") or "",
        "InjectionTime": entry.get("injection_time"),
        "WT_Area": round(wt_area, 2),
        "Mutant_Area": round(mut_area, 2),
        "Ratio": round(ratio, 4),
        "Status": status,
        "Details": details,
        "Expectation": expectation,
    }


def _flt3_control_entries(entries: list[dict]) -> list[dict]:
    return [
        entry for entry in entries
        if entry.get("group") in {"negative_control", "positive_control", "reactive_control"}
    ]


def _build_flt3_qc_trend_frames(entries: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    control_entries = _flt3_control_entries(entries)
    run_rows = []
    peak_rows = []

    for entry in control_entries:
        peak_summary = _summarize_detected_peaks(entry)
        wt_area, mut_area = _summarize_peak_areas(entry)
        peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
        interpretation = _interpret_entry(entry)

        run_rows.append(
            {
                "File": entry["fsa"].file_name,
                "ControlGroup": entry.get("group") or "",
                "Assay": entry.get("assay") or "",
                "Treatment": entry.get("analysis_type") or "",
                "DIT": entry.get("dit") or "",
                "SpecimenID": entry.get("specimen_id") or "",
                "Well": entry.get("well_id") or "",
                "RunDate": entry.get("run_date") or "",
                "RunTime": entry.get("run_time") or "",
                "RunName": entry.get("run_name") or "",
                "SourceRunDir": entry.get("source_run_dir") or "",
                "InjectionProtocol": entry.get("injection_protocol") or "",
                "InjectionTime": entry.get("injection_time"),
                "SelectedInjection": entry.get("selected_injection") or "",
                "PreferredInjection": (
                    f"{int(entry.get('preferred_injection_time', 0) or 0)}s"
                    if entry.get("preferred_injection_time") else ""
                ),
                "ProtocolInjectionTime": entry.get("protocol_injection_time"),
                "SelectionReason": entry.get("selection_reason") or "",
                "AlternateInjections": entry.get("alternate_injections_summary") or "",
                "SizingMethod": entry.get("sizing_method") or "",
                "Ladder": entry.get("ladder") or "",
                "LadderQC": entry.get("ladder_qc_status") or "",
                "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
                "PeakQC": entry.get("peak_qc_status") or "",
                "WT_bp": round(float(peak_summary["wt_bp"]), 2) if not np.isnan(peak_summary["wt_bp"]) else "",
                "WT_Area": round(wt_area, 2),
                "MutantMain_bp": round(float(peak_summary["mut_main_bp"]), 2) if not np.isnan(peak_summary["mut_main_bp"]) else "",
                "MutantMain_Area": round(float(peak_summary["mut_main_area"]), 2),
                "Mutant_bp_List": ", ".join(f"{bp:.2f}" for bp in peak_summary["mut_bps"]),
                "Mutant_Area_List": ", ".join(f"{area:.2f}" for area in peak_summary["mut_areas"]),
                "Mutant_Area_Total": round(mut_area, 2),
                "RatioNumeratorArea": round(float(entry.get("ratio_numerator_area", 0.0)), 2),
                "RatioDenominatorArea": round(float(entry.get("ratio_denominator_area", 0.0)), 2),
                "Ratio": round(float(entry.get("ratio", 0.0)), 4),
                "MutantFractionMutPlusWT": round(float(entry.get("mutant_fraction", 0.0)), 4),
                "RatioMode": peak_summary.get("ratio_mode", "auto"),
                "ManualSelectionEnabled": bool(peak_summary.get("manual_ratio_selection_valid", False)),
                "ManualSelectionReason": peak_summary.get("manual_ratio_selection_reason") or "",
                "SelectedWT_PeakID": peak_summary.get("selected_wt_peak_id") or "",
                "SelectedMutant_PeakIDs": ", ".join(str(v) for v in peak_summary.get("selected_mutant_peak_ids", [])),
                "Interpretation": interpretation,
            }
        )

        for idx, peak in enumerate(peaks.sort_values(["label", "basepairs", "peaks"], ascending=[True, True, False]).itertuples(index=False), start=1):
            peak_rows.append(
                {
                    "File": entry["fsa"].file_name,
                    "ControlGroup": entry.get("group") or "",
                    "Assay": entry.get("assay") or "",
                    "Well": entry.get("well_id") or "",
                    "RunDate": entry.get("run_date") or "",
                    "RunTime": entry.get("run_time") or "",
                    "SelectedInjection": entry.get("selected_injection") or "",
                    "LadderQC": entry.get("ladder_qc_status") or "",
                    "PeakRank": idx,
                    "PeakLabel": getattr(peak, "label", ""),
                    "PeakBP": round(float(getattr(peak, "basepairs", np.nan)), 2) if not np.isnan(getattr(peak, "basepairs", np.nan)) else "",
                    "PeakHeight": round(float(getattr(peak, "peaks", 0.0)), 2),
                    "PeakArea": round(float(getattr(peak, "area", 0.0)), 2),
                    "Keep": bool(getattr(peak, "keep", True)),
                    "PrimaryChannel": entry.get("primary_peak_channel") or "",
                }
            )

    return pd.DataFrame(run_rows), pd.DataFrame(peak_rows)


def update_flt3_qc_trends(excel_path: Path, entries: list[dict]) -> None:
    excel_path.parent.mkdir(parents=True, exist_ok=True)

    df_runs, df_peaks = _build_flt3_qc_trend_frames(entries)
    if df_runs.empty and df_peaks.empty:
        return

    if excel_path.exists():
        try:
            with pd.ExcelFile(excel_path, engine="openpyxl") as xls:
                has_runs = "Control_Runs" in xls.sheet_names
                has_peaks = "Control_Peaks" in xls.sheet_names
        except Exception:
            has_runs = False
            has_peaks = False

        old_runs = pd.read_excel(excel_path, sheet_name="Control_Runs", engine="openpyxl") if has_runs else pd.DataFrame()
        old_peaks = pd.read_excel(excel_path, sheet_name="Control_Peaks", engine="openpyxl") if has_peaks else pd.DataFrame()

        if not df_runs.empty and not old_runs.empty and "File" in old_runs.columns:
            old_runs = old_runs[~old_runs["File"].isin(df_runs["File"])]
        if not df_peaks.empty and not old_peaks.empty and "File" in old_peaks.columns:
            old_peaks = old_peaks[~old_peaks["File"].isin(df_peaks["File"])]

        all_runs = pd.concat([old_runs, df_runs], ignore_index=True)
        all_peaks = pd.concat([old_peaks, df_peaks], ignore_index=True)

        if not all_runs.empty and "File" in all_runs.columns:
            all_runs = all_runs.drop_duplicates(subset=["File"], keep="last")
        if not all_peaks.empty and {"File", "PeakRank"}.issubset(all_peaks.columns):
            all_peaks = all_peaks.drop_duplicates(subset=["File", "PeakRank"], keep="last")

        with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            all_runs.to_excel(writer, sheet_name="Control_Runs", index=False)
            all_peaks.to_excel(writer, sheet_name="Control_Peaks", index=False)
    else:
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            df_runs.to_excel(writer, sheet_name="Control_Runs", index=False)
            df_peaks.to_excel(writer, sheet_name="Control_Peaks", index=False)

    print_green(f"FLT3 QC trends updated in {excel_path}")


def generate_flt3_bp_validation_report(entries: list[dict], outdir: Path) -> None:
    rows = []
    for entry in entries:
        peak_summary = _summarize_detected_peaks(entry)
        assay_cfg = ASSAY_CONFIG.get(entry["assay"], {})
        if not np.isnan(peak_summary["wt_bp"]):
            expected_wt = float(assay_cfg.get("wt_bp", np.nan))
            rows.append(
                {
                    "DIT": entry.get("dit") or "",
                    "File": entry["fsa"].file_name,
                    "Assay": entry["assay"],
                    "Group": entry.get("group") or "",
                    "Well": entry.get("well_id") or "",
                    "InjectionTime": entry.get("injection_time"),
                    "SelectedInjection": entry.get("selected_injection") or "",
                    "SizingMethod": entry.get("sizing_method") or "",
                    "PeakType": "WT",
                    "ExpectedBP": round(expected_wt, 2) if np.isfinite(expected_wt) else "",
                    "ObservedBP": round(float(peak_summary["wt_bp"]), 2),
                    "DeltaBP": round(float(peak_summary["wt_bp"]) - expected_wt, 2) if np.isfinite(expected_wt) else "",
                    "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
                }
            )

        if entry["assay"] in {"FLT3-D835", "NPM1"} and not np.isnan(peak_summary["mut_main_bp"]):
            expected_mut = float(assay_cfg.get("mut_bp", np.nan))
            rows.append(
                {
                    "DIT": entry.get("dit") or "",
                    "File": entry["fsa"].file_name,
                    "Assay": entry["assay"],
                    "Group": entry.get("group") or "",
                    "Well": entry.get("well_id") or "",
                    "InjectionTime": entry.get("injection_time"),
                    "SelectedInjection": entry.get("selected_injection") or "",
                    "SizingMethod": entry.get("sizing_method") or "",
                    "PeakType": "MUT",
                    "ExpectedBP": round(expected_mut, 2) if np.isfinite(expected_mut) else "",
                    "ObservedBP": round(float(peak_summary["mut_main_bp"]), 2),
                    "DeltaBP": round(float(peak_summary["mut_main_bp"]) - expected_mut, 2) if np.isfinite(expected_mut) else "",
                    "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
                }
            )

    if not rows:
        return

    df = pd.DataFrame(rows)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "FLT3_BP_Validation.csv"
    df.to_csv(csv_path, index=False)
    print_green(f"FLT3 bp validation report saved to {csv_path}")


def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
    tracking_excel_path: Path | None = None,
    update_tracking_workbook: bool = True,
    progress_callback=None,
) -> list[dict] | None:
    """
    Kjor FLT3-pipeline pa alle .fsa-filer i fsa_dir.
    """
    fsa_dir, assay_dir = normalize_pipeline_paths(fsa_dir, base_outdir, assay_folder_name)

    raw_files = _scan_files(fsa_dir, mode=mode)

    if _should_use_multiprocessing() and len(raw_files) >= 2:
        from multiprocessing import Pool, cpu_count
        n_workers = max(1, cpu_count() - 1)
        try:
            with Pool(n_workers) as pool:
                meta_results = pool.map(classify_fsa, raw_files)
        except Exception:
            meta_results = [classify_fsa(p) for p in raw_files]
    else:
        meta_results = [classify_fsa(p) for p in raw_files]
    classified = [(p, m) for p, m in zip(raw_files, meta_results) if m is not None]

    if not classified:
        return [] if return_entries else None

    groups: dict[str, list[tuple[Path, dict]]] = defaultdict(list)
    for path, meta in classified:
        groups[meta["selection_key"]].append((path, meta))

    sorted_groups = sorted(groups.items())
    candidates_list = [c for _, c in sorted_groups]

    if _should_use_multiprocessing() and len(candidates_list) >= 2:
        from multiprocessing import Pool, cpu_count
        n_workers = max(1, cpu_count() - 1)
        try:
            with Pool(n_workers) as pool:
                results = pool.map(_select_best_entry, candidates_list)
        except Exception as ex:
            print_warning(f"[PARALLEL] Multiprocessing failed during FLT3 selection ({ex}), falling back to sequential.")
            results = [_select_best_entry(c) for c in candidates_list]
    else:
        results = [_select_best_entry(c) for c in candidates_list]

    entries = []
    for i, entry in enumerate(results):
        if entry is None:
            selection_key = sorted_groups[i][0]
            candidates = sorted_groups[i][1]
            first_file = candidates[0][0].name if candidates else selection_key
            print_warning(f"FLT3 selection failed for {first_file}")
            continue
        entries.append(entry)

    if not entries:
        return [] if return_entries else None

    _calculate_ratios(entries)
    generate_flt3_peak_report(entries, assay_dir)
    generate_flt3_bp_validation_report(entries, assay_dir)
    resolved_tracking_excel_path = tracking_excel_path or resolve_analysis_excel_output_path(
        "flt3",
        assay_dir,
        FLT3_QC_TRENDS_FILENAME,
    )
    if update_tracking_workbook:
        update_flt3_qc_trends(resolved_tracking_excel_path, entries)

    return finalize_pipeline_run(
        entries,
        assay_dir,
        return_entries=return_entries,
        make_dit_reports=make_dit_reports,
        mode=mode,
    )
