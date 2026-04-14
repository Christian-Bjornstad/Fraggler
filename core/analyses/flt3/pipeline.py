from __future__ import annotations

import copy
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
    LADDER_FIT_PROFILE_FLT3_GS500ROX,
    MIN_R2_QUALITY,
    _select_best_ladder_candidate,
    apply_manual_ladder_mapping,
    analyse_fsa_rox,
    compute_ladder_qc_metrics,
    estimate_running_baseline,
    get_ladder_candidates,
)
from core.analyses.flt3.classification import classify_fsa
from core.analyses.flt3.config import (
    ASSAY_CONFIG,
    BP_CORRECTION_OFFSETS,
    PREFERRED_INJECTION_TIME,
    ROX_LADDER as FLT3_ROX_LADDER,
)
from core.analyses.shared_pipeline import finalize_pipeline_run, normalize_pipeline_paths, scan_fsa_files
from core.html_reports import extract_dit_from_name
from core.area import compute_peak_area_gaussian


FLT3_LADDER_QC_THRESHOLD = 0.99
RELEVANT_PEAK_LABELS = {"WT", "MUT", "ITD"}
FLT3_QC_TRENDS_FILENAME = "FLT3_QC_TRENDS.xlsx"
FLT3_MANUAL_RATIO_VERSION = 2
MANUAL_RATIO_ASSAYS = {"FLT3-ITD", "FLT3-D835"}
FLT3_REVIEW_MAX_RESIDUAL_BP = 3.0
FLT3_LOW_SIGNAL_LADDER_PROMINENCE = 100.0
FLT3_LOW_END_LADDER_PROMINENCE = 120.0
FLT3_LATE_TEMPLATE_TOLERANCE = 950.0

FLT3_TEMPLATE_TIMES: dict[tuple[str, str], tuple[float, ...]] = {
    ("FLT3-ITD", "10x_diluted"): (
        1575.0,
        1647.0,
        1798.5,
        1942.0,
        2174.5,
        2233.5,
        2292.5,
        2537.5,
        2834.5,
        3166.0,
        3411.5,
        3475.5,
        3797.5,
        4093.5,
        4335.0,
        4385.0,
    ),
    ("FLT3-ITD", "25x_diluted"): (
        1618.0,
        1686.0,
        1842.5,
        1986.5,
        2221.5,
        2281.5,
        2340.0,
        2588.0,
        2883.5,
        3220.5,
        3465.5,
        3532.0,
        3859.5,
        4162.0,
        4410.0,
        4461.0,
    ),
    ("FLT3-ITD", "25x_diluted_compact"): (
        1510.0,
        1578.0,
        1719.0,
        1854.0,
        2070.0,
        2125.0,
        2180.0,
        2408.0,
        2686.0,
        2990.0,
        3215.0,
        3273.0,
        3568.0,
        3841.0,
        4064.0,
        4110.0,
    ),
    ("FLT3-D835", "standard"): (
        1578.0,
        1652.5,
        1799.5,
        1944.0,
        2174.5,
        2233.5,
        2293.0,
        2538.0,
        2839.0,
        3163.5,
        3408.0,
        3469.0,
        3786.5,
        4083.5,
        4325.0,
        4376.5,
    ),
    ("FLT3-D835", "standard_compact"): (
        1489.0,
        1560.0,
        1695.0,
        1827.0,
        2037.0,
        2091.0,
        2146.0,
        2366.0,
        2637.0,
        2928.0,
        3147.0,
        3201.0,
        3486.0,
        3752.0,
        3969.0,
        4015.0,
    ),
    ("FLT3-ITD", "standard"): (
        1613.5,
        1685.5,
        1837.0,
        1980.5,
        2213.0,
        2272.0,
        2331.0,
        2576.0,
        2873.0,
        3204.5,
        3450.0,
        3514.0,
        3836.0,
        4132.0,
        4373.5,
        4423.5,
    ),
}

FLT3_TEMPLATE_STEPS: tuple[float, ...] = (
    35.0,
    50.0,
    75.0,
    100.0,
    139.0,
    150.0,
    160.0,
    200.0,
    250.0,
    300.0,
    340.0,
    350.0,
    400.0,
    450.0,
    490.0,
    500.0,
)
FLT3_TEMPLATE_STEP_INDEX = {int(round(bp)): idx for idx, bp in enumerate(FLT3_TEMPLATE_STEPS)}
FLT3_SHAPE_GAP_RULES: tuple[tuple[int, int, float], ...] = (
    (490, 500, 4.2),
    (450, 490, 2.0),
    (400, 450, 1.6),
    (340, 350, 2.8),
    (300, 350, 1.2),
    (139, 150, 2.6),
    (150, 160, 2.6),
    (100, 139, 0.9),
)
FLT3_SHAPE_TIME_WEIGHTS: tuple[tuple[int, float], ...] = (
    (35, 0.12),
    (50, 0.08),
    (340, 0.06),
    (350, 0.06),
    (500, 0.30),
    (490, 0.18),
    (450, 0.10),
    (400, 0.05),
)
FLT3_SHAPE_INTENSITY_STEPS: tuple[int, ...] = (139, 150, 160, 340, 350, 400, 450, 490, 500)
FLT3_INVALID_SHAPE_PENALTY = 250.0
FLT3_ACCEPTABLE_RESCUE_SHAPE_PENALTY = 120.0
FLT3_MIN_LADDER_PEAK_INTENSITY = 50.0
FLT3_EMPIRICAL_GAP_PROFILE: dict[tuple[int, int], dict[str, float]] = {
    (35, 50): {"median": 73.0, "p10": 69.0, "p90": 87.0, "weight": 0.5},
    (50, 75): {"median": 143.5, "p10": 136.0, "p90": 165.0, "weight": 0.7},
    (75, 100): {"median": 139.0, "p10": 134.0, "p90": 158.0, "weight": 0.7},
    (100, 139): {"median": 226.0, "p10": 217.0, "p90": 256.0, "weight": 1.0},
    (139, 150): {"median": 57.0, "p10": 55.0, "p90": 65.0, "weight": 2.8},
    (150, 160): {"median": 58.0, "p10": 56.0, "p90": 66.0, "weight": 2.8},
    (160, 200): {"median": 239.0, "p10": 229.0, "p90": 274.0, "weight": 1.0},
    (200, 250): {"median": 291.5, "p10": 281.0, "p90": 338.0, "weight": 0.9},
    (250, 300): {"median": 323.0, "p10": 305.0, "p90": 370.0, "weight": 0.8},
    (300, 340): {"median": 239.5, "p10": 228.0, "p90": 277.0, "weight": 0.9},
    (340, 350): {"median": 61.0, "p10": 56.0, "p90": 72.0, "weight": 2.4},
    (350, 400): {"median": 311.5, "p10": 296.0, "p90": 362.0, "weight": 0.8},
    (400, 450): {"median": 283.5, "p10": 272.0, "p90": 332.0, "weight": 1.0},
    (450, 490): {"median": 230.0, "p10": 221.0, "p90": 270.0, "weight": 2.2},
    (490, 500): {"median": 47.0, "p10": 45.0, "p90": 58.0, "weight": 3.4},
}
FLT3_EMPIRICAL_INTENSITY_PROFILE: dict[int, dict[str, float]] = {
    35: {"target": 0.80, "low": 0.45, "high": 1.20, "weight": 0.3},
    50: {"target": 0.78, "low": 0.45, "high": 1.20, "weight": 0.3},
    75: {"target": 0.90, "low": 0.55, "high": 1.30, "weight": 0.3},
    100: {"target": 1.03, "low": 0.65, "high": 1.40, "weight": 0.4},
    139: {"target": 1.08, "low": 0.70, "high": 1.45, "weight": 0.8},
    150: {"target": 1.10, "low": 0.72, "high": 1.45, "weight": 0.9},
    160: {"target": 1.12, "low": 0.74, "high": 1.50, "weight": 0.9},
    200: {"target": 1.11, "low": 0.74, "high": 1.50, "weight": 0.7},
    250: {"target": 1.11, "low": 0.72, "high": 1.45, "weight": 0.6},
    300: {"target": 1.07, "low": 0.70, "high": 1.40, "weight": 0.5},
    340: {"target": 0.96, "low": 0.62, "high": 1.30, "weight": 0.6},
    350: {"target": 0.99, "low": 0.64, "high": 1.30, "weight": 0.7},
    400: {"target": 0.91, "low": 0.58, "high": 1.20, "weight": 0.7},
    450: {"target": 0.84, "low": 0.50, "high": 1.10, "weight": 1.0},
    490: {"target": 0.81, "low": 0.45, "high": 1.05, "weight": 1.3},
    500: {"target": 0.78, "low": 0.40, "high": 1.00, "weight": 1.7},
}
FLT3_EMPIRICAL_WIDTH_PROFILE: dict[int, dict[str, float]] = {
    139: {"target": 4.0, "low": 2.5, "high": 7.0, "weight": 0.2},
    150: {"target": 4.0, "low": 2.5, "high": 7.0, "weight": 0.2},
    160: {"target": 4.0, "low": 2.5, "high": 7.0, "weight": 0.2},
    340: {"target": 5.0, "low": 3.0, "high": 8.0, "weight": 0.3},
    350: {"target": 5.0, "low": 3.0, "high": 8.0, "weight": 0.3},
    450: {"target": 5.5, "low": 3.5, "high": 9.5, "weight": 0.5},
    490: {"target": 6.0, "low": 3.5, "high": 10.5, "weight": 0.6},
    500: {"target": 6.0, "low": 3.5, "high": 10.5, "weight": 0.7},
}


def _flt3_expected_ladder_steps(fsa: FsaFile) -> np.ndarray:
    expected_steps = np.asarray(
        getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])),
        dtype=float,
    )
    if expected_steps.size == 0:
        return np.asarray(FLT3_TEMPLATE_STEPS, dtype=float)
    return expected_steps


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
    return compute_peak_area_gaussian(trace, time_all, bp_all, center_bp, half_width_bp)


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


def _flt3_template_key(assay: str, analysis_type: str | None) -> tuple[str, str] | None:
    exact = (str(assay or ""), str(analysis_type or ""))
    if exact in FLT3_TEMPLATE_TIMES:
        return exact
    if assay == "FLT3-ITD" and str(analysis_type or "") == "standard":
        return ("FLT3-ITD", "standard")
    if assay == "FLT3-ITD":
        return ("FLT3-ITD", "10x_diluted")
    assay_only = (str(assay or ""), "standard")
    if assay_only in FLT3_TEMPLATE_TIMES:
        return assay_only
    return None


def _candidate_flt3_template_keys(assay: str, analysis_type: str | None) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []

    def add(key: tuple[str, str] | None) -> None:
        if key is None or key in candidates or key not in FLT3_TEMPLATE_TIMES:
            return
        candidates.append(key)

    exact = (str(assay or ""), str(analysis_type or ""))
    add(exact)
    add(_flt3_template_key(assay, analysis_type))

    if assay == "FLT3-ITD":
        for key in (
            ("FLT3-ITD", "10x_diluted"),
            ("FLT3-ITD", "25x_diluted"),
            ("FLT3-ITD", "25x_diluted_compact"),
            ("FLT3-ITD", "standard"),
        ):
            add(key)
    elif assay == "FLT3-D835":
        for key in (
            ("FLT3-D835", "standard"),
            ("FLT3-D835", "standard_compact"),
        ):
            add(key)

    return candidates


def _flt3_template_key_allowed_for_fsa(fsa: FsaFile, template_key: tuple[str, str]) -> bool:
    if "compact" not in str(template_key[1]):
        return True

    source_text = " ".join(
        str(getattr(fsa, attr, "") or "")
        for attr in ("file", "file_name", "path", "source_file")
    )
    if not source_text:
        return True

    compact_tokens = ("310725", "H9C0ZIZJ", "rerun_all")
    if any(token in source_text for token in compact_tokens):
        return True

    noncompact_tokens = ("130326", "180326", "H9H1DI0C", "H9H1DHZL")
    if any(token in source_text for token in noncompact_tokens):
        return False

    return True


def _rank_flt3_template_keys_for_fsa(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> list[tuple[str, str]]:
    template_keys = [
        key
        for key in _candidate_flt3_template_keys(assay, analysis_type)
        if _flt3_template_key_allowed_for_fsa(fsa, key)
    ]
    if len(template_keys) <= 1:
        return template_keys

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return template_keys

    candidate_df = _flt3_candidate_pool(fsa).sort_values("time").reset_index(drop=True)
    raw_trace = np.asarray(getattr(getattr(fsa, "fsa", {}), "get", lambda *_: [])("DATA4"), dtype=float)
    if raw_trace.shape != trace.shape:
        raw_trace = trace
    anchor_candidates = _merged_anchor_candidates(trace, candidate_df, reference_trace=raw_trace)
    if anchor_candidates.empty:
        return template_keys

    latest_peak = float(anchor_candidates["time"].astype(float).max())
    ranked: list[tuple[int, float, int, tuple[str, str]]] = []
    for index, template_key in enumerate(template_keys):
        template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
        combos = _rank_flt3_high_end_anchor_combos(
            anchor_candidates,
            float(template_times[-3]),
            float(template_times[-2]),
            float(template_times[-1]),
            limit=3,
        )
        if combos:
            anchor_450, anchor_490, anchor_500 = combos[0]
            score = (
                abs((anchor_500 - anchor_490) - float(template_times[-1] - template_times[-2])) * 3.5
                + abs((anchor_490 - anchor_450) - float(template_times[-2] - template_times[-3])) * 1.7
                + abs(anchor_500 - float(template_times[-1])) * 0.18
                + abs(anchor_490 - float(template_times[-2])) * 0.10
                + abs(anchor_450 - float(template_times[-3])) * 0.06
            )
            ranked.append((0, float(score), index, template_key))
            continue

        fallback_score = abs(latest_peak - float(template_times[-1])) + 600.0
        ranked.append((1, float(fallback_score), index, template_key))

    return [template_key for _, _, _, template_key in sorted(ranked)]


def _resolved_flt3_template_key(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> tuple[str, str] | None:
    ranked = _rank_flt3_template_keys_for_fsa(fsa, assay, analysis_type)
    if ranked:
        return ranked[0]
    return _flt3_template_key(assay, analysis_type)


def _mark_flt3_short_trace_if_needed(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> bool:
    missing_steps, trace_last_index = _flt3_short_trace_missing_steps(fsa, assay, analysis_type)
    if not missing_steps:
        return False

    setattr(fsa, "ladder_fit_strategy", "short_trace")
    setattr(
        fsa,
        "ladder_fit_note",
        (
            f"ROX DATA4 trace ends at scan {trace_last_index:.0f}, before the expected "
            "GS500ROX high-end ladder steps: "
            + ", ".join(f"{step:.0f}" for step in missing_steps)
            + ". Full ladder assignment is not reliable."
        ),
    )
    setattr(fsa, "ladder_missing_expected_steps", missing_steps)
    setattr(fsa, "ladder_review_required", True)
    return True


def _flt3_short_trace_missing_steps(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> tuple[list[float], float]:
    if getattr(fsa, "ladder_fit_strategy", "") == "manual_adjustment":
        return [], float("nan")

    template_key = _resolved_flt3_template_key(fsa, assay, analysis_type)
    if template_key is None:
        return [], float("nan")

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return [], float("nan")

    template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
    expected_steps = _flt3_expected_ladder_steps(fsa)
    if template_times.size != expected_steps.size:
        return [], float("nan")

    trace_last_index = float(trace.size - 1)
    short_margin = 50.0
    missing_steps = [
        float(step_bp)
        for step_bp, template_time in zip(expected_steps, template_times, strict=False)
        if float(template_time) > trace_last_index - short_margin
    ]
    if not missing_steps:
        return [], trace_last_index

    if min(missing_steps) < 400.0:
        return [], trace_last_index

    return missing_steps, trace_last_index


def _mapped_peak_time_for_step(fsa: FsaFile, target_bp: float) -> float | None:
    ladder_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    peak_times = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    if ladder_steps.size == 0 or peak_times.size == 0 or ladder_steps.size != peak_times.size:
        return None
    matches = np.where(np.isclose(ladder_steps, float(target_bp), atol=1e-6))[0]
    if matches.size == 0:
        return None
    return float(peak_times[int(matches[0])])


def _snap_trace_peak(
    trace: np.ndarray,
    center: float,
    *,
    search_radius: int,
    lower_bound: float | None = None,
    upper_bound: float | None = None,
) -> float | None:
    if trace.size == 0:
        return None
    lo = max(0, int(round(float(center))) - int(search_radius))
    hi = min(trace.size - 1, int(round(float(center))) + int(search_radius))
    if lower_bound is not None:
        lo = max(lo, int(np.ceil(float(lower_bound))))
    if upper_bound is not None:
        hi = min(hi, int(np.floor(float(upper_bound))))
    if hi <= lo:
        return None
    window = np.asarray(trace[lo : hi + 1], dtype=float)
    if window.size == 0 or not np.any(np.isfinite(window)):
        return None
    local_index = int(np.nanargmax(window))
    return float(lo + local_index)


def _choose_template_trace_peak(
    trace: np.ndarray,
    target_time: float,
    previous_time: float | None,
    target_gap: float | None,
    intensity_reference: float | None,
    *,
    search_radius: int,
    lower_bound: float | None = None,
    upper_bound: float | None = None,
    min_intensity: float = FLT3_MIN_LADDER_PEAK_INTENSITY,
    min_prominence: float = 0.0,
) -> float | None:
    if trace.size == 0:
        return None
    lo = max(0, int(round(float(target_time))) - int(search_radius))
    hi = min(trace.size - 1, int(round(float(target_time))) + int(search_radius))
    if lower_bound is not None:
        lo = max(lo, int(np.ceil(float(lower_bound))))
    if upper_bound is not None:
        hi = min(hi, int(np.floor(float(upper_bound))))
    if hi <= lo:
        return None

    from scipy.signal import find_peaks

    window = np.asarray(trace[lo : hi + 1], dtype=float)
    if window.size == 0 or not np.any(np.isfinite(window)):
        return None
    peak_kwargs = {"distance": 6}
    if float(min_prominence) > 0.0:
        peak_kwargs["prominence"] = float(min_prominence)
    peak_idx, _ = find_peaks(window, **peak_kwargs)
    if peak_idx.size == 0:
        if float(min_prominence) > 0.0:
            return None
        snapped = _snap_trace_peak(
            trace,
            target_time,
            search_radius=search_radius,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )
        intensity = _trace_intensity_at_time(trace, snapped)
        if intensity is None or float(intensity) < float(min_intensity):
            return None
        return snapped

    peak_times = pd.Series((peak_idx + lo).astype(float))
    peak_heights = pd.Series(trace[(peak_idx + lo).astype(int)].astype(float), index=peak_times.index)
    strong = peak_heights.astype(float) >= float(min_intensity)
    if strong.any():
        peak_times = peak_times[strong]
        peak_heights = peak_heights[strong]
    else:
        return None
    if previous_time is not None and target_gap is not None and target_gap > 0:
        actual_gap = float(previous_time) - peak_times.astype(float)
        keep = (
            (actual_gap >= max(8.0, float(target_gap) * 0.35))
            & (actual_gap <= float(target_gap) * 1.85)
        )
        peak_times = peak_times[keep]
        peak_heights = peak_heights[keep]
        if peak_times.empty:
            return None

    score = _score_flt3_template_peak_choice(
        peak_times,
        peak_heights,
        target_time,
        previous_time,
        target_gap,
        intensity_reference,
    )
    best_idx = score.idxmin()
    return float(peak_times.loc[best_idx])


def _late_trace_peak_candidates(
    trace: np.ndarray,
    *,
    min_time: float = 3400.0,
    min_intensity: float = FLT3_MIN_LADDER_PEAK_INTENSITY,
    reference_trace: np.ndarray | None = None,
) -> pd.DataFrame:
    if trace.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    from scipy.signal import find_peaks

    reference = np.asarray(reference_trace, dtype=float) if reference_trace is not None else trace
    if reference.shape != trace.shape:
        reference = trace

    peak_idx, props = find_peaks(
        reference,
        distance=15,
        prominence=FLT3_LOW_SIGNAL_LADDER_PROMINENCE,
    )
    if peak_idx.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    rows = pd.DataFrame(
        {
            "time": peak_idx.astype(float),
            "intensity": trace[peak_idx].astype(float),
            "prominence": np.asarray(props.get("prominences", []), dtype=float),
            "source": "trace",
        }
    )
    low_signal = rows["intensity"].astype(float) < float(min_intensity)
    rows.loc[low_signal, "source"] = "late_prominent"
    rows = rows[
        (rows["time"].astype(float) >= float(min_time))
        & (
            (rows["intensity"].astype(float) >= float(min_intensity))
            | (rows["prominence"].astype(float) >= FLT3_LOW_SIGNAL_LADDER_PROMINENCE)
        )
    ].copy()
    if rows.empty:
        return pd.DataFrame(columns=["time", "intensity", "source"])
    return rows.sort_values(["time", "intensity"], ascending=[True, False]).reset_index(drop=True)


def _low_end_trace_peak_candidates(
    trace: np.ndarray,
    *,
    min_time: float,
    max_time: float,
) -> pd.DataFrame:
    if trace.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    from scipy.signal import find_peaks

    peak_idx, props = find_peaks(
        trace,
        distance=5,
        prominence=FLT3_LOW_END_LADDER_PROMINENCE,
    )
    if peak_idx.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    rows = pd.DataFrame(
        {
            "time": peak_idx.astype(float),
            "intensity": trace[peak_idx].astype(float),
            "prominence": np.asarray(props.get("prominences", []), dtype=float),
            "source": "low_prominent",
        }
    )
    rows = rows[
        (rows["time"].astype(float) >= float(min_time))
        & (rows["time"].astype(float) <= float(max_time))
        & (rows["prominence"].astype(float) >= FLT3_LOW_END_LADDER_PROMINENCE)
    ].copy()
    if rows.empty:
        return pd.DataFrame(columns=["time", "intensity", "source"])
    return rows.sort_values(["time", "intensity"], ascending=[True, False]).reset_index(drop=True)


def _bootstrap_trace_peak_candidates(
    trace: np.ndarray,
    *,
    min_time: float = 1400.0,
    max_time: float = 4700.0,
) -> pd.DataFrame:
    if trace.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    from scipy.signal import find_peaks

    peak_idx, props = find_peaks(
        trace,
        distance=10,
        prominence=30.0,
        height=35.0,
    )
    if peak_idx.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    rows = pd.DataFrame(
        {
            "time": peak_idx.astype(float),
            "intensity": trace[peak_idx].astype(float),
            "prominence": np.asarray(props.get("prominences", []), dtype=float),
            "source": "trace_bootstrap",
        }
    )
    rows = rows[
        (rows["time"].astype(float) >= float(min_time))
        & (rows["time"].astype(float) <= float(max_time))
    ].copy()
    if rows.empty:
        return pd.DataFrame(columns=["time", "intensity", "source"])
    rows = rows.sort_values(["prominence", "intensity"], ascending=[False, False]).head(80)
    return rows.sort_values(["time", "intensity"], ascending=[True, False]).reset_index(drop=True)


def _flt3_candidate_pool(
    fsa: FsaFile,
    *,
    reference_trace: np.ndarray | None = None,
) -> pd.DataFrame:
    candidate_df = get_ladder_candidates(fsa)
    if "source" not in candidate_df.columns:
        candidate_df["source"] = "auto"

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    bootstrap_df = _bootstrap_trace_peak_candidates(trace)
    frames: list[pd.DataFrame] = []
    if not candidate_df.empty:
        frames.append(candidate_df.loc[:, ["time", "intensity", "source"]].copy())
    if not bootstrap_df.empty:
        frames.append(bootstrap_df)
    if not frames:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    merged = pd.concat(frames, ignore_index=True)
    if merged.empty:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    merged["time_key"] = merged["time"].astype(float).round().astype(int)
    bootstrap_keys = set(
        merged.loc[merged["source"].astype(str) == "trace_bootstrap", "time_key"].astype(int).tolist()
    )
    merged = (
        merged.sort_values(["intensity", "time"], ascending=[False, False])
        .drop_duplicates(subset=["time_key"], keep="first")
        .sort_values("time")
        .reset_index(drop=True)
    )
    if bootstrap_keys:
        bootstrap_mask = merged["time_key"].astype(int).isin(bootstrap_keys) & (
            merged["source"].astype(str) != "auto"
        )
        merged.loc[bootstrap_mask, "source"] = "trace_bootstrap"
    return merged.loc[:, ["time", "intensity", "source"]]


def _merged_anchor_candidates(
    trace: np.ndarray,
    candidate_df: pd.DataFrame,
    *,
    reference_trace: np.ndarray | None = None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    if not candidate_df.empty:
        frames.append(candidate_df.loc[:, ["time", "intensity", "source"]].copy())
    late_trace = _late_trace_peak_candidates(trace, reference_trace=reference_trace)
    if not late_trace.empty:
        frames.append(late_trace)
    if not frames:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    merged = pd.concat(frames, ignore_index=True)
    if merged.empty:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    merged["time_key"] = merged["time"].astype(float).round().astype(int)
    late_prominent_keys = set(
        merged.loc[merged["source"].astype(str) == "late_prominent", "time_key"].astype(int).tolist()
    )
    merged = (
        merged.sort_values(["intensity", "time"], ascending=[False, False])
        .drop_duplicates(subset=["time_key"], keep="first")
        .sort_values("time")
        .reset_index(drop=True)
    )
    if late_prominent_keys:
        late_mask = merged["time_key"].astype(int).isin(late_prominent_keys)
        merged.loc[late_mask, "source"] = "late_prominent"
    return merged.loc[:, ["time", "intensity", "source"]]


def _candidate_index_for_time(candidate_times: list[float], peak_time: float, tolerance: float = 1.0) -> int | None:
    for idx, candidate_time in enumerate(candidate_times):
        if abs(float(candidate_time) - float(peak_time)) <= float(tolerance):
            return idx
    return None


def _mapping_times_from_fsa(
    fsa: FsaFile,
    expected_steps: np.ndarray,
) -> dict[int, float]:
    ladder_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    peak_times = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    if ladder_steps.size == 0 or peak_times.size == 0 or ladder_steps.size != peak_times.size:
        return {}

    mapping_times: dict[int, float] = {}
    for index, step_bp in enumerate(expected_steps):
        matches = np.where(np.isclose(ladder_steps, float(step_bp), atol=1e-6))[0]
        if matches.size == 0:
            continue
        mapping_times[int(index)] = float(peak_times[int(matches[0])])
    return mapping_times




def _trace_intensity_at_time(trace: np.ndarray, peak_time: float | None) -> float | None:
    if peak_time is None or trace.size == 0:
        return None
    idx = int(round(float(peak_time)))
    if idx < 0 or idx >= trace.size:
        return None
    return float(trace[idx])


def _trace_peak_width_points(trace: np.ndarray, peak_time: float | None) -> float | None:
    if peak_time is None or trace.size == 0:
        return None
    idx = int(round(float(peak_time)))
    if idx < 0 or idx >= trace.size:
        return None

    lo = max(0, idx - 30)
    hi = min(trace.size - 1, idx + 30)
    peak_height = float(trace[idx])
    baseline = float(np.min(trace[lo : hi + 1]))
    half_level = baseline + (peak_height - baseline) * 0.5

    left = idx
    while left > lo and float(trace[left]) > half_level:
        left -= 1
    right = idx
    while right < hi and float(trace[right]) > half_level:
        right += 1
    if right <= left:
        return None
    return float(right - left)


def _flt3_ladder_intensity_reference(
    trace: np.ndarray,
    mapping_times: dict[int, float],
) -> float | None:
    if trace.size == 0 or not mapping_times:
        return None

    intensities: list[float] = []
    for step_index, peak_time in mapping_times.items():
        step_bp = FLT3_TEMPLATE_STEPS[int(step_index)] if int(step_index) < len(FLT3_TEMPLATE_STEPS) else 0.0
        if step_bp < 139.0:
            continue
        intensity = _trace_intensity_at_time(trace, peak_time)
        if intensity is None or not np.isfinite(float(intensity)) or float(intensity) <= 0.0:
            continue
        intensities.append(float(intensity))

    if not intensities:
        return None
    return float(np.median(np.asarray(intensities, dtype=float)))


def _flt3_ladder_intensity_penalty(
    intensity: float,
    intensity_reference: float | None,
) -> float:
    if not np.isfinite(float(intensity)):
        return 50.0

    intensity = float(intensity)
    if intensity_reference is not None and np.isfinite(float(intensity_reference)) and float(intensity_reference) > 0.0:
        ratio = max(float(intensity) / float(intensity_reference), 0.05)
        penalty = abs(float(np.log2(ratio))) * 4.0
        if ratio < 0.40:
            penalty += (0.40 - ratio) * 45.0
        elif ratio > 2.40:
            penalty += (ratio - 2.40) * 60.0
        return float(penalty)

    if intensity < 300.0:
        return float((300.0 - intensity) * 0.08)
    if intensity > 2000.0:
        return float((intensity - 2000.0) * 0.02)
    return 0.0


def _score_flt3_template_peak_choice(
    peak_times: pd.Series | np.ndarray,
    intensities: pd.Series | np.ndarray,
    target_time: float,
    previous_time: float | None,
    target_gap: float | None,
    intensity_reference: float | None,
) -> pd.Series:
    times = pd.Series(peak_times, dtype=float)
    heights = pd.Series(intensities, index=times.index, dtype=float)
    score = (times - float(target_time)).abs()

    if previous_time is not None and target_gap is not None and target_gap > 0:
        gap_penalty = (float(previous_time) - times - float(target_gap)).abs()
        gap_weight = 1.5
        if float(target_gap) < 90.0:
            gap_weight = 2.0
        elif float(target_gap) < 150.0:
            gap_weight = 1.25
        score = score + gap_penalty * gap_weight

    intensity_penalty = heights.apply(
        lambda value: _flt3_ladder_intensity_penalty(float(value), intensity_reference)
    )
    score = score + intensity_penalty - heights.clip(lower=0.0, upper=2000.0) * 0.0015
    return score.astype(float)


def _fit_flt3_template_affine_alignment(
    mapping_times: dict[int, float],
    template_times: np.ndarray,
) -> tuple[float, float]:
    if not mapping_times:
        return 0.0, 1.0

    x_vals: list[float] = []
    y_vals: list[float] = []
    for step_index, actual_time in mapping_times.items():
        idx = int(step_index)
        if idx < 0 or idx >= int(template_times.size):
            continue
        x_vals.append(float(template_times[idx]))
        y_vals.append(float(actual_time))

    if len(x_vals) < 2:
        return 0.0, 1.0

    x = np.asarray(x_vals, dtype=float)
    y = np.asarray(y_vals, dtype=float)
    try:
        slope, intercept = np.polyfit(x, y, 1)
    except Exception:
        return 0.0, 1.0
    if not np.isfinite(float(slope)) or not np.isfinite(float(intercept)):
        return 0.0, 1.0
    slope = float(np.clip(float(slope), 0.45, 2.5))
    intercept = float(intercept)
    return intercept, slope


def _flt3_mapping_shape_penalty(
    mapping_times: dict[int, float],
    template_times: np.ndarray,
    trace: np.ndarray | None = None,
) -> float:
    if not mapping_times:
        return float("inf")

    penalty = 0.0
    intercept, slope = _fit_flt3_template_affine_alignment(mapping_times, template_times)
    aligned_template_times = intercept + slope * np.asarray(template_times, dtype=float)

    for step_bp, weight in FLT3_SHAPE_TIME_WEIGHTS:
        step_index = FLT3_TEMPLATE_STEP_INDEX[int(step_bp)]
        actual = mapping_times.get(step_index)
        if actual is None:
            penalty += 80.0 * float(weight)
            continue
        penalty += abs(float(actual) - float(aligned_template_times[step_index])) * float(weight)

    for left_bp, right_bp, weight in FLT3_SHAPE_GAP_RULES:
        left_index = FLT3_TEMPLATE_STEP_INDEX[int(left_bp)]
        right_index = FLT3_TEMPLATE_STEP_INDEX[int(right_bp)]
        left_time = mapping_times.get(left_index)
        right_time = mapping_times.get(right_index)
        if left_time is None or right_time is None:
            penalty += 120.0 * float(weight)
            continue

        target_gap = float(aligned_template_times[right_index] - aligned_template_times[left_index])
        actual_gap = float(right_time - left_time)
        penalty += abs(actual_gap - target_gap) * float(weight)

        min_allowed = max(12.0, target_gap * 0.45)
        max_allowed = max(min_allowed + 12.0, target_gap * 1.65)
        if actual_gap < min_allowed:
            penalty += (float(min_allowed) - actual_gap) * float(weight) * 2.5
        elif actual_gap > max_allowed:
            penalty += (actual_gap - float(max_allowed)) * float(weight) * 2.5

    for pair, profile in FLT3_EMPIRICAL_GAP_PROFILE.items():
        left_bp, right_bp = pair
        left_index = FLT3_TEMPLATE_STEP_INDEX[int(left_bp)]
        right_index = FLT3_TEMPLATE_STEP_INDEX[int(right_bp)]
        left_time = mapping_times.get(left_index)
        right_time = mapping_times.get(right_index)
        if left_time is None or right_time is None:
            continue
        actual_gap = float(right_time - left_time)
        median_gap = float(profile["median"]) * float(slope)
        low_gap = float(profile["p10"]) * float(slope)
        high_gap = float(profile["p90"]) * float(slope)
        weight = float(profile["weight"])
        penalty += abs(actual_gap - median_gap) * weight * 0.25
        if actual_gap < low_gap:
            penalty += (low_gap - actual_gap) * weight * 1.8
        elif actual_gap > high_gap:
            penalty += (actual_gap - high_gap) * weight * 1.8

    gap_139_150 = None
    gap_150_160 = None
    gap_340_350 = None
    gap_450_490 = None
    gap_490_500 = None
    if FLT3_TEMPLATE_STEP_INDEX[139] in mapping_times and FLT3_TEMPLATE_STEP_INDEX[150] in mapping_times:
        gap_139_150 = float(mapping_times[FLT3_TEMPLATE_STEP_INDEX[150]] - mapping_times[FLT3_TEMPLATE_STEP_INDEX[139]])
    if FLT3_TEMPLATE_STEP_INDEX[150] in mapping_times and FLT3_TEMPLATE_STEP_INDEX[160] in mapping_times:
        gap_150_160 = float(mapping_times[FLT3_TEMPLATE_STEP_INDEX[160]] - mapping_times[FLT3_TEMPLATE_STEP_INDEX[150]])
    if FLT3_TEMPLATE_STEP_INDEX[340] in mapping_times and FLT3_TEMPLATE_STEP_INDEX[350] in mapping_times:
        gap_340_350 = float(mapping_times[FLT3_TEMPLATE_STEP_INDEX[350]] - mapping_times[FLT3_TEMPLATE_STEP_INDEX[340]])
    if FLT3_TEMPLATE_STEP_INDEX[450] in mapping_times and FLT3_TEMPLATE_STEP_INDEX[490] in mapping_times:
        gap_450_490 = float(mapping_times[FLT3_TEMPLATE_STEP_INDEX[490]] - mapping_times[FLT3_TEMPLATE_STEP_INDEX[450]])
    if FLT3_TEMPLATE_STEP_INDEX[490] in mapping_times and FLT3_TEMPLATE_STEP_INDEX[500] in mapping_times:
        gap_490_500 = float(mapping_times[FLT3_TEMPLATE_STEP_INDEX[500]] - mapping_times[FLT3_TEMPLATE_STEP_INDEX[490]])

    if gap_139_150 is not None and gap_150_160 is not None:
        family_spread = abs(gap_139_150 - gap_150_160)
        penalty += family_spread * 0.9
        if family_spread > 12.0:
            penalty += (family_spread - 12.0) * 2.2
    if gap_340_350 is not None:
        low_340_350 = 50.0 * float(slope)
        high_340_350 = 78.0 * float(slope)
        if gap_340_350 < low_340_350:
            penalty += (low_340_350 - gap_340_350) * 2.4
        elif gap_340_350 > high_340_350:
            penalty += (gap_340_350 - high_340_350) * 2.4
    if gap_450_490 is not None and gap_490_500 is not None:
        ratio = gap_490_500 / max(gap_450_490, 1.0)
        if ratio > 0.34:
            penalty += (ratio - 0.34) * 170.0
        if gap_490_500 >= gap_450_490:
            penalty += (gap_490_500 - gap_450_490) * 2.5

    if trace is not None and trace.size:
        intensity_reference = _flt3_ladder_intensity_reference(trace, mapping_times)
        for step_bp in FLT3_SHAPE_INTENSITY_STEPS:
            step_index = FLT3_TEMPLATE_STEP_INDEX[int(step_bp)]
            peak_time = mapping_times.get(step_index)
            intensity = _trace_intensity_at_time(trace, peak_time)
            if intensity is None:
                continue
            if intensity < 300.0:
                penalty += (300.0 - float(intensity)) * 0.05
            elif intensity > 2000.0:
                penalty += (float(intensity) - 2000.0) * 0.005

            profile = FLT3_EMPIRICAL_INTENSITY_PROFILE.get(int(step_bp))
            if (
                profile is not None
                and intensity_reference is not None
                and np.isfinite(float(intensity_reference))
                and float(intensity_reference) > 0.0
            ):
                ratio = float(intensity) / float(intensity_reference)
                target = float(profile["target"])
                low = float(profile["low"])
                high = float(profile["high"])
                weight = float(profile["weight"])
                penalty += abs(ratio - target) * weight * 8.0
                if ratio < low:
                    penalty += (low - ratio) * weight * 20.0
                elif ratio > high:
                    penalty += (ratio - high) * weight * 16.0

            width_profile = FLT3_EMPIRICAL_WIDTH_PROFILE.get(int(step_bp))
            if width_profile is not None:
                peak_width = _trace_peak_width_points(trace, peak_time)
                if peak_width is not None and np.isfinite(float(peak_width)):
                    target = float(width_profile["target"])
                    low = float(width_profile["low"])
                    high = float(width_profile["high"])
                    weight = float(width_profile["weight"])
                    penalty += abs(float(peak_width) - target) * weight
                    if float(peak_width) < low:
                        penalty += (low - float(peak_width)) * weight * 3.0
                    elif float(peak_width) > high:
                        penalty += (float(peak_width) - high) * weight * 2.2

    return float(penalty)


def _flt3_fit_is_geometrically_invalid(
    mapping_times: dict[int, float],
    expected_steps: np.ndarray,
    template_times: np.ndarray,
    trace: np.ndarray,
) -> bool:
    if len(mapping_times) < int(expected_steps.size):
        return True
    penalty = _flt3_mapping_shape_penalty(mapping_times, template_times, trace)
    return not np.isfinite(penalty) or float(penalty) > FLT3_INVALID_SHAPE_PENALTY


def _accept_lenient_raw_flt3_fit(
    qc: dict[str, float | int],
    assay: str,
    analysis_type: str | None,
) -> bool:
    max_bp = float(qc.get("max_abs_error_bp", float("inf")))
    mean_bp = float(qc.get("mean_abs_error_bp", float("inf")))
    r2 = float(qc.get("r2", float("-inf")))
    normalized_analysis_type = str(analysis_type or "").strip().lower()

    if not np.isfinite(max_bp) or not np.isfinite(mean_bp) or not np.isfinite(r2):
        return False

    # Keep the new recovery path, but reject obviously bad "rescues" that only
    # convert fit_failed into huge-residual auto_full fits.
    if normalized_analysis_type in {"ratio_quant", "10x_diluted", "25x_diluted", "undiluted"}:
        return bool(max_bp <= 8.0 and mean_bp <= 2.8 and r2 >= 0.9990)
    if assay == "FLT3-ITD":
        return bool(max_bp <= 12.0 and mean_bp <= 4.0 and r2 >= 0.9985)
    return bool(max_bp <= 10.0 and mean_bp <= 3.5 and r2 >= 0.9988)


def _flt3_high_end_anchors_are_plausible(
    anchor_490: float | None,
    anchor_500: float | None,
    template_490: float,
    template_500: float,
) -> bool:
    if anchor_490 is None or anchor_500 is None:
        return False
    gap = float(anchor_500) - float(anchor_490)
    if gap < 35.0 or gap > 80.0:
        return False
    if abs(float(anchor_500) - float(template_500)) > FLT3_LATE_TEMPLATE_TOLERANCE:
        return False
    if abs(float(anchor_490) - float(template_490)) > FLT3_LATE_TEMPLATE_TOLERANCE:
        return False
    return True


def _flt3_peak_meets_min_intensity(trace: np.ndarray, peak_time: float | None) -> bool:
    intensity = _trace_intensity_at_time(trace, peak_time)
    return intensity is not None and float(intensity) >= FLT3_MIN_LADDER_PEAK_INTENSITY


def _flt3_template_rescue_trace_min_intensity(step_bp: float) -> float:
    return -220.0


def _flt3_template_rescue_trace_min_prominence(step_bp: float) -> float:
    if step_bp <= 75.0:
        return 18.0
    if step_bp <= 160.0:
        return 20.0
    if step_bp >= 450.0:
        return 25.0
    if step_bp >= 340.0:
        return 20.0
    return 18.0


def _flt3_template_window(template_time: float, step_bp: float) -> tuple[float, float]:
    if step_bp >= 490.0:
        return float(template_time - 110.0), float(template_time + 1200.0)
    if step_bp >= 450.0:
        return float(template_time - 120.0), float(template_time + 320.0)
    if step_bp >= 400.0:
        return float(template_time - 110.0), float(template_time + 100.0)
    return float(template_time - 120.0), float(template_time + 120.0)


def _select_flt3_high_end_anchor_combo(
    anchor_candidates: pd.DataFrame,
    template_450: float,
    template_490: float,
    template_500: float,
) -> tuple[float, float, float] | None:
    ranked = _rank_flt3_high_end_anchor_combos(
        anchor_candidates,
        template_450,
        template_490,
        template_500,
    )
    return ranked[0] if ranked else None


def _rank_flt3_high_end_anchor_combos(
    anchor_candidates: pd.DataFrame,
    template_450: float,
    template_490: float,
    template_500: float,
    *,
    limit: int = 8,
) -> list[tuple[float, float, float]]:
    if anchor_candidates.empty:
        return []

    rows = anchor_candidates.copy()
    rows["time"] = rows["time"].astype(float)
    rows["intensity"] = rows["intensity"].astype(float)
    if "source" not in rows.columns:
        rows["source"] = ""
    rows = rows[
        (rows["intensity"] >= FLT3_MIN_LADDER_PEAK_INTENSITY)
        | (rows["source"].astype(str) == "late_prominent")
    ].copy()
    if rows.empty:
        return []

    target_gap_490_500 = float(template_500 - template_490)
    target_gap_450_490 = float(template_490 - template_450)
    target_gap_450_500 = float(template_500 - template_450)

    # Compact FLT3 runs can place the true 450/490/500 anchors below the old
    # absolute floors (for example 3722/3933/3977), so keep the window tied to
    # the selected template rather than to long-run scan positions.
    candidates_500 = rows[
        (rows["time"] >= template_500 - 300.0)
        & (rows["time"] <= template_500 + 1200.0)
    ].copy()
    if candidates_500.empty:
        return []

    scored_combos: list[tuple[float, tuple[float, float, float]]] = []

    for peak_500 in candidates_500.itertuples(index=False):
        candidates_490 = rows[
            (rows["time"] < peak_500.time - 20.0)
            & (rows["time"] >= max(template_490 - 190.0, peak_500.time - 95.0))
            & (rows["time"] <= peak_500.time - 30.0)
        ].copy()
        if candidates_490.empty:
            continue

        for peak_490 in candidates_490.itertuples(index=False):
            candidates_450 = rows[
                (rows["time"] < peak_490.time - 120.0)
                & (rows["time"] >= max(template_450 - 230.0, peak_490.time - 340.0))
                & (rows["time"] <= peak_490.time - 170.0)
            ].copy()
            if candidates_450.empty:
                continue

            for peak_450 in candidates_450.itertuples(index=False):
                gap_490_500 = float(peak_500.time - peak_490.time)
                gap_450_490 = float(peak_490.time - peak_450.time)
                gap_450_500 = float(peak_500.time - peak_450.time)
                score = (
                    abs(gap_490_500 - target_gap_490_500) * 3.5
                    + abs(gap_450_490 - target_gap_450_490) * 1.7
                    + abs(gap_450_500 - target_gap_450_500) * 0.8
                    + abs(float(peak_500.time) - float(template_500)) * 0.18
                    + abs(float(peak_490.time) - float(template_490)) * 0.10
                    + abs(float(peak_450.time) - float(template_450)) * 0.06
                    - (
                        float(peak_500.intensity)
                        + float(peak_490.intensity)
                        + float(peak_450.intensity)
                    )
                    * 0.02
                )
                scored_combos.append(
                    (
                        float(score),
                        (
                            float(peak_450.time),
                            float(peak_490.time),
                            float(peak_500.time),
                        ),
                    )
                )

    if not scored_combos:
        return []

    ranked: list[tuple[float, float, float]] = []
    seen: set[tuple[int, int, int]] = set()
    for _, combo in sorted(scored_combos, key=lambda item: item[0]):
        combo_key = tuple(int(round(value)) for value in combo)
        if combo_key in seen:
            continue
        seen.add(combo_key)
        ranked.append(combo)
        if len(ranked) >= int(limit):
            break
    return ranked


def _choose_template_candidate_time(
    candidate_df: pd.DataFrame,
    target_time: float,
    previous_time: float | None,
    target_gap: float | None,
    intensity_reference: float | None = None,
) -> float | None:
    if candidate_df.empty:
        return None

    candidates = candidate_df.copy()
    candidates["time"] = candidates["time"].astype(float)
    candidates["intensity"] = candidates["intensity"].astype(float)

    lower_bound = float(target_time - 65.0)
    upper_bound = float(target_time + 55.0)
    if previous_time is not None:
        upper_bound = min(upper_bound, float(previous_time) - 8.0)

    candidates = candidates[
        (candidates["time"] >= lower_bound)
        & (candidates["time"] <= upper_bound)
    ].copy()
    if candidates.empty:
        return None
    candidates = candidates[
        candidates["intensity"].astype(float) >= FLT3_MIN_LADDER_PEAK_INTENSITY
    ].copy()
    if candidates.empty:
        return None

    if previous_time is not None and target_gap is not None and target_gap > 0:
        actual_gap = float(previous_time) - candidates["time"].astype(float)
        candidates = candidates[
            (actual_gap >= max(10.0, float(target_gap) * 0.40))
            & (actual_gap <= float(target_gap) * 1.75)
        ].copy()
        if candidates.empty:
            return None

    score = _score_flt3_template_peak_choice(
        candidates["time"].astype(float),
        candidates["intensity"].astype(float),
        target_time,
        previous_time,
        target_gap,
        intensity_reference,
    )
    best_idx = score.astype(float).idxmin()
    max_score = 95.0
    if target_gap is not None and target_gap > 0:
        if float(target_gap) < 90.0:
            max_score = 80.0
        elif float(target_gap) < 150.0:
            max_score = 70.0
    if float(score.loc[best_idx]) > max_score:
        return None
    return float(candidates.loc[best_idx, "time"])


def _choose_flt3_forward_candidate_time(
    candidate_df: pd.DataFrame,
    target_time: float,
    previous_time: float | None,
    target_gap: float | None,
    intensity_reference: float | None = None,
) -> float | None:
    if candidate_df.empty:
        return None

    candidates = candidate_df.copy()
    candidates["time"] = candidates["time"].astype(float)
    candidates["intensity"] = candidates["intensity"].astype(float)

    lower_bound = float(target_time - 70.0)
    upper_bound = float(target_time + 70.0)
    if previous_time is not None:
        lower_bound = max(lower_bound, float(previous_time) + 8.0)

    candidates = candidates[
        (candidates["time"] >= lower_bound)
        & (candidates["time"] <= upper_bound)
        & (candidates["intensity"] >= FLT3_MIN_LADDER_PEAK_INTENSITY)
    ].copy()
    if candidates.empty:
        return None

    score = (candidates["time"].astype(float) - float(target_time)).abs()
    if previous_time is not None and target_gap is not None and target_gap > 0:
        actual_gap = candidates["time"].astype(float) - float(previous_time)
        candidates = candidates[
            (actual_gap >= max(10.0, float(target_gap) * 0.45))
            & (actual_gap <= float(target_gap) * 1.75)
        ].copy()
        if candidates.empty:
            return None
        score = (candidates["time"].astype(float) - float(target_time)).abs()
        score = score + (actual_gap.loc[candidates.index] - float(target_gap)).abs() * 1.8

    intensity_penalty = candidates["intensity"].astype(float).apply(
        lambda value: _flt3_ladder_intensity_penalty(float(value), intensity_reference)
    )
    score = score + intensity_penalty - candidates["intensity"].astype(float).clip(lower=0.0, upper=2000.0) * 0.0015
    best_idx = score.astype(float).idxmin()
    return float(candidates.loc[best_idx, "time"])


def _rank_flt3_short_trace_triads(
    candidate_df: pd.DataFrame,
    template_times: np.ndarray,
    trace_last_index: float,
    *,
    limit: int = 10,
) -> list[tuple[float, float, float]]:
    if candidate_df.empty:
        return []

    idx_139 = FLT3_TEMPLATE_STEP_INDEX[139]
    idx_150 = FLT3_TEMPLATE_STEP_INDEX[150]
    idx_160 = FLT3_TEMPLATE_STEP_INDEX[160]
    template_gap_139_150 = float(template_times[idx_150] - template_times[idx_139])
    template_gap_150_160 = float(template_times[idx_160] - template_times[idx_150])

    rows = candidate_df.copy()
    rows["time"] = rows["time"].astype(float)
    rows["intensity"] = rows["intensity"].astype(float)
    rows = rows[
        (rows["time"] >= float(template_times[idx_139]) - 140.0)
        & (rows["time"] <= min(float(trace_last_index) - 120.0, float(template_times[idx_160]) + 950.0))
        & (rows["intensity"] >= FLT3_MIN_LADDER_PEAK_INTENSITY)
    ].copy()
    if rows.empty:
        return []

    scored: list[tuple[float, tuple[float, float, float]]] = []
    for peak_139 in rows.itertuples(index=False):
        candidates_150 = rows[
            (rows["time"] > float(peak_139.time) + 25.0)
            & (rows["time"] < float(peak_139.time) + 90.0)
        ].copy()
        if candidates_150.empty:
            continue
        for peak_150 in candidates_150.itertuples(index=False):
            gap_139_150 = float(peak_150.time - peak_139.time)
            if gap_139_150 < 40.0 or gap_139_150 > 80.0:
                continue
            candidates_160 = rows[
                (rows["time"] > float(peak_150.time) + 25.0)
                & (rows["time"] < float(peak_150.time) + 90.0)
            ].copy()
            if candidates_160.empty:
                continue
            for peak_160 in candidates_160.itertuples(index=False):
                gap_150_160 = float(peak_160.time - peak_150.time)
                span_139_160 = float(peak_160.time - peak_139.time)
                if gap_150_160 < 40.0 or gap_150_160 > 80.0:
                    continue
                if span_139_160 < 95.0 or span_139_160 > 145.0:
                    continue

                score = (
                    abs(gap_139_150 - template_gap_139_150) * 3.0
                    + abs(gap_150_160 - template_gap_150_160) * 3.0
                    + abs((gap_139_150 + gap_150_160) - (template_gap_139_150 + template_gap_150_160)) * 1.2
                    - (
                        float(peak_139.intensity)
                        + float(peak_150.intensity)
                        + float(peak_160.intensity)
                    )
                    * 0.01
                )
                scored.append(
                    (
                        float(score),
                        (
                            float(peak_139.time),
                            float(peak_150.time),
                            float(peak_160.time),
                        ),
                    )
                )

    ranked: list[tuple[float, float, float]] = []
    seen: set[tuple[int, int, int]] = set()
    for _, combo in sorted(scored, key=lambda item: item[0]):
        combo_key = tuple(int(round(value)) for value in combo)
        if combo_key in seen:
            continue
        seen.add(combo_key)
        ranked.append(combo)
        if len(ranked) >= int(limit):
            break
    return ranked


def _attempt_flt3_short_trace_partial_fit(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
    missing_steps: list[float],
    trace_last_index: float,
) -> FsaFile | None:
    template_key = _resolved_flt3_template_key(fsa, assay, analysis_type)
    if template_key is None:
        return None

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return None

    full_expected_steps = _flt3_expected_ladder_steps(fsa)
    template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
    if template_times.size != full_expected_steps.size:
        return None

    idx_139 = FLT3_TEMPLATE_STEP_INDEX[139]
    idx_300 = FLT3_TEMPLATE_STEP_INDEX[300]
    last_missing_idx = FLT3_TEMPLATE_STEP_INDEX[int(round(min(missing_steps)))]
    last_available_idx = max(0, last_missing_idx - 1)
    if last_available_idx < idx_300:
        return None

    candidate_df = get_ladder_candidates(fsa).copy()
    if candidate_df.empty:
        return None

    candidate_df["time"] = candidate_df["time"].astype(float)
    candidate_df["intensity"] = candidate_df["intensity"].astype(float)
    candidate_df = candidate_df.sort_values("time").reset_index(drop=True)

    triads = _rank_flt3_short_trace_triads(candidate_df, template_times, trace_last_index)
    if not triads:
        return None

    best_trial: FsaFile | None = None
    best_qc: dict[str, float | int] | None = None
    best_subset_indices: list[int] | None = None

    template_139 = float(template_times[idx_139])
    template_gap_139_150 = float(template_times[idx_139 + 1] - template_times[idx_139])
    template_gap_150_160 = float(template_times[idx_139 + 2] - template_times[idx_139 + 1])

    for peak_139, peak_150, peak_160 in triads:
        scale = (
            ((float(peak_150) - float(peak_139)) / template_gap_139_150)
            + ((float(peak_160) - float(peak_150)) / template_gap_150_160)
        ) / 2.0
        if scale < 0.90 or scale > 1.25:
            continue

        aligned_times = (template_times - template_139) * scale + float(peak_139)
        mapping_times: dict[int, float] = {
            idx_139: float(peak_139),
            idx_139 + 1: float(peak_150),
            idx_139 + 2: float(peak_160),
        }
        previous_time = float(peak_160)
        previous_target_time = float(aligned_times[idx_139 + 2])

        for step_idx in range(idx_139 + 3, last_available_idx + 1):
            target_time = float(aligned_times[step_idx])
            if target_time >= float(trace_last_index) - 10.0:
                break
            target_gap = float(target_time - previous_target_time)
            intensity_reference = _flt3_ladder_intensity_reference(trace, mapping_times)
            peak_time = _choose_flt3_forward_candidate_time(
                candidate_df,
                target_time,
                previous_time,
                target_gap,
                intensity_reference,
            )
            if peak_time is None:
                break
            mapping_times[int(step_idx)] = float(peak_time)
            previous_time = float(peak_time)
            previous_target_time = float(target_time)

        last_kept_idx = max(mapping_times)
        if last_kept_idx < idx_300:
            continue
        subset_indices = list(range(idx_139, last_kept_idx + 1))
        if any(index not in mapping_times for index in subset_indices):
            continue

        try:
            trial = copy.deepcopy(fsa)
            subset_steps = full_expected_steps[subset_indices].copy()
            subset_peak_times = np.asarray(
                [float(mapping_times[index]) for index in subset_indices],
                dtype=float,
            )
            if np.any(np.diff(subset_peak_times) <= 0):
                continue
            trial.expected_ladder_steps = subset_steps.copy()
            trial.ladder_steps = subset_steps.copy()
            trial.best_size_standard = subset_peak_times.copy()
            trial = fit_size_standard_to_ladder(trial)
        except Exception:
            continue

        if not getattr(trial, "fitted_to_model", False):
            continue

        qc = compute_ladder_qc_metrics(trial)
        if not np.isfinite(float(qc.get("r2", float("nan")))):
            continue
        if float(qc.get("max_abs_error_bp", float("inf"))) > FLT3_REVIEW_MAX_RESIDUAL_BP:
            continue

        if best_qc is None:
            best_trial = trial
            best_qc = qc
            best_subset_indices = subset_indices
            continue

        best_max = float(best_qc.get("max_abs_error_bp", float("inf")))
        best_r2 = float(best_qc.get("r2", float("-inf")))
        current_max = float(qc.get("max_abs_error_bp", float("inf")))
        current_r2 = float(qc.get("r2", float("-inf")))
        if (
            max(subset_indices) > max(best_subset_indices or [])
            and current_max <= best_max + 0.35
        ) or current_max + 0.05 < best_max or (
            abs(current_max - best_max) <= 0.05 and current_r2 > best_r2 + 1e-6
        ):
            best_trial = trial
            best_qc = qc
            best_subset_indices = subset_indices

    if best_trial is None or best_subset_indices is None:
        return None

    kept_step_values = [float(full_expected_steps[index]) for index in best_subset_indices]
    omitted_steps = [
        float(full_expected_steps[index])
        for index in range(max(best_subset_indices) + 1, int(full_expected_steps.size))
    ]
    best_trial.expected_ladder_steps = full_expected_steps.copy()
    best_trial.ladder_missing_expected_steps = omitted_steps
    best_trial.ladder_fit_strategy = "short_trace_partial"
    best_trial.ladder_review_required = False
    best_trial.ladder_fit_note = (
        f"ROX DATA4 trace ends at scan {trace_last_index:.0f}. "
        "Partial FLT3 ladder fit anchored from 139/150/160 and trimmed to the recorded tail "
        f"(kept {', '.join(f'{step:.0f}' for step in kept_step_values)} bp; "
        f"omitted {', '.join(f'{step:.0f}' for step in omitted_steps)} bp)."
    )
    return best_trial


def _template_mapping_payload_for_reference_times(
    fsa: FsaFile,
    template_times: np.ndarray,
    expected_steps: np.ndarray,
    candidate_df: pd.DataFrame,
    anchor_450: float,
    anchor_490: float,
    anchor_500: float,
    *,
    template_label: object,
) -> dict[str, object] | None:
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    candidate_times = [float(value) for value in candidate_df.get("time", pd.Series(dtype=float)).tolist()]
    template_500 = float(template_times[-1])
    template_450 = float(template_times[-3])

    if anchor_500 <= anchor_490 or anchor_490 <= anchor_450:
        return None

    template_gap = template_500 - template_450
    anchor_gap = anchor_500 - anchor_450
    if template_gap <= 0 or anchor_gap <= 0:
        return None

    # The high-end anchors are excellent for selecting the run family, but
    # extrapolating that local 450-500 scale down to 35 bp over-stretches
    # some GS500ROX runs.  A simple run offset better matches the manual fits.
    aligned_times = template_times + (float(anchor_500) - template_500)
    expected_490 = float(aligned_times[-2])
    if abs(float(anchor_490) - expected_490) > 80.0:
        return None

    mapping_times: dict[int, float] = {}
    manual_candidates: list[float] = []
    previous_time: float | None = None
    previous_target_time: float | None = None

    for step_idx in range(len(expected_steps) - 1, -1, -1):
        target_time = float(aligned_times[step_idx])
        if int(step_idx) == len(expected_steps) - 1:
            peak_time = float(anchor_500)
        elif int(step_idx) == len(expected_steps) - 2:
            peak_time = float(anchor_490)
        elif int(step_idx) == len(expected_steps) - 3:
            peak_time = float(anchor_450)
        else:
            target_gap = None
            if previous_target_time is not None:
                target_gap = float(previous_target_time - target_time)
            intensity_reference = _flt3_ladder_intensity_reference(trace, mapping_times)
            peak_time = _choose_template_candidate_time(
                candidate_df,
                target_time,
                previous_time,
                target_gap,
                intensity_reference,
            )
            if peak_time is None:
                step_bp = float(expected_steps[step_idx])
                if step_bp <= 100.0:
                    radius = 75
                else:
                    radius = 65 if step_bp >= 400.0 else (42 if step_idx >= len(expected_steps) - 6 else 95)
                lower_bound = None
                if previous_time is not None and target_gap is not None:
                    lower_bound = float(previous_time) - float(target_gap) * 1.75
                upper_bound = None if previous_time is None else previous_time - 8.0
                peak_time = _choose_template_trace_peak(
                    trace,
                    target_time,
                    previous_time,
                    target_gap,
                    intensity_reference,
                    search_radius=radius,
                    lower_bound=lower_bound,
                    upper_bound=upper_bound,
                    min_intensity=_flt3_template_rescue_trace_min_intensity(step_bp),
                    min_prominence=_flt3_template_rescue_trace_min_prominence(step_bp),
                )
        if peak_time is None:
            return None
        if previous_time is not None and peak_time >= previous_time:
            return None

        mapping_times[step_idx] = float(peak_time)
        if _candidate_index_for_time(candidate_times, peak_time, tolerance=1.5) is None:
            manual_candidates.append(float(peak_time))
        previous_time = float(peak_time)
        previous_target_time = float(target_time)

    return {
        "mapping": {},
        "mapping_times": {int(k): float(v) for k, v in mapping_times.items()},
        "manual_candidates": sorted({float(value) for value in manual_candidates}),
        "template_label": template_label,
    }


def _template_mapping_payload_for_anchors(
    fsa: FsaFile,
    template_key: tuple[str, str],
    template_times: np.ndarray,
    expected_steps: np.ndarray,
    candidate_df: pd.DataFrame,
    anchor_450: float,
    anchor_490: float,
    anchor_500: float,
) -> dict[str, object] | None:
    return _template_mapping_payload_for_reference_times(
        fsa,
        template_times,
        expected_steps,
        candidate_df,
        anchor_450,
        anchor_490,
        anchor_500,
        template_label=template_key,
    )


def _template_mapping_payloads_for_scaled_endpoints(
    fsa: FsaFile,
    template_key: tuple[str, str],
    template_times: np.ndarray,
    expected_steps: np.ndarray,
    candidate_df: pd.DataFrame,
    anchor_candidates: pd.DataFrame,
    anchor_450: float,
    anchor_490: float,
    anchor_500: float,
    *,
    limit: int = 8,
) -> list[dict[str, object]]:
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0 or candidate_df.empty or anchor_candidates.empty:
        return []

    template_450 = float(template_times[-3])
    template_490 = float(template_times[-2])
    template_500 = float(template_times[-1])

    candidates = anchor_candidates.copy()
    candidates["time"] = candidates["time"].astype(float)
    candidates["intensity"] = candidates["intensity"].astype(float)

    candidate_times = [
        float(value) for value in candidate_df.get("time", pd.Series(dtype=float)).tolist()
    ]
    payloads: list[tuple[float, dict[str, object]]] = []
    seen: set[tuple[int, ...]] = set()

    for low_anchor_idx in (0, 1, 2, 4):
        low_anchor_bp = float(expected_steps[low_anchor_idx])
        template_low = float(template_times[low_anchor_idx])
        template_span = template_500 - template_low
        if template_span <= 0.0:
            continue

        low_min = float(anchor_500) - template_span * 1.18
        low_max = float(anchor_500) - template_span * 0.90
        low_candidates = candidates[
            (candidates["time"] >= low_min)
            & (candidates["time"] <= low_max)
            & (candidates["time"] < float(anchor_450) - 1500.0)
            & (
                (candidates["intensity"] >= FLT3_MIN_LADDER_PEAK_INTENSITY)
                | (candidates.get("source", "").astype(str) == "low_prominent")
            )
        ].copy()
        low_trace_candidates = _low_end_trace_peak_candidates(
            trace,
            min_time=low_min,
            max_time=min(low_max, float(anchor_450) - 1500.0),
        )
        if not low_trace_candidates.empty:
            low_candidates = pd.concat([low_candidates, low_trace_candidates], ignore_index=True)
        if low_candidates.empty:
            continue
        low_candidates["time_key"] = low_candidates["time"].astype(float).round().astype(int)
        low_candidates = (
            low_candidates.sort_values(["intensity", "time"], ascending=[False, False])
            .drop_duplicates(subset=["time_key"], keep="first")
            .sort_values("time")
            .reset_index(drop=True)
        )

        for low_peak in low_candidates.itertuples(index=False):
            anchor_low = float(low_peak.time)
            scale = (float(anchor_500) - anchor_low) / template_span
            if scale < 0.90 or scale > 1.18:
                continue

            aligned_times = (template_times - template_low) * scale + anchor_low
            if abs(float(anchor_450) - float(aligned_times[-3])) > 90.0:
                continue
            if abs(float(anchor_490) - float(aligned_times[-2])) > 90.0:
                continue

            mapping_times: dict[int, float] = {
                int(low_anchor_idx): anchor_low,
                int(len(expected_steps) - 3): float(anchor_450),
                int(len(expected_steps) - 2): float(anchor_490),
                int(len(expected_steps) - 1): float(anchor_500),
            }
            if int(low_anchor_idx) == FLT3_TEMPLATE_STEP_INDEX[139]:
                peak_160 = _choose_template_candidate_time(
                    candidate_df,
                    float(aligned_times[6]),
                    None,
                    None,
                    None,
                )
                if peak_160 is None:
                    peak_160 = _choose_template_trace_peak(
                        trace,
                        float(aligned_times[6]),
                        None,
                        None,
                        None,
                        search_radius=55,
                        lower_bound=float(anchor_low) + 20.0,
                        min_intensity=_flt3_template_rescue_trace_min_intensity(float(expected_steps[6])),
                        min_prominence=_flt3_template_rescue_trace_min_prominence(float(expected_steps[6])),
                    )
                if peak_160 is None or peak_160 <= float(anchor_low) + 15.0:
                    continue

                target_gap_150_160 = float(aligned_times[6] - aligned_times[5])
                peak_150 = _choose_template_candidate_time(
                    candidate_df,
                    float(aligned_times[5]),
                    float(peak_160),
                    target_gap_150_160,
                    None,
                )
                if peak_150 is None:
                    peak_150 = _choose_template_trace_peak(
                        trace,
                        float(aligned_times[5]),
                        float(peak_160),
                        target_gap_150_160,
                        None,
                        search_radius=45,
                        lower_bound=float(anchor_low) + 10.0,
                        upper_bound=float(peak_160) - 8.0,
                        min_intensity=_flt3_template_rescue_trace_min_intensity(float(expected_steps[5])),
                        min_prominence=_flt3_template_rescue_trace_min_prominence(float(expected_steps[5])),
                    )
                if peak_150 is None or peak_150 <= float(anchor_low) + 8.0 or peak_150 >= float(peak_160) - 8.0:
                    continue

                mapping_times[5] = float(peak_150)
                mapping_times[6] = float(peak_160)
            previous_time = float(anchor_500)
            previous_target_time = float(aligned_times[-1])
            ok = True

            for step_idx in range(len(expected_steps) - 2, -1, -1):
                if step_idx in mapping_times:
                    previous_time = float(mapping_times[step_idx])
                    previous_target_time = float(aligned_times[step_idx])
                    continue

                target_time = float(aligned_times[step_idx])
                target_gap = float(previous_target_time - target_time)
                intensity_reference = _flt3_ladder_intensity_reference(trace, mapping_times)
                peak_time = _choose_template_candidate_time(
                    candidate_df,
                    target_time,
                    previous_time,
                    target_gap,
                    intensity_reference,
                )
                if peak_time is None:
                    step_bp = float(expected_steps[step_idx])
                    radius = 75 if step_bp <= 100.0 else 42
                    if step_bp >= 400.0:
                        radius = 65
                    elif step_idx >= len(expected_steps) - 6:
                        radius = 42
                    else:
                        radius = 95
                    lower_bound = float(previous_time) - float(target_gap) * 1.75
                    upper_bound = previous_time - 8.0
                    peak_time = _choose_template_trace_peak(
                        trace,
                        target_time,
                        previous_time,
                        target_gap,
                        intensity_reference,
                        search_radius=radius,
                        lower_bound=lower_bound,
                        upper_bound=upper_bound,
                        min_intensity=_flt3_template_rescue_trace_min_intensity(step_bp),
                        min_prominence=_flt3_template_rescue_trace_min_prominence(step_bp),
                    )
                if peak_time is None or peak_time >= previous_time:
                    ok = False
                    break

                mapping_times[int(step_idx)] = float(peak_time)
                previous_time = float(peak_time)
                previous_target_time = target_time

            if not ok or len(mapping_times) != int(expected_steps.size):
                continue

            mapping_key = tuple(int(round(mapping_times[index])) for index in range(len(expected_steps)))
            if mapping_key in seen:
                continue
            seen.add(mapping_key)

            manual_candidates: list[float] = []
            for peak_time in mapping_times.values():
                if _candidate_index_for_time(candidate_times, peak_time, tolerance=1.5) is None:
                    manual_candidates.append(float(peak_time))

            shape_penalty = _flt3_mapping_shape_penalty(mapping_times, template_times, trace)
            payloads.append(
                (
                    float(shape_penalty),
                    {
                        "mapping": {},
                        "mapping_times": {int(k): float(v) for k, v in mapping_times.items()},
                        "manual_candidates": sorted({float(value) for value in manual_candidates}),
                        "template_key": template_key,
                        "template_anchor_mode": f"scaled_{int(low_anchor_bp)}_500",
                    },
                )
            )

    payloads = sorted(payloads, key=lambda item: item[0])
    return [payload for _, payload in payloads[: int(limit)]]


def _template_mapping_payload(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
    *,
    template_key: tuple[str, str] | None = None,
) -> dict[str, object] | None:
    template_key = template_key or _resolved_flt3_template_key(fsa, assay, analysis_type)
    if template_key is None:
        return None

    template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
    expected_steps = _flt3_expected_ladder_steps(fsa)
    if expected_steps.size == 0 or template_times.size != expected_steps.size:
        return None

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return None

    candidate_df = _flt3_candidate_pool(fsa).sort_values("time").reset_index(drop=True)
    raw_trace = np.asarray(getattr(getattr(fsa, "fsa", {}), "get", lambda *_: [])("DATA4"), dtype=float)
    if raw_trace.shape != trace.shape:
        raw_trace = trace
    anchor_candidates = _merged_anchor_candidates(trace, candidate_df, reference_trace=raw_trace)

    anchor_500 = _mapped_peak_time_for_step(fsa, 500.0)
    anchor_490 = _mapped_peak_time_for_step(fsa, 490.0)
    anchor_450 = _mapped_peak_time_for_step(fsa, 450.0)
    template_500 = float(template_times[-1])
    template_490 = float(template_times[-2])
    template_450 = float(template_times[-3])

    if not _flt3_high_end_anchors_are_plausible(anchor_490, anchor_500, template_490, template_500):
        anchor_500 = None
        anchor_490 = None
        anchor_450 = None

    if anchor_500 is None or anchor_490 is None or anchor_450 is None:
        anchor_combo = _select_flt3_high_end_anchor_combo(
            anchor_candidates,
            template_450,
            template_490,
            template_500,
        )
        if anchor_combo is not None:
            anchor_450, anchor_490, anchor_500 = anchor_combo

    if anchor_500 is None:
        low_500, high_500 = _flt3_template_window(template_500, 500.0)
        anchor_500 = _snap_trace_peak(
            trace,
            template_500,
            search_radius=180,
            lower_bound=low_500,
            upper_bound=high_500,
        )
    if anchor_490 is None:
        guess_490 = template_490 if anchor_500 is None else anchor_500 - (template_500 - template_490)
        low_490, high_490 = _flt3_template_window(template_490, 490.0)
        upper_bound_490 = min(high_490, float(anchor_500) - 20.0) if anchor_500 is not None else high_490
        anchor_490 = _snap_trace_peak(
            trace,
            guess_490,
            search_radius=180,
            lower_bound=low_490,
            upper_bound=upper_bound_490,
        )
    if anchor_450 is None:
        guess_450 = template_450 if anchor_490 is None else anchor_490 - (template_490 - template_450)
        low_450, high_450 = _flt3_template_window(template_450, 450.0)
        upper_bound_450 = min(high_450, float(anchor_490) - 140.0) if anchor_490 is not None else high_450
        anchor_450 = _snap_trace_peak(
            trace,
            guess_450,
            search_radius=180,
            lower_bound=low_450,
            upper_bound=upper_bound_450,
        )

    if anchor_500 is None or anchor_490 is None or anchor_450 is None:
        return None
    return _template_mapping_payload_for_anchors(
        fsa,
        template_key,
        template_times,
        expected_steps,
        candidate_df,
        float(anchor_450),
        float(anchor_490),
        float(anchor_500),
    )


def _template_review_scaffold_payload(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
    *,
    template_key: tuple[str, str] | None = None,
) -> dict[str, object] | None:
    payload = _template_mapping_payload(
        fsa,
        assay,
        analysis_type,
        template_key=template_key,
    )
    if payload is not None:
        return payload

    template_key = template_key or _resolved_flt3_template_key(fsa, assay, analysis_type)
    if template_key is None:
        return None

    template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
    expected_steps = _flt3_expected_ladder_steps(fsa)
    if expected_steps.size == 0 or template_times.size != expected_steps.size:
        return None

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return None

    candidate_df = _flt3_candidate_pool(fsa).sort_values("time").reset_index(drop=True)
    candidate_times = candidate_df["time"].astype(float).tolist() if not candidate_df.empty else []

    mapping_times: dict[int, float] = {}
    manual_candidates: list[float] = []
    previous_time: float | None = None
    for index, (step_bp, template_time) in enumerate(zip(expected_steps, template_times, strict=False)):
        low, high = _flt3_template_window(float(template_time), float(step_bp))
        search_radius = 220 if float(step_bp) >= 400.0 else 170
        snapped = _snap_trace_peak(
            trace,
            float(template_time),
            search_radius=search_radius,
            lower_bound=max(0.0, float(low) - 120.0),
            upper_bound=min(float(trace.size - 1), float(high) + 220.0),
        )
        peak_time = float(snapped) if snapped is not None else float(template_time)
        if previous_time is not None:
            peak_time = max(float(previous_time) + 6.0, peak_time)
        mapping_times[int(index)] = peak_time
        previous_time = peak_time
        if _candidate_index_for_time(candidate_times, peak_time, tolerance=1.5) is None:
            manual_candidates.append(float(peak_time))

    return {
        "mapping": {},
        "mapping_times": mapping_times,
        "manual_candidates": sorted({float(value) for value in manual_candidates}),
        "template_key": template_key,
        "template_anchor_mode": "review_scaffold",
    }


def _flt3_trace_peak_candidates(trace: np.ndarray) -> pd.DataFrame:
    from scipy.signal import find_peaks

    if trace.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "prominence", "source"])

    try:
        baseline = estimate_running_baseline(trace, bin_size=200, quantile=0.10)
        corrected = np.clip(np.asarray(trace, dtype=float) - baseline, a_min=0, a_max=None)
    except Exception:
        corrected = np.asarray(trace, dtype=float)

    peak_idx, props = find_peaks(
        corrected,
        distance=12,
        prominence=35.0,
        height=40.0,
    )
    if peak_idx.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "prominence", "source"])

    rows = pd.DataFrame(
        {
            "time": peak_idx.astype(float),
            "intensity": corrected[peak_idx].astype(float),
            "prominence": np.asarray(props.get("prominences", []), dtype=float),
            "source": "trace_bootstrap",
        }
    )
    rows = rows[
        (rows["time"].astype(float) >= 1350.0)
        & (rows["time"].astype(float) <= 4400.0)
    ].copy()
    if rows.empty:
        return pd.DataFrame(columns=["time", "intensity", "prominence", "source"])
    rows["time_key"] = rows["time"].round().astype(int)
    rows = (
        rows.sort_values(["prominence", "intensity"], ascending=[False, False])
        .drop_duplicates(subset=["time_key"], keep="first")
        .sort_values("time")
        .reset_index(drop=True)
    )
    return rows.loc[:, ["time", "intensity", "prominence", "source"]]


def _attempt_flt3_d835_family_bootstrap_fit(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> FsaFile | None:
    if assay != "FLT3-D835" or str(analysis_type or "").strip().lower() not in {"", "standard"}:
        return None

    expected_steps = _flt3_expected_ladder_steps(fsa)
    if expected_steps.size == 0:
        return None

    template_key = _resolved_flt3_template_key(fsa, assay, analysis_type)
    if template_key is None:
        return None
    template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
    if template_times.size != expected_steps.size:
        return None

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return None

    candidate_df = _flt3_trace_peak_candidates(trace)
    if candidate_df.empty or len(candidate_df) < int(expected_steps.size):
        return None

    candidate_times = candidate_df["time"].astype(float).tolist()
    candidate_intensities = candidate_df["intensity"].astype(float).tolist()
    ref_gaps = np.diff(template_times).astype(float)
    start_min = float(template_times[0] - 140.0)
    start_max = float(template_times[0] + 120.0)

    best_trial: FsaFile | None = None
    best_qc: dict[str, float | int] | None = None
    best_score = float("inf")

    for start_idx, start_time in enumerate(candidate_times):
        if start_time < start_min or start_time > start_max:
            continue

        mapping_times = [float(start_time)]
        picked_indices = [int(start_idx)]
        prev_idx = int(start_idx)
        prev_time = float(start_time)
        ok = True
        score = 0.0

        for step_offset, ref_gap in enumerate(ref_gaps, start=1):
            target_time = float(start_time + (template_times[step_offset] - template_times[0]))
            low_gap = max(18.0, float(ref_gap) * 0.45)
            high_gap = max(low_gap + 18.0, float(ref_gap) * 1.65)
            chosen_idx: int | None = None
            chosen_score = float("inf")
            for candidate_idx in range(prev_idx + 1, len(candidate_times)):
                candidate_time = float(candidate_times[candidate_idx])
                gap = candidate_time - prev_time
                if gap < low_gap:
                    continue
                if gap > high_gap:
                    break
                local_score = (
                    abs(gap - float(ref_gap)) * 1.8
                    + abs(candidate_time - target_time) * 0.35
                    - min(float(candidate_intensities[candidate_idx]), 2500.0) * 0.002
                )
                if step_offset == 1 and candidate_time <= prev_time + 35.0:
                    local_score += 200.0
                if local_score < chosen_score:
                    chosen_idx = int(candidate_idx)
                    chosen_score = float(local_score)
            if chosen_idx is None:
                ok = False
                break
            prev_idx = int(chosen_idx)
            prev_time = float(candidate_times[chosen_idx])
            picked_indices.append(prev_idx)
            mapping_times.append(prev_time)
            score += chosen_score

        if not ok or len(mapping_times) != int(expected_steps.size):
            continue

        try:
            trial = copy.deepcopy(fsa)
            trial.expected_ladder_steps = expected_steps.copy()
            trial.ladder_steps = expected_steps.copy()
            trial = apply_manual_ladder_mapping(
                trial,
                {
                    "mapping": {},
                    "mapping_times": {idx: float(value) for idx, value in enumerate(mapping_times)},
                    "manual_candidates": [],
                },
            )
        except Exception:
            continue

        qc = compute_ladder_qc_metrics(trial)
        trial_r2 = float(qc.get("r2", float("-inf")))
        trial_max = float(qc.get("max_abs_error_bp", float("inf")))
        if not np.isfinite(trial_r2) or not np.isfinite(trial_max):
            continue
        if trial_r2 < 0.9992 or trial_max > 3.0:
            continue

        score += _flt3_mapping_shape_penalty(
            {idx: float(value) for idx, value in enumerate(mapping_times)},
            template_times,
            trace,
        )
        if (
            best_qc is None
            or trial_max + 0.05 < float(best_qc.get("max_abs_error_bp", float("inf")))
            or (
                abs(trial_max - float(best_qc.get("max_abs_error_bp", float("inf")))) <= 0.05
                and trial_r2 > float(best_qc.get("r2", float("-inf"))) + 1e-6
            )
            or (
                abs(trial_max - float(best_qc.get("max_abs_error_bp", float("inf")))) <= 0.1
                and abs(trial_r2 - float(best_qc.get("r2", float("-inf")))) <= 1e-6
                and score < best_score
            )
        ):
            best_trial = trial
            best_qc = qc
            best_score = float(score)

    if best_trial is None:
        return None

    setattr(best_trial, "ladder_fit_strategy", "flt3_template_rescue")
    setattr(
        best_trial,
        "ladder_fit_note",
        "FLT3 D835 family bootstrap rescue selected a full GS500ROX ladder sequence from the raw trace.",
    )
    setattr(best_trial, "ladder_review_required", False)
    setattr(best_trial, "ladder_missing_expected_steps", [])
    return best_trial


def _attempt_flt3_template_fit(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> FsaFile | None:
    template_keys = _rank_flt3_template_keys_for_fsa(fsa, assay, analysis_type)
    if not template_keys:
        return None

    expected_steps = _flt3_expected_ladder_steps(fsa)
    if expected_steps.size == 0:
        return None

    current_qc = compute_ladder_qc_metrics(fsa)
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return None

    candidate_df = _flt3_candidate_pool(fsa).sort_values("time").reset_index(drop=True)
    raw_trace = np.asarray(getattr(getattr(fsa, "fsa", {}), "get", lambda *_: [])("DATA4"), dtype=float)
    if raw_trace.shape != trace.shape:
        raw_trace = trace
    anchor_candidates = _merged_anchor_candidates(trace, candidate_df, reference_trace=raw_trace)
    auto_500 = _mapped_peak_time_for_step(fsa, 500.0)
    auto_490 = _mapped_peak_time_for_step(fsa, 490.0)
    auto_450 = _mapped_peak_time_for_step(fsa, 450.0)
    current_mapping_times = _mapping_times_from_fsa(fsa, expected_steps)

    current_max = float(current_qc.get("max_abs_error_bp", float("inf")))
    current_mean = float(current_qc.get("mean_abs_error_bp", float("inf")))
    current_r2 = float(current_qc.get("r2", float("-inf")))
    if not np.isfinite(current_max):
        current_max = float("inf")
    if not np.isfinite(current_mean):
        current_mean = float("inf")
    if not np.isfinite(current_r2):
        current_r2 = float("-inf")
    best_trial: FsaFile | None = None
    best_payload: dict[str, object] | None = None
    best_qc: dict[str, float | int] | None = None

    seen_combos: set[tuple[int, int, int]] = set()
    seen_payloads: set[tuple[int, ...]] = set()
    for template_key in template_keys:
        template_times = np.asarray(FLT3_TEMPLATE_TIMES[template_key], dtype=float)
        if template_times.size != expected_steps.size:
            continue
        current_shape = _flt3_mapping_shape_penalty(current_mapping_times, template_times, trace)
        current_invalid = _flt3_fit_is_geometrically_invalid(
            current_mapping_times,
            expected_steps,
            template_times,
            trace,
        )

        template_500 = float(template_times[-1])
        template_490 = float(template_times[-2])
        template_450 = float(template_times[-3])

        combos: list[tuple[float, float, float]] = []
        if (
            auto_450 is not None
            and _flt3_high_end_anchors_are_plausible(auto_490, auto_500, template_490, template_500)
            and _flt3_peak_meets_min_intensity(trace, auto_450)
            and _flt3_peak_meets_min_intensity(trace, auto_490)
            and _flt3_peak_meets_min_intensity(trace, auto_500)
        ):
            combos.append((float(auto_450), float(auto_490), float(auto_500)))
        combos.extend(
            _rank_flt3_high_end_anchor_combos(
                anchor_candidates,
                template_450,
                template_490,
                template_500,
                limit=20,
            )
        )

        fallback_payload = _template_mapping_payload(
            fsa,
            assay,
            analysis_type,
            template_key=template_key,
        )
        if fallback_payload is not None:
            combos.append(
                (
                    float(fallback_payload["mapping_times"][len(expected_steps) - 3]),
                    float(fallback_payload["mapping_times"][len(expected_steps) - 2]),
                    float(fallback_payload["mapping_times"][len(expected_steps) - 1]),
                )
            )

        for anchor_450, anchor_490, anchor_500 in combos:
            combo_key = (
                int(round(anchor_450)),
                int(round(anchor_490)),
                int(round(anchor_500)),
            )
            if combo_key in seen_combos:
                continue
            seen_combos.add(combo_key)
            candidate_payloads: list[dict[str, object]] = []
            payload = _template_mapping_payload_for_anchors(
                fsa,
                template_key,
                template_times,
                expected_steps,
                candidate_df,
                float(anchor_450),
                float(anchor_490),
                float(anchor_500),
            )
            if payload is not None:
                candidate_payloads.append(payload)
            candidate_payloads.extend(
                _template_mapping_payloads_for_scaled_endpoints(
                    fsa,
                    template_key,
                    template_times,
                    expected_steps,
                    candidate_df,
                    anchor_candidates,
                    float(anchor_450),
                    float(anchor_490),
                    float(anchor_500),
                )
            )

            for payload in candidate_payloads:
                payload_mapping = payload["mapping_times"]
                payload_key = tuple(
                    int(round(float(payload_mapping[index]))) for index in range(len(expected_steps))
                )
                if payload_key in seen_payloads:
                    continue
                seen_payloads.add(payload_key)
                try:
                    trial = copy.deepcopy(fsa)
                    trial.expected_ladder_steps = expected_steps.copy()
                    trial.ladder_steps = trial.expected_ladder_steps.copy()
                    trial = apply_manual_ladder_mapping(
                        trial,
                        {
                            "mapping": payload["mapping"],
                            "mapping_times": payload["mapping_times"],
                            "manual_candidates": payload["manual_candidates"],
                        },
                    )
                except Exception:
                    continue

                rescued_qc = compute_ladder_qc_metrics(trial)
                rescued_max = float(rescued_qc.get("max_abs_error_bp", float("inf")))
                rescued_mean = float(rescued_qc.get("mean_abs_error_bp", float("inf")))
                rescued_r2 = float(rescued_qc.get("r2", float("-inf")))
                rescued_shape = _flt3_mapping_shape_penalty(
                    payload["mapping_times"],
                    template_times,
                    trace,
                )

                improved = (
                    rescued_max + 0.2 < current_max
                    or rescued_mean + 0.1 < current_mean
                    or rescued_r2 > current_r2 + 1e-5
                    or (
                        current_invalid
                        and rescued_r2 >= 0.999
                        and rescued_shape <= FLT3_ACCEPTABLE_RESCUE_SHAPE_PENALTY
                        and rescued_shape + 100.0 < current_shape
                    )
                    or (
                        rescued_r2 >= current_r2 - 2e-5
                        and rescued_max <= current_max + 0.15
                        and rescued_mean <= current_mean + 0.05
                        and rescued_shape + 25.0 < current_shape
                    )
                )
                if not improved:
                    continue

                if best_qc is None:
                    best_trial = trial
                    best_payload = payload
                    best_qc = dict(rescued_qc)
                    best_qc["_shape_penalty"] = rescued_shape
                    continue

                best_r2 = float(best_qc.get("r2", float("-inf")))
                best_max = float(best_qc.get("max_abs_error_bp", float("inf")))
                best_mean = float(best_qc.get("mean_abs_error_bp", float("inf")))
                best_shape = float(best_qc.get("_shape_penalty", float("inf")))
                both_high_quality = rescued_r2 >= 0.999 and best_r2 >= 0.999
                if (
                    (
                        both_high_quality
                        and rescued_shape + 20.0 < best_shape
                        and rescued_max <= best_max + 0.50
                        and rescued_mean <= best_mean + 0.25
                    )
                    or (
                        rescued_r2 > best_r2 + 1e-6
                        and (not both_high_quality or rescued_shape <= best_shape + 40.0)
                    )
                    or (abs(rescued_r2 - best_r2) <= 1e-6 and rescued_max + 0.05 < best_max)
                    or (
                        abs(rescued_r2 - best_r2) <= 1e-6
                        and abs(rescued_max - best_max) <= 0.05
                        and rescued_mean + 0.02 < best_mean
                    )
                    or (
                        abs(rescued_r2 - best_r2) <= 2e-5
                        and abs(rescued_max - best_max) <= 0.2
                        and rescued_shape + 15.0 < best_shape
                    )
                ):
                    best_trial = trial
                    best_payload = payload
                    best_qc = dict(rescued_qc)
                    best_qc["_shape_penalty"] = rescued_shape

    if best_trial is None or best_payload is None or best_qc is None:
        return None

    rescued_max = float(best_qc.get("max_abs_error_bp", float("inf")))
    template_label = best_payload.get("template_label", best_payload.get("template_key", ("FLT3", "template")))
    if isinstance(template_label, tuple) and len(template_label) >= 2:
        template_label_text = f"{template_label[0]} / {template_label[1]}"
    else:
        template_label_text = str(template_label)
    setattr(best_trial, "ladder_fit_strategy", "flt3_template_rescue")
    setattr(
        best_trial,
        "ladder_fit_note",
        (
            "FLT3 GS500ROX template rescue applied from a high-end anchor pattern "
            f"({template_label_text})."
        ),
    )
    setattr(best_trial, "ladder_review_required", bool(rescued_max > FLT3_REVIEW_MAX_RESIDUAL_BP))
    setattr(best_trial, "ladder_missing_expected_steps", [])
    return best_trial


def _attempt_flt3_bootstrap_template_fit(
    fsa: FsaFile,
    assay: str,
    analysis_type: str | None,
) -> FsaFile | None:
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    ss_peaks = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    if trace.size == 0 or ss_peaks.size < 3:
        return None

    rescued = _attempt_flt3_template_fit(fsa, assay, analysis_type)
    if rescued is None:
        return None

    sizing_method = str(getattr(rescued, "_flt3_sizing_method", "") or "")
    if not sizing_method:
        setattr(rescued, "_flt3_sizing_method", "template_bootstrap")
    else:
        setattr(rescued, "_flt3_sizing_method", f"{sizing_method}+template_bootstrap")
    return rescued


def _attempt_lenient_rox_fit(
    fsa_path: Path,
    sample_channel: str,
    assay: str,
    analysis_type: str | None,
) -> FsaFile | None:
    normalized_analysis_type = str(analysis_type or "").strip().lower()
    configs = [{"min_h": 5, "min_d": 3}]
    if assay == "FLT3-ITD":
        configs.extend(
            [
                {"min_h": 4, "min_d": 3},
                {"min_h": 4, "min_d": 2},
                {"min_h": 3, "min_d": 2},
            ]
        )
    if normalized_analysis_type in {"ratio_quant", "10x_diluted", "25x_diluted", "undiluted"}:
        configs = [
            {"min_h": 4, "min_d": 2},
            {"min_h": 3, "min_d": 2},
            {"min_h": 3, "min_d": 1},
            {"min_h": 5, "min_d": 3},
        ]

    seen_configs: set[tuple[int, int]] = set()
    best_seed_fsa: FsaFile | None = None
    best_seed_peak_count = -1
    for cfg in configs:
        cfg_key = (int(cfg["min_h"]), int(cfg["min_d"]))
        if cfg_key in seen_configs:
            continue
        seen_configs.add(cfg_key)
        try:
            fsa = FsaFile(
                file=str(fsa_path),
                ladder=FLT3_ROX_LADDER,
                sample_channel=sample_channel,
                min_distance_between_peaks=cfg["min_d"],
                min_size_standard_height=cfg["min_h"],
                size_standard_channel="DATA4",
            )
            fsa = find_size_standard_peaks(fsa)
            ss_peaks = getattr(fsa, "size_standard_peaks", None)
            ss_peak_count = 0 if ss_peaks is None else int(getattr(ss_peaks, "shape", [0])[0])
            if ss_peak_count < 3:
                rescued = _attempt_flt3_d835_family_bootstrap_fit(fsa, assay, analysis_type)
                if rescued is not None:
                    setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                    return rescued
                trace_peak_count = len(_flt3_trace_peak_candidates(np.asarray(getattr(fsa, "size_standard", []), dtype=float)))
                if trace_peak_count > best_seed_peak_count:
                    best_seed_fsa = copy.deepcopy(fsa)
                    best_seed_peak_count = trace_peak_count
                continue
            seed_peak_count = ss_peak_count
            if seed_peak_count > best_seed_peak_count:
                best_seed_fsa = copy.deepcopy(fsa)
                best_seed_peak_count = seed_peak_count
            fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=1.5)
            for _ in range(20):
                fsa = generate_combinations(fsa)
                best = getattr(fsa, "best_size_standard_combinations", None)
                if best is not None and best.shape[0] > 0:
                    break
                fsa.maxium_allowed_distance_between_size_standard_peaks += 10

            best = getattr(fsa, "best_size_standard_combinations", None)
            if best is None or best.shape[0] == 0:
                rescued = _attempt_flt3_bootstrap_template_fit(fsa, assay, analysis_type)
                if rescued is not None:
                    return rescued
                rescued = _attempt_flt3_d835_family_bootstrap_fit(fsa, assay, analysis_type)
                if rescued is not None:
                    setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                    return rescued
                continue

            selected_fit = _select_best_ladder_candidate(fsa)
            if selected_fit is not None:
                fsa = selected_fit
            else:
                fsa = calculate_best_combination_of_size_standard_peaks(fsa)
                if not getattr(fsa, "fitted_to_model", False):
                    rescued = _attempt_flt3_bootstrap_template_fit(fsa, assay, analysis_type)
                    if rescued is not None:
                        return rescued
                    rescued = _attempt_flt3_d835_family_bootstrap_fit(fsa, assay, analysis_type)
                    if rescued is not None:
                        setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                        return rescued
                    fsa = fit_size_standard_to_ladder(fsa)

            if not getattr(fsa, "fitted_to_model", False):
                rescued = _attempt_flt3_bootstrap_template_fit(fsa, assay, analysis_type)
                if rescued is not None:
                    return rescued
                rescued = _attempt_flt3_d835_family_bootstrap_fit(fsa, assay, analysis_type)
                if rescued is not None:
                    setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                    return rescued
                continue

            qc = compute_ladder_qc_metrics(fsa)
            if qc["r2"] >= MIN_R2_QUALITY:
                rescued = _attempt_flt3_template_fit(fsa, assay, analysis_type)
                if rescued is not None:
                    setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                    return rescued
                if _accept_lenient_raw_flt3_fit(qc, assay, analysis_type):
                    setattr(fsa, "_flt3_sizing_method", "spline_lenient")
                    return fsa
                continue

            rescued = _attempt_flt3_template_fit(fsa, assay, analysis_type)
            if rescued is not None:
                setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                return rescued
            short_trace_missing_steps, trace_last_index = _flt3_short_trace_missing_steps(
                fsa,
                assay,
                analysis_type,
            )
            if short_trace_missing_steps:
                short_trace_rescue = _attempt_flt3_short_trace_partial_fit(
                    fsa,
                    assay,
                    analysis_type,
                    short_trace_missing_steps,
                    trace_last_index,
                )
                if short_trace_rescue is not None:
                    setattr(short_trace_rescue, "_flt3_sizing_method", _infer_sizing_method(short_trace_rescue))
                    return short_trace_rescue
            if _accept_lenient_raw_flt3_fit(qc, assay, analysis_type):
                setattr(fsa, "_flt3_sizing_method", "spline_lenient")
                return fsa
        except Exception:
            continue

    if best_seed_fsa is not None:
        rescued = _attempt_flt3_d835_family_bootstrap_fit(best_seed_fsa, assay, analysis_type)
        if rescued is not None:
            setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
            return rescued
        expected_steps = _flt3_expected_ladder_steps(best_seed_fsa)
        review_payload = _template_review_scaffold_payload(
            best_seed_fsa,
            assay,
            analysis_type,
        )
        setattr(best_seed_fsa, "ladder_fit_strategy", "trace_bootstrap_review")
        setattr(
            best_seed_fsa,
            "ladder_fit_note",
            (
                "No valid automatic GS500ROX ladder combination was found, but raw ladder trace peaks were detected. "
                "Manual ladder review is required."
            ),
        )
        setattr(
            best_seed_fsa,
            "ladder_missing_expected_steps",
            [float(value) for value in expected_steps.tolist()] if expected_steps.size else [],
        )
        setattr(best_seed_fsa, "ladder_review_required", True)
        setattr(best_seed_fsa, "_flt3_sizing_method", "trace_bootstrap_review")
        if review_payload is not None:
            setattr(best_seed_fsa, "ladder_review_mapping_times", review_payload["mapping_times"])
            setattr(best_seed_fsa, "ladder_review_manual_candidates", review_payload["manual_candidates"])
            setattr(best_seed_fsa, "ladder_review_template_key", review_payload.get("template_key"))
        return best_seed_fsa

    return None


def _analyse_fsa_candidate(
    fsa_path: Path,
    sample_channel: str,
    assay: str,
    analysis_type: str | None = None,
) -> FsaFile | None:
    fsa = analyse_fsa_rox(
        fsa_path,
        sample_channel,
        ladder_name=FLT3_ROX_LADDER,
        ladder_fit_profile=LADDER_FIT_PROFILE_FLT3_GS500ROX,
    )
    if fsa is not None:
        short_trace_missing_steps, trace_last_index = _flt3_short_trace_missing_steps(
            fsa,
            assay,
            analysis_type,
        )
        if short_trace_missing_steps:
            rescued = _attempt_flt3_short_trace_partial_fit(
                fsa,
                assay,
                analysis_type,
                short_trace_missing_steps,
                trace_last_index,
            )
            if rescued is not None:
                setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
                return rescued
        if _mark_flt3_short_trace_if_needed(fsa, assay, analysis_type):
            setattr(fsa, "_flt3_sizing_method", _infer_sizing_method(fsa))
            return fsa
        rescued = _attempt_flt3_template_fit(fsa, assay, analysis_type)
        if rescued is not None:
            setattr(rescued, "_flt3_sizing_method", _infer_sizing_method(rescued))
            return rescued
        setattr(fsa, "_flt3_sizing_method", _infer_sizing_method(fsa))
        return fsa
    return _attempt_lenient_rox_fit(
        fsa_path,
        sample_channel,
        assay,
        analysis_type,
    )


def _build_entry_from_candidate(fsa_path: Path, meta: dict) -> dict | None:
    fsa = _analyse_fsa_candidate(
        fsa_path,
        meta["primary_peak_channel"],
        meta["assay"],
        meta.get("analysis_type"),
    )
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
    expected_ladder_steps = list(
        map(float, getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])))
    )
    fitted_ladder_steps = list(map(float, getattr(fsa, "ladder_steps", [])))
    ladder_fit_strategy = str(getattr(fsa, "ladder_fit_strategy", "auto_full"))
    ladder_missing_expected_steps = list(
        map(float, getattr(fsa, "ladder_missing_expected_steps", []))
    )
    ladder_fit_note = str(
        getattr(
            fsa,
            "ladder_fit_note",
            "All expected ladder steps were fitted." if not ladder_missing_expected_steps else "Manual ladder review recommended.",
        )
    )
    ladder_max_residual_bp = float(metrics.get("max_abs_error_bp", float("inf")))
    ladder_review_required = bool(
        getattr(fsa, "ladder_review_required", bool(ladder_missing_expected_steps))
        or ladder_max_residual_bp > FLT3_REVIEW_MAX_RESIDUAL_BP
    )
    if ladder_fit_strategy == "manual_adjustment":
        ladder_qc_status = "manual_adjustment"
    elif ladder_review_required:
        ladder_qc_status = "review_required"
    elif float(metrics.get("r2", float("nan"))) > FLT3_LADDER_QC_THRESHOLD:
        ladder_qc_status = "ok"
    else:
        ladder_qc_status = "ladder_qc_failed"
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
        "ladder_fit_strategy": ladder_fit_strategy,
        "ladder_missing_expected_steps": ladder_missing_expected_steps,
        "ladder_fit_note": ladder_fit_note,
        "ladder_review_required": ladder_review_required,
        "ladder_expected_step_count": len(expected_ladder_steps),
        "ladder_fitted_step_count": len(fitted_ladder_steps),
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


def _entry_ranking_key(entry: dict, preferred_injection: int) -> tuple[int, int, float, int, str, str]:
    selected_injection = int(entry.get("selected_injection_time", entry.get("injection_time", 0)) or 0)
    ladder_r2 = float(entry.get("ladder_r2", float("nan")))
    if not np.isfinite(ladder_r2):
        ladder_r2 = float("-inf")
    ladder_steps = int(entry.get("n_ladder_steps", 0) or 0)
    return (
        0 if selected_injection == preferred_injection else 1,
        abs(selected_injection - preferred_injection),
        -ladder_r2,
        -ladder_steps,
        str(entry.get("source_run_dir", "")),
        str(getattr(entry.get("fsa"), "file_name", "")),
    )


def _fallback_entry_ranking_key(entry: dict, preferred_injection: int) -> tuple[int, int, int, int, float, int, str, str]:
    selected_injection = int(entry.get("selected_injection_time", entry.get("injection_time", 0)) or 0)
    ladder_r2 = float(entry.get("ladder_r2", float("nan")))
    if not np.isfinite(ladder_r2):
        ladder_r2 = float("-inf")
    ladder_steps = int(entry.get("n_ladder_steps", 0) or 0)
    return (
        0 if selected_injection == preferred_injection else 1,
        abs(selected_injection - preferred_injection),
        0 if entry.get("ladder_qc_status") == "ok" else 1,
        0 if entry.get("peak_qc_pass") else 1,
        -ladder_r2,
        -ladder_steps,
        str(entry.get("source_run_dir", "")),
        str(getattr(entry.get("fsa"), "file_name", "")),
    )


def _select_best_entry(candidates: list[tuple[Path, dict]]) -> dict | None:
    if not candidates:
        return None

    preferred_injection = _preferred_injection_time(candidates[0][1])
    preferred_only = [
        item for item in candidates
        if int(item[1].get("injection_time", 0) or 0) == preferred_injection
    ]
    pool = preferred_only if preferred_only else candidates
    ordered = sorted(pool, key=lambda item: _candidate_sort_key(item, preferred_injection))

    audit_records: list[dict] = []
    acceptable_entries: list[tuple[dict, str, bool]] = []
    fallback_entries: list[tuple[dict, str, bool]] = []
    preferred_candidates_present = bool(preferred_only)
    preferred_reason = "preferred injection unavailable"

    for path, meta in ordered:
        same_as_preferred = int(meta.get("injection_time", 0) or 0) == preferred_injection
        entry = _build_entry_from_candidate(path, meta)
        if entry is None:
            reason = "ladder_fit_failed"
            if same_as_preferred:
                preferred_reason = reason
            audit_records.append(_candidate_audit_record(path, meta, "rejected", reason))
            continue

        acceptable = entry["ladder_qc_status"] in {"ok", "manual_adjustment"} and entry["peak_qc_pass"]
        candidate_reason = "qc_pass"
        if not acceptable:
            if entry["ladder_qc_status"] != "ok":
                candidate_reason = entry["ladder_qc_status"]
            else:
                candidate_reason = entry["peak_qc_status"]
            if same_as_preferred:
                preferred_reason = candidate_reason
            fallback_entries.append((entry, candidate_reason, same_as_preferred))
            audit_records.append(_candidate_audit_record(path, meta, "rejected", candidate_reason))
            continue

        acceptable_entries.append((entry, candidate_reason, same_as_preferred))
        audit_records.append(_candidate_audit_record(path, meta, "not_selected", "better candidate selected"))

    if acceptable_entries:
        entry, _candidate_reason, selected_is_preferred = min(
            acceptable_entries,
            key=lambda item: _entry_ranking_key(item[0], preferred_injection),
        )
        if selected_is_preferred:
            selection_reason = f"Preferred {preferred_injection}s injection selected"
        elif preferred_candidates_present:
            selection_reason = f"Preferred {preferred_injection}s failed ({preferred_reason}); selected {entry['selected_injection']} fallback"
        else:
            selection_reason = f"Preferred {preferred_injection}s unavailable; selected {entry['selected_injection']}"

        selected_file = getattr(entry["fsa"], "file_name", "")
        entry["selection_reason"] = selection_reason
        entry["preferred_injection_time"] = preferred_injection
        entry["alternate_injections"] = [
            record for record in audit_records
            if record["file"] != selected_file
        ]
        entry["alternate_injections_summary"] = "; ".join(
            f"{alt['selected_injection']} {alt['file']} ({alt['status']}: {alt['reason']})"
            for alt in entry["alternate_injections"]
        )
        return entry

    if not fallback_entries:
        return None

    entry, best_reason, _selected_is_preferred = min(
        fallback_entries,
        key=lambda item: _fallback_entry_ranking_key(item[0], preferred_injection),
    )
    selected_file = getattr(entry["fsa"], "file_name", "")
    entry["selection_reason"] = f"No candidate passed QC; kept {entry['selected_injection']} ({best_reason})"
    entry["preferred_injection_time"] = preferred_injection
    entry["alternate_injections"] = [record for record in audit_records if record["file"] != selected_file]
    entry["alternate_injections_summary"] = "; ".join(
        f"{alt['selected_injection']} {alt['file']} ({alt['status']}: {alt['reason']})"
        for alt in entry["alternate_injections"]
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
