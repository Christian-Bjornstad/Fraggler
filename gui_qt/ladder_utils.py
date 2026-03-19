from __future__ import annotations

import importlib
from pathlib import Path

from core.analyses.clonality.classification import classify_fsa as classify_clonality_fsa
from core.analyses.flt3.classification import classify_fsa as classify_flt3_fsa
from core.analysis import analyse_fsa_liz, analyse_fsa_rox
from core.assay_config import (
    LIZ_LADDER,
    MIN_DISTANCE_BETWEEN_PEAKS_LIZ,
    MIN_DISTANCE_BETWEEN_PEAKS_ROX,
    MIN_SIZE_STANDARD_HEIGHT_LIZ,
    MIN_SIZE_STANDARD_HEIGHT_ROX,
    ROX_LADDER,
)
from config import (
    GENERAL_DEFAULT_BP_MAX,
    GENERAL_DEFAULT_BP_MIN,
    GENERAL_DEFAULT_LADDER,
    GENERAL_DEFAULT_PRIMARY_CHANNEL,
    GENERAL_DEFAULT_TRACE_CHANNELS,
    GENERAL_LADDERS,
    GENERAL_TRACE_CHANNELS,
    get_analysis_settings,
)
from fraggler.fraggler import FsaFile

GENERAL_LADDER_RUNTIME = {
    "LIZ500_250": (
        float(MIN_DISTANCE_BETWEEN_PEAKS_LIZ),
        float(MIN_SIZE_STANDARD_HEIGHT_LIZ),
    ),
    "ROX400HD": (
        float(MIN_DISTANCE_BETWEEN_PEAKS_ROX),
        float(MIN_SIZE_STANDARD_HEIGHT_ROX),
    ),
    "GS500ROX": (
        float(MIN_DISTANCE_BETWEEN_PEAKS_ROX),
        float(MIN_SIZE_STANDARD_HEIGHT_ROX),
    ),
}


def _analysis_config_module(analysis_id: str):
    if analysis_id == "general":
        return None
    try:
        return importlib.import_module(f"core.analyses.{analysis_id}.config")
    except ImportError:
        return importlib.import_module("core.analyses.clonality.config")


def _resolve_exact_ladder_name(analysis_id: str, ladder_name: str) -> str:
    if analysis_id == "general":
        upper = ladder_name.upper()
        if upper == "LIZ500":
            return "LIZ500_250"
        if upper == "LIZ":
            return "LIZ500_250"
        if upper == "ROX":
            return GENERAL_DEFAULT_LADDER
        if upper in GENERAL_LADDERS:
            return upper
        return ladder_name
    config_mod = _analysis_config_module(analysis_id)
    upper = ladder_name.upper()
    if upper == "LIZ":
        return getattr(config_mod, "LIZ_LADDER", LIZ_LADDER)
    if upper == "ROX":
        return getattr(config_mod, "ROX_LADDER", ROX_LADDER)
    return ladder_name


def _resolve_ladder_runtime(analysis_id: str, ladder_name: str) -> tuple[str, float, float]:
    if analysis_id == "general":
        exact_ladder = _resolve_exact_ladder_name(analysis_id, ladder_name)
        min_distance, min_height = GENERAL_LADDER_RUNTIME.get(
            exact_ladder,
            GENERAL_LADDER_RUNTIME[GENERAL_DEFAULT_LADDER],
        )
        return exact_ladder, float(min_distance), float(min_height)
    config_mod = _analysis_config_module(analysis_id)
    exact_ladder = _resolve_exact_ladder_name(analysis_id, ladder_name)
    if exact_ladder.upper().startswith("LIZ"):
        return (
            exact_ladder,
            float(getattr(config_mod, "MIN_DISTANCE_BETWEEN_PEAKS_LIZ", MIN_DISTANCE_BETWEEN_PEAKS_LIZ)),
            float(getattr(config_mod, "MIN_SIZE_STANDARD_HEIGHT_LIZ", MIN_SIZE_STANDARD_HEIGHT_LIZ)),
        )
    return (
        exact_ladder,
        float(getattr(config_mod, "MIN_DISTANCE_BETWEEN_PEAKS_ROX", MIN_DISTANCE_BETWEEN_PEAKS_ROX)),
        float(getattr(config_mod, "MIN_SIZE_STANDARD_HEIGHT_ROX", MIN_SIZE_STANDARD_HEIGHT_ROX)),
    )


def _normalize_clonality_classification(fsa_path: Path, classified: tuple) -> dict:
    assay, group, ladder, trace_channels, peak_channels, primary_peak_channel, bp_min, bp_max = classified
    exact_ladder = _resolve_exact_ladder_name("clonality", ladder)
    return {
        "analysis": "clonality",
        "assay": assay,
        "group": group,
        "ladder": exact_ladder,
        "trace_channels": trace_channels,
        "peak_channels": peak_channels,
        "primary_peak_channel": primary_peak_channel,
        "bp_min": float(bp_min),
        "bp_max": float(bp_max),
        "sample_channel": trace_channels[0] if trace_channels else None,
        "raw": classified,
        "file_path": fsa_path,
    }


def _normalize_flt3_classification(fsa_path: Path, classified: dict) -> dict:
    exact_ladder = _resolve_exact_ladder_name("flt3", classified["ladder"])
    return {
        "analysis": "flt3",
        "assay": classified["assay"],
        "group": classified["group"],
        "ladder": exact_ladder,
        "trace_channels": classified["trace_channels"],
        "peak_channels": classified["peak_channels"],
        "primary_peak_channel": classified["primary_peak_channel"],
        "bp_min": float(classified["bp_min"]),
        "bp_max": float(classified["bp_max"]),
        "sample_channel": classified["trace_channels"][0] if classified.get("trace_channels") else None,
        "raw": classified,
        "file_path": fsa_path,
    }


def _general_runtime_metadata(fsa_path: Path) -> dict:
    profile = get_analysis_settings("general")
    pipeline = profile.get("pipeline", {})
    ladder = _resolve_exact_ladder_name("general", pipeline.get("ladder", GENERAL_DEFAULT_LADDER))

    trace_channels = pipeline.get("trace_channels", list(GENERAL_DEFAULT_TRACE_CHANNELS))
    if not isinstance(trace_channels, list):
        trace_channels = list(GENERAL_DEFAULT_TRACE_CHANNELS)
    trace_channels = [ch for ch in trace_channels if ch in GENERAL_TRACE_CHANNELS]
    if not trace_channels:
        trace_channels = list(GENERAL_DEFAULT_TRACE_CHANNELS)

    primary_peak_channel = str(pipeline.get("primary_peak_channel", GENERAL_DEFAULT_PRIMARY_CHANNEL))
    if primary_peak_channel not in trace_channels:
        primary_peak_channel = trace_channels[0]

    try:
        bp_min = float(pipeline.get("bp_min", GENERAL_DEFAULT_BP_MIN))
    except (TypeError, ValueError):
        bp_min = GENERAL_DEFAULT_BP_MIN
    try:
        bp_max = float(pipeline.get("bp_max", GENERAL_DEFAULT_BP_MAX))
    except (TypeError, ValueError):
        bp_max = GENERAL_DEFAULT_BP_MAX

    return {
        "analysis": "general",
        "assay": "General",
        "group": "sample",
        "ladder": ladder,
        "trace_channels": trace_channels,
        "peak_channels": list(trace_channels),
        "primary_peak_channel": primary_peak_channel,
        "sample_channel": primary_peak_channel,
        "bp_min": bp_min,
        "bp_max": bp_max,
        "raw": pipeline,
        "file_path": fsa_path,
    }


def _classify_for_analysis(fsa_path: Path, analysis_id: str) -> dict | None:
    if analysis_id == "clonality":
        classified = classify_clonality_fsa(fsa_path)
        if classified:
            return _normalize_clonality_classification(fsa_path, classified)
        return None
    if analysis_id == "flt3":
        classified = classify_flt3_fsa(fsa_path)
        if classified:
            return _normalize_flt3_classification(fsa_path, classified)
        return None
    if analysis_id == "general":
        return _general_runtime_metadata(fsa_path)
    return None


def detect_fsa_for_ladder(fsa_path: Path, preferred_analysis: str | None = None) -> dict | None:
    if preferred_analysis:
        preferred = _classify_for_analysis(fsa_path, preferred_analysis)
        if preferred:
            return preferred

    for analysis_id in ("clonality", "flt3"):
        if analysis_id == preferred_analysis:
            continue
        detected = _classify_for_analysis(fsa_path, analysis_id)
        if detected:
            return detected

    return None


def load_adjustable_fsa(
    fsa_path: Path,
    preferred_analysis: str | None = None,
    metadata: dict | None = None,
) -> tuple[FsaFile, dict]:
    metadata = metadata or detect_fsa_for_ladder(fsa_path, preferred_analysis=preferred_analysis)
    if not metadata:
        raise ValueError(f"Could not classify {fsa_path.name}")

    sample_channel = metadata.get("sample_channel")
    if not sample_channel:
        raise ValueError(f"No sample channel found for {fsa_path.name}")

    ladder, min_distance, min_height = _resolve_ladder_runtime(metadata["analysis"], metadata["ladder"])
    metadata["ladder"] = ladder

    if ladder.upper().startswith("LIZ"):
        fsa = analyse_fsa_liz(
            fsa_path,
            sample_channel,
            ladder_name=ladder,
            min_distance_between_peaks=min_distance,
            min_size_standard_height=min_height,
        )
    else:
        fsa = analyse_fsa_rox(
            fsa_path,
            sample_channel,
            ladder_name=ladder,
            min_distance_between_peaks=min_distance,
            min_size_standard_height=min_height,
        )

    if fsa is None:
        if ladder.upper().startswith("LIZ"):
            fsa = FsaFile(
                str(fsa_path),
                ladder,
                sample_channel,
                min_distance,
                min_height,
                size_standard_channel="DATA105",
            )
        else:
            fsa = FsaFile(
                str(fsa_path),
                ladder,
                sample_channel,
                min_distance,
                min_height,
                size_standard_channel="DATA4",
            )

    return fsa, metadata
