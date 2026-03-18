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
from fraggler.fraggler import FsaFile


def _analysis_config_module(analysis_id: str):
    try:
        return importlib.import_module(f"core.analyses.{analysis_id}.config")
    except ImportError:
        return importlib.import_module("core.analyses.clonality.config")


def _resolve_exact_ladder_name(analysis_id: str, ladder_name: str) -> str:
    config_mod = _analysis_config_module(analysis_id)
    upper = ladder_name.upper()
    if upper == "LIZ":
        return getattr(config_mod, "LIZ_LADDER", LIZ_LADDER)
    if upper == "ROX":
        return getattr(config_mod, "ROX_LADDER", ROX_LADDER)
    return ladder_name


def _resolve_ladder_runtime(analysis_id: str, ladder_name: str) -> tuple[str, float, float]:
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
