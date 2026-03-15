from __future__ import annotations

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


def _normalize_clonality_classification(fsa_path: Path, classified: tuple) -> dict:
    assay, group, ladder, trace_channels, peak_channels, primary_peak_channel, bp_min, bp_max = classified
    return {
        "analysis": "clonality",
        "assay": assay,
        "group": group,
        "ladder": ladder,
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
    return {
        "analysis": "flt3",
        "assay": classified["assay"],
        "group": classified["group"],
        "ladder": classified["ladder"],
        "trace_channels": classified["trace_channels"],
        "peak_channels": classified["peak_channels"],
        "primary_peak_channel": classified["primary_peak_channel"],
        "bp_min": float(classified["bp_min"]),
        "bp_max": float(classified["bp_max"]),
        "sample_channel": classified["trace_channels"][0] if classified.get("trace_channels") else None,
        "raw": classified,
        "file_path": fsa_path,
    }


def detect_fsa_for_ladder(fsa_path: Path, preferred_analysis: str | None = None) -> dict | None:
    normalized: list[dict] = []

    clonality = classify_clonality_fsa(fsa_path)
    if clonality:
        normalized.append(_normalize_clonality_classification(fsa_path, clonality))

    flt3 = classify_flt3_fsa(fsa_path)
    if flt3:
        normalized.append(_normalize_flt3_classification(fsa_path, flt3))

    if not normalized:
        return None

    if preferred_analysis:
        for item in normalized:
            if item["analysis"] == preferred_analysis:
                return item

    return normalized[0]


def load_adjustable_fsa(fsa_path: Path, preferred_analysis: str | None = None) -> tuple[FsaFile, dict]:
    metadata = detect_fsa_for_ladder(fsa_path, preferred_analysis=preferred_analysis)
    if not metadata:
        raise ValueError(f"Could not classify {fsa_path.name}")

    sample_channel = metadata.get("sample_channel")
    if not sample_channel:
        raise ValueError(f"No sample channel found for {fsa_path.name}")

    ladder = metadata["ladder"]
    if ladder == "LIZ":
        fsa = analyse_fsa_liz(fsa_path, sample_channel)
    else:
        fsa = analyse_fsa_rox(fsa_path, sample_channel)

    if fsa is None:
        if ladder == "LIZ":
            fsa = FsaFile(
                str(fsa_path),
                LIZ_LADDER,
                sample_channel,
                MIN_DISTANCE_BETWEEN_PEAKS_LIZ,
                MIN_SIZE_STANDARD_HEIGHT_LIZ,
                size_standard_channel="DATA105",
            )
        else:
            fsa = FsaFile(
                str(fsa_path),
                ROX_LADDER,
                sample_channel,
                MIN_DISTANCE_BETWEEN_PEAKS_ROX,
                MIN_SIZE_STANDARD_HEIGHT_ROX,
                size_standard_channel="DATA4",
            )

    return fsa, metadata
