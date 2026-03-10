"""
FLT3 Analysis — FSA Classification (Skeleton).
"""
from __future__ import annotations
from pathlib import Path
from core.analyses.flt3.config import ASSAY_CONFIG

def detect_assay(name: str) -> str:
    if "flt3" in name.lower():
        return "FLT3-ITD"
    return "UNKNOWN"

def classify_fsa(fsa_path: Path) -> tuple[str, str, str, list[str], list[str], str, float, float] | None:
    assay = detect_assay(fsa_path.name)
    if assay not in ASSAY_CONFIG:
        return None
    
    cfg = ASSAY_CONFIG[assay]
    return (
        assay,
        "unknown",
        "ROX",
        cfg["trace_channels"],
        cfg["peak_channels"],
        cfg["peak_channels"][0],
        cfg["bp_min"],
        cfg["bp_max"],
    )
