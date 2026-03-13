"""
FLT3 Analysis — FSA Classification (Skeleton).
"""
from __future__ import annotations
from pathlib import Path
from core.analyses.flt3.config import ASSAY_CONFIG

from Bio import SeqIO
from fraggler.fraggler import print_warning
from core.analyses.flt3.config import ASSAY_CONFIG

def get_injection_metadata(fsa_path: Path) -> dict:
    """Extracts injection time from FSA metadata."""
    try:
        record = SeqIO.read(str(fsa_path), "abi")
        tags = record.annotations.get("abif_raw", {})
        return {
            "injection_time": tags.get("InSc1", 0),
            "injection_voltage": tags.get("InVt1", 0),
        }
    except Exception as e:
        print_warning(f"Could not read metadata for {fsa_path.name}: {e}")
        return {"injection_time": 0, "injection_voltage": 0}

def detect_assay(name: str) -> str:
    """Detects FLT3/NPM1 assay from filename."""
    lower = name.lower()
    if "itd" in lower or "ratio" in lower:
        return "FLT3-ITD"
    if "d835" in lower or "tkd" in lower or "d8365" in lower:
        return "FLT3-D835"
    if "npm1" in lower:
        return "NPM1"
    return "UNKNOWN"

def classify_fsa(fsa_path: Path) -> dict | None:
    """
    Classifies an FSA file for FLT3 analysis.
    Returns a dictionary with all relevant metadata.
    """
    name = fsa_path.name
    assay = detect_assay(name)
    
    if assay not in ASSAY_CONFIG:
        return None
    
    cfg = ASSAY_CONFIG[assay]
    meta = get_injection_metadata(fsa_path)
    
    # Determine group and specific type
    lower = name.lower()
    group = "sample"
    if "ntc" in lower:
        group = "negative_control"
    elif "ivs-0000" in lower:
        group = "reactive_control"
    elif "ivs-p001" in lower:
        group = "positive_control"
        
    analysis_type = "standard"
    if "10x" in lower or "1-10" in lower:
        analysis_type = "10x_diluted"
    elif "25x" in lower or "1-25" in lower:
        analysis_type = "25x_diluted"
    elif "ratio" in lower:
        analysis_type = "ratio_quant"
    elif "ufort" in lower:
        analysis_type = "undiluted"
    elif "tkd" in lower or "kutting" in lower:
        analysis_type = "TKD_digested"

    parallel = None
    if "p1" in lower:
        parallel = "p1"
    elif "p2" in lower:
        parallel = "p2"

    return {
        "assay": assay,
        "group": group,
        "analysis_type": analysis_type,
        "parallel": parallel,
        "ladder": "ROX",
        "trace_channels": cfg["trace_channels"],
        "peak_channels": cfg["peak_channels"],
        "primary_peak_channel": cfg["peak_channels"][0],
        "bp_min": cfg["bp_min"],
        "bp_max": cfg["bp_max"],
        "wt_bp": cfg.get("wt_bp"),
        "mut_bp": cfg.get("mut_bp"),
        "injection_time": meta["injection_time"],
        "injection_voltage": meta["injection_voltage"],
    }
