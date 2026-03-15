"""FLT3 / NPM1 FSA classification."""
from __future__ import annotations
import re
from pathlib import Path

from Bio import SeqIO
from fraggler.fraggler import print_warning
from core.analyses.flt3.config import ASSAY_CONFIG, PREFERRED_INJECTION_TIME

PARALLEL_RE = re.compile(r"(^|[_-])(p[12])([_-]|$)", re.IGNORECASE)

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
    if any(token in lower for token in ("itd", "ratio", "itdr")):
        return "FLT3-ITD"
    if any(token in lower for token in ("d835", "tkd", "d8365", "cutting", "kutting")):
        return "FLT3-D835"
    if any(token in lower for token in ("npm1", "npm-1", "npm_1")):
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

    protocol_injection_time = PREFERRED_INJECTION_TIME.get(analysis_type)
    if protocol_injection_time is None:
        protocol_injection_time = PREFERRED_INJECTION_TIME.get(assay, meta["injection_time"])

    parallel = None
    match = PARALLEL_RE.search(lower)
    if match:
        parallel = match.group(2).lower()

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
        "protocol_injection_time": protocol_injection_time,
    }
