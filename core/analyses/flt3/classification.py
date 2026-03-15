"""FLT3 / NPM1 FSA classification."""
from __future__ import annotations
import re
from pathlib import Path

from Bio import SeqIO
from fraggler.fraggler import print_warning
from core.analyses.flt3.config import ASSAY_CONFIG, PREFERRED_INJECTION_TIME

PARALLEL_RE = re.compile(r"(^|[_-])(p[12])([_-]|$)", re.IGNORECASE)
DIT_RE = re.compile(r"(\d{2}OUM\d{5})", re.IGNORECASE)


def _decode_abi_value(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore").strip()
    if value is None:
        return ""
    return str(value).strip()


def _normalize_well_id(value: str) -> str | None:
    if not value:
        return None
    match = re.fullmatch(r"([A-Ha-h])0*([1-9]|1[0-2])", value.replace("-", "").replace("_", ""))
    if not match:
        return value.upper()
    row, col = match.groups()
    return f"{row.upper()}{int(col):02d}"


def _specimen_id_from_name(name: str) -> str:
    dit_match = DIT_RE.search(name)
    if dit_match:
        return dit_match.group(1).upper()

    upper = name.upper()
    if "IVS-0000" in upper:
        return "IVS-0000"
    if "IVS-P001" in upper:
        return "IVS-P001"
    if "NTC" in upper:
        return "NTC"
    if upper.startswith("V_") or upper.startswith("V__") or upper == "V":
        return "V"
    prefix = name.split("__", 1)[0].strip("_-")
    return prefix or name


def _treatment_family(analysis_type: str) -> str:
    mapping = {
        "10x_diluted": "10x_diluted",
        "25x_diluted": "25x_diluted",
        "ratio_quant": "ratio_quant",
        "undiluted": "undiluted",
        "TKD_digested": "TKD_digested",
        "standard": "standard",
    }
    return mapping.get(analysis_type, analysis_type or "standard")


def _build_selection_key(
    specimen_id: str,
    assay: str,
    analysis_type: str,
    well_id: str | None,
    parallel: str | None,
    stem: str,
) -> str:
    identity = well_id or parallel or stem
    return "::".join(
        [
            specimen_id or "unknown",
            assay or "unknown",
            _treatment_family(analysis_type),
            identity or "unknown",
        ]
    )

def get_injection_metadata(fsa_path: Path) -> dict:
    """Extracts injection time from FSA metadata."""
    try:
        record = SeqIO.read(str(fsa_path), "abi")
        tags = record.annotations.get("abif_raw", {})
        return {
            "injection_time": tags.get("InSc1", 0),
            "injection_voltage": tags.get("InVt1", 0),
            "well_id": _normalize_well_id(_decode_abi_value(tags.get("TUBE1"))),
            "run_name": _decode_abi_value(tags.get("RunN1")),
            "run_date": _decode_abi_value(tags.get("RUND1")),
            "run_time": _decode_abi_value(tags.get("RUNT1")),
            "injection_protocol": _decode_abi_value(tags.get("RPrN1")),
        }
    except Exception as e:
        print_warning(f"Could not read metadata for {fsa_path.name}: {e}")
        return {
            "injection_time": 0,
            "injection_voltage": 0,
            "well_id": None,
            "run_name": "",
            "run_date": "",
            "run_time": "",
            "injection_protocol": "",
        }

def detect_assay(name: str, *, default_to_d835: bool = False) -> str:
    """Detects FLT3/NPM1 assay from filename."""
    lower = name.lower()
    if any(token in lower for token in ("itd", "ratio", "itdr")):
        return "FLT3-ITD"
    if any(token in lower for token in ("d835", "tkd", "d8365", "cutting", "kutting")):
        return "FLT3-D835"
    if any(token in lower for token in ("npm1", "npm-1", "npm_1")):
        return "NPM1"
    if default_to_d835:
        return "FLT3-D835"
    return "UNKNOWN"

def classify_fsa(fsa_path: Path) -> dict | None:
    """
    Classifies an FSA file for FLT3 analysis.
    Returns a dictionary with all relevant metadata.
    """
    name = fsa_path.name
    assay = detect_assay(name, default_to_d835=True)
    
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
    if "10x" in lower or "1-10" in lower or "x10" in lower:
        analysis_type = "10x_diluted"
    elif "25x" in lower or "1-25" in lower or "x25" in lower:
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

    specimen_id = _specimen_id_from_name(name)
    selection_key = _build_selection_key(
        specimen_id=specimen_id,
        assay=assay,
        analysis_type=analysis_type,
        well_id=meta.get("well_id"),
        parallel=parallel,
        stem=fsa_path.stem,
    )

    return {
        "assay": assay,
        "group": group,
        "analysis_type": analysis_type,
        "parallel": parallel,
        "well_id": meta.get("well_id"),
        "specimen_id": specimen_id,
        "selection_key": selection_key,
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
        "injection_protocol": meta.get("injection_protocol", ""),
        "run_name": meta.get("run_name", ""),
        "run_date": meta.get("run_date", ""),
        "run_time": meta.get("run_time", ""),
        "source_run_dir": fsa_path.parent.name,
    }
