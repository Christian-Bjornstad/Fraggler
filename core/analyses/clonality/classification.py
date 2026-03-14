"""
Clonality Analysis — FSA Classification.
"""
from __future__ import annotations
import re
from pathlib import Path

from fraggler.fraggler import print_warning
from core.analyses.clonality.config import ASSAY_CONFIG
from core.utils import strip_stage_prefix

def detect_assay(name: str) -> str:
    """
    Returnerer assay-navn slik at det matcher nøkkelen i ASSAY_CONFIG.
    """
    if not name:
        return "UNKNOWN"

    s = name.strip()
    lower = s.lower()

    # 1) SIZE LADDER synonymer
    if any(m in lower for m in ["_sl_", "sizeladder", "size_ladder", "size-ladder", "sladder"]) or lower.endswith("_sl.fsa"):
        return "SL"

    # 2) TCRγ (LIZ) - try specific patterns first
    if any(m in lower for m in ["tcrga", "tcrg_a", "tcrg-a", "trga", "trg_a", "trg-a", "tcrg_mix_a", "tcrgmixa"]):
        return "TCRgA"
    if any(m in lower for m in ["tcrgb", "tcrg_b", "tcrg-b", "trgb", "trg_b", "trg-b", "tcrg_mix_b", "tcrgmixb"]):
        return "TCRgB"

    # 3) TCRβ (ROX)
    if any(m in lower for m in ["tcrba", "tcrb_a", "tcrb-a", "trba", "trb_a", "trb-a", "tcrb_mix_a", "tcrbmixa"]):
        return "TCRbA"
    if any(m in lower for m in ["tcrbb", "tcrb_b", "tcrb-b", "trbb", "trb_b", "trb-b", "tcrb_mix_b", "tcrbmixb"]):
        return "TCRbB"
    if any(m in lower for m in ["tcrbc", "tcrb_c", "tcrb-c", "trbc", "trb_c", "trb-c", "tcrb_mix_c", "tcrbmixc"]):
        return "TCRbC"

    # 4) IgH-regionene
    if "fr1" in lower: return "FR1"
    if "fr2" in lower: return "FR2"
    if "fr3" in lower: return "FR3"

    # 5) DHJH-mikser
    if any(m in lower for m in ["dhjh_d", "dhjhd", "dhjh_mixd", "dhjh_mix_d", "dhjhmixd"]):
        return "DHJH_D"
    if any(m in lower for m in ["dhjh_e", "dhjhe", "dhjh_mixe", "dhjh_mix_e", "dhjhmixe"]):
        return "DHJH_E"

    # 6) LIZ IgK / KDE
    if "igk" in lower: return "IGK"
    if "kde" in lower: return "KDE"

    return "UNKNOWN"

def classify_fsa(fsa_path: Path) -> tuple[str, str, str, list[str], list[str], str, float, float] | None:
    """
    Returnerer klassifisering for Clonality.
    """
    name = fsa_path.name
    clean_name = strip_stage_prefix(name)
    assay = detect_assay(clean_name)

    if assay not in ASSAY_CONFIG:
        print_warning(
            f"[CLASSIFY] {name}: assay '{assay}' ikke i ASSAY_CONFIG – hopper over."
        )
        return None

    cfg = ASSAY_CONFIG[assay]

    parts = clean_name.split("_")
    prefix = parts[0].lower() if parts else ""
    if prefix.startswith("pk"):
        group = "positive"
    elif prefix.startswith("nk"):
        group = "negative"
    elif prefix.startswith("rk"):
        group = "reactive"
    else:
        group = "unknown"

    dye = cfg["dye"]
    ladder = "LIZ" if dye.upper() == "LIZ" else "ROX"

    trace_channels = cfg["trace_channels"]
    peak_channels = cfg["peak_channels"]
    if not trace_channels or not peak_channels:
        return None

    primary_peak_channel = peak_channels[0]
    bp_min = float(cfg["bp_min"])
    bp_max = float(cfg["bp_max"])

    return (
        assay,
        group,
        ladder,
        trace_channels,
        peak_channels,
        primary_peak_channel,
        bp_min,
        bp_max,
    )
