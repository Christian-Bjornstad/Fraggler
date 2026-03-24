"""
Clonality Analysis — FSA Classification.
"""
from __future__ import annotations
import re
from pathlib import Path

from fraggler.fraggler import print_warning
from core.analyses.clonality.config import ASSAY_CONFIG
from core.utils import strip_stage_prefix


def _assay_pattern(*roots: str, suffix: str) -> re.Pattern[str]:
    joined_roots = "|".join(re.escape(root) for root in roots)
    return re.compile(
        rf"(?<![a-z0-9])(?:{joined_roots})(?:[\s_-]*mix)?[\s_-]*{re.escape(suffix)}(?![a-z0-9])",
        re.IGNORECASE,
    )


TCRG_PATTERNS = {
    "TCRgA": _assay_pattern("tcrg", "trg", "tcrgamma", "trgamma", suffix="a"),
    "TCRgB": _assay_pattern("tcrg", "trg", "tcrgamma", "trgamma", suffix="b"),
}

TCRB_PATTERNS = {
    "TCRbA": _assay_pattern("tcrb", "trb", "tcrbeta", "trbeta", suffix="a"),
    "TCRbB": _assay_pattern("tcrb", "trb", "tcrbeta", "trbeta", suffix="b"),
    "TCRbC": _assay_pattern("tcrb", "trb", "tcrbeta", "trbeta", suffix="c"),
}

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

    # 2) IKZF1 / IKAROS and albumin control
    if (
        "ktralbumin" in lower.replace("_", "").replace("-", "").replace(" ", "")
        or re.search(r"(?<![a-z0-9])albumin(?![a-z0-9])", lower)
        or re.search(r"(?<![a-z0-9])ktr(?![a-z0-9])", lower)
    ):
        return "Ktr-albumin"
    if "ikzf1" in lower.replace("_", "").replace("-", "").replace(" ", "") or re.search(
        r"(?<![a-z0-9])ikz(?:f1)?(?:[\s_-]*mix[\s_-]*\d+)?(?![a-z0-9])",
        lower,
    ):
        return "IKZF1"

    # 3) TCRγ (LIZ) - try specific patterns first
    for assay, pattern in TCRG_PATTERNS.items():
        if pattern.search(lower):
            return assay

    # 4) TCRβ (ROX)
    for assay, pattern in TCRB_PATTERNS.items():
        if pattern.search(lower):
            return assay

    # 5) IgH-regionene
    if "fr1" in lower: return "FR1"
    if "fr2" in lower: return "FR2"
    if "fr3" in lower: return "FR3"

    # 6) DHJH-mikser
    if any(m in lower for m in ["dhjh_d", "dhjhd", "dhjh_mixd", "dhjh_mix_d", "dhjhmixd"]):
        return "DHJH_D"
    if any(m in lower for m in ["dhjh_e", "dhjhe", "dhjh_mixe", "dhjh_mix_e", "dhjhmixe"]):
        return "DHJH_E"

    # 7) LIZ IgK / KDE
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
