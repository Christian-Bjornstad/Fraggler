"""
Fraggler Diagnostics — FSA Classification.

``detect_assay`` normalises file names to canonical assay keys.
``classify_fsa`` returns the full classification tuple for a single FSA.
"""
from __future__ import annotations

import re
from pathlib import Path

from fraggler.fraggler import print_warning

from core.assay_config import ASSAY_CONFIG


# ==================================================================
# ===================== KLASSIFISERING =============================
# ==================================================================

def detect_assay(name: str) -> str:
    """
    Returnerer assay-navn slik at det matcher nøkkelen i ASSAY_CONFIG.
    Støtter:
    - Size Ladder-synonymer
    - TRB/TRG-synonymer
    - tcrb_a/tcrg_b osv.
    - MIX-varianter: mixa, mix_a, mixA, mix-a (TCRbeta/TCRb/TRB og TCRg/TRG)
    """
    if not name:
        return "UNKNOWN"

    s = name.strip()
    lower = s.lower()

    # 1) SIZE LADDER synonymer
    if (
        "_sl_" in lower
        or lower.endswith("_sl.fsa")
        or "sizeladder" in lower
        or "size_ladder" in lower
        or "size-ladder" in lower
        or "sladder" in lower
    ):
        return "SL"

    # 2) MIX-deteksjon (robust)
    mix_pat_beta = re.compile(
        r"(tcrbeta|tcrb|trb)[\s\-_]*mix[\s\-_]*([abc])",
        re.IGNORECASE,
    )
    mix_pat_gamma = re.compile(
        r"(tcrg|trg)[\s\-_]*mix[\s\-_]*([ab])",
        re.IGNORECASE,
    )

    m_beta = mix_pat_beta.search(s)
    if m_beta:
        mix_letter = m_beta.group(2).lower()
        if mix_letter == "a":
            return "TCRbA"
        if mix_letter == "b":
            return "TCRbB"
        if mix_letter == "c":
            return "TCRbC"

    m_gamma = mix_pat_gamma.search(s)
    if m_gamma:
        mix_letter = m_gamma.group(2).lower()
        if mix_letter == "a":
            return "TCRgA"
        if mix_letter == "b":
            return "TCRgB"

    # 3) TCRγ (LIZ) synonymer
    if (
        "tcrga" in lower or "tcrg_a" in lower or "tcrg-a" in lower
        or "trga" in lower or "trg_a" in lower or "trg-a" in lower
    ):
        return "TCRgA"

    if (
        "tcrgb" in lower or "tcrg_b" in lower or "tcrg-b" in lower
        or "trgb" in lower or "trg_b" in lower or "trg-b" in lower
    ):
        return "TCRgB"

    # 4) TCRβ (ROX) synonymer
    if (
        "tcrba" in lower or "tcrb_a" in lower or "tcrb-a" in lower
        or "trba" in lower or "trb_a" in lower or "trb-a" in lower
    ):
        return "TCRbA"

    if (
        "tcrbb" in lower or "tcrb_b" in lower or "tcrb-b" in lower
        or "trbb" in lower or "trb_b" in lower or "trb-b" in lower
    ):
        return "TCRbB"

    if (
        "tcrbc" in lower or "tcrb_c" in lower or "tcrb-c" in lower
        or "trbc" in lower or "trb_c" in lower or "trb-c" in lower
    ):
        return "TCRbC"

    # 5) IgH-regionene
    if "fr1" in lower:
        return "FR1"
    if "fr2" in lower:
        return "FR2"
    if "fr3" in lower:
        return "FR3"

    # 6) DHJH-mikser
    if "dhjh_d" in lower or "dhjhd" in lower or "dhjh_mixd" in lower or "dhjh_mix_d" in lower:
        return "DHJH_D"
    if "dhjh_e" in lower or "dhjhe" in lower or "dhjh_mixe" in lower or "dhjh_mix_e" in lower:
        return "DHJH_E"

    # 7) LIZ IgK / KDE
    if "igk" in lower:
        return "IGK"
    if "kde" in lower:
        return "KDE"

    # 8) Fallback
    return "UNKNOWN"


def strip_stage_prefix(name: str) -> str:
    return re.sub(r"^\d{5}_[a-f0-9]{8}_", "", name, flags=re.IGNORECASE)

def classify_fsa(fsa_path: Path):
    """
    Returnerer:
      assay (str),
      group (positive/negative/reactive/unknown),
      ladder ('LIZ'/'ROX'),
      trace_channels (list[str]),
      peak_channels (list[str]),
      primary_peak_channel (str),
      bp_min, bp_max
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

    # Gruppe (PK/NK/RK) – kun label
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
    if not trace_channels:
        print_warning(f"[CLASSIFY] {name}: trace_channels tom – hopper.")
        return None
    if not peak_channels:
        print_warning(f"[CLASSIFY] {name}: peak_channels tom – hopper.")
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
