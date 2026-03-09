"""
Fraggler QC — Rules, aliases, and control regex.
"""
from __future__ import annotations

import re
from dataclasses import dataclass


# --------------------------------------------------------------
ASSAY_ALIASES_QC = {
    # ---- Size Ladder (SL) synonymer ----
    "SIZELADDER": "SL",
    "SIZE_LADDER": "SL",
    "SIZE-LADDER": "SL",
    "SLADDER": "SL",

    # ---- TCRγ (LIZ) synonymer ----
    "TCRGA": "TCRgA",
    "TCRG_A": "TCRgA",
    "TCRG-A": "TCRgA",
    "TRGA": "TCRgA",
    "TRG_A": "TCRgA",
    "TRG-A": "TCRgA",

    "TCRGB": "TCRgB",
    "TCRG_B": "TCRgB",
    "TCRG-B": "TCRgB",
    "TRGB": "TCRgB",
    "TRG_B": "TCRgB",
    "TRG-B": "TCRgB",

    # ---- TCRβ (ROX) synonymer ----
    "TCRBA": "TCRbA",
    "TCRB_A": "TCRbA",
    "TCRB-A": "TCRbA",
    "TRBA": "TCRbA",
    "TRB_A": "TCRbA",
    "TRB-A": "TCRbA",

    "TCRBB": "TCRbB",
    "TCRB_B": "TCRbB",
    "TCRB-B": "TCRbB",
    "TRBB": "TCRbB",
    "TRB_B": "TCRbB",
    "TRB-B": "TCRbB",

    "TCRBC": "TCRbC",
    "TCRB_C": "TCRbC",
    "TCRB-C": "TCRbC",
    "TRBC": "TCRbC",
    "TRB_C": "TCRbC",
    "TRB-C": "TCRbC",
}

def normalize_assay_qc(name: str) -> str:
    if not name:
        return name
    up = name.upper()
    return ASSAY_ALIASES_QC.get(up, name)

# ----------------------------
# QC-regler (justerbare)
# ----------------------------
@dataclass
class QCRules:
    # Ladder-fit
    min_r2_ok: float = 0.999
    min_r2_warn: float = 0.995


    # NK: auto-scale, men ymax skal ikke bli mindre enn dette
    nk_ymax_floor: float = 250.


    # Peak search vinduer (bp)
    sample_peak_window_bp: float = 2.0
    ladder_peak_window_bp: float = 2.0

    # SL signal
    min_sl_total_area: float = 1e4


CONTROL_PREFIX_RE = re.compile(r"^(PK1|PK2|PK|NK|RK)_", re.IGNORECASE)

RUN_CODE_RE = re.compile(r"_([A-Za-z0-9]{6,})\\.fsa$", re.IGNORECASE)

from pathlib import Path
