"""
FLT3 Analysis Configuration (Skeleton).
"""
from __future__ import annotations

ASSAY_CONFIG = {
    "FLT3-ITD": {
        "dye": "ROX",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 50.0,
        "bp_max": 1000.0,
        "wt_bp": 330.0,
    },
    "FLT3-D835": {
        "dye": "ROX",
        "trace_channels": ["DATA3"],
        "peak_channels": ["DATA3"],
        "bp_min": 50.0,
        "bp_max": 200.0,
        "wt_bp": 77.0,
        "mut_bp": 126.0,
    },
    "NPM1": {
        "dye": "ROX",
        "trace_channels": ["DATA3"],
        "peak_channels": ["DATA3"],
        "bp_min": 50.0,
        "bp_max": 1000.0,
        "wt_bp": 300.0,
        "mut_bp": 304.0,
    },
}

ROX_LADDER = "GS500ROX"

ASSAY_DISPLAY_ORDER = ["FLT3-ITD", "FLT3-D835", "NPM1"]
NONSPECIFIC_PEAKS = {}
ASSAY_REFERENCE_RANGES = {
    "FLT3-ITD": [(320.0, 1000.0)],
    "FLT3-D835": [(70.0, 150.0)],
    "NPM1": [(290.0, 310.0)],
}
ASSAY_REFERENCE_LABEL = {
    "FLT3-ITD": "Villtype: ~330 bp, Mutert: >335 bp",
    "FLT3-D835": "Villtype: ~80 bp, Mutert: ~129 bp",
    "NPM1": "Villtype: ~300 bp, Mutert: ~304 bp",
}

# Injection time preference (seconds)
PREFERRED_INJECTION_TIME = {
    "FLT3-D835": 3,
    "TKD_digested": 3,
    "undiluted": 3,
    "ratio_quant": 1,
    "10x_diluted": 1,
    "25x_diluted": 1,
    "FLT3-ITD": 1,
    "NPM1": 1,
}
