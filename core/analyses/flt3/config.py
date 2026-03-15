"""FLT3 / NPM1 analysis configuration."""
from __future__ import annotations

ASSAY_CONFIG = {
    "FLT3-ITD": {
        "dye": "ROX",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 50.0,
        "bp_max": 1000.0,
        "wt_bp": 330.0,
        "itd_min_bp": 335.0,
        "positive_ratio": 0.02,
        "control_wt_min_area": 10000.0,
    },
    "FLT3-D835": {
        "dye": "ROX",
        "trace_channels": ["DATA3"],
        "peak_channels": ["DATA3"],
        "bp_min": 50.0,
        "bp_max": 250.0,
        "wt_bp": 80.0,
        "mut_bp": 129.0,
        "wt_range": (76.0, 83.5),
        "mut_ranges": [(121.0, 130.5)],
        "peak_height_min": 50.0,
        "peak_distance": 8,
        "positive_ratio": 0.05,
        "control_wt_min_area": 1000.0,
    },
    "NPM1": {
        "dye": "ROX",
        "trace_channels": ["DATA3"],
        "peak_channels": ["DATA3"],
        "bp_min": 50.0,
        "bp_max": 1000.0,
        "wt_bp": 300.0,
        "mut_bp": 304.0,
        "positive_ratio": 0.01,
    },
}

ROX_LADDER = "GS500ROX"
BP_CORRECTION_OFFSETS = {
    "FLT3-ITD": 0.0,
    "FLT3-D835": 0.0,
    "NPM1": 0.0,
}

ASSAY_DISPLAY_ORDER = ["FLT3-ITD", "FLT3-D835", "NPM1"]
NONSPECIFIC_PEAKS = {}
ASSAY_REFERENCE_RANGES = {
    "FLT3-ITD": [(300.0, 1000.0)],
    "FLT3-D835": [(50.0, 250.0)],
    "NPM1": [(299.0, 301.0), (303.0, 305.0)],
}
ASSAY_REFERENCE_LABEL = {
    "FLT3-ITD": "Analysevindu: 300-1000 bp. Villtype forventet rundt 330 bp, mutert >335 bp.",
    "FLT3-D835": "Analysevindu: 50-250 bp. Villtype: 80 bp, Mutert >129 bp.",
    "NPM1": "Villtype: 299/300-301 bp, Mutert: 303-305 bp",
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
    "NPM1": 3,
}
