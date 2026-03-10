"""
FLT3 Analysis Configuration (Skeleton).
"""
from __future__ import annotations

ASSAY_CONFIG = {
    "FLT3-ITD": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 300.0,
        "bp_max": 600.0,
    }
}

ASSAY_DISPLAY_ORDER = ["FLT3-ITD"]
NONSPECIFIC_PEAKS = {}
ASSAY_REFERENCE_RANGES = {}
ASSAY_REFERENCE_LABEL = {}
