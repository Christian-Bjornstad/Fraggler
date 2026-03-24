"""
Clonality Analysis Configuration.
"""
from __future__ import annotations
from pathlib import Path

# ============================================================
# ======================= KONFIG =============================
# ============================================================

LIZ_LADDER = "LIZ500_250"
ROX_LADDER = "ROX400HD"

MIN_DISTANCE_BETWEEN_PEAKS_LIZ = 30
MIN_SIZE_STANDARD_HEIGHT_LIZ = 300

MIN_DISTANCE_BETWEEN_PEAKS_ROX = 15
MIN_SIZE_STANDARD_HEIGHT_ROX = 200

# --------------------- Peak-parametre ------------------------
ABSOLUTE_MIN_PEAK_HEIGHT = 400.0
RELATIVE_PEAK_HEIGHT_MIN = 0.40
MAX_PEAKS_PER_CHANNEL = 12
MIN_INTERPEAK_DISTANCE_BP = 3.0

LOCAL_BACKGROUND_WINDOW_BP = 10.0
MIN_PEAK_TO_LOCAL_BACKGROUND_RATIO = 2.5

# --------------------- Klonalitetsregler ---------------------
CLONAL_MAX_LABELLED_PEAKS = 3
CLONAL_CLUSTER_WINDOW_BP = 12.0
CLONAL_DOMINANCE_RATIO = 1.7

# --------------------- Polyklonal sjekk ----------------------
POLY_LOCAL_WINDOW_BP = 12.0
POLY_LOCAL_REL_HEIGHT = 0.40
POLY_LOCAL_MAX_PEAKS = 4

# --------------------- SL-spesifikke terskler ----------------
ABSOLUTE_MIN_PEAK_HEIGHT_SL = 500.0
RELATIVE_PEAK_HEIGHT_MIN_SL = 0.10
MAX_PEAKS_PER_CHANNEL_SL = 10
MIN_PEAK_TO_LOCAL_BACKGROUND_RATIO_SL = 1.5

SL_TARGET_FRAGMENTS_BP = [100.0, 200.0, 300.0, 400.0, 600.0]
SL_WINDOW_BP = 20.0

# ============================================================
# Rekkefølge på assays i pasientrapport
# ============================================================
ASSAY_DISPLAY_ORDER = [
    "FR1", "FR2", "FR3",
    "IKZF1", "Ktr-albumin",
    "IGK", "KDE",
    "DHJH_D", "DHJH_E",
    "TCRbA", "TCRbB", "TCRbC",
    "TCRgA", "TCRgB",
]

# ============================================================
# ========= ASSAY-SPESIFIKK KONFIGURASJON ====================
# ============================================================

ASSAY_CONFIG = {
    "FR1": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 280.0,
        "bp_max": 420.0,
    },
    "FR2": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 200.0,
        "bp_max": 400.0,
    },
    "FR3": {
        "dye": "ROX",
        "trace_channels": ["DATA2"],
        "peak_channels": ["DATA2"],
        "bp_min": 60.0,
        "bp_max": 220.0,
    },
    "IKZF1": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 50.0,
        "bp_max": 400.0,
    },
    "Ktr-albumin": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 50.0,
        "bp_max": 400.0,
    },
    "TCRbA": {
        "dye": "ROX",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 210.0,
        "bp_max": 310.0,
    },
    "TCRbB": {
        "dye": "ROX",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 210.0,
        "bp_max": 310.0,
    },
    "TCRbC": {
        "dye": "ROX",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 140.0,
        "bp_max": 360.0,
    },
    "SL": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 80.0,
        "bp_max": 700.0,
    },
    "DHJH_D": {
        "dye": "ROX",
        "trace_channels": ["DATA2"],
        "peak_channels": ["DATA2"],
        "bp_min": 90.0,
        "bp_max": 440.0,
    },
    "DHJH_E": {
        "dye": "ROX",
        "trace_channels": ["DATA1"],
        "peak_channels": ["DATA1"],
        "bp_min": 65.0,
        "bp_max": 160.0,
    },
    "IGK": {
        "dye": "LIZ",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 90.0,
        "bp_max": 330.0,
    },
    "KDE": {
        "dye": "LIZ",
        "trace_channels": ["DATA3"],
        "peak_channels": ["DATA3"],
        "bp_min": 190.0,
        "bp_max": 410.0,
    },
    "TCRgA": {
        "dye": "LIZ",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 110.0,
        "bp_max": 290.0,
    },
    "TCRgB": {
        "dye": "LIZ",
        "trace_channels": ["DATA1", "DATA2"],
        "peak_channels": ["DATA1", "DATA2"],
        "bp_min": 60.0,
        "bp_max": 250.0,
    },
}

# ============================================================
# ========= ASSAY-SPESIFIKK REFERANSE-SHADING =================
# ============================================================

REFERENCE_SHADE_COLOR = "#ebe8cb"

ASSAY_REFERENCE_RANGES: dict[str, list[tuple[float, float]]] = {
    "FR1": [(310.0, 360.0)],
    "FR2": [(250.0, 295.0)],
    "FR3": [(100.0, 170.0)],
    "IKZF1": [(100.0, 300.0)],
    "Ktr-albumin": [(100.0, 300.0)],

    "IGK": [(120.0, 300.0)],
    "KDE": [(210.0, 390.0)],

    "DHJH_D": [(110.0, 290.0), (390.0, 420.0)],
    "DHJH_E": [(100.0, 130.0)],

    "TCRgA": [(145.0, 255.0)],
    "TCRgB": [(80.0, 220.0)],

    "TCRbA": [(240.0, 285.0)],
    "TCRbB": [(240.0, 285.0)],
    "TCRbC": [(170.0, 210.0), (285.0, 325.0)],
}

ASSAY_REFERENCE_LABEL: dict[str, str] = {
    "FR1": "FR1 (IgH): 310–360 bp (VH–JH)",
    "FR2": "FR2 (IgH): 250–295 bp (VH–JH)",
    "FR3": "FR3 (IgH): 100–170 bp (VH–JH)",
    "IKZF1": "IKZF1 (IKAROS): 100–300 bp",
    "Ktr-albumin": "Ktr-albumin kontroll: 100–300 bp",

    "IGK": "IgK: 120–160, 190–210, 260–300 bp (Vκ–Jκ)",
    "KDE": "Kde: 210–250, 270–300, 350–390 bp (Kde–involveringer)",

    "DHJH_D": "DHJH mix D: 110–290 og 390–420 bp (DH–JH)",
    "DHJH_E": "DHJH mix E: 100–130 bp (DH7–JH)",

    "TCRgA": "TCRγ mix A (Vγ1–8/10): 145–255 bp",
    "TCRgB": "TCRγ mix B (Vγ9/11): 80–220 bp",

    "TCRbA": "TCRβ mix A: 240–285 bp (Vβ–Jβ1/2)",
    "TCRbB": "TCRβ mix B: 240–285 bp (Vβ–Jβ2)",
    "TCRbC": "TCRβ mix C: 170–210 og 285–325 bp (Dβ–Jβ)",
}

# --------------------------------------------------
# Non-specific peaks (bp) per assay
# --------------------------------------------------
NONSPECIFIC_PEAKS: dict[str, list[float]] = {
    "FR1": [60, 85, 98, 203, 566],
    "FR2": [199, 226, 228, 800],
    "FR3": [211, 213, 286],
    "DHJH_D": [76, 94, 96, 158, 161, 176, 179, 196, 200, 202, 345, 350, 421, 459, 501, 678, 694, 707, 748, 753, 796],
    "DHJH_E": [53, 79, 93, 123, 161, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 211, 390, 415, 416, 419, 476, 599, 602, 718, 783, 1031, 1404, 1804, 2420],
    "IGK": [217],  # ~217
    "KDE": [401, 403, 404],
    "TCRbA": [213, 273],  # ~213, ~273
    "TCRbB": [93, 126, 127, 150, 221],  # ~93, ~126, 127, 150, ~221
    "TCRbC": [128, 123],  # ~128, ~123
    # TCRg: NO non-specific peaks
}
