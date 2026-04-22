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
MIN_SIZE_STANDARD_HEIGHT_LIZ = 150

MIN_DISTANCE_BETWEEN_PEAKS_ROX = 15
MIN_SIZE_STANDARD_HEIGHT_ROX = 100

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

REFERENCE_SHADE_COLOR = "#ded7a6"

ASSAY_REFERENCE_RANGES: dict[str, list[tuple[float, float]]] = {
    "FR1": [(310.0, 360.0)],
    "FR2": [(250.0, 295.0)],
    "FR3": [(100.0, 170.0)],
    "IKZF1": [(100.0, 300.0)],
    "Ktr-albumin": [(100.0, 300.0)],

    "IGK": [(120.0, 160.0), (190.0, 300.0)],
    "KDE": [(210.0, 390.0)],

    "DHJH_D": [(110.0, 290.0), (390.0, 420.0)],
    "DHJH_E": [(100.0, 130.0)],

    "TCRgA": [(145.0, 255.0)],
    "TCRgB": [(80.0, 220.0)],

    "TCRbA": [(240.0, 285.0)],
    "TCRbB": [(240.0, 285.0)],
    "TCRbC": [(170.0, 210.0), (285.0, 325.0)],
}

# Channel → text color mapping for rearrangement info tables
CHANNEL_TEXT_COLORS: dict[str, str] = {
    "DATA1": "#2563eb",   # blue
    "DATA2": "#16a34a",   # green
    "DATA3": "#1e293b",   # black (orange trace, dark text)
}

# ============================================================
# Rearrangement info per assay (for DIT report tables)
# - "title": assay header line
# - "rows": list of rearrangement entries
#   - "name": rearrangement name
#   - "range": bp range string
#   - "channel": single channel (DATA1/DATA2/DATA3)
#   - "channels": list of channels (for dual-channel assays)
# - "prefix_parts": for IGK-style split-color prefix
# ============================================================
ASSAY_REARRANGEMENT_INFO: dict[str, dict] = {
    "FR1": {
        "title": "FR1 (IgH): 310–360 bp",
        "rows": [
            {"name": "FR1-JH", "range": "310–360", "channel": "DATA1"},
        ],
    },
    "FR2": {
        "title": "FR2 (IgH): 250–295 bp",
        "rows": [
            {"name": "FR2-JH", "range": "250–295", "channel": "DATA1"},
        ],
    },
    "FR3": {
        "title": "FR3 (IgH): 100–170 bp",
        "rows": [
            {"name": "FR3-JH", "range": "100–170", "channel": "DATA2"},
        ],
    },
    "DHJH_D": {
        "title": "DHJH mix D: 110–290, 390–420 bp",
        "rows": [
            {"name": "IGHD1-IGHJ", "range": "260–290", "channel": "DATA2"},
            {"name": "IGHD2-IGHJ", "range": "230–260", "channel": "DATA2"},
            {"name": "IGHD3-IGHJ", "range": "390–420", "channel": "DATA2"},
            {"name": "IGHD4-IGHJ", "range": "175–205", "channel": "DATA2"},
            {"name": "IGHD5-IGHJ", "range": "225–255", "channel": "DATA2"},
            {"name": "IGHD6-IGHJ", "range": "110–150", "channel": "DATA2"},
        ],
    },
    "DHJH_E": {
        "title": "DHJH mix E: 100–130 bp",
        "rows": [
            {"name": "IGHD7-IGHJ", "range": "100–130", "channel": "DATA1"},
        ],
    },
    "IGK": {
        "title": "IgK: 120–160, 190–300 bp",
        "prefix_parts": [("Jk1-4", "DATA2"), ("Jk5", "DATA1")],
        "rows": [
            {"name": "Vk1f/6", "range": "140–160"},
            {"name": "Vk2f", "range": "280–300"},
            {"name": "Vk3f", "range": "190–210"},
            {"name": "Vk4", "range": "275–295"},
            {"name": "Vk5", "range": "260–280"},
            {"name": "Vk7", "range": "120–140"},
        ],
    },
    "KDE": {
        "title": "Kde: 210–390 bp",
        "rows": [
            {"name": "Vk1f/6-Kde", "range": "225–245", "channel": "DATA3"},
            {"name": "Vk2f-Kde", "range": "360–390", "channel": "DATA3"},
            {"name": "Vk3f-Kde", "range": "279–300", "channel": "DATA3"},
            {"name": "Vk4-Kde", "range": "255–385", "channel": "DATA3"},
            {"name": "Vk5-Kde", "range": "350–380", "channel": "DATA3"},
            {"name": "Vk7-Kde", "range": "210–230", "channel": "DATA3"},
            {"name": "IntronRSS-Kde", "range": "270–300", "channel": "DATA3"},
        ],
    },
    "TCRgA": {
        "title": "TCRγ mix A: 145–255 bp",
        "rows": [
            {"name": "Vg1-8 + Jg1,1/2,1", "range": "230–255", "channel": "DATA1"},
            {"name": "Vg1-8 + Jg1,3/2,3", "range": "195–230", "channel": "DATA2"},
            {"name": "Vg10 + Jg1,1/2,1", "range": "175–195", "channel": "DATA1"},
            {"name": "Vg10 + Jg1,3/2,3", "range": "145–175", "channel": "DATA2"},
        ],
    },
    "TCRgB": {
        "title": "TCRγ mix B: 80–220 bp",
        "rows": [
            {"name": "Vg9 + Jg1,1/2,1", "range": "195–220", "channel": "DATA1"},
            {"name": "Vg9 + Jg1,3/2,3", "range": "160–195", "channel": "DATA2"},
            {"name": "Vg11 + Jg1,1/2,1", "range": "110–140", "channel": "DATA1"},
            {"name": "Vg11 + Jg1,3/2,3", "range": "80–110", "channel": "DATA2"},
        ],
    },
    "TCRbA": {
        "title": "TCRβ mix A: 240–285 bp",
        "prefix_parts": [("Jβ1.X", "DATA2"), ("Jβ2.X", "DATA1")],
        "rows": [
            {"name": "Vβ+Jβ1/2", "range": "240–285"},
        ],
    },
    "TCRbB": {
        "title": "TCRβ mix B: 240–285 bp",
        "prefix_parts": [("Jβ1.X", "DATA2"), ("Jβ2.X", "DATA1")],
        "rows": [
            {"name": "Vβ+Jβ2", "range": "240–285"},
        ],
    },
    "TCRbC": {
        "title": "TCRβ mix C: 170–210, 285–325 bp",
        "prefix_parts": [("Jβ1.X", "DATA2"), ("Jβ2.X", "DATA1")],
        "rows": [
            {"name": "Dβ1+Jβ1/2", "range": "285–325"},
            {"name": "Dβ2 + Jβ2", "range": "170–210"},
        ],
    },
}

# Keep backward-compat simple label for assays without rearrangement info
ASSAY_REFERENCE_LABEL: dict[str, str] = {
    "IKZF1": "IKZF1 (IKAROS): 100–300 bp",
    "Ktr-albumin": "Ktr-albumin kontroll: 100–300 bp",
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
