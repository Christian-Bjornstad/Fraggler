"""
HemaFrag Diagnostics — Assay Configuration & Constants.

Dispatcher for analysis-specific configurations.
"""
from __future__ import annotations
from pathlib import Path
import datetime
import importlib

# ==================================================================
# ============================= OPTIONS ============================
# ==================================================================

from dataclasses import dataclass

@dataclass
class PlotOptions:
    y_min: int = 100
    show_peaks: bool = True
    show_ladder_peaks: bool = True
    ladder_output_mode: str = "in_report"
    split_html_tabs: bool = True


# ============================================================
# ======================= SHARED CONFIG ======================
# ============================================================

DEFAULT_FSA_DIR = Path.cwd() / "data"

DEFAULT_LIZ_LADDER = "LIZ500_250"
DEFAULT_ROX_LADDER = "ROX400HD"

DEFAULT_REFERENCE_SHADE_COLOR = "#ded7a6"
DEFAULT_REFERENCE_SHADE_ALPHA = 0.24

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

# --------------------- Klonalitetsregler (Keep as defaults) ---------------------
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
# DELEGATION to active analysis
# ============================================================
from core.analyses.registry import get_analysis_module
from typing import Any

def _get_analysis_attr(attr: str, default: Any = None) -> Any:
    mod = get_analysis_module("config")
    return getattr(mod, attr, default)

def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    value = str(color or DEFAULT_REFERENCE_SHADE_COLOR).strip()
    if value.startswith("#"):
        value = value[1:]
    if len(value) == 3:
        value = "".join(ch * 2 for ch in value)
    if len(value) != 6:
        value = DEFAULT_REFERENCE_SHADE_COLOR[1:]
    try:
        return int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16)
    except ValueError:
        fallback = DEFAULT_REFERENCE_SHADE_COLOR[1:]
        return int(fallback[0:2], 16), int(fallback[2:4], 16), int(fallback[4:6], 16)

def reference_shade_rgba(alpha: float | None = None) -> str:
    """Return the shared Plotly fill color for assay reference windows."""
    color = _get_analysis_attr("REFERENCE_SHADE_COLOR", DEFAULT_REFERENCE_SHADE_COLOR)
    r, g, b = _hex_to_rgb(str(color))
    a = DEFAULT_REFERENCE_SHADE_ALPHA if alpha is None else float(alpha)
    a = max(0.0, min(1.0, a))
    return f"rgba({r},{g},{b},{a:.2f})"


def merged_analysis_attr(attr: str) -> Any:
    """Return an analysis config mapping merged across known analysis modules."""
    merged: dict = {}
    for module_path in (
        "core.analyses.general.config",
        "core.analyses.clonality.config",
        "core.analyses.flt3.config",
    ):
        try:
            value = getattr(importlib.import_module(module_path), attr, None)
        except ModuleNotFoundError:
            value = None
        if isinstance(value, dict):
            merged.update(value)

    active_value = _get_analysis_attr(attr, None)
    if isinstance(active_value, dict):
        merged.update(active_value)
    return merged

def __getattr__(name: str) -> Any:
    """Module-level getattr (Python 3.7+) for dynamic delegation."""
    if name == "ASSAY_DISPLAY_ORDER":
        return _get_analysis_attr("ASSAY_DISPLAY_ORDER", [])
    if name == "ASSAY_CONFIG":
        return _get_analysis_attr("ASSAY_CONFIG", {})
    if name == "ASSAY_REFERENCE_RANGES":
        return _get_analysis_attr("ASSAY_REFERENCE_RANGES", {})
    if name == "ASSAY_REFERENCE_LABEL":
        return _get_analysis_attr("ASSAY_REFERENCE_LABEL", {})
    if name == "NONSPECIFIC_PEAKS":
        return _get_analysis_attr("NONSPECIFIC_PEAKS", {})
    if name == "ASSAY_REARRANGEMENT_INFO":
        return _get_analysis_attr("ASSAY_REARRANGEMENT_INFO", {})
    if name == "CHANNEL_TEXT_COLORS":
        return _get_analysis_attr("CHANNEL_TEXT_COLORS", {})
    if name == "REFERENCE_SHADE_COLOR":
        return _get_analysis_attr("REFERENCE_SHADE_COLOR", DEFAULT_REFERENCE_SHADE_COLOR)
    
    if name == "ROX_LADDER":
        return _get_analysis_attr("ROX_LADDER", DEFAULT_ROX_LADDER)
    if name == "LIZ_LADDER":
        return _get_analysis_attr("LIZ_LADDER", DEFAULT_LIZ_LADDER)
    
    raise AttributeError(f"module {__name__} has no attribute {name}")

# --------------------------------------------------
# Farger per kanal (brukes i både Matplotlib & Plotly)
# --------------------------------------------------
from core.utils import CHANNEL_COLORS, DEFAULT_TRACE_COLOR
 
# --------------------------------------------------
# Output directory names
# --------------------------------------------------
OUTDIR_NAME_TEMPLATE = "reports_{date}"

def get_default_outdir_name() -> str:
    """Returns a formatted output directory name using the current date."""
    now_str = datetime.datetime.now().strftime("%Y-%m-%d")
    return OUTDIR_NAME_TEMPLATE.format(date=now_str)

OUTDIR_NAME = get_default_outdir_name()
