"""
Fraggler Diagnostics — Shared Utilities.
"""
from __future__ import annotations
import re

# Centralized Regex for Control identification
CONTROL_PREFIX_RE = re.compile(r"^(PK1|PK2|PK|NK|RK|DIT|KTR)_", re.IGNORECASE)
WATER_RE = re.compile(r"^(v|water|h2o)([_\-]|(?=\.fsa)|$)", re.IGNORECASE)

def strip_stage_prefix(name: str) -> str:
    """Removes the 5-digit prefix and 8-character hash from filenames."""
    return re.sub(r"^\d{5}_[a-f0-9]{8}_", "", name, flags=re.IGNORECASE)

def is_water_file(filename: str) -> bool:
    """Returns True if the filename looks like a water/negative control."""
    clean_name = strip_stage_prefix(filename)
    return bool(WATER_RE.match(clean_name))

def is_control_file(filename: str) -> bool:
    """Returns True if the filename starts with a known control prefix (after stripping stage)."""
    clean_name = strip_stage_prefix(filename)
    return bool(CONTROL_PREFIX_RE.match(clean_name))

# Shared Color configuration (can be moved here if it helps centralize)
CHANNEL_COLORS = {
    "DATA1": "#3b82f6",     # Vibrant Blue
    "DATA2": "#10b981",     # Vibrant Emerald Green
    "DATA3": "#334155",     # Slate/Black
    "DATA4": "#ef4444",     # Vibrant Red/Rose
    "DATA105": "#f97316",   # Vibrant Orange
}
DEFAULT_TRACE_COLOR = "#3b82f6"
