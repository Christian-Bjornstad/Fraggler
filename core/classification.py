"""
Fraggler Diagnostics — FSA Classification.

Dispatcher for analysis-specific classification logic.
"""
from __future__ import annotations
from pathlib import Path

from core.analyses.registry import get_analysis_module

def detect_assay(name: str) -> str:
    """Delegates assay detection to the active analysis module."""
    mod = get_analysis_module("classification")
    return mod.detect_assay(name)

def classify_fsa(fsa_path: Path) -> tuple[str, str, str, list[str], list[str], str, float, float] | None:
    """Delegates FSA classification to the active analysis module."""
    mod = get_analysis_module("classification")
    return mod.classify_fsa(fsa_path)
