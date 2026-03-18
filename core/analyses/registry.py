"""
Analysis Registry — Manages multiple analysis types.
"""
from __future__ import annotations
import importlib
from typing import Any

from config import APP_SETTINGS

def get_active_analysis_name() -> str:
    """Returns the name of the active analysis (clonality, flt3, etc)."""
    from config import APP_SETTINGS
    return APP_SETTINGS.get("active_analysis", "clonality")

def get_analysis_module(submodule: str) -> Any:
    """
    Returns the module for the active analysis and submodule (config, classification, pipeline).
    Example: get_analysis_module("config") -> core.analyses.clonality.config
    """
    name = get_active_analysis_name()
    module_path = f"core.analyses.{name}.{submodule}"
    try:
        return importlib.import_module(module_path)
    except ModuleNotFoundError as exc:
        # Fall back only when the missing module is the requested analysis package/module,
        # not when an inner dependency import failed from inside that module.
        missing_name = exc.name or ""
        analysis_package = f"core.analyses.{name}"
        if missing_name not in {analysis_package, module_path}:
            raise
        fallback_path = f"core.analyses.clonality.{submodule}"
        return importlib.import_module(fallback_path)

def get_assay_config() -> dict:
    return get_analysis_module("config").ASSAY_CONFIG

def get_assay_display_order() -> list:
    return get_analysis_module("config").ASSAY_DISPLAY_ORDER
