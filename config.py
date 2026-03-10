"""
Fraggler Diagnostics — Configuration & Settings

Centralized settings persistence (YAML), defaults, and path helpers.
Compatible with Python 3.10+.
"""
from __future__ import annotations

import yaml
from pathlib import Path
from typing import Any, Dict

# ============================================================
# PATHS
# ============================================================

SETTINGS_PATH = Path.home() / ".fraggler_gui.yaml"

# ============================================================
# DEFAULTS
# ============================================================

DEFAULT_SETTINGS: Dict[str, Any] = {
    "theme": "default",  # "default" | "dark"
    "active_analysis": "clonality",
    "pipeline": {
        "input_dir": str(Path.home()),
        "output_base": str(Path.home()),
        "out_folder_name": "ASSAY_REPORTS",
        "mode": "all",               # "all" | "controls" | "custom"
        "assay_filter_substring": "",
    },
    "qc": {
        "input_dir": str(Path.home()),
        "output_base": "",
        "outfile_html": "QC_REPORT.html",
        "excel_name": "QC_TRENDS.xlsx",
        "min_r2_ok": 0.995,
        "min_r2_warn": 0.990,
        "max_mse_ok": 2.0,
        "max_mse_warn": 5.0,
        "nk_ymax_floor": 250.0,
        "w_sample": 3.0,
        "w_ladder": 3.0,
    },
    "batch": {
        "base_input_dir": str(Path.home()),
        "job_source": "subfolders",   # "subfolders" | "yaml"
        "yaml_path": "",
        "job_type": "pipeline",       # "pipeline" | "qc" | "dit"
        "output_base": str(Path.home()),
        "out_folder_tmpl": "ASSAY_REPORTS",
        "outfile_html_tmpl": "QC_REPORT_{name}.html",
        "excel_name_tmpl": "QC_TRENDS_{name}.xlsx",
        "mode": "all",
        "assay_filter_substring": "",
        "aggregate_by_patient": True,
        "patient_id_regex": r"\d{2}OUM\d{5}",
        "aggregate_dit_reports": True,
    },
}


# ============================================================
# LOAD / SAVE
# ============================================================

def _deep_update(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge *override* into *base*."""
    result = base.copy()
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = _deep_update(result[k], v)
        else:
            result[k] = v
    return result

import os

def _validate_settings(settings: Dict[str, Any]) -> None:
    """Basic validation for critical settings."""
    pipeline = settings.get("pipeline", {})
    if not isinstance(pipeline.get("out_folder_name"), str):
        pipeline["out_folder_name"] = "ASSAY_REPORTS"
    
    # Ensure min_r2 is within 0-1 range
    qc = settings.get("qc", {})
    if not (0 <= qc.get("min_r2_ok", 0.995) <= 1):
        qc["min_r2_ok"] = 0.995

def load_settings() -> Dict[str, Any]:
    """Load settings from YAML, merged over defaults, then env vars."""
    settings = DEFAULT_SETTINGS.copy()
    
    # 1) Load from YAML if exists
    if SETTINGS_PATH.exists():
        try:
            with open(SETTINGS_PATH, "r", encoding="utf-8") as f:
                user = yaml.safe_load(f) or {}
            settings = _deep_update(settings, user)
        except Exception:
            pass
            
    # 2) Override with Environment Variables (e.g. FRAGGLER_THEME=dark)
    # Mapping: FRAGGLER_CATEGORY_KEY
    for key, value in os.environ.items():
        if key.startswith("FRAGGLER_"):
            parts = key.split("_")[1:] # E.g. ["PIPELINE", "INPUT", "DIR"]
            if len(parts) == 1:
                settings[parts[0].lower()] = value
            elif len(parts) == 2:
                cat, k = parts[0].lower(), parts[1].lower()
                if cat in settings and isinstance(settings[cat], dict):
                    settings[cat][k] = value

    _validate_settings(settings)
    return settings


def save_settings(settings: Dict[str, Any]) -> None:
    """Persist settings to YAML."""
    try:
        # Avoid saving ephemeral env overrides back to file if desired, 
        # but here we just save the current active state.
        with open(SETTINGS_PATH, "w", encoding="utf-8") as f:
            yaml.safe_dump(settings, f, sort_keys=False, allow_unicode=True)
    except Exception:
        pass


# Singleton — imported by other modules
APP_SETTINGS = load_settings()
