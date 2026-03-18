"""
Fraggler Diagnostics — Configuration & Settings

Centralized settings persistence (YAML), defaults, and path helpers.
Compatible with Python 3.10+.
"""
from __future__ import annotations

import yaml
import copy
import logging
import os
import warnings
from pathlib import Path
from typing import Any, Dict, Mapping

# ============================================================
# PATHS
# ============================================================

SETTINGS_PATH = Path.home() / ".fraggler_gui.yaml"
_LOG = logging.getLogger(__name__)

LAST_SETTINGS_LOAD_ERROR: str | None = None
LAST_SETTINGS_SAVE_ERROR: str | None = None

# ============================================================
# DEFAULTS
# ============================================================

DEFAULT_SETTINGS: Dict[str, Any] = {
    "theme": "default",  # "default" | "dark"
    "active_analysis": "clonality",
    "general": {
        "author": "OUS",
        "default_output": "",
    },
    # Keep the legacy top-level key for backward compatibility.
    "default_output": "",
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
    "analyses": {
        "clonality": {
            "batch": {
                "base_input_dir": str(Path.home()),
                "output_base": str(Path.home()),
                "aggregate_by_patient": True,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": True,
            },
            "pipeline": {
                "mode": "all",
                "assay_filter_substring": "",
            },
        },
        "flt3": {
            "batch": {
                "base_input_dir": str(Path.home()),
                "output_base": str(Path.home()),
                "aggregate_by_patient": True,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": True,
            },
            "pipeline": {
                "mode": "all",
                "assay_filter_substring": "",
            },
        },
    },
}


# ============================================================
# LOAD / SAVE
# ============================================================

def _deep_update(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge *override* into *base*."""
    result = copy.deepcopy(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = _deep_update(result[k], v)
        else:
            result[k] = v
    return result

def _coerce_env_value(value: str) -> Any:
    lower = value.strip().lower()
    if lower in {"true", "false"}:
        return lower == "true"
    try:
        if "." in value:
            return float(value)
        return int(value)
    except ValueError:
        return value


def _apply_env_overrides(settings: Dict[str, Any], env: Mapping[str, str]) -> None:
    """Apply FRAGGLER_* environment overrides to nested settings."""
    def _assign(target: Dict[str, Any], schema: Any, parts: list[str], raw_value: str) -> None:
        if not parts:
            return

        if not isinstance(schema, dict):
            target["_".join(parts)] = _coerce_env_value(raw_value)
            return

        for width in range(len(parts), 0, -1):
            key = "_".join(parts[:width])
            if key not in schema:
                continue

            if width == len(parts) or not isinstance(schema.get(key), dict):
                target[key] = _coerce_env_value(raw_value)
            else:
                if key not in target or not isinstance(target.get(key), dict):
                    target[key] = {}
                _assign(target[key], schema[key], parts[width:], raw_value)
            return

        target["_".join(parts)] = _coerce_env_value(raw_value)

    for key, value in env.items():
        if not key.startswith("FRAGGLER_"):
            continue

        parts = [part.lower() for part in key.split("_")[1:] if part]
        if not parts:
            continue

        _assign(settings, DEFAULT_SETTINGS, parts, value)


def _migrate_legacy_settings(settings: Dict[str, Any]) -> Dict[str, Any]:
    """Fold legacy keys into the canonical settings structure."""
    general = settings.setdefault("general", {})
    batch = settings.setdefault("batch", {})
    pipeline = settings.setdefault("pipeline", {})
    analyses = settings.setdefault("analyses", {})

    legacy_default_output = settings.get("default_output", "")
    if legacy_default_output and not general.get("default_output"):
        general["default_output"] = legacy_default_output

    default_output = general.get("default_output", "") or legacy_default_output
    if default_output:
        if not batch.get("output_base") or batch.get("output_base") == DEFAULT_SETTINGS["batch"]["output_base"]:
            batch["output_base"] = default_output
        if not pipeline.get("output_base") or pipeline.get("output_base") == DEFAULT_SETTINGS["pipeline"]["output_base"]:
            pipeline["output_base"] = default_output

    for analysis_id, analysis_defaults in DEFAULT_SETTINGS.get("analyses", {}).items():
        profile = analyses.setdefault(analysis_id, {})
        profile_batch = profile.setdefault("batch", {})
        profile_pipeline = profile.setdefault("pipeline", {})

        default_batch = analysis_defaults.get("batch", {})
        default_pipeline = analysis_defaults.get("pipeline", {})

        if not profile_batch.get("base_input_dir") or profile_batch.get("base_input_dir") == default_batch.get("base_input_dir"):
            profile_batch["base_input_dir"] = batch.get("base_input_dir", default_batch.get("base_input_dir", str(Path.home())))
        if not profile_batch.get("output_base") or profile_batch.get("output_base") == default_batch.get("output_base"):
            profile_batch["output_base"] = batch.get("output_base", default_batch.get("output_base", str(Path.home())))
        if (
            "aggregate_by_patient" not in profile_batch
            or profile_batch.get("aggregate_by_patient") == default_batch.get("aggregate_by_patient")
        ):
            profile_batch["aggregate_by_patient"] = batch.get("aggregate_by_patient", default_batch.get("aggregate_by_patient", True))
        if (
            not profile_batch.get("patient_id_regex")
            or profile_batch.get("patient_id_regex") == default_batch.get("patient_id_regex")
        ):
            profile_batch["patient_id_regex"] = batch.get("patient_id_regex", default_batch.get("patient_id_regex", r"\d{2}OUM\d{5}"))
        if (
            "aggregate_dit_reports" not in profile_batch
            or profile_batch.get("aggregate_dit_reports") == default_batch.get("aggregate_dit_reports")
        ):
            profile_batch["aggregate_dit_reports"] = batch.get("aggregate_dit_reports", default_batch.get("aggregate_dit_reports", True))

        if not profile_pipeline.get("mode") or profile_pipeline.get("mode") == default_pipeline.get("mode"):
            profile_pipeline["mode"] = pipeline.get("mode", default_pipeline.get("mode", "all"))
        if not profile_pipeline.get("assay_filter_substring"):
            profile_pipeline["assay_filter_substring"] = pipeline.get("assay_filter_substring", default_pipeline.get("assay_filter_substring", ""))

    return settings

def _validate_settings(settings: Dict[str, Any]) -> None:
    """Basic validation for critical settings."""
    pipeline = settings.get("pipeline", {})
    if not isinstance(pipeline.get("out_folder_name"), str):
        pipeline["out_folder_name"] = "ASSAY_REPORTS"
    for key in ("input_dir", "output_base", "assay_filter_substring"):
        if key in pipeline and not isinstance(pipeline.get(key), str):
            pipeline[key] = str(pipeline.get(key, ""))
    if pipeline.get("mode") not in {"all", "controls", "custom"}:
        pipeline["mode"] = "all"
    
    # Ensure min_r2 is within 0-1 range
    qc = settings.get("qc", {})
    if not (0 <= qc.get("min_r2_ok", 0.995) <= 1):
        qc["min_r2_ok"] = 0.995
    if not (0 <= qc.get("min_r2_warn", 0.990) <= 1):
        qc["min_r2_warn"] = 0.990

    batch = settings.get("batch", {})
    if not isinstance(batch.get("base_input_dir"), str):
        batch["base_input_dir"] = str(Path.home())
    if not isinstance(batch.get("output_base"), str):
        batch["output_base"] = str(Path.home())
    if not isinstance(batch.get("patient_id_regex"), str):
        batch["patient_id_regex"] = r"\d{2}OUM\d{5}"

    general = settings.get("general", {})
    if not isinstance(general.get("author", "OUS"), str):
        general["author"] = "OUS"
    if not isinstance(general.get("default_output", ""), str):
        general["default_output"] = ""

    if settings.get("active_analysis") not in DEFAULT_SETTINGS["analyses"]:
        settings["active_analysis"] = DEFAULT_SETTINGS["active_analysis"]

    settings["default_output"] = general.get("default_output", "")

    analyses = settings.setdefault("analyses", {})
    for analysis_id, defaults in DEFAULT_SETTINGS.get("analyses", {}).items():
        profile = analyses.setdefault(analysis_id, {})
        profile_batch = profile.setdefault("batch", {})
        profile_pipeline = profile.setdefault("pipeline", {})

        if not isinstance(profile_batch.get("base_input_dir"), str):
            profile_batch["base_input_dir"] = defaults["batch"]["base_input_dir"]
        if not isinstance(profile_batch.get("output_base"), str):
            profile_batch["output_base"] = defaults["batch"]["output_base"]
        if not isinstance(profile_batch.get("patient_id_regex"), str):
            profile_batch["patient_id_regex"] = defaults["batch"]["patient_id_regex"]
        if not isinstance(profile_batch.get("aggregate_by_patient"), bool):
            profile_batch["aggregate_by_patient"] = defaults["batch"]["aggregate_by_patient"]
        if not isinstance(profile_batch.get("aggregate_dit_reports"), bool):
            profile_batch["aggregate_dit_reports"] = defaults["batch"]["aggregate_dit_reports"]

        if profile_pipeline.get("mode") not in {"all", "controls", "custom"}:
            profile_pipeline["mode"] = defaults["pipeline"]["mode"]
        if not isinstance(profile_pipeline.get("assay_filter_substring"), str):
            profile_pipeline["assay_filter_substring"] = defaults["pipeline"]["assay_filter_substring"]


def get_analysis_settings(analysis_id: str | None = None, settings: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """Return the merged per-analysis settings profile."""
    settings = settings or APP_SETTINGS
    analysis_id = analysis_id or settings.get("active_analysis", "clonality")
    analysis_defaults = DEFAULT_SETTINGS.get("analyses", {}).get(
        analysis_id,
        DEFAULT_SETTINGS["analyses"]["clonality"],
    )
    stored = settings.setdefault("analyses", {}).get(analysis_id, {})
    return _deep_update(analysis_defaults, stored)


def _report_settings_issue(kind: str, settings_path: Path, exc: Exception) -> str:
    message = f"{kind} settings at {settings_path}: {exc}"
    warnings.warn(message, RuntimeWarning, stacklevel=2)
    _LOG.warning(message)
    return message

def load_settings(
    settings_path: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> Dict[str, Any]:
    """Load settings from YAML, merged over defaults, then env vars."""
    global LAST_SETTINGS_LOAD_ERROR
    settings = copy.deepcopy(DEFAULT_SETTINGS)
    settings_path = settings_path or SETTINGS_PATH
    env = env or os.environ
    LAST_SETTINGS_LOAD_ERROR = None
    
    # 1) Load from YAML if exists
    if settings_path.exists():
        try:
            with open(settings_path, "r", encoding="utf-8") as f:
                user = yaml.safe_load(f) or {}
            settings = _deep_update(settings, user)
        except Exception as exc:
            LAST_SETTINGS_LOAD_ERROR = _report_settings_issue("Failed to load", settings_path, exc)

    _apply_env_overrides(settings, env)
    settings = _migrate_legacy_settings(settings)
    _validate_settings(settings)
    return settings


def save_settings(settings: Dict[str, Any], settings_path: Path | None = None) -> None:
    """Persist settings to YAML."""
    global LAST_SETTINGS_SAVE_ERROR
    settings_path = settings_path or SETTINGS_PATH
    LAST_SETTINGS_SAVE_ERROR = None
    try:
        payload = copy.deepcopy(settings)
        payload = _migrate_legacy_settings(payload)
        payload["default_output"] = payload.get("general", {}).get("default_output", "")
        settings_path.parent.mkdir(parents=True, exist_ok=True)
        with open(settings_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)
    except Exception as exc:
        LAST_SETTINGS_SAVE_ERROR = _report_settings_issue("Failed to save", settings_path, exc)


# Singleton — imported by other modules
APP_SETTINGS = load_settings()
