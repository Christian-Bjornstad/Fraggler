"""General analysis FSA classification."""
from __future__ import annotations

from pathlib import Path

from core.utils import strip_stage_prefix

from .config import GENERAL_ASSAY_NAME, resolve_runtime_config


def classify_fsa(fsa_path: Path) -> dict | None:
    """Return a generic metadata payload for arbitrary .fsa inputs."""
    if not fsa_path.name.lower().endswith(".fsa"):
        return None

    runtime = resolve_runtime_config()
    clean_name = strip_stage_prefix(fsa_path.name)
    sample_id = fsa_path.stem
    return {
        "analysis": "general",
        "assay": GENERAL_ASSAY_NAME,
        "group": "sample",
        "ladder": runtime["ladder"],
        "trace_channels": runtime["trace_channels"],
        "peak_channels": runtime["peak_channels"],
        "primary_peak_channel": runtime["primary_peak_channel"],
        "sample_channel": runtime["sample_channel"],
        "bp_min": runtime["bp_min"],
        "bp_max": runtime["bp_max"],
        "file_path": fsa_path,
        "source_run_dir": fsa_path.parent.name,
        "sample_id": sample_id,
        "clean_name": clean_name,
    }
