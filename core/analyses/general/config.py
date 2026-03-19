"""General analysis configuration and runtime helpers."""
from __future__ import annotations

from config import APP_SETTINGS
from core.assay_config import DEFAULT_LIZ_LADDER, DEFAULT_ROX_LADDER

GENERAL_ASSAY_NAME = "GENERAL"
LIZ_LADDER = DEFAULT_LIZ_LADDER
ROX_LADDER = DEFAULT_ROX_LADDER
ALLOWED_LADDERS = ("LIZ500_250", "ROX400HD", "GS500ROX")
ALLOWED_TRACE_CHANNELS = ("DATA1", "DATA2", "DATA3")
DEFAULT_TRACE_CHANNELS = ("DATA1",)
DEFAULT_BP_MIN = 50.0
DEFAULT_BP_MAX = 1000.0

ASSAY_CONFIG = {
    GENERAL_ASSAY_NAME: {
        "dye": "ROX",
        "trace_channels": list(ALLOWED_TRACE_CHANNELS),
        "peak_channels": list(ALLOWED_TRACE_CHANNELS),
        "bp_min": DEFAULT_BP_MIN,
        "bp_max": DEFAULT_BP_MAX,
    }
}
ASSAY_DISPLAY_ORDER = [GENERAL_ASSAY_NAME]
ASSAY_REFERENCE_RANGES: dict[str, list[tuple[float, float]]] = {}
ASSAY_REFERENCE_LABEL: dict[str, str] = {}
NONSPECIFIC_PEAKS: dict[str, list[float]] = {}
REFERENCE_SHADE_COLOR = "#ebe8cb"


def _general_profile(settings: dict | None = None) -> dict:
    settings = settings or APP_SETTINGS
    analyses = settings.get("analyses", {})
    profile = analyses.get("general", {})
    return profile if isinstance(profile, dict) else {}


def get_general_pipeline_settings(settings: dict | None = None) -> dict:
    profile = _general_profile(settings)
    pipeline = profile.get("pipeline", {})
    return pipeline if isinstance(pipeline, dict) else {}


def normalize_ladder_name(ladder_name: str | None) -> str:
    if not ladder_name:
        return ROX_LADDER
    cleaned = str(ladder_name).strip().upper().replace("-", "").replace(" ", "")
    mapping = {
        "LIZ500": LIZ_LADDER,
        "LIZ500_250": LIZ_LADDER,
        "ROX400HD": ROX_LADDER,
        "GS500ROX": "GS500ROX",
    }
    return mapping.get(cleaned, ROX_LADDER)


def normalize_trace_channels(trace_channels: object | None) -> list[str]:
    if trace_channels is None:
        return list(DEFAULT_TRACE_CHANNELS)
    if isinstance(trace_channels, str):
        trace_channels = [trace_channels]

    cleaned: list[str] = []
    for channel in trace_channels:
        value = str(channel).strip().upper()
        if value in ALLOWED_TRACE_CHANNELS and value not in cleaned:
            cleaned.append(value)
    return cleaned or list(DEFAULT_TRACE_CHANNELS)


def choose_primary_channel(trace_channels: list[str], preferred: str | None = None) -> str:
    preferred = (preferred or "").strip().upper()
    if preferred and preferred in trace_channels:
        return preferred
    return trace_channels[0] if trace_channels else DEFAULT_TRACE_CHANNELS[0]


def resolve_runtime_config(settings: dict | None = None) -> dict:
    pipeline = get_general_pipeline_settings(settings)
    trace_channels = normalize_trace_channels(pipeline.get("trace_channels"))
    primary_channel = choose_primary_channel(
        trace_channels,
        pipeline.get("primary_peak_channel") or pipeline.get("sample_channel"),
    )
    ladder = normalize_ladder_name(pipeline.get("ladder"))
    return {
        "ladder": ladder,
        "trace_channels": trace_channels,
        "peak_channels": list(trace_channels),
        "primary_peak_channel": primary_channel,
        "sample_channel": primary_channel,
        "bp_min": float(pipeline.get("bp_min", DEFAULT_BP_MIN) or DEFAULT_BP_MIN),
        "bp_max": float(pipeline.get("bp_max", DEFAULT_BP_MAX) or DEFAULT_BP_MAX),
    }
