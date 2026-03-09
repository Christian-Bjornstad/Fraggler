"""
Core package — Fraggler Diagnostics engine.

Re-exports key symbols from submodules for convenient access.
"""
from core.assay_config import ASSAY_CONFIG, ASSAY_DISPLAY_ORDER  # noqa: F401
from core.classification import detect_assay, classify_fsa  # noqa: F401
from core.pipeline import run_pipeline  # noqa: F401
from core.plotly_offline import local_plotly_tag, plotly_inline_script_tag  # noqa: F401
