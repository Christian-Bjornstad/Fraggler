"""
Backward-compatible facade — imports from core.* submodules.

This file is kept for backward compatibility so that existing code
that does ``import fraggler_master_assay_channels`` or
``from fraggler_master_assay_channels import run_pipeline`` still works.
"""
from __future__ import annotations

# --------------- assay_config ---------------
from core.assay_config import *  # noqa: F401,F403

# --------------- classification ---------------
from core.classification import detect_assay, classify_fsa  # noqa: F401

# --------------- analysis ---------------
from core.analysis import (  # noqa: F401
    analyse_fsa_liz,
    analyse_fsa_rox,
    auto_detect_sl_peaks,
    compute_ladder_qc_metrics,
    compute_sl_area_metrics,
    _find_local_maxima,
    estimate_running_baseline,
)

# --------------- plotting (matplotlib) ---------------
from core.plotting_mpl import (  # noqa: F401
    compute_zoom_ymax,
    draw_multi_channel_zoom_on_ax,
)

# --------------- plotting (plotly interactive) ---------------
from core.plotting_plotly import (  # noqa: F401
    compute_group_ymax,
    compute_group_ymax_all_channels,
    compute_group_ymax_for_entries,
    build_interactive_peak_plot_for_entry,
    build_interactive_assay_batch_plot_html,
)

# --------------- html reports ---------------
from core.html_reports import (  # noqa: F401
    DIT_PATTERN,
    extract_dit_from_name,
    dit_to_year,
    build_dit_html_reports,
    interpret_sl_quality,
)

# --------------- plotly offline ---------------
from core.plotly_offline import local_plotly_tag  # noqa: F401

# --------------- pipeline ---------------
from core.pipeline import run_pipeline  # noqa: F401

from core.assay_config import DEFAULT_FSA_DIR


def main():
    """CLI-inngangspunkt: bruker DEFAULT_FSA_DIR."""
    run_pipeline(DEFAULT_FSA_DIR)


if __name__ == "__main__":
    main()
