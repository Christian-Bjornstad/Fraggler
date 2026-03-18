"""
Fraggler Diagnostics — Main Pipeline.

``run_pipeline`` processes all .fsa files in a directory, classifies them,
fits ladders, detects peaks, builds DIT reports, and orchestrates the full
analysis flow.
"""
from __future__ import annotations

import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd

from fraggler.fraggler import print_green, print_warning

from core.analyses.clonality.config import (
    ASSAY_CONFIG,
    SL_TARGET_FRAGMENTS_BP,
    SL_WINDOW_BP,
)
from core.analyses.clonality.classification import classify_fsa
from core.analysis import (
    analyse_fsa_liz,
    analyse_fsa_rox,
    auto_detect_sl_peaks,
    compute_ladder_qc_metrics,
    compute_sl_area_metrics,
)
from core.plotting_plotly import (
    compute_group_ymax_for_entries,
    build_interactive_assay_batch_plot_html,
)
from core.plotting_mpl import compute_zoom_ymax
from core.html_reports import (
    extract_dit_from_name,
)
from core.analyses.shared_pipeline import (
    finalize_pipeline_run,
    normalize_pipeline_paths,
    scan_fsa_files,
)
from core.utils import strip_stage_prefix


def _scan_files(fsa_dir: Path, mode: str = "all") -> list[Path]:
    """Scans for .fsa files, filtering out water files and optionally non-controls."""
    fsa_files = scan_fsa_files(fsa_dir, mode=mode)
    if fsa_files:
        print_green(f"Fant {len(fsa_files)} .fsa-filer: {[p.name for p in fsa_files]}")
    return fsa_files


def _analyze_files(fsa_files: list[Path]) -> tuple[list[dict], int]:
    """Performs analysis (ladder fitting, peak detection) on a list of FSA files."""
    entries = []
    skipped = 0

    for fsa_path in fsa_files:
        classified = classify_fsa(fsa_path)
        if classified is None:
            skipped += 1
            continue

        (
            assay,
            group,
            ladder,
            trace_channels,
            peak_channels,
            primary_peak_channel,
            bp_min,
            bp_max,
        ) = classified

        sample_channel = trace_channels[0]

        if ladder == "LIZ":
            fsa = analyse_fsa_liz(fsa_path, sample_channel)
        else:
            fsa = analyse_fsa_rox(fsa_path, sample_channel)

        if fsa is None:
            skipped += 1
            continue

        peaks_by_channel: dict[str, pd.DataFrame | None] = {}
        if assay == "SL":
            peaks_by_channel = auto_detect_sl_peaks(
                fsa,
                peak_channels=peak_channels,
                targets_bp=SL_TARGET_FRAGMENTS_BP,
                window_bp=SL_WINDOW_BP,
                min_height=800.0,
            )
        else:
            for ch in peak_channels:
                peaks_by_channel[ch] = pd.DataFrame(columns=["basepairs", "peaks", "keep"])

        ymax = compute_zoom_ymax(fsa, bp_min, bp_max, trace_channels, assay_name=assay)

        # Ladder QC
        expected_ladder_steps = list(
            map(float, getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])))
        )
        fitted_ladder_steps = list(map(float, getattr(fsa, "ladder_steps", [])))
        ladder_fit_strategy = str(getattr(fsa, "ladder_fit_strategy", "auto_full"))
        ladder_missing_expected_steps = list(
            map(float, getattr(fsa, "ladder_missing_expected_steps", []))
        )
        ladder_fit_note = str(
            getattr(
                fsa,
                "ladder_fit_note",
                "All expected ladder steps were fitted." if not ladder_missing_expected_steps else "Manual ladder review recommended.",
            )
        )
        ladder_review_required = bool(
            getattr(fsa, "ladder_review_required", bool(ladder_missing_expected_steps))
        )
        ladder_qc_status = "ok"
        ladder_r2, n_ladder_steps, n_size_standard_peaks = np.nan, np.nan, np.nan
        try:
            metrics = compute_ladder_qc_metrics(fsa)
            ladder_r2 = metrics["r2"]
            n_ladder_steps = metrics["n_ladder_steps"]
            n_size_standard_peaks = metrics["n_size_standard_peaks"]
            if ladder_fit_strategy == "manual_adjustment":
                ladder_qc_status = "manual_adjustment"
            elif ladder_review_required:
                ladder_qc_status = "review_required"
        except Exception as ex:
            print_warning(f"[LADDER_QC] Klarte ikke beregne QC for {fsa.file_name}: {ex}")
            ladder_qc_status = "ladder_qc_failed"

        # SL-area
        sl_metrics = None
        if assay == "SL":
            try:
                sl_metrics = compute_sl_area_metrics(
                    fsa,
                    trace_channel=primary_peak_channel,
                    targets_bp=SL_TARGET_FRAGMENTS_BP,
                    window_bp=SL_WINDOW_BP,
                )
            except Exception as ex:
                print_warning(f"[SL] Klarte ikke beregne SL-area for {fsa.file_name}: {ex}")

        entries.append({
            "fsa": fsa,
            "peaks_by_channel": peaks_by_channel,
            "trace_channels": trace_channels,
            "primary_peak_channel": primary_peak_channel,
            "ymax": ymax,
            "assay": assay,
            "group": group,
            "ladder": ladder,
            "bp_min": bp_min,
            "bp_max": bp_max,
            "dit": extract_dit_from_name(fsa.file_name),
            "ladder_qc_status": ladder_qc_status,
            "ladder_r2": ladder_r2,
            "n_ladder_steps": n_ladder_steps,
            "n_size_standard_peaks": n_size_standard_peaks,
            "ladder_fit_strategy": ladder_fit_strategy,
            "ladder_missing_expected_steps": ladder_missing_expected_steps,
            "ladder_fit_note": ladder_fit_note,
            "ladder_review_required": ladder_review_required,
            "ladder_expected_step_count": len(expected_ladder_steps),
            "ladder_fitted_step_count": len(fitted_ladder_steps),
            "sl_metrics": sl_metrics,
        })

    print_green(f"[MASTER] Totalt {len(entries)} filer analysert. {skipped} skippet.")
    return entries, skipped


def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
) -> list[dict] | None:

    """
    Kjør full Fraggler-pipeline på alle .fsa-filer i fsa_dir.
    """
    fsa_dir, assay_dir = normalize_pipeline_paths(fsa_dir, base_outdir, assay_folder_name)

    # 1) Scan
    fsa_files = _scan_files(fsa_dir, mode)
    if not fsa_files:
        return [] if return_entries else None

    # 2) Analyze
    entries, _ = _analyze_files(fsa_files)
    if not entries:
        print_warning("Ingen gyldige entries etter analyse – avslutter.")
        return [] if return_entries else None

    return finalize_pipeline_run(
        entries,
        assay_dir,
        return_entries=return_entries,
        make_dit_reports=make_dit_reports,
        mode=mode,
    )
