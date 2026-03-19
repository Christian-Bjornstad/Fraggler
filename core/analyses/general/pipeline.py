from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from core.analysis import analyse_fsa_liz, analyse_fsa_rox, compute_ladder_qc_metrics
from core.analyses.general.classification import classify_fsa
from core.analyses.general.config import DEFAULT_BP_MAX, DEFAULT_BP_MIN, LIZ_LADDER, ROX_LADDER
from core.analyses.general.reporting import build_general_html_report
from core.analyses.shared_pipeline import finalize_pipeline_run, normalize_pipeline_paths, scan_fsa_files
from fraggler.fraggler import print_green, print_warning


def _scan_files(fsa_dir: Path, mode: str = "all") -> list[Path]:
    return scan_fsa_files(fsa_dir, mode=mode, recursive=False)


def _load_fsa(fsa_path: Path, meta: dict):
    ladder = str(meta.get("ladder") or ROX_LADDER).upper()
    sample_channel = str(meta.get("sample_channel") or meta.get("primary_peak_channel") or "DATA1")
    if ladder == LIZ_LADDER:
        return analyse_fsa_liz(
            fsa_path,
            sample_channel,
            ladder_name=ladder,
        )
    return analyse_fsa_rox(
        fsa_path,
        sample_channel,
        ladder_name=ladder,
    )


def _analyze_files(fsa_files: list[Path]) -> tuple[list[dict], int]:
    entries: list[dict] = []
    skipped = 0

    for fsa_path in fsa_files:
        classified = classify_fsa(fsa_path)
        if classified is None:
            skipped += 1
            continue

        fsa = _load_fsa(fsa_path, classified)
        if fsa is None:
            skipped += 1
            continue

        trace_channels = list(classified.get("trace_channels") or ["DATA1"])
        primary_peak_channel = str(classified.get("primary_peak_channel") or trace_channels[0])
        peaks_by_channel = {
            ch: pd.DataFrame(columns=["basepairs", "peaks", "area", "keep"])
            for ch in trace_channels
        }

        try:
            qc_metrics = compute_ladder_qc_metrics(fsa)
            ladder_r2 = float(qc_metrics.get("r2", np.nan))
        except Exception as exc:
            print_warning(f"[GENERAL] Ladder QC failed for {fsa.file_name}: {exc}")
            qc_metrics = {}
            ladder_r2 = float("nan")

        ladder_strategy = str(getattr(fsa, "ladder_fit_strategy", "auto_full"))
        ladder_missing = list(map(float, getattr(fsa, "ladder_missing_expected_steps", [])))
        ladder_status = "ok"
        if ladder_strategy == "manual_adjustment":
            ladder_status = "manual_adjustment"
        elif ladder_missing:
            ladder_status = "review_required"

        entries.append(
            {
                "analysis": "general",
                "assay": classified.get("assay", "GENERAL"),
                "group": classified.get("group", "sample"),
                "fsa": fsa,
                "trace_channels": trace_channels,
                "peak_channels": list(classified.get("peak_channels") or trace_channels),
                "primary_peak_channel": primary_peak_channel,
                "sample_channel": classified.get("sample_channel") or primary_peak_channel,
                "bp_min": float(classified.get("bp_min", DEFAULT_BP_MIN)),
                "bp_max": float(classified.get("bp_max", DEFAULT_BP_MAX)),
                "ladder": classified.get("ladder") or ROX_LADDER,
                "ladder_qc_status": ladder_status,
                "ladder_r2": ladder_r2,
                "ladder_fit_strategy": ladder_strategy,
                "ladder_missing_expected_steps": ladder_missing,
                "ladder_fit_note": str(getattr(fsa, "ladder_fit_note", "")),
                "ladder_review_required": bool(getattr(fsa, "ladder_review_required", False)),
                "ladder_expected_step_count": int(getattr(fsa, "ladder_expected_step_count", 0)),
                "ladder_fitted_step_count": int(getattr(fsa, "ladder_fitted_step_count", 0)),
                "n_ladder_steps": int(len(getattr(fsa, "ladder_steps", []))),
                "n_size_standard_peaks": int(len(getattr(fsa, "size_standard_peaks", []))),
                "peaks_by_channel": peaks_by_channel,
                "source_run_dir": classified.get("source_run_dir") or fsa_path.parent.name,
                "file_name": fsa.file_name,
            }
        )

    print_green(f"[GENERAL] Totalt {len(entries)} filer analysert. {skipped} skippet.")
    return entries, skipped


def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
) -> list[dict] | None:
    fsa_dir, assay_dir = normalize_pipeline_paths(fsa_dir, base_outdir, assay_folder_name)

    fsa_files = _scan_files(fsa_dir, mode)
    if not fsa_files:
        return [] if return_entries else None

    entries, _ = _analyze_files(fsa_files)
    if not entries:
        print_warning("Ingen gyldige entries etter analyse – avslutter.")
        return [] if return_entries else None

    if make_dit_reports:
        build_general_html_report(entries, assay_dir, run_label=assay_dir.name or fsa_dir.name)

    return finalize_pipeline_run(
        entries,
        assay_dir,
        return_entries=return_entries,
        make_dit_reports=False,
        mode=mode,
    )
