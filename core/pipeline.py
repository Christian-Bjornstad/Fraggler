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

from core.assay_config import (
    ASSAY_CONFIG,
    SL_TARGET_FRAGMENTS_BP,
    SL_WINDOW_BP,
    OUTDIR_NAME,
)
from core.classification import classify_fsa
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
    build_dit_html_reports,
)


def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
) -> None:

    """
    Kjør full Fraggler-pipeline på alle .fsa-filer i fsa_dir.

    Parameters
    ----------
    fsa_dir : Path
        Katalog med .fsa-filer.
    base_outdir : Path | None
        Rotkatalog for rapporter. Hvis None brukes fsa_dir.
        Assay-figurer og kontrollrapport havner i base_outdir/OUTDIR_NAME,
        DIT-HTML-rapporter + tilhørende figurer havner i base_outdir/DIT_HTML_SUBDIR.
    """


    CONTROL_PREFIX_RE = re.compile(r"^(PK1|PK2|PK|RK|NK)_", re.IGNORECASE)

    WATER_RE = re.compile(r"^(v|water|h2o)[_\-]?", re.IGNORECASE)
    
    STAGE_PREFIX_RE = re.compile(r"^\d{5}_[a-f0-9]{8}_", re.IGNORECASE)
    
    def strip_stage_prefix(filename: str) -> str:
        return STAGE_PREFIX_RE.sub("", filename)

    def is_water_file(filename: str) -> bool:
        return WATER_RE.search(strip_stage_prefix(filename)) is not None

    def is_control_file(filename: str) -> bool:
        return CONTROL_PREFIX_RE.search(strip_stage_prefix(filename)) is not None


    if not fsa_dir.exists():
        print_warning(f"FSA-katalog finnes ikke: {fsa_dir}")
        return

    fsa_dir = Path(fsa_dir).expanduser()

    if base_outdir is None:
        base_outdir = fsa_dir
    else:
        base_outdir = Path(base_outdir).expanduser()


    assay_folder_name = (assay_folder_name or OUTDIR_NAME).strip()

    if not assay_folder_name:
        assay_folder_name = OUTDIR_NAME

    assay_dir = base_outdir / assay_folder_name
    dit_dir = assay_dir

    # Hent alle FSA-filer unntatt vann
    fsa_files = [
        p for p in sorted(fsa_dir.glob("*.fsa"))
        if not is_water_file(p.name)
    ]

    if mode == "controls":
        fsa_files = [p for p in fsa_files if is_control_file(p.name)]
        print_green(f"[INFO] Controls mode: {len(fsa_files)} control files selected.")


    if not fsa_files:
        print_warning(f"Fant ingen .fsa-filer i {fsa_dir}.")
        return

    print_green(f"Fant {len(fsa_files)} .fsa-filer: {[p.name for p in fsa_files]}")

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

        sample_channel = trace_channels[0]  # brukes i FsaFile-konstruktør

        if ladder == "LIZ":
            fsa = analyse_fsa_liz(fsa_path, sample_channel)
        else:
            fsa = analyse_fsa_rox(fsa_path, sample_channel)

        if fsa is None:
            skipped += 1
            continue

        # ---------------------------------------------------------
        # Peaks:
        #   - SL: automatisk peaks rundt definerte fragment-størrelser
        #   - alle andre assays: starter med tomme dataframes (kun manuell plukking)
        # ---------------------------------------------------------
        peaks_by_channel: dict[str, pd.DataFrame | None] = {}

        if assay == "SL":
            # Automatisk peak-detection kun for SL
            peaks_by_channel = auto_detect_sl_peaks(
                fsa,
                peak_channels=peak_channels,
                targets_bp=SL_TARGET_FRAGMENTS_BP,
                window_bp=SL_WINDOW_BP,
                min_height=800.0,  # juster terskel hvis du vil
            )
        else:
            # Alle andre assays: tomme peaks → kun manuell plukking
            for ch in peak_channels:
                peaks_by_channel[ch] = pd.DataFrame(
                    columns=["basepairs", "peaks", "keep"]
                )



        # Y-maks basert på trace-kanalene for denne assayen
        ymax = compute_zoom_ymax(
            fsa,
            bp_min,
            bp_max,
            trace_channels,
            assay_name=assay,
        )


        # Ladder QC-metrikker
        ladder_qc_status = "ok"
        ladder_r2 = np.nan
        n_ladder_steps = np.nan
        n_size_standard_peaks = np.nan
        try:
            metrics = compute_ladder_qc_metrics(fsa)
            ladder_r2 = metrics["r2"]
            n_ladder_steps = metrics["n_ladder_steps"]
            n_size_standard_peaks = metrics["n_size_standard_peaks"]
        except Exception as ex:
            print_warning(
                f"[LADDER_QC] Klarte ikke beregne QC for {fsa.file_name}: {ex}"
            )
            ladder_qc_status = "ladder_qc_failed"

        # SL-area-metrikker (gjelder kun SL-assay)
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
                print_warning(
                    f"[SL] Klarte ikke beregne SL-area for {fsa.file_name}: {ex}"
                )
                sl_metrics = None

        # Hent DIT-nummer fra filnavn, f.eks. 25OUM10166
        dit = extract_dit_from_name(fsa.file_name)

        entries.append(
            {
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
                "dit": dit,
                "ladder_qc_status": ladder_qc_status,
                "ladder_r2": ladder_r2,
                "n_ladder_steps": n_ladder_steps,
                "n_size_standard_peaks": n_size_standard_peaks,
                # Nytt felt for SL-spesifikke metrikker
                "sl_metrics": sl_metrics,
            }
        )

    print_green(
        f"[MASTER] Totalt {len(entries)} filer analysert. "
        f"{skipped} filer ble skippet (mangler assay-konfig eller ladderproblem)."
    )

    if not entries:
        print_warning("Ingen gyldige entries etter analyse – avslutter.")
        return

    # Her lagrer vi Plotly-HTML-fragmenter for DIT-rapporter
    dit_assay_html: dict[str, dict[str, str]] = defaultdict(dict)
    dit_combo_html: dict[str, dict[str, str]] = defaultdict(dict)


    # --------------------------------------------------------------
    # 1) Per-DIT grupperinger
    # --------------------------------------------------------------
    per_dit_entries: dict[str, list[dict]] = defaultdict(list)
    for e in entries:
        dit = e.get("dit")
        if dit:
            per_dit_entries[dit].append(e)

    dit_assay_pngs: dict[str, dict[str, str]] = defaultdict(dict)
    dit_combo_pngs: dict[str, dict[str, str]] = defaultdict(dict)

    for dit, dit_entries in per_dit_entries.items():
        # Gruppér innenfor DIT per assay
        per_assay: dict[str, list[dict]] = defaultdict(list)
        for e in dit_entries:
            per_assay[e["assay"]].append(e)

        # ---- Per-DIT-per-assay interaktive figurer (Plotly) ----
        for assay_name, a_entries in per_assay.items():
            a_sorted = sorted(a_entries, key=lambda x: x["fsa"].file_name)
            bp_min = a_sorted[0]["bp_min"]
            bp_max = a_sorted[0]["bp_max"]
            title = (
                f"{assay_name} – DIT {dit} (n={len(a_sorted)}) "
                f"[{bp_min:.0f}–{bp_max:.0f} bp]"
            )
            html_frag = build_interactive_assay_batch_plot_html(
                a_sorted,
                title=title,
                assay_name=assay_name,
            )
            dit_assay_html[dit][assay_name] = html_frag


        # ---- Kombinasjonsfigurer for TCRb A+B+C ----
        tcrb_assays = ["TCRbA", "TCRbB", "TCRbC"]
        have_all_tcrb = all(a in per_assay for a in tcrb_assays)

        if have_all_tcrb:
            # Sorter per assay, så vi har en stabil rekkefølge (rep1, rep2)
            per_tcrb_sorted = {
                a: sorted(per_assay[a], key=lambda x: x["fsa"].file_name)
                for a in tcrb_assays
            }

            for rep_idx in (1, 2):  # replikat 1 og 2
                zero_idx = rep_idx - 1
                rep_entries: list[dict] = []

                # Hent A/B/C for denne replikaten (hvis de finnes)
                for a in tcrb_assays:
                    lst = per_tcrb_sorted[a]
                    if len(lst) > zero_idx:
                        rep_entries.append(lst[zero_idx])

                # Vi trenger A + B + C for å lage kombinasjon
                if len(rep_entries) == 3:
                    bp_min = min(e["bp_min"] for e in rep_entries)
                    bp_max = max(e["bp_max"] for e in rep_entries)

                    # Lag nye entries med felles bp_min/bp_max
                    combo_entries: list[dict] = []
                    for e in rep_entries:
                        e2 = dict(e)
                        e2["bp_min"] = bp_min
                        e2["bp_max"] = bp_max
                        combo_entries.append(e2)

                    # Felles y-maks for denne replikaten, basert på *disse* tre
                    group_ymax = compute_group_ymax_for_entries(combo_entries)

                    title = (
                        f"TCRb A + B + C – replikat {rep_idx} "
                        f"(DIT {dit}) [{bp_min:.0f}–{bp_max:.0f} bp]"
                    )
                    html_frag = build_interactive_assay_batch_plot_html(
                        combo_entries,
                        title=title,
                        assay_name="TCRb_combo",
                        ymax_override=group_ymax,
                    )
                    dit_combo_html[dit][f"TCRb_combo_rep{rep_idx}"] = html_frag





        # ---- Kombinasjonsfigur for TCRg A1/A2/B1/B2 ----
        if "TCRgA" in per_assay and "TCRgB" in per_assay:
            tcrgA_list = sorted(per_assay["TCRgA"], key=lambda x: x["fsa"].file_name)
            tcrgB_list = sorted(per_assay["TCRgB"], key=lambda x: x["fsa"].file_name)

            if len(tcrgA_list) >= 2 and len(tcrgB_list) >= 2:
                all_entries = [
                    tcrgA_list[0], tcrgA_list[1],
                    tcrgB_list[0], tcrgB_list[1],
                ]

                bp_min = min(e["bp_min"] for e in all_entries)
                bp_max = max(e["bp_max"] for e in all_entries)

                combo_entries: list[dict] = []
                for e in all_entries:
                    e2 = dict(e)
                    e2["bp_min"] = bp_min
                    e2["bp_max"] = bp_max
                    combo_entries.append(e2)

                # Felles y-maks for alle 4 TCRg-plott
                group_ymax = compute_group_ymax_for_entries(combo_entries)

                title = (
                    f"TCRg A1/A2/B1/B2 (DIT {dit}) "
                    f"[{bp_min:.0f}–{bp_max:.0f} bp]"
                )

                html_frag = build_interactive_assay_batch_plot_html(
                    combo_entries,
                    title=title,
                    assay_name="TCRg_combo",
                    ymax_override=group_ymax,
                )
                dit_combo_html[dit]["TCRg_combo"] = html_frag



    # --------------------------------------------------------------
    # 2) Bygg pasient-HTML per DIT (i DIT-MAPPEN, kun Plotly-HTML)
    # --------------------------------------------------------------
    if make_dit_reports and mode != "controls":
        build_dit_html_reports(entries, assay_dir)

    # Returnér entries til QC-systemet hvis ønsket
    if return_entries:
        return entries

    return None
