"""
Fraggler Diagnostics — Matplotlib Plot Helpers.

Static matplotlib plotting for zoom views, multi-channel traces, and peak labels.
"""
from __future__ import annotations

import numpy as np

from fraggler.fraggler import print_warning
from core.utils import CHANNEL_COLORS

from core.assay_config import (
    ASSAY_REFERENCE_RANGES,
    ASSAY_REFERENCE_LABEL,
    REFERENCE_SHADE_COLOR,
)


def compute_zoom_ymax(
    fsa,
    bp_min: float,
    bp_max: float,
    trace_channels,
    assay_name: str | None = None,
):
    """
    Beregn Y-maks for zoom basert på kanalene vi faktisk tegner som trace.

    Viktig:
      - Hvis assay_name finnes i ASSAY_REFERENCE_RANGES bruker vi UNION av disse
        vinduene som grunnlag for y-aksen (dvs. kun det gule referanseområdet).
      - Ellers bruker vi bp_min–bp_max fra ASSAY_CONFIG.
    """
    raw_df = fsa.sample_data_with_basepairs
    if raw_df is None:
        print_warning(
            "sample_data_with_basepairs er None – kan ikke beregne zoom-ymax."
        )
        return 0.0

    if "basepairs" not in raw_df.columns or "time" not in raw_df.columns:
        print_warning(
            f"sample_data_with_basepairs mangler 'basepairs'/'time' for {fsa.file_name} – kan ikke beregne zoom-ymax."
        )
        return 0.0

    bp_vals = raw_df["basepairs"].to_numpy()
    time_vals = raw_df["time"].astype(int).to_numpy()

    # -----------------------------
    # Velg bp-vindu for zoom
    # -----------------------------
    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        # Union av alle referansevinduer (slik vi fargelegger gult)
        mask_bp = np.zeros(bp_vals.shape, dtype=bool)
        for a, b in ASSAY_REFERENCE_RANGES[assay_name]:
            mask_bp |= (bp_vals >= float(a)) & (bp_vals <= float(b))
    else:
        # Fallback: bp_min–bp_max for assayen
        mask_bp = (bp_vals >= bp_min) & (bp_vals <= bp_max)

    if not np.any(mask_bp):
        print_warning(
            f"Ingen data i referanse-/bp-intervallet for {fsa.file_name} "
            f"(assay={assay_name}, {bp_min}–{bp_max} bp)."
        )
        return 0.0

    time_zoom = time_vals[mask_bp]

    # -----------------------------
    # Finn maks over alle relevante kanaler
    # -----------------------------
    channels_to_plot = [ch for ch in trace_channels if ch in fsa.fsa]
    if not channels_to_plot:
        print_warning(
            f"Ingen av trace-kanalene {trace_channels} funnet i {fsa.file_name}."
        )
        return 0.0

    ymax = 0.0
    for ch in channels_to_plot:
        trace = np.asarray(fsa.fsa[ch])
        mask_t = (time_zoom >= 0) & (time_zoom < len(trace))
        if not np.any(mask_t):
            continue
        y = trace[time_zoom[mask_t]]
        if y.size > 0 and np.any(np.isfinite(y)):
            ymax = max(ymax, float(np.nanmax(y)))

    return ymax


def draw_multi_channel_zoom_on_ax(
    ax,
    fsa,
    peaks_by_channel: dict,
    trace_channels,
    primary_peak_channel: str,
    bp_min: float,
    bp_max: float,
    assay_name: str | None = None,
):

    """
    - Plotter trace for relevante kanaler (trace_channels).
    - Overlay peaks for alle peak_channels (fra peaks_by_channel).
    - bp-labels på ALLE peak-kanaler, med enkel "repel" per kanal.
    - Ladder-peaks plottes som kryss.
    """


    # --------------------------------------------------
    # 0) Shade referanseområder (lys gul bakgrunn)
    # --------------------------------------------------
    shade_ranges = None

    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        shade_ranges = ASSAY_REFERENCE_RANGES[assay_name]
    else:
        # fallback: hele bp-vinduet som referanse
        shade_ranges = [(bp_min, bp_max)]

    for start_bp, end_bp in shade_ranges:
        # klipp mot plottvinduet så vi ikke tegner utenfor
        s = max(bp_min, float(start_bp))
        e = min(bp_max, float(end_bp))
        if e <= s:
            continue
        ax.axvspan(
            s,
            e,
            color=REFERENCE_SHADE_COLOR,
            alpha=0.25,      # «svak» gulfarge
            zorder=0,
        )


    raw_df = fsa.sample_data_with_basepairs
    if raw_df is None:
        print_warning(
            "sample_data_with_basepairs er None – kan ikke lage multi-channel bp-zoom."
        )
        return

    zoom_df = raw_df[(raw_df["basepairs"] >= bp_min) & (raw_df["basepairs"] <= bp_max)]
    if zoom_df.empty:
        print_warning(
            f"Ingen data i bp-intervallet {bp_min}–{bp_max} for {fsa.file_name}"
        )
        return

    time_zoom = zoom_df["time"].astype(int).to_numpy()
    bp_zoom = zoom_df["basepairs"].to_numpy()

    # Traces
    available = [k for k in fsa.fsa.keys() if k.startswith("DATA")]
    channels_to_plot = [ch for ch in trace_channels if ch in available]

    colors = CHANNEL_COLORS

    for ch in channels_to_plot:
        trace = np.asarray(fsa.fsa[ch])
        mask = (time_zoom >= 0) & (time_zoom < len(trace))
        if not np.any(mask):
            continue

        bp = bp_zoom[mask]
        y = trace[time_zoom[mask]]

        ax.plot(
            bp,
            y,
            linewidth=0.6,
            alpha=0.85,
            label=f"{ch} trace",
            color=colors.get(ch, None),
        )

    # Peaks
    peak_color_map = {
        "DATA1": "k",
        "DATA2": "tab:purple",
        "DATA3": "tab:orange",
    }

    for ch, df in peaks_by_channel.items():
        if df is None:
            continue
        if df.empty or "basepairs" not in df.columns or "peaks" not in df.columns:
            continue

        # Plot selve peak-punktene
        ax.scatter(
            df["basepairs"],
            df["peaks"],
            color=peak_color_map.get(ch, "k"),
            s=20,
            label=f"Detected peaks ({ch})",
        )

        # === NYTT: sett bp-labels for ALLE kanaler ===
        if df["peaks"].empty:
            continue

        # Litt x-offset slik at teksten ikke står rett over toppen
        x_offset = (bp_max - bp_min) * 0.004
        base_factor = 1.04   # litt luft over toppen
        extra_factor = 0.08  # ekstra hvis to topper er veldig nærme
        last_x = None

        for idx, (_, row) in enumerate(df.iterrows()):
            x = float(row["basepairs"]) + x_offset
            y = float(row["peaks"]) * base_factor

            # Enkel repel: hvis to labels står veldig tett i x, flytt den ene
            if last_x is not None and abs(x - last_x) < 2.5:
                y = float(row["peaks"]) * (base_factor + extra_factor)
                x += x_offset * (1 + (idx % 2))

            ax.text(
                x,
                y,
                f"{row['basepairs']:.1f}",
                fontsize=7,
                rotation=0,
                ha="left",
                va="bottom",
            )

            last_x = x

    # Ladder-peaks
    ladder_df = getattr(fsa, "size_standard_peaks", None)
    ladder_zoom = None
    if ladder_df is not None and hasattr(ladder_df, "columns"):
        ladder_bp_col = None
        for cand in ["basepairs", "bp", "size", "ladder_bp"]:
            if cand in ladder_df.columns:
                ladder_bp_col = cand
                break

        if ladder_bp_col is not None:
            ladder_zoom = ladder_df[
                (ladder_df[ladder_bp_col] >= bp_min)
                & (ladder_df[ladder_bp_col] <= bp_max)
            ].copy()
            ladder_zoom.rename(columns={ladder_bp_col: "bp"}, inplace=True)

    if ladder_zoom is not None and not ladder_zoom.empty:
        ladder_height_col = None
        for cand in ["peaks", "height", "max_peak"]:
            if cand in ladder_zoom.columns:
                ladder_height_col = cand
                break

        if ladder_height_col is not None:
            ax.scatter(
                ladder_zoom["bp"],
                ladder_zoom[ladder_height_col],
                marker="x",
                s=35,
                color="tab:gray",
                label="Ladder peaks",
            )

    ax.set_xlim(bp_min, bp_max)
    ax.grid(True, linestyle=":", linewidth=0.3)
    ax.set_title(fsa.file_name, fontsize=9)

    # Liten referanse-tekst i hjørnet av plottet
    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        ranges_str = ", ".join(
            f"{int(a)}–{int(b)} bp" for (a, b) in ASSAY_REFERENCE_RANGES[assay_name]
        )
        label_txt = ASSAY_REFERENCE_LABEL.get(assay_name, ranges_str)

        ax.text(
            0.01,
            0.98,
            label_txt,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            bbox=dict(
                boxstyle="round,pad=0.15",
                facecolor="white",
                alpha=0.7,
                edgecolor="none",
            ),
        )
