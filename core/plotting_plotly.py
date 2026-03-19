"""
Fraggler Diagnostics — Interactive Plotly Plot Builders.

Interactive peak editors and assay batch plots using Plotly.
"""
from __future__ import annotations

import json
import uuid
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from html import escape

from core.assay_config import (
    ASSAY_REFERENCE_RANGES,
    CHANNEL_COLORS,
    DEFAULT_TRACE_COLOR,
    NONSPECIFIC_PEAKS,
)
from core.analysis import (
    estimate_running_baseline,
    BASELINE_BIN_SIZE,
    BASELINE_QUANTILE,
    YMAX_PADDING_FACTOR,
)
from core.plotly_offline import local_plotly_tag as _local_plotly_tag

FLT3_NEGATIVE_CONTROL_YMIN = 250.0


def _flt3_peak_id(row: pd.Series, index: int) -> str:
    """Build a stable peak id for persisted FLT3 manual selections."""
    bp = float(row.get("basepairs", np.nan))
    height = float(row.get("peaks", np.nan))
    area = float(row.get("area", np.nan))
    bp_key = int(round(bp * 10)) if np.isfinite(bp) else 0
    height_key = int(round(height)) if np.isfinite(height) else 0
    area_key = int(round(area)) if np.isfinite(area) else 0
    return f"flt3_pk_{bp_key}_{height_key}_{area_key}_{index}"


def _flt3_candidate_peaks_for_entry(entry: dict) -> list[dict]:
    """Return FLT3 candidate peaks with channel metadata for the HTML editor."""
    peaks_by_channel = entry.get("peaks_by_channel", {})
    primary_channel = entry.get("primary_peak_channel")
    peaks = peaks_by_channel.get(primary_channel, pd.DataFrame())
    if peaks.empty:
        return []

    candidates: list[dict] = []
    for idx, (_, row) in enumerate(peaks.sort_values(["basepairs", "peaks"], ascending=[True, False]).iterrows()):
        bp = float(row.get("basepairs", np.nan))
        height = float(row.get("peaks", np.nan))
        area = float(row.get("area", np.nan))
        blue_area = float(row.get("area_DATA1", np.nan))
        green_area = float(row.get("area_DATA2", np.nan))
        dominant_channel = "DATA1" if np.nan_to_num(blue_area, nan=-1.0) >= np.nan_to_num(green_area, nan=-1.0) else "DATA2"
        if np.isfinite(blue_area) and np.isfinite(green_area) and abs(blue_area - green_area) < 1e-6:
            dominant_channel = "BOTH"
        peak_id = str(row.get("peak_id") or _flt3_peak_id(row, idx))
        candidates.append(
            {
                "peak_id": peak_id,
                "x": bp,
                "y": height,
                "area": area,
                "active": bool(row.get("keep", True)),
                "label": str(row.get("label", "")),
                "blue_area": blue_area,
                "green_area": green_area,
                "dominant_channel": dominant_channel,
            }
        )
    return candidates


def _flt3_preferred_channel(peak: dict, role: str = "mutant") -> str:
    blue_area = float(peak.get("blue_area", 0.0) or 0.0)
    green_area = float(peak.get("green_area", 0.0) or 0.0)
    if role == "wt":
        return "DATA1" if blue_area >= green_area else "DATA2"
    return "DATA2" if green_area >= blue_area else "DATA1"


def _flt3_manual_selection_defaults(candidate_peaks: list[dict]) -> dict:
    """Build an auto-prefilled FLT3 manual selection payload."""
    wt_candidates = [
        peak for peak in candidate_peaks
        if str(peak.get("label", "")) == "WT" and peak.get("active", True)
    ]
    wt_candidates.sort(key=lambda peak: (float(peak.get("area", 0.0)) * -1.0, float(peak.get("y", 0.0)) * -1.0, float(peak.get("x", 0.0))))
    mutant_ids = [
        {
            "peak_id": str(peak.get("peak_id")),
            "channel": _flt3_preferred_channel(peak, role="mutant"),
        }
        for peak in candidate_peaks
        if str(peak.get("label", "")) in {"MUT", "ITD"} and peak.get("active", True)
    ]
    mutant_ids = [item for item in mutant_ids if item.get("peak_id")]
    return {
        "enabled": False,
        "version": 1,
        "wt": {
            "peak_id": wt_candidates[0]["peak_id"] if wt_candidates else None,
            "channel": _flt3_preferred_channel(wt_candidates[0], role="wt") if wt_candidates else None,
        },
        "mutants": mutant_ids,
    }


def compute_group_ymax(entries: list[dict]) -> float:
    """
    Finn maksimal topp (RFU) i primary_peak_channel for alle entries i gruppen,
    innenfor hvert entry sitt bp-min/bp-max. Brukes for å gi felles y-akse.
    """
    group_max = 0.0

    for e in entries:
        fsa = e["fsa"]
        primary_ch = e["primary_peak_channel"]
        bp_min = float(e["bp_min"])
        bp_max = float(e["bp_max"])

        raw_df = getattr(fsa, "sample_data_with_basepairs", None)
        if raw_df is None or raw_df.empty:
            continue
        if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
            continue

        if primary_ch not in fsa.fsa:
            continue

        time_all = raw_df["time"].astype(int).to_numpy()
        bp_all = raw_df["basepairs"].to_numpy()
        trace_full = np.asarray(fsa.fsa[primary_ch])

        mask = (time_all >= 0) & (time_all < len(trace_full))
        if not np.any(mask):
            continue

        x_bp = bp_all[mask]
        y = trace_full[time_all[mask]]

        # Begrens til bp-vindu for denne analysen
        win = (x_bp >= bp_min) & (x_bp <= bp_max)
        if np.any(win):
            y_win = y[win]
        else:
            y_win = y

        if y_win.size == 0 or not np.any(np.isfinite(y_win)):
            continue

        local_max = float(np.nanmax(y_win))
        group_max = max(group_max, local_max)

    return group_max


def compute_group_ymax_all_channels(entries: list[dict]) -> float:
    """
    Finn maksimal RFU i alle aktuelle kanaler for en gruppe entries.

    Brukes for å gi felles y-akse i kombinasjonsfigurer (TCRb/TCRg).
    Vi tar maks over alle trace-kanaler som faktisk finnes i FSA-fila.
    """
    group_max = 0.0

    for e in entries:
        fsa = e["fsa"]
        primary_ch = e.get("primary_peak_channel")
        trace_channels = e.get("trace_channels") or [primary_ch]

        for ch in trace_channels:
            if not ch or ch not in fsa.fsa:
                continue
            arr = np.asarray(fsa.fsa[ch])
            if arr.size == 0 or not np.any(np.isfinite(arr)):
                continue
            local_max = float(np.nanmax(arr))
            if local_max > group_max:
                group_max = local_max

    return group_max


def compute_group_ymax_for_entries(entries: list[dict]) -> float:
    """
    Beregn felles Y-maks for en gruppe entries, basert på:

      - Alle trace-kanaler i entry["trace_channels"] (f.eks. DATA1 + DATA2).
      - Referansevindu (ASSAY_REFERENCE_RANGES[assay]) hvis det finnes.
      - Ellers bp_min–bp_max for entryet.

    Brukes for å gi kombinasjonsplott (TCRb/TCRg) en felles og
    *reelt* nødvendig y-akse.
    """
    group_ymax = 0.0

    for e in entries:
        fsa = e["fsa"]
        assay = e.get("assay")
        primary_ch = e.get("primary_peak_channel")
        bp_min = float(e["bp_min"])
        bp_max = float(e["bp_max"])

        raw_df = getattr(fsa, "sample_data_with_basepairs", None)
        if raw_df is None or raw_df.empty:
            continue
        if "basepairs" not in raw_df.columns or "time" not in raw_df.columns:
            continue

        bp_all = raw_df["basepairs"].to_numpy()
        time_all = raw_df["time"].astype(int).to_numpy()

        # Hvilke kanaler skal vi se på? (DATA1 + DATA2 …)
        trace_channels = e.get("trace_channels") or []
        available = [k for k in fsa.fsa.keys() if k.startswith("DATA")]
        channels_to_use = [ch for ch in trace_channels if ch in available]

        # Fallback: bruk primærkanal hvis trace_channels er tomme eller feil
        if not channels_to_use:
            if primary_ch and primary_ch in fsa.fsa:
                channels_to_use = [primary_ch]
            else:
                continue

        # --- Velg bp-vindu: referanse(r) hvis definert, ellers bp_min–bp_max ---
        if assay and assay in ASSAY_REFERENCE_RANGES:
            mask_bp = np.zeros_like(bp_all, dtype=bool)
            for a, b in ASSAY_REFERENCE_RANGES[assay]:
                mask_bp |= (bp_all >= float(a)) & (bp_all <= float(b))
        else:
            mask_bp = (bp_all >= bp_min) & (bp_all <= bp_max)

        if not np.any(mask_bp):
            continue

        time_win = time_all[mask_bp]

        for ch in channels_to_use:
            trace = np.asarray(fsa.fsa[ch])
            mask_t = (time_win >= 0) & (time_win < len(trace))
            if not np.any(mask_t):
                continue

            y = trace[time_win[mask_t]]
            if y.size == 0 or not np.any(np.isfinite(y)):
                continue

            local_max = float(np.nanmax(y))
            if np.isfinite(local_max) and local_max > group_ymax:
                group_ymax = local_max

    if group_ymax <= 0 or not np.isfinite(group_ymax):
        group_ymax = 1000.0

    return group_ymax


def _prepare_plot_data(entry: dict) -> dict | None:
    """Extracts and prepares data for plotting from an entry dict."""
    fsa = entry["fsa"]
    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return None
    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return None

    time_all = raw_df["time"].astype(int).to_numpy()
    bp_all = raw_df["basepairs"].to_numpy()
    
    primary_ch = entry["primary_peak_channel"]
    trace_channels = entry.get("trace_channels", [primary_ch])
    available = [k for k in fsa.fsa.keys() if k.startswith("DATA")]
    channels_to_plot = [ch for ch in trace_channels if ch in available]
    if not channels_to_plot:
        channels_to_plot = [primary_ch] if primary_ch in fsa.fsa else []

    if not channels_to_plot:
        return None

    # Common x-axis (bp) based on first channel
    first_ch = channels_to_plot[0]
    trace_first = np.asarray(fsa.fsa[first_ch])
    mask = (time_all >= 0) & (time_all < len(trace_first))
    if not np.any(mask):
        return None

    return {
        "fsa": fsa,
        "time_all": time_all,
        "bp_all": bp_all,
        "mask": mask,
        "bp_trace": bp_all[mask],
        "channels_to_plot": channels_to_plot,
        "primary_ch": primary_ch,
        "bp_min": float(entry["bp_min"]),
        "bp_max": float(entry["bp_max"]),
        "assay_name": entry.get("assay"),
        "group": entry.get("group"),
        "forced_ymax": entry.get("forced_ymax") or entry.get("force_ymax"),
        "forced_xmin": entry.get("forced_xmin"),
        "forced_xmax": entry.get("forced_xmax"),
        "peaks_by_channel": entry["peaks_by_channel"],
        "wt_bp": entry.get("wt_bp"),
        "mut_bp": entry.get("mut_bp"),
        "file_name": fsa.file_name,
        "sample_id": f"{fsa.file_name}_{primary_ch}"
    }


def _create_plotly_figure(data: dict) -> tuple[go.Figure, float, int]:
    """Constructs the Plotly figure and calculates y-axis limits."""
    fsa, mask, bp_trace = data["fsa"], data["mask"], data["bp_trace"]
    channels_to_plot, primary_ch = data["channels_to_plot"], data["primary_ch"]
    bp_min, bp_max, assay_name = data["bp_min"], data["bp_max"], data["assay_name"]
    time_all = data["time_all"]

    fig = go.Figure()
    ymax_auto_primary = 0.0
    ymax_auto_all = 0.0

    # Window for auto-y
    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        win_bp = np.zeros_like(bp_trace, dtype=bool)
        for a, b in ASSAY_REFERENCE_RANGES[assay_name]:
            win_bp |= (bp_trace >= float(a)) & (bp_trace <= float(b))
    else:
        win_bp = (bp_trace >= bp_min) & (bp_trace <= bp_max)

    for ch in channels_to_plot:
        full_trace = np.asarray(fsa.fsa[ch])
        baseline = estimate_running_baseline(full_trace, bin_size=BASELINE_BIN_SIZE, quantile=BASELINE_QUANTILE)
        full_corr = np.maximum(full_trace - baseline, 0.0)
        y_corr = full_corr[time_all[mask]]
        
        if y_corr.size == 0: continue

        color = CHANNEL_COLORS.get(ch, DEFAULT_TRACE_COLOR)
        fig.add_trace(go.Scatter(x=bp_trace, y=y_corr, mode="lines", name=f"{ch} trace", line=dict(width=1, color=color), hoverinfo="x+y"))

        y_win = y_corr[win_bp] if np.any(win_bp) else y_corr
        if y_win.size > 0 and np.any(np.isfinite(y_win)):
            local_max = float(np.nanmax(y_win))
            ymax_auto_all = max(ymax_auto_all, local_max)
            if ch == primary_ch: ymax_auto_primary = max(ymax_auto_primary, local_max)

    # 4) Select final ymax
    forced_ymax = data["forced_ymax"]
    if forced_ymax and float(forced_ymax) > 0:
        ymax = float(forced_ymax)
    else:
        multi_channel_assays = {"TCRgA", "TCRgB", "TCRg", "TCRγA", "TCRγB", "TCRγ", "TCRbA", "TCRbB", "TCRbC", "TCRβA", "TCRβB", "TCRβC"}
        use_all_channel_ymax = assay_name in multi_channel_assays or len(channels_to_plot) > 1
        base = ymax_auto_all if use_all_channel_ymax else (ymax_auto_primary or ymax_auto_all)
        ymax = base if base > 0 else 1000.0

    if assay_name in {"FLT3-ITD", "FLT3-D835", "NPM1"} and data.get("group") == "negative_control":
        ymax = max(ymax, FLT3_NEGATIVE_CONTROL_YMIN)

    # Shapes
    shapes = []
    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        for (a, b) in ASSAY_REFERENCE_RANGES[assay_name]:
            shapes.append(dict(type="rect", x0=float(a), x1=float(b), y0=0, y1=1, xref="x", yref="paper", fillcolor="rgba(235,232,203,0.5)", line_width=0))
    else:
        shapes.append(dict(type="rect", x0=float(bp_min), x1=float(bp_max), y0=0, y1=1, xref="x", yref="paper", fillcolor="rgba(235,232,203,0.15)", line_width=0))

    # NS Peaks
    if assay_name in NONSPECIFIC_PEAKS:
        trace_data = np.asarray(fsa.fsa[primary_ch]).astype(float)
        baseline = estimate_running_baseline(trace_data, bin_size=200, quantile=0.10)
        corr_trace = np.maximum(trace_data - baseline, 0.0)
        ns_x, ns_y, ns_text = [], [], []
        for ns_bp in NONSPECIFIC_PEAKS[assay_name]:
            shapes.append(dict(type="line", x0=float(ns_bp), x1=float(ns_bp), y0=0, y1=1, xref="x", yref="paper", line=dict(color="rgba(100, 116, 139, 0.7)", width=1.5, dash="dashdot")))
            mask_ns = (bp_trace >= (ns_bp - 3)) & (bp_trace <= (ns_bp + 3))
            if np.any(mask_ns):
                y_win = corr_trace[time_all[mask_ns]]
                if y_win.size > 0:
                    best_idx = np.argmax(y_win)
                    if float(y_win[best_idx]) > 100:
                        ns_x.append(float(bp_trace[np.where(mask_ns)[0][best_idx]]))
                        ns_y.append(float(y_win[best_idx]))
                        ns_text.append(f"Potensiell uspesifikk peak ({ns_bp}bp)<br>Høyde: {float(y_win[best_idx]):.0f}")
        if ns_x:
            fig.add_trace(go.Scatter(x=ns_x, y=ns_y, mode="markers", name="Uspesifikke peaks", marker=dict(symbol="x", size=8, color="#64748b", line=dict(color="white", width=0.5)), hovertext=ns_text, hoverinfo="text"))

    fig.update_layout(shapes=shapes)
    return fig, ymax, len(channels_to_plot)


def build_interactive_peak_plot_for_entry(entry: dict) -> str | None:
    """Main entry point for building an interactive peak plot."""
    data = _prepare_plot_data(entry)
    if not data: return None

    initial_peaks = []
    if data["assay_name"] == "SL":
        df0 = data["peaks_by_channel"].get(data["primary_ch"])
        if df0 is not None and not df0.empty:
            for _, row in df0.iterrows():
                x, y = float(row.get("basepairs", np.nan)), float(row.get("peaks", np.nan))
                if np.isfinite(x) and np.isfinite(y):
                    peak = {"x": x, "y": y, "active": True}
                    area = float(row.get("area", np.nan))
                    if np.isfinite(area):
                        peak["area"] = area
                    initial_peaks.append(peak)

    fig, ymax, peaks_trace_index = _create_plotly_figure(data)
    
    nice_title = f"{data['assay_name']} – {data['sample_id']}" if data["assay_name"] else data["sample_id"]
    fig.update_layout(
        title=nice_title, xaxis_title="Basepairs (bp)", yaxis_title="RFU", height=420,
        margin=dict(l=60, r=30, t=40, b=40), paper_bgcolor="white", plot_bgcolor="white",
        template="simple_white", font=dict(family="Inter, sans-serif", color="#0f172a"),
        clickmode="event", showlegend=True,
        hoverlabel=dict(bgcolor="white", font_size=12, font_family="Inter, sans-serif"),
        hoverdistance=20, spikedistance=20
    )
    # Select final x-range
    forced_xmin = data.get("forced_xmin")
    forced_xmax = data.get("forced_xmax")
    x_range = [
        float(forced_xmin) if forced_xmin is not None else data["bp_min"],
        float(forced_xmax) if forced_xmax is not None else data["bp_max"]
    ]

    fig.update_yaxes(range=[0.0, ymax * YMAX_PADDING_FACTOR], gridcolor="#f1f5f9", zerolinecolor="#cbd5e1")
    fig.update_xaxes(range=x_range, gridcolor="#f1f5f9", zerolinecolor="#cbd5e1")

    # Empty peaks trace for JS
    fig.add_trace(go.Scatter(x=[], y=[], mode="markers", name="Peaks", marker=dict(size=8, color="red", line=dict(color="black", width=1)), hovertemplate="bp=%{x:.2f}<br>height=%{y:.0f}<extra></extra>"))
    
    # The Peaks trace is always the LAST trace added so far.
    final_peaks_trace_index = len(fig.data) - 1
    primary_trace_index = data["channels_to_plot"].index(data["primary_ch"]) if data["primary_ch"] in data["channels_to_plot"] else 0

    div_id = f"peakplot_{data['sample_id'].replace('.','_')}_{uuid.uuid4().hex}"
    entry["_report_plot_id"] = div_id
    fig_json = json.dumps(fig.to_plotly_json())
    initial_peaks_json = json.dumps(initial_peaks)
    manual_ratio_assays = {"FLT3-ITD", "FLT3-D835"}
    is_manual_ratio_assay = data["assay_name"] in manual_ratio_assays
    manual_trace_channels = [ch for ch in data["channels_to_plot"] if str(ch).startswith("DATA")]
    flt3_initial_selection = entry.get("manual_ratio_selection") if is_manual_ratio_assay else {}
    flt3_manual_selection_json = json.dumps(
        flt3_initial_selection or {"enabled": False, "version": 2, "mutant_peak_ids": [], "wt_peak_ids": []}
    )
    manual_trace_channels_json = json.dumps(manual_trace_channels)

    manual_panel_html = ""
    if is_manual_ratio_assay:
        channel_button_html = ""
        if len(manual_trace_channels) > 1:
            button_parts = [
                '<button type="button" class="comment-toggle-btn" data-channel-choice="AUTO" '
                'style="width:auto; padding:6px 10px; border:1px solid #cbd5e1; border-radius:999px;">Auto</button>'
            ]
            for channel in manual_trace_channels:
                label = "Blue" if channel == "DATA1" else "Green" if channel == "DATA2" else channel
                button_parts.append(
                    f'<button type="button" class="comment-toggle-btn" data-channel-choice="{escape(channel)}" '
                    'style="width:auto; padding:6px 10px; border:1px solid #cbd5e1; border-radius:999px;">'
                    f'{escape(label)}</button>'
                )
            channel_button_html = (
                '<div id="{div_id}_flt3_channel_picker" style="display:flex; gap:0.4rem; flex-wrap:wrap; margin-top:0.4rem;">'
                '<span class="small" style="align-self:center;">Neste klikk:</span>'
                + "".join(button_parts)
                + "</div>"
            )

        if data["assay_name"] == "FLT3-ITD":
            table_head_html = (
                "<tr>"
                "<th>Kanal</th><th>WT</th><th>Mut</th><th>bp</th><th>RFU</th>"
                "<th>Blue area</th><th>Green area</th><th>Status</th><th>Fjern</th>"
                "</tr>"
            )
            helper_text = "Legg til peaks manuelt. Velg WT med radio-knappen, velg muterte peaks med avkrysning."
        else:
            table_head_html = (
                "<tr>"
                "<th>Kanal</th><th>WT</th><th>Mut</th><th>bp</th><th>RFU</th>"
                "<th>Area</th><th>Status</th><th>Fjern</th>"
                "</tr>"
            )
            helper_text = "Legg til peaks manuelt. Velg WT med radio-knappen, velg muterte peaks med avkrysning."

        manual_panel_html = f"""
<div id="{div_id}_flt3_panel" class="peak-table-container" style="display:none; margin-top:0.75rem;">
    <div style="display:flex; align-items:center; gap:0.5rem; flex-wrap:wrap; margin-bottom:0.5rem;">
        <span id="{div_id}_flt3_mode_badge" class="status-badge warning">Manuell ratio</span>
        <span id="{div_id}_flt3_mode_text" class="small">{helper_text}</span>
        <button type="button" id="{div_id}_flt3_reset" class="comment-toggle-btn" style="width:auto; padding:6px 10px;">Nullstill valg</button>
    </div>
    {channel_button_html.format(div_id=div_id) if channel_button_html else ""}
    <div id="{div_id}_flt3_ratio_cards" style="display:grid; grid-template-columns:repeat(auto-fit, minmax(130px, 1fr)); gap:0.5rem; margin:0.75rem 0;">
        <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:0.6rem 0.75rem;">
            <div class="small">WT area</div>
            <div id="{div_id}_flt3_denominator" style="font-weight:700;">&mdash;</div>
        </div>
        <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:0.6rem 0.75rem;">
            <div class="small">Mut area</div>
            <div id="{div_id}_flt3_numerator" style="font-weight:700;">&mdash;</div>
        </div>
        <div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:10px; padding:0.6rem 0.75rem;">
            <div class="small">Ratio</div>
            <div id="{div_id}_flt3_ratio_value" style="font-weight:700;">&mdash;</div>
        </div>
    </div>
    <table id="{div_id}_flt3_table">
        <thead>{table_head_html}</thead>
        <tbody></tbody>
    </table>
</div>
"""

    html_fragment = f"""
<div id="{div_id}" class="peak-editor-block"></div>
<div id="{div_id}_table_container" class="peak-table-container" style="display:none;">
    <table id="{div_id}_table">
        <thead>
            <tr><th>Peak Størrelse (bp)</th><th>Høyde (RFU)</th><th>Area</th></tr>
        </thead>
        <tbody></tbody>
    </table>
</div>
{manual_panel_html}
<script type="text/javascript">
(function() {{
  var fig = {fig_json};
  var initialPeaks = {initial_peaks_json};
  var flt3InitialSelection = {flt3_manual_selection_json};
  var divId = "{div_id}";
  var gd = document.getElementById(divId);
  if (!gd) return;

  var areaWindowBp = 5.0;
  var peaksTraceIndex = {final_peaks_trace_index};
  var primaryTraceIndex = {primary_trace_index};
  var assayName = {json.dumps(data["assay_name"])};
  var isFlt3ItD = assayName === "FLT3-ITD";
  var isManualRatioAssay = assayName === "FLT3-ITD" || assayName === "FLT3-D835";
  var manualTraceChannels = {manual_trace_channels_json};
  var manualClickChannel = manualTraceChannels.length > 1 ? "AUTO" : (manualTraceChannels[0] || "AUTO");
  var expectedWtBp = {json.dumps(data.get("wt_bp"))};
  var expectedMutBp = {json.dumps(data.get("mut_bp"))};
  var plotFileName = {json.dumps(data.get("file_name", ""))};
  var overviewIdPrefix = "overview_" + plotFileName.replace(/\./g, "_").replace(/ /g, "_");
  var initialPlotState = (window.ReportPlotManager && window.ReportPlotManager.getInitialStateForPlot)
    ? window.ReportPlotManager.getInitialStateForPlot(divId)
    : null;

  if (initialPlotState && typeof initialPlotState === "object") {{
    fig.layout = fig.layout || {{}};
    fig.layout.xaxis = fig.layout.xaxis || {{}};
    fig.layout.yaxis = fig.layout.yaxis || {{}};
    if (Array.isArray(initialPlotState.xaxis_range) && initialPlotState.xaxis_range.length === 2) {{
      fig.layout.xaxis.range = initialPlotState.xaxis_range;
      fig.layout.xaxis.autorange = false;
    }}
    if (Array.isArray(initialPlotState.yaxis_range) && initialPlotState.yaxis_range.length === 2) {{
      fig.layout.yaxis.range = initialPlotState.yaxis_range;
      fig.layout.yaxis.autorange = false;
    }}
  }}

  Plotly.newPlot(gd, fig.data, fig.layout, {{ responsive: true, displaylogo: false }}).then(function(g) {{
    if (window.ReportPlotManager) {{ window.ReportPlotManager.register(g); }}
    var baseShapes = (g.layout.shapes || []).slice();
    var baseAnnots = (g.layout.annotations || []).slice();
    var primaryTrace = g.data[primaryTraceIndex] || {{}};

    function decodePlotlyArray(val) {{
      if (Array.isArray(val)) return val;
      if (ArrayBuffer.isView(val)) return Array.from(val);
      if (!val || typeof val !== "object") return [];
      if (typeof val.length === "number") {{
        try {{ return Array.from(val); }} catch (e) {{}}
      }}
      if (typeof val.bdata === "string" && typeof val.dtype === "string") {{
        var binary = atob(val.bdata);
        var bytes = new Uint8Array(binary.length);
        for (var i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
        var buf = bytes.buffer;
        switch (val.dtype) {{
          case "f8": return Array.from(new Float64Array(buf));
          case "f4": return Array.from(new Float32Array(buf));
          case "i1": return Array.from(new Int8Array(buf));
          case "u1": return Array.from(new Uint8Array(buf));
          case "i2": return Array.from(new Int16Array(buf));
          case "u2": return Array.from(new Uint16Array(buf));
          case "i4": return Array.from(new Int32Array(buf));
          case "u4": return Array.from(new Uint32Array(buf));
          default: return [];
        }}
      }}
      return [];
    }}

    var traceXYCache = g.data.map(function(trace) {{
      return {{
        x: decodePlotlyArray(trace && trace.x),
        y: decodePlotlyArray(trace && trace.y)
      }};
    }});

    function peakHalfWidthBp(xCenter) {{
      if (assayName === "FLT3-D835") {{
        if (Number.isFinite(expectedMutBp) && Math.abs(xCenter - Number(expectedMutBp)) <= 3.0) return 0.5;
        if (Number.isFinite(expectedWtBp) && Math.abs(xCenter - Number(expectedWtBp)) <= 6.0) return 1.2;
        if (Math.abs(xCenter - 150.0) <= 6.0) return 0.8;
        return 0.8;
      }}
      if (assayName === "FLT3-ITD") {{
        if (Number.isFinite(expectedWtBp) && Math.abs(xCenter - Number(expectedWtBp)) <= 8.0) return 2.0;
        if (xCenter >= 335.0) return 1.0;
        return 2.0;
      }}
      return areaWindowBp;
    }}

    function computePeakArea(xCenter, traceIndex) {{
      var traceData = traceXYCache[Number.isFinite(traceIndex) ? traceIndex : primaryTraceIndex] || traceXYCache[primaryTraceIndex] || {{}};
      var traceX = Array.isArray(traceData.x) ? traceData.x : [];
      var traceY = Array.isArray(traceData.y) ? traceData.y : [];
      var halfWidth = peakHalfWidthBp(xCenter);
      var total = 0.0;
      for (var i = 0; i < traceX.length; i++) {{
        var x = Number(traceX[i]);
        var y = Number(traceY[i]);
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        if (Math.abs(x - xCenter) <= halfWidth) total += y;
      }}
      return total;
    }}

    function channelLabel(channel) {{
      if (channel === "DATA1") return "blue";
      if (channel === "DATA2") return "green";
      if (channel === "DATA3") return "orange";
      return "auto";
    }}

    function channelAccent(channel) {{
      if (channel === "DATA1") return "#2563eb";
      if (channel === "DATA2") return "#16a34a";
      if (channel === "DATA3") return "#ea580c";
      return "#dc2626";
    }}

    function findPeakById(peakId) {{
      for (var i = 0; i < peaks.length; i++) {{
        if (peaks[i].peak_id === peakId) return peaks[i];
      }}
      return null;
    }}

    function normalizeSourceChannel(raw, fallbackCurveNumber) {{
      var channel = raw ? String(raw).toUpperCase() : "";
      if (channel === "DATA1" || channel === "DATA2" || channel === "DATA3") return channel;
      if (Number.isFinite(fallbackCurveNumber)) {{
        var traceName = String((g.data[fallbackCurveNumber] && g.data[fallbackCurveNumber].name) || "");
        if (traceName.indexOf("DATA1") === 0) return "DATA1";
        if (traceName.indexOf("DATA2") === 0) return "DATA2";
        if (traceName.indexOf("DATA3") === 0) return "DATA3";
      }}
      return null;
    }}

    function traceIndexForChannel(channel) {{
      for (var i = 0; i < g.data.length; i++) {{
        var traceName = String((g.data[i] && g.data[i].name) || "");
        if (traceName.indexOf(channel) === 0) return i;
      }}
      return primaryTraceIndex;
    }}

    function computePeakAreaForChannel(xCenter, channel) {{
      if (!channel) return 0.0;
      return computePeakArea(xCenter, traceIndexForChannel(channel));
    }}

    function nearestTracePointForChannel(xClick, channel) {{
      var traceIndex = traceIndexForChannel(channel);
      var traceData = traceXYCache[traceIndex] || {{ x: [], y: [] }};
      var traceX = Array.isArray(traceData.x) ? traceData.x : [];
      var traceY = Array.isArray(traceData.y) ? traceData.y : [];
      if (!traceX.length) return null;
      var bestIdx = -1;
      var bestDist = Infinity;
      for (var i = 0; i < traceX.length; i++) {{
        var x = Number(traceX[i]);
        var y = Number(traceY[i]);
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        var dist = Math.abs(x - xClick);
        if (dist < bestDist) {{
          bestDist = dist;
          bestIdx = i;
        }}
      }}
      if (bestIdx < 0) return null;
      return {{
        x: Number(traceX[bestIdx]),
        y: Number(traceY[bestIdx]),
        curveNumber: traceIndex
      }};
    }}

    function makePeakId(p, idx) {{
      if (p && p.peak_id) return String(p.peak_id);
      var xKey = Number.isFinite(Number(p && p.x)) ? Math.round(Number(p.x) * 10) : 0;
      var yKey = Number.isFinite(Number(p && p.y)) ? Math.round(Number(p.y)) : 0;
      var areaKey = Number.isFinite(Number(p && p.area)) ? Math.round(Number(p.area)) : 0;
      var channelKey = normalizeSourceChannel(p && p.source_channel, Number(p && p.curve_number)) || "AUTO";
      return "flt3_pk_" + xKey + "_" + yKey + "_" + areaKey + "_" + channelKey + "_" + idx;
    }}

    function peakIdFor(peak, idx) {{
      if (!peak) return "";
      if (!peak.peak_id) {{
        peak.peak_id = makePeakId(peak, idx);
      }}
      return String(peak.peak_id);
    }}

    function ensurePeakIds() {{
      for (var i = 0; i < peaks.length; i++) {{
        peakIdFor(peaks[i], i);
      }}
    }}

    function normalizePeak(p, idx) {{
      var x = Number(p && p.x);
      var y = Number(p && p.y);
      var sourceChannel = normalizeSourceChannel(p && p.source_channel, Number(p && p.curve_number));
      if (!sourceChannel && manualTraceChannels.length === 1) {{
        sourceChannel = manualTraceChannels[0];
      }}
      var blueArea = Number(p && p.blue_area);
      var greenArea = Number(p && p.green_area);
      if (!Number.isFinite(blueArea)) blueArea = computePeakAreaForChannel(x, "DATA1");
      if (!Number.isFinite(greenArea)) greenArea = computePeakAreaForChannel(x, "DATA2");
      var area = Number(p && p.area);
      if (!Number.isFinite(area)) {{
        if (sourceChannel === "DATA2") area = greenArea;
        else if (sourceChannel === "DATA1") area = blueArea;
        else area = computePeakAreaForChannel(x, sourceChannel || manualTraceChannels[0] || null);
      }}
      return {{
        x: x,
        y: y,
        area: area,
        active: !(p && p.active === false),
        peak_id: makePeakId(p, idx),
        blue_area: blueArea,
        green_area: greenArea,
        source_channel: sourceChannel,
        curve_number: Number(p && p.curve_number)
      }};
    }}

    function cloneSelection(selection) {{
      var out = {{
        enabled: false,
        version: 2,
        mutant_peak_ids: [],
        wt_peak_ids: []
      }};
      if (selection && typeof selection === "object") {{
        out.enabled = !!selection.enabled;
        out.version = Number.isFinite(Number(selection.version)) ? Number(selection.version) : 2;
        if (Array.isArray(selection.mutant_peak_ids)) {{
          out.mutant_peak_ids = selection.mutant_peak_ids.filter(function(v) {{ return !!v; }}).map(String);
        }} else if (Array.isArray(selection.mutants)) {{
          out.mutant_peak_ids = selection.mutants
            .filter(function(item) {{ return item && item.peak_id; }})
            .map(function(item) {{ return String(item.peak_id); }});
        }}
        if (Array.isArray(selection.wt_peak_ids)) {{
          out.wt_peak_ids = selection.wt_peak_ids.filter(function(v) {{ return !!v; }}).map(String);
        }} else if (selection.wt_peak_id) {{
          out.wt_peak_ids = [String(selection.wt_peak_id)];
        }} else if (selection.wt && selection.wt.peak_id) {{
          out.wt_peak_ids = [String(selection.wt.peak_id)];
        }}
      }}
      return out;
    }}

    var peaks = [];
    var initialPeakData = (window.PeakManager && window.PeakManager.getInitialPeakDataForPlot)
      ? window.PeakManager.getInitialPeakDataForPlot(divId)
      : null;
    if (initialPeakData && Array.isArray(initialPeakData.peaks) && initialPeakData.peaks.length) {{
      peaks = initialPeakData.peaks.slice();
    }} else if (window.PeakManager) {{
      peaks = window.PeakManager.getInitialPeaksForPlot(divId);
    }}
    if (!peaks || peaks.length === 0) {{ peaks = (initialPeaks && Array.isArray(initialPeaks)) ? initialPeaks.slice() : []; }}
    peaks = peaks.map(function(p, idx) {{ return normalizePeak(p, idx); }}).filter(function(p) {{ return Number.isFinite(p.x) && Number.isFinite(p.y); }});
    ensurePeakIds();

    var manualSelection = cloneSelection(initialPeakData && initialPeakData.flt3_manual_ratio_selection ? initialPeakData.flt3_manual_ratio_selection : flt3InitialSelection);

    function updatePeakManagerRegistration() {{
      if (!window.PeakManager) return;
      ensurePeakIds();
      window.PeakManager.registerPlot(divId, {{
        getPeaks: function() {{ return peaks; }},
        getPeakData: function() {{
          return {{
            peaks: peaks,
            flt3_manual_ratio_selection: {{
              enabled: !!manualSelection.enabled && manualSelection.mutant_peak_ids.length > 0,
              version: 2,
              mutant_peak_ids: manualSelection.mutant_peak_ids.slice(),
              wt_peak_ids: manualSelection.wt_peak_ids.slice()
            }}
          }};
        }}
      }});
    }}

    function peakAreaForSourceChannel(peak) {{
      if (!peak) return 0.0;
      if (peak.source_channel === "DATA2") return Number(peak.green_area) || 0.0;
      if (peak.source_channel === "DATA1") return Number(peak.blue_area) || 0.0;
      return Number(peak.area) || 0.0;
    }}

    function wtToleranceBp() {{
      if (assayName === "FLT3-D835") return 4.0;
      if (assayName === "FLT3-ITD") return 8.0;
      return 5.0;
    }}

    function inferredWtByChannel() {{
      var wtMap = {{}};
      for (var i = 0; i < peaks.length; i++) {{
        var peak = peaks[i];
        if (!peak.active || !peak.source_channel) continue;
        var distance = Math.abs(Number(peak.x) - Number(expectedWtBp));
        if (!Number.isFinite(distance) || distance > wtToleranceBp()) continue;
        if (!wtMap[peak.source_channel]) {{
          wtMap[peak.source_channel] = peak;
          continue;
        }}
        var current = wtMap[peak.source_channel];
        var currentDistance = Math.abs(Number(current.x) - Number(expectedWtBp));
        if (distance < currentDistance || (distance === currentDistance && peakAreaForSourceChannel(peak) > peakAreaForSourceChannel(current))) {{
          wtMap[peak.source_channel] = peak;
        }}
      }}
      return wtMap;
    }}

    function isWtPeak(peak, wtMap) {{
      if (!peak || !peak.source_channel) return false;
      return !!(wtMap[peak.source_channel] && wtMap[peak.source_channel].peak_id === peak.peak_id);
    }}

    function isManualWt(peakId) {{
      return manualSelection.wt_peak_ids.indexOf(peakId) >= 0;
    }}

    function selectionSummary() {{
      var wtMap;
      if (manualSelection.wt_peak_ids.length > 0) {{
        wtMap = {{}};
        for (var w = 0; w < manualSelection.wt_peak_ids.length; w++) {{
          var manualWtPeak = findPeakById(manualSelection.wt_peak_ids[w]);
          if (manualWtPeak && manualWtPeak.active && manualWtPeak.source_channel) {{
            wtMap[manualWtPeak.source_channel] = manualWtPeak;
          }}
        }}
      }} else {{
        wtMap = inferredWtByChannel();
      }}
      var selectedMutants = peaks.filter(function(peak, idx) {{
        return peak.active && manualSelection.mutant_peak_ids.indexOf(peakIdFor(peak, idx)) >= 0;
      }});
      var activeChannels = [];
      if (assayName === "FLT3-ITD") {{
        for (var i = 0; i < selectedMutants.length; i++) {{
          if (selectedMutants[i].source_channel && activeChannels.indexOf(selectedMutants[i].source_channel) < 0) {{
            activeChannels.push(selectedMutants[i].source_channel);
          }}
        }}
      }} else if (selectedMutants.length) {{
        activeChannels = [manualTraceChannels[0] || "DATA3"];
      }}
      var numerator = 0.0;
      var denominator = 0.0;
      var wtParts = [];
      var mutParts = [];
      var missingWtChannels = [];
      for (var c = 0; c < activeChannels.length; c++) {{
        var channel = activeChannels[c];
        var wtPeak = wtMap[channel];
        if (!wtPeak) {{
          missingWtChannels.push(channelLabel(channel));
          continue;
        }}
        var wtArea = peakAreaForSourceChannel(wtPeak);
        denominator += wtArea;
        wtParts.push(channelLabel(channel) + ": " + wtPeak.x.toFixed(1) + " bp");
      }}
      for (var m = 0; m < selectedMutants.length; m++) {{
        var mutPeak = selectedMutants[m];
        var mutArea = peakAreaForSourceChannel(mutPeak);
        numerator += mutArea;
        mutParts.push(mutPeak.x.toFixed(1) + " bp (" + channelLabel(mutPeak.source_channel) + ")");
      }}
      return {{
        wtMap: wtMap,
        wtText: wtParts.length ? wtParts.join(", ") : "Ingen WT valgt ennå",
        mutText: mutParts.length ? mutParts.join(", ") : "Ingen mutant valgt",
        selectedMutants: selectedMutants,
        missingWtChannels: missingWtChannels,
        numerator: numerator,
        denominator: denominator,
        ratio: denominator > 0 ? numerator / denominator : 0.0,
        valid: selectedMutants.length > 0 && missingWtChannels.length === 0 && denominator > 0
      }};
    }}

    function updateRatioCards(summary) {{
      var numeratorEl = document.getElementById(divId + "_flt3_numerator");
      var denominatorEl = document.getElementById(divId + "_flt3_denominator");
      var ratioEl = document.getElementById(divId + "_flt3_ratio_value");
      if (numeratorEl) numeratorEl.textContent = summary.selectedMutants.length ? summary.numerator.toFixed(0) : "—";
      if (denominatorEl) denominatorEl.textContent = summary.denominator > 0 ? summary.denominator.toFixed(0) : "—";
      if (ratioEl) ratioEl.textContent = summary.valid ? summary.ratio.toFixed(4) : "—";
    }}

    function applyChannelChoiceStyles() {{
      var picker = document.getElementById(divId + "_flt3_channel_picker");
      if (!picker) return;
      var buttons = picker.querySelectorAll("[data-channel-choice]");
      for (var i = 0; i < buttons.length; i++) {{
        var isActive = buttons[i].dataset.channelChoice === manualClickChannel;
        buttons[i].style.background = isActive ? "#dbeafe" : "transparent";
        buttons[i].style.color = isActive ? "#1d4ed8" : "#475569";
        buttons[i].style.borderColor = isActive ? "#93c5fd" : "#cbd5e1";
      }}
    }}

    function redrawSelectionControls() {{
      if (!isManualRatioAssay) return;
      var badge = document.getElementById(divId + "_flt3_mode_badge");
      var modeText = document.getElementById(divId + "_flt3_mode_text");
      var summaryStatus = selectionSummary();
      if (badge) {{
        if (summaryStatus.valid) {{
          badge.textContent = "Ratio klar";
          badge.className = "status-badge manual";
        }} else if (summaryStatus.selectedMutants.length > 0 && summaryStatus.denominator === 0) {{
          badge.textContent = "Velg WT";
          badge.className = "status-badge warning";
        }} else if (summaryStatus.selectedMutants.length > 0) {{
          badge.textContent = "Mangler WT";
          badge.className = "status-badge warning";
        }} else {{
          badge.textContent = "Velg WT og mutant";
          badge.className = "status-badge warning";
        }}
      }}
      if (modeText) {{
        var summary = summaryStatus;
        var baseText = "WT: " + summary.wtText + " - Mut: " + summary.mutText;
        if (summary.valid) {{
          baseText += " - Ratio: " + summary.ratio.toFixed(4) + " (" + summary.numerator.toFixed(0) + "/" + summary.denominator.toFixed(0) + ")";
        }} else if (summary.missingWtChannels.length) {{
          baseText += " - Mangler WT i " + summary.missingWtChannels.join(", ");
        }}
        modeText.textContent = baseText;
      }}
      updateRatioCards(summaryStatus);
      applyChannelChoiceStyles();

      // Update the overview summary table at the top of the report
      if (plotFileName) {{
        var wtEl = document.getElementById("overview_wt_" + overviewIdPrefix.replace("overview_", ""));
        var mutEl = document.getElementById("overview_mut_" + overviewIdPrefix.replace("overview_", ""));
        var ratEl = document.getElementById("overview_ratio_" + overviewIdPrefix.replace("overview_", ""));
        if (wtEl) {{
          if (summaryStatus.wtText !== "Ingen WT valgt ennå") {{
            wtEl.innerHTML = summaryStatus.wtText;
          }} else {{
            wtEl.innerHTML = "<span class='small'>Manuell WT</span>";
          }}
        }}
        if (mutEl) {{
          if (summaryStatus.mutText !== "Ingen mutant valgt") {{
            mutEl.innerHTML = summaryStatus.mutText;
          }} else {{
            mutEl.innerHTML = "<span class='small'>Velg mutantpeaks manuelt</span>";
          }}
        }}
        if (ratEl) {{
          if (summaryStatus.valid) {{
            ratEl.innerHTML = summaryStatus.ratio.toFixed(4) + " <span class='status-badge manual'>Manual</span>";
          }} else {{
            ratEl.innerHTML = "\u2014";
          }}
        }}
      }}

      var table = document.getElementById(divId + "_flt3_table");
      if (!table) return;
      var summaryForTable = summaryStatus;
      var checks = table.querySelectorAll('input[type="checkbox"][data-peak-id]');
      for (var j = 0; j < checks.length; j++) {{
        var role = checks[j].dataset.role || "";
        if (role === "wt") {{
          checks[j].checked = manualSelection.wt_peak_ids.indexOf(checks[j].dataset.peakId) >= 0;
        }} else {{
          checks[j].checked = manualSelection.mutant_peak_ids.indexOf(checks[j].dataset.peakId) >= 0;
        }}
      }}
      var rows = table.querySelectorAll('tbody tr[data-peak-id]');
      for (var k = 0; k < rows.length; k++) {{
        var row = rows[k];
        var peakId = row.dataset.peakId;
        var peak = findPeakById(peakId);
        row.classList.toggle("selected-wt", isWtPeak(peak, summaryForTable.wtMap));
        row.classList.toggle("selected-mut", manualSelection.mutant_peak_ids.indexOf(peakId) >= 0);
      }}
      updatePeakManagerRegistration();
    }}

    function clearMutantSelection() {{
      manualSelection.mutant_peak_ids = [];
      manualSelection.wt_peak_ids = [];
      manualSelection.enabled = false;
      redrawSelectionControls();
    }}

    function toggleMutantSelection(peakId, checked) {{
      if (checked) {{
        if (manualSelection.mutant_peak_ids.indexOf(peakId) < 0) {{
          manualSelection.mutant_peak_ids.push(peakId);
        }}
      }} else {{
        manualSelection.mutant_peak_ids = manualSelection.mutant_peak_ids.filter(function(id) {{ return id !== peakId; }});
      }}
      manualSelection.enabled = manualSelection.mutant_peak_ids.length > 0;
    }}

    function removePeakById(peakId) {{
      peaks = peaks.filter(function(peak) {{ return peak.peak_id !== peakId; }});
      manualSelection.mutant_peak_ids = manualSelection.mutant_peak_ids.filter(function(id) {{ return id !== peakId; }});
      manualSelection.wt_peak_ids = manualSelection.wt_peak_ids.filter(function(id) {{ return id !== peakId; }});
      manualSelection.enabled = manualSelection.mutant_peak_ids.length > 0;
    }}

    function renderFlt3CandidateTable() {{
      if (!isManualRatioAssay) return;
      var panel = document.getElementById(divId + "_flt3_panel");
      var tbody = document.querySelector("#" + divId + "_flt3_table tbody");
      if (!panel || !tbody) return;
      panel.style.display = "block";
      ensurePeakIds();
      var summary = selectionSummary();
      var wtMap = summary.wtMap;
      var sortedPeaks = peaks.slice().sort(function(a, b) {{ return a.x - b.x; }});
      var html = "";
      for (var i = 0; i < sortedPeaks.length; i++) {{
        var peak = sortedPeaks[i];
        var peakId = peakIdFor(peak, peaks.indexOf(peak));
        var blueArea = Number.isFinite(Number(peak.blue_area)) ? Number(peak.blue_area) : 0.0;
        var greenArea = Number.isFinite(Number(peak.green_area)) ? Number(peak.green_area) : 0.0;
        var wtHere = isManualWt(peakId) || isWtPeak(peak, wtMap);
        var isMut = manualSelection.mutant_peak_ids.indexOf(peakId) >= 0;
        var sourceText = channelLabel(peak.source_channel);
        html += "<tr data-peak-id='" + peakId + "'>";
        html += "<td>" + sourceText + "</td>";
        html += "<td><input type='checkbox' data-role='wt' data-peak-id='" + peakId + "'" + (wtHere ? " checked" : "") + (isMut ? " disabled" : "") + "></td>";
        html += "<td><input type='checkbox' data-role='mutant' data-peak-id='" + peakId + "'" + (wtHere ? " disabled" : "") + "></td>";
        html += "<td>" + peak.x.toFixed(1) + "</td>";
        html += "<td>" + peak.y.toFixed(0) + "</td>";
        if (isFlt3ItD) {{
          html += "<td>" + blueArea.toFixed(0) + "</td>";
          html += "<td>" + greenArea.toFixed(0) + "</td>";
        }} else {{
          html += "<td>" + peakAreaForSourceChannel(peak).toFixed(0) + "</td>";
        }}
        html += "<td>" + (wtHere ? "WT" : isMut ? "Mutant" : "Peak") + "</td>";
        html += "<td><button type='button' class='comment-toggle-btn' data-role='delete-peak' data-peak-id='" + peakId + "' style='width:auto; padding:4px 8px;'>Slett</button></td>";
        html += "</tr>";
      }}
      tbody.innerHTML = html;
      redrawSelectionControls();
    }}

    function nearestPeakIdx(xClick, preferredChannel) {{
      if (!peaks.length) return -1;
      var bestIdx = -1;
      var bestScore = Infinity;
      for (var i = 0; i < peaks.length; i++) {{
        var peak = peaks[i];
        if (preferredChannel && peak.source_channel !== preferredChannel) continue;
        var d = Math.abs(peak.x - xClick);
        var score = d;
        if (score < bestScore) {{
          bestScore = score;
          bestIdx = i;
        }}
      }}
      return bestIdx;
    }}

    function displayXForPeak(peak) {{
      if (!peak) return 0.0;
      if (peak.source_channel === "DATA1") return peak.x - 0.16;
      if (peak.source_channel === "DATA2") return peak.x + 0.16;
      return peak.x;
    }}

    function redrawPeaks() {{
      ensurePeakIds();
      var summary = isManualRatioAssay ? selectionSummary() : null;
      var xs = peaks.map(function(p) {{ return displayXForPeak(p); }});
      var ys = peaks.map(function(p) {{ return p.y; }});
      var op = peaks.map(function(p) {{ return p.active ? 1.0 : 0.3; }});
      var col = peaks.map(function(p) {{ return p.active ? channelAccent(p.source_channel) : "#94a3b8"; }});
      var texts = peaks.map(function(p) {{
        if (!p.active) return "";
        var shortChannel = p.source_channel ? channelLabel(p.source_channel).charAt(0).toUpperCase() : "";
        return shortChannel ? (p.x.toFixed(1) + " " + shortChannel) : p.x.toFixed(1);
      }});
      var symbols = peaks.map(function(p) {{
        if (!p.active) return "circle-open";
        if (summary && isWtPeak(p, summary.wtMap)) return "diamond";
        if (manualSelection.mutant_peak_ids.indexOf(p.peak_id) >= 0) return "square";
        return "circle";
      }});
      var sizes = peaks.map(function(p) {{
        if (!p.active) return 8;
        if (summary && isWtPeak(p, summary.wtMap)) return 11;
        if (manualSelection.mutant_peak_ids.indexOf(p.peak_id) >= 0) return 10;
        return 8;
      }});

      Plotly.restyle(
        g,
        {{
          x: [xs],
          y: [ys],
          "marker.opacity": [op],
          "marker.color": [col],
          "marker.symbol": [symbols],
          "marker.size": [sizes],
          text: [texts]
        }},
        [peaksTraceIndex]
      );

      var ann = [];
      var tbody = document.querySelector("#{div_id}_table tbody");
      var tableHtml = "";
      
      // Sort peaks by size (X) before rendering the table
      var sortedPeaks = peaks.slice().sort(function(a, b) {{ return a.x - b.x; }});

      for (var i = 0; i < sortedPeaks.length; i++) {{
        var p = sortedPeaks[i];
        if (!p.active) continue;
        
        // Add annotation
        ann.push({{
          x: displayXForPeak(p),
          y: p.y * 1.03,
          xref: "x",
          yref: "y",
          text: texts[peaks.indexOf(p)],
          showarrow: false,
          font: {{ size: 9, color: "#222" }},
          xanchor: "left",
          yanchor: "bottom"
        }});
        
        // Add table row
        tableHtml += "<tr><td>" + p.x.toFixed(1) + "</td><td>" + p.y.toFixed(0) + "</td><td>" + p.area.toFixed(0) + "</td></tr>";
      }}

      Plotly.relayout(g, {{ shapes: baseShapes, annotations: baseAnnots.concat(ann) }});
      
      // Update Table
      if (tbody) tbody.innerHTML = tableHtml;
      var tCont = document.getElementById("{div_id}_table_container");
      if (tCont) tCont.style.display = isManualRatioAssay ? "none" : ((tableHtml !== "") ? "block" : "none");
      if (isManualRatioAssay) {{
        renderFlt3CandidateTable();
      }}
    }}

    if (peaks.length) {{ redrawPeaks(); }}
    if (isManualRatioAssay) {{
      renderFlt3CandidateTable();
      redrawSelectionControls();
    }}

    gd.on("plotly_click", function(ev) {{
      if (!ev.points || !ev.points.length) return;
      var pt = ev.points[0];
      var xVal = pt.x;
      var yVal = pt.y;
      var isShift = !!(ev.event && ev.event.shiftKey);
      var requestedChannel = manualClickChannel === "AUTO" ? normalizeSourceChannel(null, pt.curveNumber) : manualClickChannel;
      if (!requestedChannel && manualTraceChannels.length === 1) {{
        requestedChannel = manualTraceChannels[0];
      }}
      if (isManualRatioAssay && requestedChannel) {{
        var channelPoint = nearestTracePointForChannel(xVal, requestedChannel);
        if (channelPoint) {{
          xVal = channelPoint.x;
          yVal = channelPoint.y;
          pt.curveNumber = channelPoint.curveNumber;
        }}
      }}

      if (isShift) {{
        var idxDel = nearestPeakIdx(xVal, requestedChannel);
        if (idxDel >= 0) {{
          var removedPeakId = peaks[idxDel].peak_id;
          peaks.splice(idxDel, 1);
          manualSelection.mutant_peak_ids = manualSelection.mutant_peak_ids.filter(function(id) {{ return id !== removedPeakId; }});
          manualSelection.wt_peak_ids = manualSelection.wt_peak_ids.filter(function(id) {{ return id !== removedPeakId; }});
          manualSelection.enabled = manualSelection.mutant_peak_ids.length > 0;
          redrawPeaks();
        }}
        return;
      }}

      var idx = nearestPeakIdx(xVal, requestedChannel);
      if (idx >= 0 && Math.abs(peaks[idx].x - xVal) < 0.4) {{
        peaks[idx].active = !peaks[idx].active;
        if (!peaks[idx].active) {{
          manualSelection.mutant_peak_ids = manualSelection.mutant_peak_ids.filter(function(id) {{ return id !== peaks[idx].peak_id; }});
          manualSelection.enabled = manualSelection.mutant_peak_ids.length > 0;
        }}
        redrawPeaks();
        return;
      }}

      var sourceChannel = requestedChannel || normalizeSourceChannel(null, pt.curveNumber);
      var blueArea = computePeakAreaForChannel(xVal, "DATA1");
      var greenArea = computePeakAreaForChannel(xVal, "DATA2");
      var peakArea = sourceChannel ? computePeakAreaForChannel(xVal, sourceChannel) : computePeakArea(xVal, pt.curveNumber);
      var newPeak = {{
        x: xVal,
        y: yVal,
        area: peakArea,
        blue_area: blueArea,
        green_area: greenArea,
        source_channel: sourceChannel,
        curve_number: pt.curveNumber,
        active: true
      }};
      newPeak.peak_id = makePeakId(newPeak, peaks.length);
      peaks.push(newPeak);
      redrawPeaks();
    }});

    if (isManualRatioAssay) {{
      var table = document.getElementById(divId + "_flt3_table");
      var resetBtn = document.getElementById(divId + "_flt3_reset");
      var picker = document.getElementById(divId + "_flt3_channel_picker");
      if (table) {{
        table.addEventListener("change", function(ev) {{
          var target = ev.target || {{}};
          var peakId = target.dataset ? target.dataset.peakId : "";
          if (!peakId) return;
          if (target.type === "checkbox" && target.dataset.role === "wt") {{
            if (target.checked) {{
              if (manualSelection.wt_peak_ids.indexOf(peakId) < 0) {{
                manualSelection.wt_peak_ids.push(peakId);
              }}
              manualSelection.mutant_peak_ids = manualSelection.mutant_peak_ids.filter(function(id) {{ return id !== peakId; }});
            }} else {{
              manualSelection.wt_peak_ids = manualSelection.wt_peak_ids.filter(function(id) {{ return id !== peakId; }});
            }}
            renderFlt3CandidateTable();
            redrawPeaks();
            return;
          }}
          if (target.type === "checkbox" && target.dataset.role === "mutant") {{
            toggleMutantSelection(peakId, !!target.checked);
            if (!!target.checked) {{
              manualSelection.wt_peak_ids = manualSelection.wt_peak_ids.filter(function(id) {{ return id !== peakId; }});
            }}
          }}
          renderFlt3CandidateTable();
          redrawPeaks();
        }});
        table.addEventListener("click", function(ev) {{
          var target = ev.target || {{}};
          if (!target.dataset) return;
          if (target.dataset.role === "delete-peak" && target.dataset.peakId) {{
            removePeakById(target.dataset.peakId);
            redrawPeaks();
          }}
        }});
      }}
      if (resetBtn) {{
        resetBtn.addEventListener("click", function() {{
          clearMutantSelection();
          redrawPeaks();
        }});
      }}
      if (picker) {{
        picker.addEventListener("click", function(ev) {{
          var target = ev.target || {{}};
          if (!target.dataset || !target.dataset.channelChoice) return;
          manualClickChannel = target.dataset.channelChoice;
          applyChannelChoiceStyles();
        }});
      }}
      updatePeakManagerRegistration();
      applyChannelChoiceStyles();
    }}
  }});
}}).call(this);
</script>
"""
    return html_fragment

def build_interactive_assay_batch_plot_html(
    entries: list[dict],
    title: str,
    assay_name: str | None = None,
    ymax_override: float | None = None,   # <--- NY PARAM
) -> str:
    """
    Interaktiv visning for én assay med én liten Plotly-editor per entry.

    - Vanlig modus:
        y-aksen auto-tilpasses per entry (som før).
    - Kombinasjonsmodus (TCRb A+B+C, alle TCRg):
        pass inn ymax_override = max(e["ymax"] for e in entries),
        slik at alle småplott får samme y-akse og kan sammenliknes direkte.
    """

    if not entries:
        return "<p><em>Ingen entries for denne assayen.</em></p>"

    html_parts: list[str] = []
    html_parts.append(f"<h2>{escape(title)}</h2>")

    # Vi inkluderer Plotly-script første gang
    plotly_script_included = False

    for idx, e in enumerate(entries, start=1):
        fsa = e["fsa"]
        primary_ch = e["primary_peak_channel"]
        bp_min = float(e["bp_min"])
        bp_max = float(e["bp_max"])
        assay = e.get("assay", assay_name)

        raw_df = getattr(fsa, "sample_data_with_basepairs", None)
        if raw_df is None or raw_df.empty:
            html_parts.append(
                f"<p><strong>{escape(fsa.file_name)}</strong>: "
                f"<em>mangler sample_data_with_basepairs – kan ikke lage interaktiv figur.</em></p>"
            )
            continue

        if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
            html_parts.append(
                f"<p><strong>{escape(fsa.file_name)}</strong>: "
                f"<em>sample_data_with_basepairs mangler 'time'/'basepairs'.</em></p>"
            )
            continue

        time_all = raw_df["time"].astype(int).to_numpy()
        bp_all = raw_df["basepairs"].to_numpy()

        if primary_ch not in fsa.fsa:
            html_parts.append(
                f"<p><strong>{escape(fsa.file_name)}</strong>: "
                f"<em>fant ikke kanal {escape(primary_ch)} i FSA-filen.</em></p>"
            )
            continue

        trace_full = np.asarray(fsa.fsa[primary_ch])

        mask = (time_all >= 0) & (time_all < len(trace_full))
        if not np.any(mask):
            html_parts.append(
                f"<p><strong>{escape(fsa.file_name)}</strong>: "
                f"<em>ingen gyldige punkter i trace.</em></p>"
            )
            continue

        bp_trace = bp_all[mask]
        y_trace = trace_full[time_all[mask]]

        # --- Beregn auto-ymax ut fra referanse-/bp-vindu ---
        # Hvis assay har egne referansevinduer bruker vi UNION av disse
        # (slik som shadingen i plottet). Ellers bruker vi bp_min–bp_max.
        if assay and assay in ASSAY_REFERENCE_RANGES:
            mask_win = np.zeros_like(bp_trace, dtype=bool)
            for a, b in ASSAY_REFERENCE_RANGES[assay]:
                mask_win |= (bp_trace >= float(a)) & (bp_trace <= float(b))
        else:
            mask_win = (bp_trace >= bp_min) & (bp_trace <= bp_max)

        if np.any(mask_win):
            y_window = y_trace[mask_win]
        else:
            y_window = y_trace


        if y_window.size == 0 or np.all(np.isnan(y_window)):
            auto_ymax = 1000.0
        else:
            auto_ymax = float(np.nanmax(y_window))
            if auto_ymax <= 0:
                auto_ymax = 1000.0

        # Hvis vi har fått inn et globalt maksimum for gruppa (kombinasjon),
        # så bruker vi det, ellers bruker vi auto_ymax som før.
        if ymax_override is not None:
            ymax = float(ymax_override)
        else:
            ymax = auto_ymax

        # --- Bygg Plotly-figur for denne ene fila ---
        fig = go.Figure()

        # velg farge basert på kanal
        color = CHANNEL_COLORS.get(primary_ch, DEFAULT_TRACE_COLOR)

        # Linje-trace
        fig.add_trace(
            go.Scatter(
                x=bp_trace,
                y=y_trace,
                mode="lines",
                name=f"{primary_ch} trace",
                line=dict(width=1, color=color),
                hoverinfo="x+y",
            )
        )

        # Tom peaks-trace
        fig.add_trace(
            go.Scatter(
                x=[],
                y=[],
                mode="markers+text",
                name="Manuelle peaks",
                marker=dict(size=9, color="red", line=dict(color="black", width=1)),
                text=[],
                textposition="top center",
                textfont=dict(size=9),
                hovertemplate="bp=%{x:.2f}<br>height=%{y:.0f}<extra></extra>",
            )
        )

        # Referanse-shapes
        shapes = []
        if assay and assay in ASSAY_REFERENCE_RANGES:
            for (a, b) in ASSAY_REFERENCE_RANGES[assay]:
                shapes.append(
                    dict(
                        type="rect",
                        x0=float(a),
                        x1=float(b),
                        y0=0,
                        y1=1,
                        xref="x",
                        yref="paper",
                        fillcolor="rgba(235,232,203,0.25)",
                        line_width=0,
                    )
                )
        # Uspesifikke topper (vertikale dotted linjer + detection) i batch-plott
        if assay in NONSPECIFIC_PEAKS:
            ns_x, ns_y, ns_text = [], [], []
            trace_data = np.asarray(fsa.fsa[primary_ch]).astype(float)
            try:
                baseline = estimate_running_baseline(trace_data, bin_size=BASELINE_BIN_SIZE, quantile=BASELINE_QUANTILE)
                corr_trace = trace_data - baseline
                corr_trace[corr_trace < 0] = 0.0
            except Exception:
                corr_trace = trace_data

            for ns_bp in NONSPECIFIC_PEAKS[assay]:
                shapes.append(dict(
                    type="line",
                    x0=float(ns_bp), x1=float(ns_bp),
                    y0=0, y1=1, xref="x", yref="paper",
                    line=dict(color="rgba(100, 116, 139, 0.6)", width=1.2, dash="dashdot"),
                    name=f"NS_{ns_bp}"
                ))
                
                # Peak deteksjon i ±3 bp vindu
                mask = (bp_trace >= (ns_bp - 3)) & (bp_trace <= (ns_bp + 3))
                if np.any(mask):
                    idx_in_mask = np.where(mask)[0]
                    y_win = corr_trace[time_all[mask]]
                    if y_win.size > 0:
                        best_local_idx = np.argmax(y_win)
                        peak_h = float(y_win[best_local_idx])
                        peak_bp = float(bp_trace[idx_in_mask[best_local_idx]])
                        if peak_h > 150: # Litt strengere i batch
                            ns_x.append(peak_bp)
                            ns_y.append(peak_h)
                            ns_text.append(f"Uspesifikk peak ({ns_bp}bp)")

            if ns_x:
                fig.add_trace(go.Scatter(
                    x=ns_x, y=ns_y, mode="markers",
                    name="Uspesifikke peaks",
                    marker=dict(symbol="x", size=7, color="#64748b", opacity=0.8),
                    hovertext=ns_text, hoverinfo="text",
                    showlegend=False
                ))

        fig.update_layout(
            title=f"{fsa.file_name} – {primary_ch}",
            xaxis_title="Basepairs (bp)",
            yaxis_title="RFU",
            height=420,
            margin=dict(l=60, r=30, t=60, b=50),
            shapes=shapes,
            paper_bgcolor="white",
            plot_bgcolor="white",
            clickmode="event",
            showlegend=True,
            template="simple_white",
            font=dict(family="Inter, -apple-system, sans-serif", color="#0f172a"),
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="Inter, -apple-system, sans-serif"),
        )

        # Y-akse: felles ymax hvis angitt, ellers per-entry auto
        fig.update_yaxes(
            rangemode="tozero",
            range=[0.0, ymax * 1.15],
            gridcolor="#f1f5f9",
            zerolinecolor="#cbd5e1"
        )

        # X-akse: start-zoom til assay-vindu
        fig.update_xaxes(
            range=[bp_min, bp_max],
            gridcolor="#f1f5f9",
            zerolinecolor="#cbd5e1"
        )

        fig_json = json.dumps(fig.to_plotly_json())

        # Unik div-id per entry
        safe_name = (
            f"{fsa.file_name}_{primary_ch}_{idx}"
            .replace(" ", "_")
            .replace(".", "_")
            .replace("/", "_")
            .replace("\\", "_")
            .replace(":", "_")
        )
        div_id = f"assay_peak_editor_{safe_name}"

        # --- HTML for denne editoren ---
        html_parts.append("<div class='assay-block'>")
        html_parts.append(f"<h3>{escape(fsa.file_name)} – {escape(primary_ch)}</h3>")
        html_parts.append(f"<div id='{div_id}'></div>")
        html_parts.append(
            f"<div id='{div_id}_table_container' class='peak-table-container' style='display:none;'>"
            f"<table id='{div_id}_table'>"
            "<thead><tr><th>Peak Størrelse (bp)</th><th>Høyde (RFU)</th><th>Area</th></tr></thead>"
            "<tbody></tbody></table></div>"
        )
        html_parts.append(
            "<p class='small'>Klikk på tracen for å legge til peaks. "
            "Shift+klikk for å slette nærmeste peak.</p>"
        )
        # Skjult JSON-buffer – ikke synlig, men lar vi stå for evt. senere bruk
        html_parts.append(
            f"<pre id='{div_id}_peaks_json' class='small' style='display:none;'>[]</pre>"
        )
        html_parts.append("</div>")

        # Inkluder Plotly-script én gang
        if not plotly_script_included:
            html_parts.append(_local_plotly_tag(Path("."), version="2.35.2"))
            plotly_script_included = True

        # JS for akkurat denne editoren – synka med PeakManager
    html_parts.append(f"""
<script type="text/javascript">
(function() {{
  var fig = {fig_json};
  var divId = "{div_id}";
  var gd = document.getElementById(divId);
  if (!gd) return;
  var areaWindowBp = 5.0;
  var assayName = {json.dumps(e.get("assay"))};
  var expectedWtBp = {json.dumps(e.get("wt_bp"))};
  var expectedMutBp = {json.dumps(e.get("mut_bp"))};
  var initialPlotState = (window.ReportPlotManager && window.ReportPlotManager.getInitialStateForPlot)
    ? window.ReportPlotManager.getInitialStateForPlot(divId)
    : null;

  if (initialPlotState && typeof initialPlotState === "object") {{
    fig.layout = fig.layout || {{}};
    fig.layout.xaxis = fig.layout.xaxis || {{}};
    fig.layout.yaxis = fig.layout.yaxis || {{}};
    if (Array.isArray(initialPlotState.xaxis_range) && initialPlotState.xaxis_range.length === 2) {{
      fig.layout.xaxis.range = initialPlotState.xaxis_range;
      fig.layout.xaxis.autorange = false;
    }}
    if (Array.isArray(initialPlotState.yaxis_range) && initialPlotState.yaxis_range.length === 2) {{
      fig.layout.yaxis.range = initialPlotState.yaxis_range;
      fig.layout.yaxis.autorange = false;
    }}
  }}

  Plotly.newPlot(gd, fig.data, fig.layout).then(function(g) {{
    if (window.ReportPlotManager) {{ window.ReportPlotManager.register(g); }}
    var primaryTrace = g.data[0] || {{}};

    function decodePlotlyArray(val) {{
      if (Array.isArray(val)) return val;
      if (ArrayBuffer.isView(val)) return Array.from(val);
      if (!val || typeof val !== "object") return [];
      if (typeof val.length === "number") {{
        try {{ return Array.from(val); }} catch (e) {{}}
      }}
      if (typeof val.bdata === "string" && typeof val.dtype === "string") {{
        var binary = atob(val.bdata);
        var bytes = new Uint8Array(binary.length);
        for (var i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
        var buf = bytes.buffer;
        switch (val.dtype) {{
          case "f8": return Array.from(new Float64Array(buf));
          case "f4": return Array.from(new Float32Array(buf));
          case "i1": return Array.from(new Int8Array(buf));
          case "u1": return Array.from(new Uint8Array(buf));
          case "i2": return Array.from(new Int16Array(buf));
          case "u2": return Array.from(new Uint16Array(buf));
          case "i4": return Array.from(new Int32Array(buf));
          case "u4": return Array.from(new Uint32Array(buf));
          default: return [];
        }}
      }}
      return [];
    }}

    var traceXYCache = g.data.map(function(trace) {{
      return {{
        x: decodePlotlyArray(trace && trace.x),
        y: decodePlotlyArray(trace && trace.y)
      }};
    }});

    function peakHalfWidthBp(xCenter) {{
      if (assayName === "FLT3-D835") {{
        if (Number.isFinite(expectedMutBp) && Math.abs(xCenter - Number(expectedMutBp)) <= 3.0) return 0.5;
        if (Number.isFinite(expectedWtBp) && Math.abs(xCenter - Number(expectedWtBp)) <= 6.0) return 1.2;
        if (Math.abs(xCenter - 150.0) <= 6.0) return 0.8;
        return 0.8;
      }}
      if (assayName === "FLT3-ITD") {{
        if (Number.isFinite(expectedWtBp) && Math.abs(xCenter - Number(expectedWtBp)) <= 8.0) return 2.0;
        if (xCenter >= 335.0) return 1.0;
        return 2.0;
      }}
      return areaWindowBp;
    }}

    function computePeakArea(xCenter, traceIndex) {{
      var traceData = traceXYCache[Number.isFinite(traceIndex) ? traceIndex : primaryTraceIndex] || traceXYCache[primaryTraceIndex] || {{}};
      var traceX = Array.isArray(traceData.x) ? traceData.x : [];
      var traceY = Array.isArray(traceData.y) ? traceData.y : [];
      var halfWidth = peakHalfWidthBp(xCenter);
      var total = 0.0;
      for (var i = 0; i < traceX.length; i++) {{
        var x = Number(traceX[i]);
        var y = Number(traceY[i]);
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        if (Math.abs(x - xCenter) <= halfWidth) total += y;
      }}
      return total;
    }}

    function normalizePeak(p) {{
      var x = Number(p && p.x);
      var y = Number(p && p.y);
      var area = Number(p && p.area);
      return {{
        x: x,
        y: y,
        area: Number.isFinite(area) ? area : computePeakArea(x, primaryTraceIndex),
        active: !(p && p.active === false)
      }};
    }}

    var peaks = [];
    if (window.PeakManager) {{
        peaks = window.PeakManager.getInitialPeaksForPlot(divId);
    }}
    peaks = (Array.isArray(peaks) ? peaks : []).map(normalizePeak).filter(function(p) {{
      return Number.isFinite(p.x) && Number.isFinite(p.y);
    }});

    if (window.PeakManager) {{
        window.PeakManager.registerPlot(divId, {{
            getPeaks: function() {{ return peaks; }}
        }});
    }}

    function redrawPeaks() {{
      var xs = peaks.map(function(p) {{ return p.x; }});
      var ys = peaks.map(function(p) {{ return p.y; }});
      var op = peaks.map(function(p) {{ return p.active ? 1.0 : 0.3; }});
      var col = peaks.map(function(p) {{ return p.active ? "red" : "gray"; }});
      var texts = peaks.map(function(p) {{ return p.active ? p.x.toFixed(1) : ""; }});

      Plotly.restyle(g, {{
        x: [xs],
        y: [ys],
        "marker.opacity": [op],
        "marker.color": [col],
        text: [texts]
      }}, [1]); // peaks-trace er index 1

      var arr = peaks.map(function(p) {{
        return {{ x: p.x, y: p.y, area: p.area, active: p.active }};
      }});
      var pre = document.getElementById(divId + "_peaks_json");
      if (pre) {{
        pre.textContent = JSON.stringify(arr, null, 2);
      }}

      // --- Table Rendering ---
      var tbody = document.querySelector("#" + divId + "_table tbody");
      var tableHtml = "";
      
      var sortedPeaks = peaks.slice().sort(function(a, b) {{ return a.x - b.x; }});
      for (var i = 0; i < sortedPeaks.length; i++) {{
        var p = sortedPeaks[i];
        if (!p.active) continue;
        tableHtml += "<tr><td>" + p.x.toFixed(1) + "</td><td>" + p.y.toFixed(0) + "</td><td>" + p.area.toFixed(0) + "</td></tr>";
      }}
      
      if (tbody) tbody.innerHTML = tableHtml;
      var tCont = document.getElementById(divId + "_table_container");
      if (tCont) tCont.style.display = (tableHtml !== "") ? "block" : "none";
    }}

    function findNearestPeakIdx(xClick) {{
      if (!peaks.length) return -1;
      var bestIdx = 0;
      var bestDist = Math.abs(peaks[0].x - xClick);
      for (var i = 1; i < peaks.length; i++) {{
        var d = Math.abs(peaks[i].x - xClick);
        if (d < bestDist) {{
          bestDist = d;
          bestIdx = i;
        }}
      }}
      return bestIdx;
    }}

    if (peaks.length) {{
        redrawPeaks();
    }}

    gd.on("plotly_click", function(ev) {{
      if (!ev || !ev.points || !ev.points.length) return;
      var pt = ev.points[0];
      var xVal = pt.x;
      var yVal = pt.y;
      var isShift = !!(ev.event && ev.event.shiftKey);

      if (isShift) {{
        var idx = findNearestPeakIdx(xVal);
        if (idx >= 0) {{
          peaks.splice(idx, 1);
          redrawPeaks();
        }}
        return;
      }}

      var idx = findNearestPeakIdx(xVal);
      if (idx >= 0 && Math.abs(peaks[idx].x - xVal) < 0.4) {{
        peaks[idx].active = !peaks[idx].active;
        redrawPeaks();
        return;
      }}

      peaks.push({{ x: xVal, y: yVal, area: computePeakArea(xVal, pt.curveNumber), active: true }});
      redrawPeaks();
    }});
  }});
}})();
</script>
""")

    return "\n".join(html_parts)
