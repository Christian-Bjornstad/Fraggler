"""
Fraggler Diagnostics — Interactive Plotly Plot Builders.

Interactive peak editors and assay batch plots using Plotly.
"""
from __future__ import annotations

import json
import uuid
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
from html import escape

from core.assay_config import (
    ASSAY_REFERENCE_RANGES,
    CHANNEL_COLORS,
    DEFAULT_TRACE_COLOR,
)
from core.analysis import estimate_running_baseline
from core.plotly_offline import local_plotly_tag as _local_plotly_tag


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


def build_interactive_peak_plot_for_entry(entry: dict) -> str | None:
    """
    Bygger interaktiv Plotly-figur for én entry (én FSA) med ALLE trace-kanaler:

      - Tegner alle kanaler i entry["trace_channels"] (f.eks. DATA1 + DATA2).
      - Start-visning: bp_min–bp_max.
      - Shade referanseområder for assay (ASSAY_REFERENCE_RANGES).
      - Manuelle peaks: klikk (legg til), klikk på peak (toggle aktiv/grå),
        Shift+klikk (slett nærmeste peak).

      - Hvis entry inneholder 'forced_ymax' (eller legacy 'force_ymax') > 0, brukes denne som
        felles y-maks (nyttig for kombinasjonsfigurer).
      - Ellers:
          * standard: auto-y basert på primary_peak_channel
          * TCRgA/TCRgB: auto-y basert på alle trace-kanaler (DATA1 + DATA2)

      - For SL:
          * peaks_by_channel[primary_ch] brukes til å forhåndsfylle peaks i editoren.
    """
    fsa = entry["fsa"]
    peaks_by_channel = entry["peaks_by_channel"]
    primary_ch = entry["primary_peak_channel"]
    trace_channels = entry.get("trace_channels", [primary_ch])
    bp_min = float(entry["bp_min"])
    bp_max = float(entry["bp_max"])
    assay_name = entry.get("assay")

        # --- Initielle peaks for SL (kun primary channel) ---
    initial_peaks = []
    if assay_name == "SL":
        try:
            df0 = peaks_by_channel.get(primary_ch)
            if df0 is not None and not df0.empty:
                for _, row in df0.iterrows():
                    x = float(row.get("basepairs", np.nan))
                    y = float(row.get("peaks", np.nan))
                    if not (np.isfinite(x) and np.isfinite(y)):
                        continue
                    initial_peaks.append({"x": x, "y": y, "active": True})
        except Exception as ex:
            print_warning(
                f"[SL_PEAKS_INIT] Klarte ikke hente initielle peaks for "
                f"{fsa.file_name} ({primary_ch}): {ex}"
            )
            initial_peaks = []

    initial_peaks_json = json.dumps(initial_peaks)


    # --- 1) Valgfri tvungen y-maks (for kombinasjonsplott) ---
    forced_ymax = entry.get("forced_ymax", None)
    if forced_ymax is None:
        # støtt også gammelt navn 'force_ymax'
        forced_ymax = entry.get("force_ymax", None)
    try:
        if forced_ymax is not None:
            forced_ymax = float(forced_ymax)
        else:
            forced_ymax = None
    except Exception:
        forced_ymax = None

    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return None

    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return None

    time_all = raw_df["time"].astype(int).to_numpy()
    bp_all = raw_df["basepairs"].to_numpy()

    # --- 2) Hvilke kanaler finnes faktisk i FSA-filen? ---
    available = [k for k in fsa.fsa.keys() if k.startswith("DATA")]
    channels_to_plot = [ch for ch in trace_channels if ch in available]

    # Fallback: hvis trace_channels er tomme/ikke finnes, prøv primary_ch
    if not channels_to_plot:
        if primary_ch in fsa.fsa:
            channels_to_plot = [primary_ch]
        else:
            return None

    # Felles x-akse (bp) basert på første kanal
    first_ch = channels_to_plot[0]
    trace_first = np.asarray(fsa.fsa[first_ch])
    mask = (time_all >= 0) & (time_all < len(trace_first))
    if not np.any(mask):
        return None

    bp_trace = bp_all[mask]

    fig = go.Figure()

    # --- 3) Tegn alle spor, og beregn auto-y ---
    ymax_auto_primary = 0.0
    ymax_auto_all = 0.0

    # Predefiner bp-vindu-mask for å slippe å lage den for hver kanal
    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        win_bp = np.zeros_like(bp_trace, dtype=bool)
        for a, b in ASSAY_REFERENCE_RANGES[assay_name]:
            win_bp |= (bp_trace >= float(a)) & (bp_trace <= float(b))
    else:
        win_bp = (bp_trace >= bp_min) & (bp_trace <= bp_max)

    for ch in channels_to_plot:
        full_trace = np.asarray(fsa.fsa[ch])

        # Rask baseline: blokkvis lav-percentil + interpolasjon
        baseline = estimate_running_baseline(
            full_trace,
            bin_size=200,   # juster opp/ned for fart vs. smoothness
            quantile=0.10,
        )
        full_corr = full_trace - baseline
        full_corr[full_corr < 0] = 0.0


        # Zoomet variant: samme time_all[mask] som før
        y_corr = full_corr[time_all[mask]]
        if y_corr.size == 0:
            continue

        color = CHANNEL_COLORS.get(ch, DEFAULT_TRACE_COLOR)
        fig.add_trace(
            go.Scatter(
                x=bp_trace,
                y=y_corr,
                mode="lines",
                name=f"{ch} trace",
                line=dict(width=1, color=color),
                hoverinfo="x+y",
            )
        )

        # Begrens til bp-vindu for auto-y
        if np.any(win_bp):
            y_win = y_corr[win_bp]
        else:
            y_win = y_corr

        if y_win.size > 0 and np.any(np.isfinite(y_win)):
            local_max = float(np.nanmax(y_win))

            # maks over alle kanaler i denne entryen
            if local_max > ymax_auto_all:
                ymax_auto_all = local_max

            # maks for primærkanal
            if ch == primary_ch and local_max > ymax_auto_primary:
                ymax_auto_primary = local_max



    # --- 4) Velg endelig ymax ---
    if forced_ymax is not None and forced_ymax > 0:
        # Kombinasjonsfigurer: felles y-akse bestemt utenfor
        ymax = forced_ymax
    else:
        # Vanlige figurer
        multi_channel_assays = {
            "TCRgA", "TCRgB", "TCRg", "TCRγA", "TCRγB", "TCRγ",
            "TCRbA", "TCRbB", "TCRbC", "TCRβA", "TCRβB", "TCRβC",
        }
        if assay_name in multi_channel_assays:
            base = ymax_auto_all
        else:
            base = ymax_auto_primary or ymax_auto_all  # fallback hvis primary er 0


        if base <= 0:
            ymax = 1000.0
        else:
            ymax = base

    if ymax <= 0:
        ymax = 1000.0

    # Index til peaks-trace = etter alle sportracene
    peaks_trace_index = len(channels_to_plot)

    # --- 5) Forhåndsfyll peaks for SL (kun primary channel) ---
    initial_peaks = []
    try:
        if assay_name == "SL":
            df0 = peaks_by_channel.get(primary_ch)
            if df0 is not None and not df0.empty:
                for _, row in df0.iterrows():
                    x = float(row.get("basepairs", np.nan))
                    y = float(row.get("peaks", np.nan))
                    if not (np.isfinite(x) and np.isfinite(y)):
                        continue
                    initial_peaks.append({"x": x, "y": y, "active": True})
    except Exception as ex:
        print_warning(f"[SL_PEAKS_INIT] Klarte ikke hente initielle peaks for {fsa.file_name}: {ex}")
        initial_peaks = []

    initial_peaks_json = json.dumps(initial_peaks)

    # --- 6) Tom peaks-trace i fig-data (JS fyller inn punktene) ---
    fig.add_trace(
        go.Scatter(
            x=[],
            y=[],
            mode="markers",
            name="Peaks",
            marker=dict(
                size=8,
                color="red",
                opacity=1.0,
                line=dict(color="black", width=1),
            ),
            hovertemplate="bp=%{x:.2f}<br>height=%{y:.0f}<extra></extra>",
        )
    )

    # --- 7) Referanse-shapes ---
    shapes = []
    if assay_name and assay_name in ASSAY_REFERENCE_RANGES:
        for (a, b) in ASSAY_REFERENCE_RANGES[assay_name]:
            shapes.append(
                dict(
                    type="rect",
                    x0=float(a),
                    x1=float(b),
                    y0=0,
                    y1=1,
                    xref="x",
                    yref="paper",
                    fillcolor="rgba(235,232,203,0.5)",
                    line_width=0,
                )
            )
    else:
        shapes.append(
            dict(
                type="rect",
                x0=float(bp_min),
                x1=float(bp_max),
                y0=0,
                y1=1,
                xref="x",
                yref="paper",
                fillcolor="rgba(235,232,203,0.15)",
                line_width=0,
            )
        )

    sample_id = f"{fsa.file_name}_{primary_ch}"
    nice_title = f"{assay_name} – {sample_id}" if assay_name else sample_id

    fig.update_layout(
        title=nice_title,
        xaxis_title="Basepairs (bp)",
        yaxis_title="RFU",
        height=420,
        margin=dict(l=60, r=30, t=40, b=40),
        shapes=shapes,
        paper_bgcolor="white",
        plot_bgcolor="white",
        clickmode="event",
        showlegend=True,
        template="simple_white",
        font=dict(family="Inter, -apple-system, sans-serif", color="#0f172a"),
        hoverlabel=dict(bgcolor="white", font_size=12, font_family="Inter, -apple-system, sans-serif"),
    )

    # Y-akse og X-akse
    fig.update_yaxes(
        range=[0.0, ymax * 1.1],
        gridcolor="#f1f5f9",
        zerolinecolor="#cbd5e1"
    )
    fig.update_xaxes(
        range=[bp_min, bp_max],
        gridcolor="#f1f5f9",
        zerolinecolor="#cbd5e1"
    )

    fig_json = json.dumps(fig.to_plotly_json())

    safe_id = (
        sample_id.replace(" ", "_")
        .replace(".", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
    )
    div_id = f"peakplot_{safe_id}_{uuid.uuid4().hex}"

    # Merk: dobbeltklammer {{ }} er for f-string escaping i Python
    html_fragment = f"""
<div id="{div_id}" class="peak-editor-block"></div>
<script type="text/javascript">
(function() {{
  var fig = {fig_json};
  var initialPeaks = {initial_peaks_json};
  var gd = document.getElementById("{div_id}");
  if (!gd) return;

  // Dynamisk index til peaks-trace (etter alle spor-tracene)
  var peaksTraceIndex = {peaks_trace_index};

  Plotly.newPlot(gd, fig.data, fig.layout).then(function(g) {{

    var baseShapes = (g.layout.shapes || []).slice();
    var baseAnnots = (g.layout.annotations || []).slice();

    // peaks = (x, y, active) – starter fra initialPeaks for SL, ellers tom
    var peaks = (initialPeaks && Array.isArray(initialPeaks))
        ? initialPeaks.slice()
        : [];


    function nearestPeakIdx(xClick) {{
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

    function rebuild() {{

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
      }}, [peaksTraceIndex]);

      var ann = [];
      for (var i = 0; i < peaks.length; i++) {{
        var p = peaks[i];
        if (!p.active) continue;
        ann.push({{
          x: p.x,
          y: p.y * 1.03,
          xref: "x",
          yref: "y",
          text: p.x.toFixed(1),
          showarrow: false,
          font: {{ size: 9, color: "#222" }},
          xanchor: "left",
          yanchor: "bottom"
        }});
      }}

      Plotly.relayout(g, {{
        shapes: baseShapes,
        annotations: baseAnnots.concat(ann)
      }});
    }}

    // Tegn initielle peaks hvis vi har noen (typisk SL)
    if (peaks.length) {{
      rebuild();
    }}

    gd.on("plotly_click", function(ev) {{
      if (!ev.points || !ev.points.length) return;

      var pt = ev.points[0];
      var xVal = pt.x;
      var yVal = pt.y;
      var isShift = ev.event && ev.event.shiftKey;

      if (isShift) {{
        // Slett nærmeste peak
        var idxDel = nearestPeakIdx(xVal);
        if (idxDel >= 0) {{
          peaks.splice(idxDel, 1);
          rebuild();
        }}
        return;
      }}

      // Sjekk om vi treffer en eksisterende peak (lite bp-vindu rundt)
      var idx = nearestPeakIdx(xVal);
      if (idx >= 0 && Math.abs(peaks[idx].x - xVal) < 0.4) {{
        // Toggle aktiv / inaktiv
        peaks[idx].active = !peaks[idx].active;
        rebuild();
        return;
      }}

      // Ellers: legg til ny peak (aktiv)
      peaks.push({{ x: xVal, y: yVal, active: true }});
      rebuild();
    }});
  }});
}})();
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
        else:
            shapes.append(
                dict(
                    type="rect",
                    x0=float(bp_min),
                    x1=float(bp_max),
                    y0=0,
                    y1=1,
                    xref="x",
                    yref="paper",
                    fillcolor="rgba(235,232,203,0.25)",
                    line_width=0,
                )
            )

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

        # JS for akkurat denne editoren – uendret logikk, bare limt inn som før
        html_parts.append(f"""
<script type="text/javascript">
(function() {{
  var fig = {fig_json};
  var divId = "{div_id}";
  var gd = document.getElementById(divId);
  if (!gd) {{
    console.error("Fant ikke div", divId);
    return;
  }}

  Plotly.newPlot(gd, fig.data, fig.layout).then(function(g) {{
    var peaksXs = [];
    var peaksYs = [];

    function redrawPeaks() {{
      var texts = peaksXs.map(function(x) {{
        return (typeof x === "number") ? x.toFixed(1) : String(x);
      }});

      Plotly.restyle(g, {{
        x: [peaksXs],
        y: [peaksYs],
        text: [texts]
      }}, [1]); // peaks-trace er index 1

      var arr = peaksXs.map(function(x, i) {{
        return {{
          bp: x,
          height: peaksYs[i]
        }};
      }});
      var pre = document.getElementById("{div_id}_peaks_json");
      if (pre) {{
        pre.textContent = JSON.stringify(arr, null, 2);
      }}
    }}

    function findNearestPeakIdx(xClick) {{
      if (!peaksXs.length) return -1;
      var bestIdx = 0;
      var bestDist = Math.abs(peaksXs[0] - xClick);
      for (var i = 1; i < peaksXs.length; i++) {{
        var d = Math.abs(peaksXs[i] - xClick);
        if (d < bestDist) {{
          bestDist = d;
          bestIdx = i;
        }}
      }}
      return bestIdx;
    }}

    gd.on("plotly_click", function(ev) {{
      if (!ev || !ev.points || !ev.points.length) return;
      var pt = ev.points[0];
      var isShift = !!(ev.event && ev.event.shiftKey);

      var xVal = pt.x;
      var yVal = pt.y;

      if (isShift && peaksXs.length) {{
        var idx = findNearestPeakIdx(xVal);
        if (idx >= 0) {{
          peaksXs.splice(idx, 1);
          peaksYs.splice(idx, 1);
          redrawPeaks();
        }}
        return;
      }}

      // vanlig klikk -> legg til peak
      peaksXs.push(xVal);
      peaksYs.push(yVal);
      redrawPeaks();
    }});
  }});
}})();
</script>
""")

    return "\n".join(html_parts)
