"""
Fraggler QC — Interactive Plotly QC plot builder.
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from html import escape
import json

from core.qc.qc_rules import QCRules
from core.qc.qc_markers import markers_for_entry, find_peak_near_bp, control_id_from_filename
from core.assay_config import (
    ASSAY_REFERENCE_RANGES,
    CHANNEL_COLORS,
    DEFAULT_TRACE_COLOR,
)
from core.analysis import estimate_running_baseline



def build_interactive_peak_plot_for_entry_qc(entry: dict, rules: QCRules) -> str | None:
    import plotly.graph_objects as go
    import uuid

    fsa = entry["fsa"]
    assay = entry.get("assay")
    primary_ch = entry.get("primary_peak_channel")
    trace_channels = entry.get("trace_channels", [primary_ch])
    bp_min = float(entry.get("bp_min", 0))
    bp_max = float(entry.get("bp_max", 0))

    # NK y-min krav
    ctrl = control_id_from_filename(fsa.file_name)
    is_pk = ctrl in {"PK", "PK1", "PK2"}
    # NK: ikke fast ymin, men vi gulver ymax senere
    ymin = 0.0

    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return None
    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return None

    time_all = raw_df["time"].astype(int).to_numpy()
    bp_all = raw_df["basepairs"].to_numpy()

    # Finn hvilke kanaler som faktisk finnes
    available = [k for k in fsa.fsa.keys() if k.startswith("DATA")]
    channels_to_plot = [ch for ch in trace_channels if ch in available]
    if not channels_to_plot:
        if primary_ch in fsa.fsa:
            channels_to_plot = [primary_ch]
        else:
            return None

# Vis ladder-kanalen som trace kun for PK (for å unngå rot for NK/RK)
    if is_pk:
        ladder_ch = "DATA4" if entry.get("ladder") == "ROX" else "DATA105"  
        if ladder_ch in available and ladder_ch not in channels_to_plot:
            channels_to_plot.append(ladder_ch)

    # Felles x-akse basert på første kanal
    first_ch = channels_to_plot[0]
    trace_first = np.asarray(fsa.fsa[first_ch])
    mask = (time_all >= 0) & (time_all < len(trace_first))
    if not np.any(mask):
        return None
    bp_trace = bp_all[mask]

    # Referanse-vindu (for auto-ymax)
    if assay and assay in ASSAY_REFERENCE_RANGES:
        win_bp = np.zeros_like(bp_trace, dtype=bool)
        for a, b in ASSAY_REFERENCE_RANGES[assay]:
            win_bp |= (bp_trace >= float(a)) & (bp_trace <= float(b))
    else:
        win_bp = (bp_trace >= bp_min) & (bp_trace <= bp_max)

    fig = go.Figure()

    # Tegn traces (baseline-korrigert slik som i master). [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    ymax_auto_primary = 0.0
    ymax_auto_all = 0.0

    scale_channels = [ch for ch in channels_to_plot if ch not in ("DATA4", "DATA105")]

    for ch in channels_to_plot:
        full_trace = np.asarray(fsa.fsa[ch]).astype(float)

        try:
            baseline = estimate_running_baseline(full_trace, bin_size=200, quantile=0.10)
            full_corr = full_trace - baseline
            full_corr[full_corr < 0] = 0.0
        except Exception:
            full_corr = full_trace

        y_corr = full_corr[time_all[mask]]

        color = CHANNEL_COLORS.get(ch, "#1f77b4")
        fig.add_trace(go.Scatter(
            x=bp_trace, y=y_corr, mode="lines",
            name=f"{ch} trace", line=dict(width=1, color=color),
            hoverinfo="x+y"
        ))

        if np.any(win_bp):
            y_win = y_corr[win_bp]
        else:
            y_win = y_corr

        if y_win.size > 0 and np.any(np.isfinite(y_win)):
            local_max = float(np.nanmax(y_win))

            # Oppdater y-skalering kun for "scale_channels"
            if ch in scale_channels:
                ymax_auto_all = max(ymax_auto_all, local_max)
                if ch == primary_ch:
                    ymax_auto_primary = max(ymax_auto_primary, local_max)

    multi_channel_assays = {
        "TCRgA", "TCRgB",
        "TCRbA", "TCRbB", "TCRbC",
        "IGK",  # to kanaler i config
        # legg til flere hvis du vil
    }

    # Bruk master-lik oppførsel: multi-kanal assays -> ymax fra alle kanaler
    if assay in multi_channel_assays:
        base = ymax_auto_all
    else:
        base = ymax_auto_primary if ymax_auto_primary > 0 else ymax_auto_all

    ymax = base if (base and base > 0) else 1000.0

    # Hvis forced_ymax settes (kombinasjoner): bruk den
    forced_ymax = entry.get("forced_ymax", entry.get("force_ymax", None))
    try:
        if forced_ymax is not None:
            forced_ymax = float(forced_ymax)
            if forced_ymax > 0:
                ymax = forced_ymax
    except Exception:
        pass

    # -----------------------------
    # MARKØRER: forventede sample-peaks + ladder-peaks
    # -----------------------------
    marker_specs = markers_for_entry(entry, rules)  # tom for ikke-PK
    marker_results = []

    # Bare beregn/legg til markører om vi faktisk har specs (dvs PK)
    if marker_specs:
        for m in marker_specs:
            ch = primary_ch if m["channel"] == "primary" else m["channel"]
            res = find_peak_near_bp(
                fsa=fsa,
                channel=ch,
                target_bp=float(m["expected_bp"]),
                window_bp=float(m["window_bp"]),
                baseline_correct=True
            )
            res2 = dict(m)
            res2.update(res)
            marker_results.append(res2)

    # Lagre til Excel (tom liste for NK/RK)
    entry["qc_marker_results"] = marker_results

    # --- Hvis markører finnes: legg inn marker-traces ---
    n_extra_traces = 0
    if marker_results:
        xs_sample, ys_sample, text_sample = [], [], []
        xs_ladder, ys_ladder, text_ladder = [], [], []
 
        for mr in marker_results:
            if not mr.get("ok"):
                continue
            delta = float(mr["found_bp"]) - float(mr["expected_bp"])
            txt = (
                f"{mr['name']}: exp {mr['expected_bp']:.1f} → {mr['found_bp']:.2f} "
                f"(Δ {delta:+.2f})<br>H={mr['height']:.0f}, A={mr['area']:.0f}"
            )
            if mr["kind"] == "ladder":
                xs_ladder.append(mr["found_bp"]); ys_ladder.append(mr["height"]); text_ladder.append(txt)
            else:
                xs_sample.append(mr["found_bp"]); ys_sample.append(mr["height"]); text_sample.append(txt)

        fig.add_trace(go.Scatter(
            x=xs_sample, y=ys_sample, mode="markers",
            name="QC markers (sample)",
            marker=dict(symbol="diamond", size=10, color="#7b2cbf", line=dict(color="black", width=1)),
            hovertext=text_sample, hoverinfo="text"
        ))
        fig.add_trace(go.Scatter(
            x=xs_ladder, y=ys_ladder, mode="markers",
            name="QC markers (ladder)",
            marker=dict(symbol="diamond", size=10, color="#f59f00", line=dict(color="black", width=1)),
            hovertext=text_ladder, hoverinfo="text"
        ))
        n_extra_traces = 2

    # Peaks trace kommer etter trace-kanaler + evt 2 marker-traces
    peaks_trace_index = len(channels_to_plot) + n_extra_traces

    fig.add_trace(go.Scatter(
        x=[], y=[], mode="markers",
        name="Peaks",
        marker=dict(size=8, color="red", opacity=1.0, line=dict(color="black", width=1)),
        hovertemplate="bp=%{x:.2f}<br>height=%{y:.0f}<extra></extra>",
    ))

    # -----------------------------
    # Shapes: referanse-shading + vertikale markerlinjer
    # -----------------------------
    shapes = []

    # Referanse-shading (samme referansevinduer som master). [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    if assay and assay in ASSAY_REFERENCE_RANGES:
        for (a, b) in ASSAY_REFERENCE_RANGES[assay]:
            shapes.append(dict(
                type="rect",
                x0=float(a), x1=float(b),
                y0=0, y1=1, xref="x", yref="paper",
                fillcolor="rgba(235,232,203,0.25)",
                line_width=0
            ))
    else:
        shapes.append(dict(
            type="rect",
            x0=float(bp_min), x1=float(bp_max),
            y0=0, y1=1, xref="x", yref="paper",
            fillcolor="rgba(235,232,203,0.25)",
            line_width=0
        ))

    if marker_specs:
        for ms in marker_specs:
            col = "rgba(245,159,0,0.7)" if ms["kind"] == "ladder" else "rgba(123,44,191,0.55)"
            shapes.append(dict(
                type="line",
                x0=float(ms["expected_bp"]), x1=float(ms["expected_bp"]),
                y0=0, y1=1, xref="x", yref="paper",
                line=dict(color=col, width=1, dash="dot")
            ))
 
    # Layout


    # Layout
    sample_id = f"{fsa.file_name}_{primary_ch}"
    nice_title = f"{assay} – {sample_id}"

    fig.update_layout(
        title=nice_title,
        xaxis_title="Basepairs (bp)",
        yaxis_title="RFU",
        height=450,
        margin=dict(l=60, r=30, t=45, b=40),
        shapes=shapes,
        paper_bgcolor="white",
        plot_bgcolor="white",
        clickmode="event",
        showlegend=True,
    )


    # Y-akse: NK starter på 250
    # NK: behold auto-skalering, men unngå "for mye zoom inn"
    if ctrl == "NK":
        ymax = max(float(ymax), float(rules.nk_ymax_floor))

    fig.update_yaxes(range=[ymin, ymax * 1.10])

    x_min = bp_min
    x_max = bp_max

    if marker_specs:
        exp_bps = [float(m["expected_bp"]) for m in marker_specs]
        margins = [max(float(m.get("window_bp", 0)), 8.0) for m in marker_specs]
        margin = max(margins) if margins else 8.0

        x_min = min(x_min, min(exp_bps) - margin)
        x_max = max(x_max, max(exp_bps) + margin)

    fig.update_xaxes(range=[x_min, x_max])

    fig_json = json.dumps(fig.to_plotly_json())
    safe_id = (sample_id.replace(" ", "_").replace(".", "_").replace("/", "_").replace("\\", "_").replace(":", "_"))
    div_id = f"qc_peakplot_{safe_id}_{uuid.uuid4().hex}"

    html_fragment = f"""
<div id="{div_id}" class="peak-editor-block"></div>
<script type="text/javascript">
(function() {{
  var fig = {fig_json};
  var gd = document.getElementById("{div_id}");
  if (!gd) return;

  var peaksTraceIndex = {peaks_trace_index};

  Plotly.newPlot(gd, fig.data, fig.layout).then(function(g) {{
    var baseShapes = (g.layout.shapes || []).slice();
    var baseAnnots = (g.layout.annotations || []).slice();

    var peaks = [];

    function nearestPeakIdx(xClick) {{
      if (!peaks.length) return -1;
      var bestIdx = 0;
      var bestDist = Math.abs(peaks[0].x - xClick);
      for (var i = 1; i < peaks.length; i++) {{
        var d = Math.abs(peaks[i].x - xClick);
        if (d < bestDist) {{ bestDist = d; bestIdx = i; }}
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
          x: p.x, y: p.y * 1.03, xref: "x", yref: "y",
          text: p.x.toFixed(1), showarrow: false,
          font: {{ size: 9, color: "#222" }},
          xanchor: "left", yanchor: "bottom"
        }});
      }}

      Plotly.relayout(g, {{
        shapes: baseShapes,
        annotations: baseAnnots.concat(ann)
      }});
    }}

    gd.on("plotly_click", function(ev) {{
      if (!ev.points || !ev.points.length) return;
      var pt = ev.points[0];
      var xVal = pt.x;
      var yVal = pt.y;
      var isShift = ev.event && ev.event.shiftKey;

      if (isShift) {{
        var idxDel = nearestPeakIdx(xVal);
        if (idxDel >= 0) {{
          peaks.splice(idxDel, 1);
          rebuild();
        }}
        return;
      }}

      var idx = nearestPeakIdx(xVal);
      if (idx >= 0 && Math.abs(peaks[idx].x - xVal) < 0.4) {{
        peaks[idx].active = !peaks[idx].active;
        rebuild();
        return;
      }}

      peaks.push({{ x: xVal, y: yVal, active: true }});
      rebuild();
    }});
  }});
}})();
</script>
"""
    return html_fragment
