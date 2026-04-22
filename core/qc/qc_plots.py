"""
HemaFrag QC — Interactive Plotly QC plot builder.
"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from html import escape
import json

from core.qc.qc_rules import QCRules
from core.qc.qc_markers import (
    control_id_from_filename,
    find_peak_near_bp,
    find_peak_near_bp_with_fallback,
    markers_for_entry,
)
import core.assay_config as assay_config
from core.assay_config import (
    CHANNEL_COLORS,
    DEFAULT_TRACE_COLOR,
    merged_analysis_attr,
    reference_shade_rgba,
)
from core.analysis import estimate_running_baseline


def _assay_reference_ranges() -> dict:
    return merged_analysis_attr("ASSAY_REFERENCE_RANGES")


def _qc_json_dumps_compact(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), ensure_ascii=False)


def _get_qc_entry_cache(entry: dict) -> dict:
    cache = entry.get("_qc_plot_cache")
    if not isinstance(cache, dict):
        cache = {}
        entry["_qc_plot_cache"] = cache
    return cache


def _get_qc_fsa_cache(fsa: object) -> dict:
    cache = getattr(fsa, "_qc_plot_cache", None)
    if not isinstance(cache, dict):
        cache = {}
        setattr(fsa, "_qc_plot_cache", cache)
    return cache


def _get_qc_axis_arrays(fsa: object) -> dict | None:
    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        return None
    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        return None

    cache = _get_qc_fsa_cache(fsa)
    cache_key = ("axis_arrays", id(raw_df), tuple(raw_df.columns))
    cached = cache.get("axis_arrays")
    if isinstance(cached, dict) and cached.get("key") == cache_key:
        return cached["value"]

    value = {
        "time_all": raw_df["time"].astype(int).to_numpy(),
        "bp_all": raw_df["basepairs"].to_numpy(),
        "available_channels": tuple(k for k in fsa.fsa.keys() if k.startswith("DATA")),
    }
    cache["axis_arrays"] = {"key": cache_key, "value": value}
    return value


def _get_qc_trace_array(fsa: object, channel: str) -> np.ndarray:
    cache = _get_qc_fsa_cache(fsa)
    trace_arrays = cache.setdefault("trace_arrays", {})
    cached = trace_arrays.get(channel)
    current = getattr(fsa, "fsa", {}).get(channel)
    current_id = id(current)
    if isinstance(cached, dict) and cached.get("source_id") == current_id:
        return cached["value"]

    value = np.asarray(current, dtype=float)
    trace_arrays[channel] = {"source_id": current_id, "value": value}
    return value


def _marker_rules_signature(rules: QCRules) -> tuple[float, float, float, float, float]:
    return (
        float(getattr(rules, "min_r2_ok", 0.0)),
        float(getattr(rules, "min_r2_warn", 0.0)),
        float(getattr(rules, "nk_ymax_floor", 0.0)),
        float(getattr(rules, "sample_peak_window_bp", 0.0)),
        float(getattr(rules, "sample_peak_window_bp_fallback", 0.0)),
    )


def _get_qc_display_trace(entry: dict, channel: str, assay: str | None) -> np.ndarray:
    fsa = entry["fsa"]
    cache = _get_qc_entry_cache(entry)
    display_cache = cache.setdefault("display_traces", {})
    trace = _get_qc_trace_array(fsa, channel)
    cache_key = (assay, channel, id(trace), trace.shape)
    cached = display_cache.get(channel)
    if isinstance(cached, dict) and cached.get("key") == cache_key:
        return cached["value"]

    value = _baseline_correct_trace_for_display(
        trace,
        channel=channel,
        assay=assay,
    )
    display_cache[channel] = {"key": cache_key, "value": value}
    return value


def _get_cached_marker_results(
    entry: dict,
    rules: QCRules,
    marker_specs: list[dict],
    *,
    fsa: object,
    primary_ch: str | None,
) -> list[dict]:
    cache = _get_qc_entry_cache(entry)
    specs_signature = tuple(
        (
            str(m.get("name") or ""),
            str(m.get("kind") or ""),
            str(m.get("channel") or ""),
            float(m.get("expected_bp", 0.0) or 0.0),
            float(m.get("window_bp", 0.0) or 0.0),
        )
        for m in marker_specs
    )
    cache_key = (primary_ch, _marker_rules_signature(rules), specs_signature)
    cached = cache.get("marker_results")
    if isinstance(cached, dict) and cached.get("key") == cache_key:
        results = [dict(item) for item in cached["value"]]
        entry["qc_marker_results"] = results
        return results

    results: list[dict] = []
    for m in marker_specs:
        ch = primary_ch if m["channel"] == "primary" else m["channel"]
        if m["kind"] == "sample":
            res = find_peak_near_bp_with_fallback(
                fsa=fsa,
                channel=ch,
                target_bp=float(m["expected_bp"]),
                window_bp=float(m["window_bp"]),
                fallback_window_bp=float(getattr(rules, "sample_peak_window_bp_fallback", m["window_bp"])),
                baseline_correct=True,
            )
        else:
            res = find_peak_near_bp(
                fsa=fsa,
                channel=ch,
                target_bp=float(m["expected_bp"]),
                window_bp=float(m["window_bp"]),
                baseline_correct=True,
            )
        res2 = dict(m)
        res2.update(res)
        results.append(res2)

    cache["marker_results"] = {"key": cache_key, "value": [dict(item) for item in results]}
    entry["qc_marker_results"] = results
    return results


def _is_ladder_channel(channel: str | None) -> bool:
    return channel in {"DATA4", "DATA105"}


def _smooth_signal(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1 or values.size == 0:
        return np.asarray(values, dtype=float)
    kernel = np.ones(int(window), dtype=float) / float(window)
    return np.convolve(np.asarray(values, dtype=float), kernel, mode="same")


def _baseline_correct_trace_for_display(
    trace: np.ndarray,
    *,
    channel: str | None,
    assay: str | None,
) -> np.ndarray:
    """Baseline-correct a trace for display without letting low-frequency drift dominate."""
    full_trace = np.asarray(trace, dtype=float)
    if full_trace.size == 0:
        return full_trace

    try:
        if assay == "SL":
            baseline = estimate_running_baseline(
                full_trace,
                bin_size=5000,
                quantile=0.01,
                use_arpls=False,
            )
            corrected = full_trace - baseline
            corrected[corrected < 0] = 0.0
            return corrected

        if _is_ladder_channel(channel):
            primary_baseline = estimate_running_baseline(
                full_trace,
                bin_size=max(400, min(full_trace.size, 1200)),
                quantile=0.20,
                use_arpls=False,
            )
            corrected = full_trace - primary_baseline
            corrected[corrected < 0] = 0.0

            residual_baseline = estimate_running_baseline(
                corrected,
                bin_size=max(150, min(corrected.size, 500)),
                quantile=0.05,
                use_arpls=False,
            )
            corrected = corrected - residual_baseline
            corrected[corrected < 0] = 0.0

            smooth_window = max(5, min((corrected.size // 20) * 2 + 1, 151))
            broad_trend = _smooth_signal(corrected, smooth_window)
            corrected = corrected - broad_trend
            corrected[corrected < 0] = 0.0
            return corrected

        baseline = estimate_running_baseline(full_trace, bin_size=200, quantile=0.10)
        corrected = full_trace - baseline
        corrected[corrected < 0] = 0.0
        return corrected
    except Exception:
        return full_trace



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

    axis_arrays = _get_qc_axis_arrays(fsa)
    if not axis_arrays:
        return None

    time_all = axis_arrays["time_all"]
    bp_all = axis_arrays["bp_all"]
    entry_cache = _get_qc_entry_cache(entry)

    # Finn hvilke kanaler som faktisk finnes
    available = list(axis_arrays["available_channels"])
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
    trace_first = _get_qc_trace_array(fsa, first_ch)
    mask = (time_all >= 0) & (time_all < len(trace_first))
    if not np.any(mask):
        return None
    bp_trace = bp_all[mask]

    # Referanse-vindu (for auto-ymax)
    ref_ranges_by_assay = _assay_reference_ranges()
    if assay and assay in ref_ranges_by_assay:
        win_bp = np.zeros_like(bp_trace, dtype=bool)
        for a, b in ref_ranges_by_assay[assay]:
            win_bp |= (bp_trace >= float(a)) & (bp_trace <= float(b))
    else:
        win_bp = (bp_trace >= bp_min) & (bp_trace <= bp_max)

    # Bare beregn/legg til markører om vi faktisk har specs (dvs PK)
    marker_specs = markers_for_entry(entry, rules)  # tom for ikke-PK
    fragment_cache_key = (
        primary_ch,
        tuple(channels_to_plot),
        float(bp_min),
        float(bp_max),
        str(assay or ""),
        str(entry.get("ladder") or ""),
        str(ctrl),
        _marker_rules_signature(rules),
        tuple(
            (
                str(m.get("name") or ""),
                str(m.get("kind") or ""),
                str(m.get("channel") or ""),
                float(m.get("expected_bp", 0.0) or 0.0),
                float(m.get("window_bp", 0.0) or 0.0),
            )
            for m in marker_specs
        ),
    )
    cached_fragment = entry_cache.get("html_fragment")
    if isinstance(cached_fragment, dict) and cached_fragment.get("key") == fragment_cache_key:
        marker_results = cached_fragment.get("marker_results")
        if isinstance(marker_results, list):
            entry["qc_marker_results"] = [dict(item) for item in marker_results]
        return cached_fragment.get("value")

    fig = go.Figure()

    # Tegn traces (baseline-korrigert slik som i master). [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    ymax_auto_primary = 0.0
    ymax_auto_all = 0.0

    scale_channels = [ch for ch in channels_to_plot if ch not in ("DATA4", "DATA105")]

    for ch in channels_to_plot:
        full_corr = _get_qc_display_trace(entry, ch, assay)

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
    marker_results = _get_cached_marker_results(
        entry,
        rules,
        marker_specs,
        fsa=fsa,
        primary_ch=primary_ch,
    )

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
            if mr.get("search_mode") == "fallback":
                txt += f"<br>Fallback window used: ±{float(mr.get('search_window_bp', 0.0)):.1f} bp"
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
    if assay and assay in ref_ranges_by_assay:
        for (a, b) in ref_ranges_by_assay[assay]:
            shapes.append(dict(
                type="rect",
                x0=float(a), x1=float(b),
                y0=0, y1=1, xref="x", yref="paper",
                fillcolor=reference_shade_rgba(),
                line_width=0,
                layer="above",
                opacity=1.0,
            ))
    else:
        shapes.append(dict(
            type="rect",
            x0=float(bp_min), x1=float(bp_max),
            y0=0, y1=1, xref="x", yref="paper",
            fillcolor=reference_shade_rgba(),
            line_width=0,
            layer="above",
            opacity=1.0,
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
        template="simple_white",
    )


    # Y-akse: NK starter på 250
    # NK: behold auto-skalering, men unngå "for mye zoom inn"
    if ctrl == "NK":
        ymax = max(float(ymax), float(rules.nk_ymax_floor))

    fig.update_yaxes(range=[ymin, ymax * 1.10], showgrid=False, zeroline=False)

    x_min = bp_min
    x_max = bp_max

    if marker_specs:
        exp_bps = [float(m["expected_bp"]) for m in marker_specs]
        actual_search_margins = {
            str(mr.get("name") or ""): float(mr.get("search_window_bp", 0.0) or 0.0)
            for mr in marker_results
            if mr.get("ok")
        }
        margins = [
            max(
                float(m.get("window_bp", 0)),
                actual_search_margins.get(str(m.get("name") or ""), 0.0),
                8.0,
            )
            for m in marker_specs
        ]
        margin = max(margins) if margins else 8.0

        x_min = min(x_min, min(exp_bps) - margin)
        x_max = max(x_max, max(exp_bps) + margin)

    fig.update_xaxes(range=[x_min, x_max], showgrid=False, zeroline=False)

    fig_json = _qc_json_dumps_compact(fig.to_plotly_json())
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
  var assayName = {_qc_json_dumps_compact(assay)};

  var plotConfig = {{ responsive: true, displaylogo: false }};
  var mountPlot = (window.ReportPlotManager && window.ReportPlotManager.mountPlot)
    ? window.ReportPlotManager.mountPlot(gd, fig.data, fig.layout, plotConfig)
    : Plotly.newPlot(gd, fig.data, fig.layout, plotConfig);

  mountPlot.then(function(g) {{
    function decodePlotlyArray(val) {{
      if (Array.isArray(val)) return val;
      if (ArrayBuffer.isView(val)) return Array.from(val);
      if (!val || typeof val !== "object" || typeof val.length !== "number") return [];
      try {{ return Array.from(val); }} catch(e) {{ return []; }}
    }}

    var traceXYCache = g.data.map(function(t) {{
      return {{ x: decodePlotlyArray(t.x), y: decodePlotlyArray(t.y) }};
    }});

    function solve3x3(mat, vec) {{
      var det = mat[0][0]*(mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) -
                mat[0][1]*(mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0]) +
                mat[0][2]*(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
      if (Math.abs(det) < 1e-12) return null;
      function getDet(m, colIdx, v) {{
        var nm = m.map(function(row) {{ return row.slice(); }});
        for (var i = 0; i < 3; i++) nm[i][colIdx] = v[i];
        return nm[0][0]*(nm[1][1]*nm[2][2] - nm[1][2]*nm[2][1]) -
               nm[0][1]*(nm[1][0]*nm[2][2] - nm[1][2]*nm[2][0]) +
               nm[0][2]*(nm[1][0]*nm[2][1] - nm[1][1]*nm[2][0]);
      }}
      return [getDet(mat,0,vec)/det, getDet(mat,1,vec)/det, getDet(mat,2,vec)/det];
    }}

    function computeGaussianArea(xCenter, traceIdx) {{
      var data = traceXYCache[traceIdx] || traceXYCache[0] || {{ x: [], y: [] }};
      var hw = (assayName === "SL") ? 20.0 : 5.0;
      var pts = [];
      for (var i = 0; i < data.x.length; i++) {{
        if (Math.abs(data.x[i] - xCenter) <= hw && data.y[i] > 0.01) {{
          pts.push({{ x: data.x[i], lny: Math.log(data.y[i]) }});
        }}
      }}
      if (pts.length < 3) return 0;
      var sx4=0, sx3=0, sx2=0, sx1=0, n=pts.length, sxy2=0, sxy1=0, sy=0;
      for (var i=0; i<n; i++) {{
        var xi = i;
        var yi = pts[i].lny;
        var xi2 = xi*xi;
        sx4 += xi2*xi2; sx3 += xi2*xi; sx2 += xi2; sx1 += xi;
        sxy2 += xi2*yi; sxy1 += xi*yi; sy += yi;
      }}
      var sol = solve3x3([[sx4, sx3, sx2], [sx3, sx2, sx1], [sx2, sx1, n]], [sxy2, sxy1, sy]);
      if (!sol || sol[0] >= 0) return 0;
      var sigma2 = -1.0 / (2.0 * sol[0]);
      var mu = sol[1] * sigma2;
      var amp = Math.exp(sol[2] + (mu*mu)/(2.0*sigma2));
      if (mu < -1 || mu > (n + 1)) return 0;
      return amp * Math.sqrt(sigma2) * Math.sqrt(2.0 * Math.PI);
    }}

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
      var texts = peaks.map(function(p) {{ 
        return p.active ? (p.x.toFixed(1) + (p.area ? " (A=" + p.area.toFixed(0) + ")" : "")) : ""; 
      }});

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

      var area = computeGaussianArea(xVal, pt.curveNumber || 0);
      peaks.push({{ x: xVal, y: yVal, active: true, area: area }});
      rebuild();
    }});
  }});
}})();
</script>
"""
    entry_cache["html_fragment"] = {
        "key": fragment_cache_key,
        "value": html_fragment,
        "marker_results": [dict(item) for item in marker_results],
    }
    return html_fragment
