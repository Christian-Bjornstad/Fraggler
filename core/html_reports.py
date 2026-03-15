"""
Fraggler Diagnostics — DIT HTML Reports & SL Quality Interpretation.

Builds per-patient DIT HTML reports with embedded interactive Plotly figures.
"""
from __future__ import annotations

import re
import uuid
import json
from pathlib import Path
from collections import defaultdict
from html import escape
from datetime import datetime
import pandas as pd

from core.analyses.registry import get_active_analysis_name

import numpy as np

from fraggler.fraggler import print_green, print_warning

from core.assay_config import (
    ASSAY_CONFIG,
    ASSAY_DISPLAY_ORDER,
    ASSAY_REFERENCE_RANGES,
    ASSAY_REFERENCE_LABEL,
    CHANNEL_COLORS,
    DEFAULT_TRACE_COLOR,
    OUTDIR_NAME,
)
from core.plotly_offline import local_plotly_tag as _local_plotly_tag
from core.plotting_plotly import (
    compute_group_ymax_for_entries,
    build_interactive_peak_plot_for_entry,
)
from config import APP_SETTINGS


DIT_PATTERN = re.compile(r"(\d{2}OUM\d{5})")

REPORT_STYLE = """
<style>
/* ── Base Typography ── */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

body {
    font-family: 'Inter', -apple-system, sans-serif;
    margin: 1.5rem 2rem 5rem;
    background: #f4f7fb;
    color: #0f172a;
    line-height: 1.5;
}
h1 { font-size: 1.6rem; font-weight: 800; color: #0f172a; margin-bottom: 0.2rem; }
h2 { font-size: 1.15rem; font-weight: 700; color: #1e293b; margin-top: 2rem; margin-bottom: 0.5rem; padding-bottom: 6px; border-bottom: 2px solid #e2e8f0; }
h3 { font-size: 1rem; font-weight: 700; color: #4338ca; margin-top: 1rem; margin-bottom: 0.3rem; }
p  { margin-top: 0.2rem; margin-bottom: 0.4rem; color: #334155; }

/* ── Header banner ── */
.report-header {
    background: linear-gradient(135deg, #06b6d4 0%, #4338ca 100%);
    color: white;
    padding: 24px 28px;
    border-radius: 12px;
    margin-bottom: 24px;
    box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 8px 10px -6px rgba(0, 0, 0, 0.05);
}
.report-header h1 { color: white; font-size: 1.8rem; font-weight: 800; margin: 0 0 4px; letter-spacing: -0.6px; }
.report-header .meta { font-size: 0.9rem; font-weight: 500; color: rgba(255, 255, 255, 0.85); }

/* ── Tables ── */
table { 
    border-collapse: collapse; 
    margin-bottom: 1.2rem; 
    width: 100%; 
    background: white; 
    border-radius: 8px; 
    overflow: hidden; 
    box-shadow: 0 4px 12px rgba(0,0,0,0.03); 
}
th, td { border-bottom: 1px solid #e2e8f0; padding: 12px 14px; font-size: 0.85rem; }
th { 
    background: #f8fafc; 
    font-weight: 800; 
    color: #4338ca; 
    text-transform: uppercase; 
    font-size: 0.75rem; 
    letter-spacing: 0.8px; 
    border-bottom: 2px solid #e2e8f0;
}
tr:nth-child(even) td { background: #fafbfc; }
tr:hover td { background: #f0fdfa; /* Soft teal hover */ transition: background 0.15s ease; }

/* ── Cards ── */
.assay-block {
    padding: 18px 22px;
    margin-bottom: 20px;
    background: #ffffff;
    border-radius: 12px;
    border: 1px solid rgba(226, 232, 240, 0.8);
    box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.05), 0 8px 10px -6px rgba(0, 0, 0, 0.01);
    transition: box-shadow 0.3s ease;
}
.assay-block:hover { box-shadow: 0 14px 28px -5px rgba(0, 0, 0, 0.08); }
.sample-header { font-size: 0.9rem; font-weight: 700; margin-top: 0.4rem; color: #0f172a; }
.small { font-size: 0.85rem; color: #64748b; font-weight: 500; }
.peak-editor-block { margin-top: 0.5rem; margin-bottom: 1.2rem; border-radius: 8px; overflow: hidden; }
.combo-grid { display: block; }
.combo-item { margin-bottom: 1.2rem; }
/* ── Floating Print Button ── */
.print-fab {
    position: fixed;
    bottom: 28px;
    right: 28px;
    z-index: 9999;
    display: flex;
    gap: 10px;
    align-items: center;
    flex-direction: column;
}
.print-btn {
    background: linear-gradient(135deg, #06b6d4, #4338ca);
    color: white;
    border: none;
    border-radius: 50px;
    padding: 14px 26px;
    font-size: 14px;
    font-weight: 700;
    font-family: inherit;
    cursor: pointer;
    box-shadow: 0 6px 16px rgba(67, 56, 202, 0.4);
    transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
    display: flex;
    align-items: center;
    gap: 8px;
    white-space: nowrap;
    letter-spacing: 0.3px;
}
.print-btn:hover { transform: translateY(-3px); box-shadow: 0 10px 20px rgba(67, 56, 202, 0.5); filter: brightness(1.1); }
.print-btn:active { transform: translateY(0); }

/* ── Print Media ── */
@media print {
    .print-fab { display: none !important; }
    body { background: white; color: black; margin: 0; font-size: 11pt; line-height: 1.4; }
    .report-header { background: #4338ca; -webkit-print-color-adjust: exact; print-color-adjust: exact; border-radius: 0; margin: 0 0 16px; box-shadow: none; }
    h1 { font-size: 18pt; }
    h2 { font-size: 13pt; page-break-after: avoid; }
    h3 { font-size: 11pt; page-break-after: avoid; }
    .assay-block { border: 1px solid #e2e8f0; page-break-inside: avoid; box-shadow: none; margin-bottom: 10px; }
    .peak-editor-block { page-break-inside: avoid; border: none; }
    table { page-break-inside: avoid; box-shadow: none; border: 1px solid #e2e8f0; }
    th { background: #f8fafc !important; color: #0f172a; -webkit-print-color-adjust: exact; print-color-adjust: exact; border-bottom: 1px solid #cbd5e1; }
    td { border-bottom: 1px solid #e2e8f0; }
    img { max-width: 100%; page-break-inside: avoid; }
    .modebar { display: none !important; }
    .js-plotly-plot .plotly .modebar { display: none !important; }
}

/* ── Save Peaks Button ── */
.save-peaks-btn {
    background: #0ea5e9;
    color: white;
    box-shadow: 0 6px 16px rgba(14, 165, 233, 0.4);
}
.save-peaks-btn:hover {
    box-shadow: 0 10px 20px rgba(14, 165, 233, 0.5);
}

/* ── Interactive Peak Tables ── */
.peak-table-container {
    margin-top: 10px;
    padding: 0 10px 10px 10px;
}
.peak-table-container table {
    width: auto;
    min-width: 300px;
    margin: 0;
    box-shadow: none;
    border: 1px solid #e2e8f0;
}
.peak-table-container th {
    padding: 8px 12px;
    background: #f1f5f9;
}
.peak-table-container td {
    padding: 6px 12px;
}

/* ── Comment Boxes ── */
.comment-box-container {
    margin-top: 15px;
    border-radius: 8px;
    border: 1px dashed #cbd5e1;
    background: #f8fafc;
    overflow: hidden;
    transition: border-color 0.2s ease;
}
.comment-box-container:has(.comment-body.open) {
    border-color: #0ea5e9;
}
.comment-toggle-btn {
    width: 100%;
    background: none;
    border: none;
    padding: 10px 14px;
    text-align: left;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 8px;
    font-family: inherit;
    font-size: 0.85rem;
    font-weight: 600;
    color: #64748b;
    transition: color 0.15s ease, background 0.15s ease;
    user-select: none;
}
.comment-toggle-btn:hover { color: #0ea5e9; background: #f1f5f9; }
.comment-toggle-btn .caret {
    margin-left: auto;
    transition: transform 0.2s ease;
    font-style: normal;
    font-size: 0.75rem;
    opacity: 0.6;
}
.comment-body {
    display: none;
    padding: 0 14px 14px;
}
.comment-body.open { display: block; }
.report-comment {
    width: 100%;
    min-height: 80px;
    padding: 10px;
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    font-family: inherit;
    font-size: 0.9rem;
    resize: vertical;
    box-sizing: border-box;
    background: white;
}
.report-comment:focus {
    outline: none;
    border-color: #0ea5e9;
    box-shadow: 0 0 0 2px rgba(14, 165, 233, 0.2);
}

@media print {
    .comment-box-container { border: none; background: transparent; }
    .comment-toggle-btn { display: none; }
    .comment-body { display: block !important; padding: 0; }
    .report-comment { border: none; padding: 0; resize: none; overflow: hidden; background: transparent; }
}
</style>
"""


def _build_plotly_reflow_script() -> str:
    """Global Plotly reflow helpers for embedded/hidden report viewers."""
    return """
<script>
window.ReportPlotManager = (function() {
    var plots = {};
    var initialStates = {};

    function loadInitialStates() {
        if (Object.keys(initialStates).length) return;
        try {
            var tag = document.getElementById('plot-state');
            if (!tag) return;
            var raw = JSON.parse(tag.textContent || '{}');
            if (raw && typeof raw === 'object') initialStates = raw;
        } catch (e) {
            initialStates = {};
        }
    }

    function cloneRange(range) {
        if (!Array.isArray(range) || range.length !== 2) return null;
        var a = Number(range[0]);
        var b = Number(range[1]);
        if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
        return [a, b];
    }

    function resizeOne(gd) {
        if (!gd || !window.Plotly) return;
        if (typeof gd.isConnected === 'boolean' && !gd.isConnected) return;
        var width = gd.clientWidth || (gd.parentElement && gd.parentElement.clientWidth) || 0;
        if (width <= 0) return;
        try { Plotly.Plots.resize(gd); } catch (e) {}
        try { Plotly.relayout(gd, {autosize: true}); } catch (e) {}
    }

    function refreshAll() {
        for (var id in plots) {
            if (Object.prototype.hasOwnProperty.call(plots, id)) resizeOne(plots[id]);
        }
    }

    function captureState(gd) {
        if (!gd || !gd.layout) return null;
        var xRange = cloneRange(gd.layout.xaxis && gd.layout.xaxis.range);
        var yRange = cloneRange(gd.layout.yaxis && gd.layout.yaxis.range);
        if (!xRange && !yRange) return null;
        return {
            xaxis_range: xRange,
            yaxis_range: yRange
        };
    }

    function scheduleRefresh() {
        var delays = [0, 80, 250, 750];
        for (var i = 0; i < delays.length; i++) {
            setTimeout(refreshAll, delays[i]);
        }
        if (window.requestAnimationFrame) {
            window.requestAnimationFrame(function() { refreshAll(); });
        }
    }

    function attachObservers(gd) {
        if (!gd || gd.__fragglerObserversAttached) return;
        gd.__fragglerObserversAttached = true;

        if (typeof ResizeObserver === 'function') {
            try {
                var ro = new ResizeObserver(function() { resizeOne(gd); });
                ro.observe(gd);
                if (gd.parentElement) ro.observe(gd.parentElement);
                gd.__fragglerResizeObserver = ro;
            } catch (e) {}
        }

        if (typeof IntersectionObserver === 'function') {
            try {
                var io = new IntersectionObserver(function(entries) {
                    for (var i = 0; i < entries.length; i++) {
                        if (entries[i].isIntersecting) {
                            scheduleRefresh();
                            break;
                        }
                    }
                });
                io.observe(gd);
                gd.__fragglerIntersectionObserver = io;
            } catch (e) {}
        }
    }

    window.addEventListener('load', scheduleRefresh);
    window.addEventListener('resize', scheduleRefresh);
    window.addEventListener('pageshow', scheduleRefresh);
    window.addEventListener('focus', scheduleRefresh);
    document.addEventListener('visibilitychange', function() {
        if (document.visibilityState === 'visible') scheduleRefresh();
    });

    loadInitialStates();

    return {
        register: function(gd) {
            if (!gd || !gd.id) return;
            plots[gd.id] = gd;
            attachObservers(gd);
            scheduleRefresh();
        },
        getInitialStateForPlot: function(id) {
            loadInitialStates();
            return initialStates[id] || null;
        },
        getAllStates: function() {
            var all = {};
            for (var id in plots) {
                if (!Object.prototype.hasOwnProperty.call(plots, id)) continue;
                var state = captureState(plots[id]);
                if (state) all[id] = state;
            }
            return all;
        },
        refreshAll: scheduleRefresh,
        resizeOne: resizeOne
    };
})();
</script>
"""

def extract_dit_from_name(name: str) -> str | None:
    """Finner første forekomst av 2-sifret år + 'OUM' + 5 siffer."""
    m = DIT_PATTERN.search(name)
    return m.group(1) if m else None

def dit_to_year(dit: str) -> int | None:
    """25OUM10166 -> 2025, 26OUMxxxxx -> 2026, etc."""
    if not dit or len(dit) < 2: return None
    try:
        return 2000 + int(dit[:2])
    except ValueError:
        return None

def _resolve_report_display_name(entries: list[dict] | None = None) -> str:
    analysis_name = get_active_analysis_name()
    if entries:
        assays = {e.get("assay") for e in entries}
        if assays and assays.issubset({"FLT3-ITD", "FLT3-D835", "NPM1"}):
            return "Flt3"
    return "Klonalitet" if analysis_name == "clonality" else analysis_name.capitalize()


def _create_html_header(
    dit: str,
    year: int | None,
    num_entries: int,
    dit_root: Path,
    html_lines: list[str],
    *,
    display_name: str,
):
    """Appends the HTML head and page header to html_lines."""
    html_lines.extend(["<!DOCTYPE html>", "<html lang='no'>", "<head>", "<meta charset='utf-8'>"])
    html_lines.append(f"<title>{escape(dit)}_{display_name}_Resultater</title>")
    html_lines.append(REPORT_STYLE)
    html_lines.append('<script id="peak-data" type="application/json">{}</script>')
    html_lines.append('<script id="plot-state" type="application/json">{}</script>')
    html_lines.append("""
<script>
// Toggle comment boxes
function toggleComment(btn) {
    var body = btn.nextElementSibling;
    var caret = btn.querySelector('.caret');
    var isOpen = body.classList.toggle('open');
    caret.textContent = isOpen ? '▲' : '▼';
    if (isOpen) btn.querySelector('.comment-label').textContent = 'Skjul kommentar';
    else btn.querySelector('.comment-label').textContent = 'Legg til kommentar';
}

window.PeakManager = {
    plots: {},
    registerPlot: function(id, plotObj) { this.plots[id] = plotObj; },
    getAllPeaks: function() { var all = {}; for (var id in this.plots) { all[id] = this.plots[id].getPeaks(); } return all; },
    getInitialPeaksForPlot: function(id) { try { var data = JSON.parse(document.getElementById('peak-data').textContent); return data[id] || []; } catch(e) { return []; } },
    downloadUpdatedHtml: function() {
        // Force textareas back to innerHTML so they persist
        var tas = document.querySelectorAll('textarea.report-comment');
        for (var i = 0; i < tas.length; i++) {
            var val = tas[i].value.trim();
            tas[i].innerHTML = val;
            
            var container = tas[i].closest('.comment-box-container');
            var body = tas[i].closest('.comment-body');
            
            if (val !== "") {
                // If there's content, make sure it's open/visible
                body.classList.add('open');
            } else {
                // If empty, hide it (even from print)
                body.classList.remove('open');
            }
        }
        
        var allPeaks = this.getAllPeaks();
        var allPlotStates = (window.ReportPlotManager && window.ReportPlotManager.getAllStates)
            ? window.ReportPlotManager.getAllStates()
            : {};
        var currentHtml = document.documentElement.outerHTML;
        var peakDataStr = JSON.stringify(allPeaks);
        var plotStateStr = JSON.stringify(allPlotStates);
        var pattern = /<script id="peak-data" type="application\/json">[\\s\\S]*?<\/script>/;
        var newTag = '<script id="peak-data" type="application/json">\\n' + peakDataStr + '\\n<\/script>';
        var plotPattern = /<script id="plot-state" type="application\/json">[\\s\\S]*?<\/script>/;
        var newPlotTag = '<script id="plot-state" type="application/json">\\n' + plotStateStr + '\\n<\/script>';
        var updatedHtml = currentHtml.replace(pattern, newTag).replace(plotPattern, newPlotTag);
        var blob = new Blob(['<!DOCTYPE html>\\n' + updatedHtml], {type: 'text/html'});
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a'); a.href = url; a.download = document.title + '.html'; a.click(); URL.revokeObjectURL(url);
    }
};
function printReport() { window.print(); }
</script>
""")
    html_lines.append(_local_plotly_tag(dit_root, version="2.35.2"))
    html_lines.append(_build_plotly_reflow_script())
    html_lines.extend(["</head>", "<body>"])

    gen_date = datetime.now().strftime("%Y-%m-%d %H:%M")
    meta = [f"År: {year}"] if year else []
    meta.extend([f"{num_entries} analyserte filer", f"Generert: {gen_date}"])
    meta_str = " &nbsp;&bull;&nbsp; ".join(meta)
    
    html_lines.append(f"""
<div class='report-header no-print'>
  <h1>{escape(dit)}_{display_name}_Resultater</h1>
  <div class='meta'>{meta_str}</div>
</div>
<div style='display:none' class='print-only-header'>
  <h1>{escape(dit)}_{display_name}_Resultater</h1>
  <p>{" | ".join(meta)}</p>
</div>
""")

def _render_file_summary_table(dit_entries: list[dict], html_lines: list[str]):
    """Renders the overview table of all analyzed files."""
    html_lines.append("<h2>Oversikt over analyserte filer</h2>")
    html_lines.append("<p class='small'>Tabellen viser alle filer for dette DIT-nummeret.</p>")
    is_flt3 = {e.get("assay") for e in dit_entries}.issubset({"FLT3-ITD", "FLT3-D835", "NPM1"})
    if is_flt3:
        html_lines.append(
            "<table><tr><th>Filnavn</th><th>Assay</th><th>Behandling</th><th>WT</th><th>Mutert</th><th>Ratio</th><th>Ladder QC</th><th>R²</th></tr>"
        )
        for e in sorted(dit_entries, key=lambda x: (x["assay"], x.get("well_id") or "", x["fsa"].file_name)):
            status = e.get("ladder_qc_status", "unknown")
            r2 = e.get("ladder_r2", None)
            r2_str = f"{r2:.4f}" if r2 is not None and not np.isnan(r2) else "&mdash;"
            peaks = e["peaks_by_channel"].get(e["primary_peak_channel"], pd.DataFrame())
            wt_rows = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False) if not peaks.empty else pd.DataFrame()
            mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].sort_values("area", ascending=False) if not peaks.empty else pd.DataFrame()
            wt_text = _peak_text(_dominant_peak(wt_rows))
            mut_text = _peak_text(_dominant_peak(mut_rows))
            ratio = float(e.get("ratio", 0.0))
            ratio_str = f"{ratio:.4f}" if ratio > 0 else "&mdash;"
            html_lines.append(
                f"<tr><td>{escape(e['fsa'].file_name)}</td><td>{escape(e['assay'])}</td>"
                f"<td>{escape(_format_flt3_treatment(e))}</td>"
                f"<td>{wt_text}</td><td>{mut_text}</td><td>{ratio_str}</td>"
                f"<td>{escape(status)}</td><td>{r2_str}</td></tr>"
            )
    else:
        html_lines.append("<table><tr><th>Filnavn</th><th>Assay</th><th>Ladder</th><th>bp-område</th><th>Ladder QC</th><th>R²</th></tr>")
        for e in sorted(dit_entries, key=lambda x: (x["assay"], x["fsa"].file_name)):
            status = e.get("ladder_qc_status", "unknown")
            r2 = e.get("ladder_r2", None)
            r2_str = f"{r2:.4f}" if r2 is not None and not np.isnan(r2) else "&mdash;"
            html_lines.append(
                f"<tr><td>{escape(e['fsa'].file_name)}</td><td>{escape(e['assay'])}</td>"
                f"<td>{escape(e['ladder'])}</td><td>{int(e['bp_min'])}–{int(e['bp_max'])} bp</td>"
                f"<td>{escape(status)}</td><td>{r2_str}</td></tr>"
            )
    html_lines.append("</table>")

def _format_flt3_treatment(entry: dict) -> str:
    atype = entry["analysis_type"]
    treatment = "Standard"
    if atype == "TKD_digested":
        treatment = "Digert (EcoRV)"
    elif atype == "10x_diluted":
        treatment = "Fortynnet 1:10"
    elif atype == "25x_diluted":
        treatment = "Fortynnet 1:25"
    elif atype == "ratio_quant":
        treatment = "Ratio-sett"
    elif atype == "undiluted":
        treatment = "Ufortynnet"
    protocol_inj = entry.get("protocol_injection_time", entry.get("injection_time", 0))
    return f"{treatment} - {protocol_inj}s protokoll"


def _format_flt3_selection(entry: dict) -> str:
    selected = entry.get("selected_injection") or f"{entry.get('injection_time', 0)}s"
    source = entry.get("source_run_dir") or "ukjent kjøring"
    sizing = entry.get("sizing_method") or "ukjent sizing"
    reason = entry.get("selection_reason") or ""
    return f"Valgt {selected} fra {source} ({sizing}). {reason}".strip()


def _flt3_display_priority(entry: dict) -> int:
    assay = entry.get("assay")
    analysis_type = entry.get("analysis_type")
    if assay == "FLT3-ITD" and analysis_type == "ratio_quant":
        return 0
    if assay == "FLT3-D835":
        return 1
    if assay == "FLT3-ITD":
        return 2
    return 3


def _flt3_display_sort_key(entry: dict) -> tuple[int, str, str]:
    return (
        _flt3_display_priority(entry),
        entry.get("well_id") or "",
        entry["fsa"].file_name,
    )


def _flt3_report_blocks(assays: dict[str, list[dict]]) -> list[tuple[str, str, list[dict]]]:
    blocks: list[tuple[str, str, list[dict]]] = []
    itd_entries = assays.get("FLT3-ITD", [])
    itd_ratio_entries = [e for e in itd_entries if e.get("analysis_type") == "ratio_quant"]
    itd_other_entries = [e for e in itd_entries if e.get("analysis_type") != "ratio_quant"]

    if itd_ratio_entries:
        blocks.append(("FLT3-ITD", "FLT3-ITD-ratio", itd_ratio_entries))
    if "FLT3-D835" in assays:
        blocks.append(("FLT3-D835", "FLT3-D835", assays["FLT3-D835"]))
    if itd_other_entries:
        blocks.append(("FLT3-ITD", "FLT3-ITD", itd_other_entries))
    if "NPM1" in assays:
        blocks.append(("NPM1", "NPM1", assays["NPM1"]))
    return blocks

def _format_peak_list(mut_rows: pd.DataFrame, max_peaks: int = 3) -> str:
    if mut_rows.empty:
        return "Ingen mutasjoner detektert"
    mut_rows = mut_rows.sort_values("basepairs")
    parts = []
    for idx, (_, row) in enumerate(mut_rows.iterrows(), start=1):
        if idx > max_peaks:
            remaining = len(mut_rows) - max_peaks
            parts.append(f"+ {remaining} andre topper")
            break
        parts.append(f"{row.basepairs:.1f} bp ({row.area:,.0f})")
    return "<br>".join(parts)

def _peak_text(row: pd.Series | None, area_key: str = "area") -> str:
    if row is None:
        return "&mdash;"
    return f"{float(row.basepairs):.1f} bp <span class='small'>({float(row.get(area_key, 0.0)):,.0f})</span>"

def _dominant_peak(rows: pd.DataFrame, area_key: str = "area") -> pd.Series | None:
    if rows.empty:
        return None
    return rows.sort_values(area_key, ascending=False).iloc[0]

def _find_peak_in_range(peaks: pd.DataFrame, bp_min: float, bp_max: float) -> pd.DataFrame:
    if peaks.empty:
        return pd.DataFrame()
    return peaks[(peaks.basepairs >= bp_min) & (peaks.basepairs <= bp_max)].copy()


def _reportable_itd_mut_rows_for_report(entry: dict, peaks: pd.DataFrame, wt_rows: pd.DataFrame, mut_rows: pd.DataFrame) -> pd.DataFrame:
    if entry.get("assay") != "FLT3-ITD" or peaks.empty or mut_rows.empty or wt_rows.empty:
        return mut_rows
    if entry.get("analysis_type") == "ratio_quant":
        return mut_rows

    wt_main = wt_rows.iloc[0]
    wt_bp = float(wt_main.basepairs)
    wt_area = float(wt_main.area)
    shoulder_bp_limit = wt_bp + 12.0
    shoulder_area_limit = max(4000.0, wt_area * 0.02)

    keep_mask = ~(
        (mut_rows.basepairs <= shoulder_bp_limit)
        & (mut_rows.area <= shoulder_area_limit)
    )
    return mut_rows[keep_mask].copy()

def _itd_concordance_text(wt_row: pd.Series | None, mut_rows: pd.DataFrame) -> str:
    if wt_row is None and mut_rows.empty:
        return ""
    wt_blue = float(wt_row.get("area_DATA1", 0.0)) if wt_row is not None else 0.0
    wt_green = float(wt_row.get("area_DATA2", 0.0)) if wt_row is not None else 0.0
    mut_blue = float(mut_rows.get("area_DATA1", pd.Series(0.0)).sum()) if not mut_rows.empty else 0.0
    mut_green = float(mut_rows.get("area_DATA2", pd.Series(0.0)).sum()) if not mut_rows.empty else 0.0

    seen_blue = mut_blue > max(1000.0, wt_blue * 0.02)
    seen_green = mut_green > max(1000.0, wt_green * 0.02)
    if seen_blue and seen_green:
        return "Mutant signal i begge kanaler"
    if seen_blue:
        return "Mutant signal mest tydelig i bla kanal"
    if seen_green:
        return "Mutant signal mest tydelig i gronn kanal"
    return ""

D835_DIGEST_HEIGHT_MIN = 100.0
D835_DIGEST_AREA_MIN = 500.0


def _d835_digest_status(peaks: pd.DataFrame, wt_row: pd.Series | None, mut_row: pd.Series | None) -> tuple[str, pd.Series | None]:
    digest_rows = _find_peak_in_range(peaks, 145.0, 155.5)
    digest_row = _dominant_peak(digest_rows)
    digest_area = float(digest_row.area) if digest_row is not None else 0.0
    digest_height = float(digest_row.peaks) if digest_row is not None else 0.0
    wt_area = float(wt_row.area) if wt_row is not None else 0.0
    mut_area = float(mut_row.area) if mut_row is not None else 0.0

    if digest_row is None or digest_height < D835_DIGEST_HEIGHT_MIN or digest_area < D835_DIGEST_AREA_MIN:
        return "", None
    if digest_area >= max(wt_area, mut_area) * 0.60:
        return "Mulig ufullstendig kutting", digest_row
    return "", digest_row

def _build_flt3_summary_table(e: dict) -> str:
    """Validation-oriented FLT3/NPM1 table below each figure."""
    assay = e["assay"]
    ratio = float(e.get("ratio", 0.0))
    ratio_str = f"{ratio:.4f}" if ratio > 0 else "&mdash;"
    positive_ratio = float(ASSAY_CONFIG.get(assay, {}).get("positive_ratio", 0.01))
    peaks = e["peaks_by_channel"].get(e["primary_peak_channel"], pd.DataFrame())

    wt_row = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False) if not peaks.empty else pd.DataFrame()
    mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].sort_values(["peaks", "basepairs"], ascending=[False, True]) if not peaks.empty else pd.DataFrame()
    wt_main = _dominant_peak(wt_row)
    mut_main = _dominant_peak(mut_rows)

    if assay == "FLT3-ITD":
        wt_blue = wt_green = mut_blue = mut_green = 0.0
        reportable_mut_rows = peaks[0:0].copy() if peaks.empty else peaks.iloc[0:0].copy()
        if wt_main is not None:
            r = wt_main
            wt_blue = r.get("area_DATA1", 0.0)
            wt_green = r.get("area_DATA2", 0.0)
        if not mut_rows.empty:
            mut_blue = mut_rows.get("area_DATA1", pd.Series(0.0)).sum()
            mut_green = mut_rows.get("area_DATA2", pd.Series(0.0)).sum()
            reportable_mut_rows = _reportable_itd_mut_rows_for_report(e, peaks, wt_row, mut_rows)
        ratio_num = float(e.get("ratio_numerator_area", 0.0))
        ratio_den = float(e.get("ratio_denominator_area", 0.0))
        mut_prop = (ratio_num / (ratio_num + ratio_den)) if (ratio_num + ratio_den) > 0 else 0.0
        label = "Positiv" if ratio >= positive_ratio else "Negativ"
        if ratio < positive_ratio and not reportable_mut_rows.empty:
            label = "Negativ, dokumentert"
        concordance = _itd_concordance_text(wt_main, reportable_mut_rows)
        validation_text = f"<strong>{label}</strong>"
        if concordance:
            validation_text += f"<br><span class='small'>{concordance}</span>"
        return (
            "<div style='margin-top:10px; margin-bottom:24px;'>"
            "<table style='width:100%; border:1px solid #e2e8f0; table-layout:fixed;'>"
            "<tr><th>WT-topp</th><th>Muterte topper</th><th>Bla kanal</th><th>Gronn kanal</th><th>Ratioer</th><th>Validering</th></tr>"
            f"<tr><td>{_peak_text(wt_main)}</td>"
            f"<td>{_format_peak_list(mut_rows, max_peaks=6)}</td>"
            f"<td>WT: {wt_blue:,.0f}<br>Mut: {mut_blue:,.0f}</td>"
            f"<td>WT: {wt_green:,.0f}<br>Mut: {mut_green:,.0f}</td>"
            f"<td>ITD-ratio: <strong>{ratio_str}</strong><br>"
            f"<span class='small'>Mut/WT: {float(e.get('ratio_numerator_area', 0.0)):,.0f} / {float(e.get('ratio_denominator_area', 0.0)):,.0f}<br>"
            f"Mut/(Mut+WT): {mut_prop:.4f}<br>Positiv grense > {positive_ratio:.2f}</span></td>"
            f"<td>{validation_text}</td></tr></table></div>"
        )

    if assay == "FLT3-D835":
        digest_status, digest_row = _d835_digest_status(peaks, wt_main, mut_main)
        label = "Positiv" if ratio >= positive_ratio else "Negativ" if mut_main is None else "Under positiv grense"
        digest_text = "&mdash;"
        if digest_row is not None:
            digest_text = _peak_text(digest_row)
            if digest_status:
                digest_text += f"<br><span class='small'>{digest_status}</span>"
        return (
            "<div style='margin-top:10px; margin-bottom:24px;'>"
            "<table style='width:100%; border:1px solid #e2e8f0; table-layout:fixed;'>"
            "<tr><th>WT-topp</th><th>Mutert topp</th><th>150 bp kontroll</th><th>TKD-ratio</th><th>Validering</th></tr>"
            f"<tr><td>{_peak_text(wt_main)}</td>"
            f"<td>{_peak_text(mut_main)}<br><span class='small'>{_format_peak_list(mut_rows, max_peaks=4)}</span></td>"
            f"<td>{digest_text}</td>"
            f"<td><strong>{ratio_str}</strong><br><span class='small'>Mut/WT: {float(e.get('ratio_numerator_area', 0.0)):,.0f} / {float(e.get('ratio_denominator_area', 0.0)):,.0f}<br>Positiv grense > {positive_ratio:.2f}</span></td>"
            f"<td><strong>{label}</strong></td></tr></table></div>"
        )

    if assay == "NPM1":
        label = "Positiv" if ratio >= positive_ratio else "Negativ" if mut_main is None else "Manuell vurdering"
        return (
            "<div style='margin-top:10px; margin-bottom:24px;'>"
            "<table style='width:100%; border:1px solid #e2e8f0; table-layout:fixed;'>"
            "<tr><th>Villtype</th><th>Mutert</th><th>Ratio</th><th>Validering</th></tr>"
            f"<tr><td>{_peak_text(wt_main)}</td><td>{_format_peak_list(mut_rows, max_peaks=4)}</td>"
            f"<td><strong>{ratio_str}</strong></td><td><strong>{label}</strong></td></tr></table></div>"
        )

    return ""

def _render_assay_block(assay_name: str, assay_entries: list[dict], html_lines: list[str]):
    """Renders a single assay block with plots for each file."""
    display_name = assay_name
    reference_assay = assay_name
    if assay_name == "FLT3-ITD-ratio":
        display_name = "FLT3-ITD-ratio"
        reference_assay = "FLT3-ITD"

    html_lines.append("<div class='assay-block'>")
    html_lines.append(f"<h3>{escape(display_name)}</h3>")
    ref_ranges = ASSAY_REFERENCE_RANGES.get(reference_assay)
    if ref_ranges:
        ranges_str = ", ".join(f"{int(a)}–{int(b)} bp" for (a, b) in ref_ranges)
        label_txt = ASSAY_REFERENCE_LABEL.get(reference_assay, ranges_str)
        html_lines.append(f"<p class='small'><strong>Referanseområde:</strong> {escape(ranges_str)}<br>{escape(label_txt)}</p>")

    sort_key = _flt3_display_sort_key if reference_assay in {"FLT3-ITD", "FLT3-D835", "NPM1"} else (lambda x: x["fsa"].file_name)
    for e in sorted(assay_entries, key=sort_key):
        fsa, primary_ch = e["fsa"], e["primary_peak_channel"]
        html_lines.append(f"<p class='sample-header'>{escape(fsa.file_name)} ({escape(primary_ch)})</p>")
        if reference_assay in {"FLT3-ITD", "FLT3-D835", "NPM1"}:
            sub = [
                f"Well: {e.get('well_id') or '&mdash;'}",
                f"Injeksjon: {e.get('selected_injection') or ''}",
            ]
            html_lines.append(f"<p class='small'>{escape(' | '.join(sub))}</p>")
        try:
            frag = build_interactive_peak_plot_for_entry(e)
            html_lines.append(frag if frag else "<p class='small'><em>Ingen data å vise.</em></p>")
        except Exception as ex:
            html_lines.append(f"<p class='small'><em>Kunne ikke lage plott: {escape(str(ex))}</em></p>")
            
        if reference_assay in {"FLT3-ITD", "FLT3-D835", "NPM1"}:
            html_lines.append(_build_flt3_summary_table(e))

    # Add collapsible Comment Box for the overall assay
    html_lines.append(
        "<div class='comment-box-container'>"
        "<button class='comment-toggle-btn' onclick='toggleComment(this)'>"
        "💬 <span class='comment-label'>Legg til kommentar</span>"
        f" <em style='font-weight:400;opacity:0.7;'>({escape(display_name)})</em>"
        "<i class='caret'>&#x25BC;</i>"
        "</button>"
        "<div class='comment-body'>"
        "<textarea class='report-comment' placeholder='Skriv inn eventuelle kommentarer her...'></textarea>"
        "</div>"
        "</div>"
    )

    html_lines.append("</div>")

def _render_tcrb_rep_block(entries: list[dict], replicate_num: str, html_lines: list[str]):
    """Renders a combination block for TCRb replicates."""
    if not entries: return
    html_lines.append("<div class='assay-block'>")
    html_lines.append(f"<h3>TCRβ – Parallell {replicate_num}</h3>")
    html_lines.append("<p class='small'>TCRβ mix A, B og C med felles y-akse for enkel sammenligning.</p>")
    html_lines.append("<div class='combo-grid'>")
    group_y = compute_group_ymax_for_entries(entries)
    
    # Calculate global X range
    forced_xmin = min((float(e["bp_min"]) for e in entries), default=0)
    forced_xmax = max((float(e["bp_max"]) for e in entries), default=1000)

    for e in sorted(entries, key=lambda x: x["assay"]):
        fsa, primary_ch = e["fsa"], e["primary_peak_channel"]
        e_combo = dict(e)
        e_combo["forced_ymax"] = group_y
        e_combo["forced_xmin"] = forced_xmin
        e_combo["forced_xmax"] = forced_xmax
        
        html_lines.append("<div class='combo-item'>")
        html_lines.append(f"<p class='sample-header'>{escape(e_combo['assay'])} – {escape(fsa.file_name)}</p>")
        try:
            frag = build_interactive_peak_plot_for_entry(e_combo)
            html_lines.append(frag if frag else "<p class='small'><em>Ingen data å vise.</em></p>")
        except Exception as ex:
            html_lines.append(f"<p class='small'><em>Kunne ikke lage plott: {escape(str(ex))}</em></p>")
        html_lines.append("</div>")
    html_lines.append("</div></div>")

def _render_tcrg_combo_block(tcrg_entries: list[dict], html_lines: list[str]):
    """Renders a combined block for TCRg assays."""
    if not tcrg_entries: return
    group_y = compute_group_ymax_for_entries(tcrg_entries)
    
    # Calculate global X range
    forced_xmin = min((float(e["bp_min"]) for e in tcrg_entries), default=0)
    forced_xmax = max((float(e["bp_max"]) for e in tcrg_entries), default=1000)

    html_lines.append("<h2>Kombinasjonsfigur – TCRγ</h2><div class='assay-block'>")
    html_lines.append("<p class='small'>TCRγ mix A og mix B (begge paralleller) med felles x- og y-akse.</p>")
    html_lines.append("<div class='combo-grid'>")
    for e in tcrg_entries:
        fsa, primary_ch = e["fsa"], e["primary_peak_channel"]
        e_combo = dict(e)
        e_combo["forced_ymax"] = group_y
        e_combo["forced_xmin"] = forced_xmin
        e_combo["forced_xmax"] = forced_xmax
        
        html_lines.append("<div class='combo-item'>")
        html_lines.append(f"<p class='sample-header'>{escape(e_combo['assay'])} – {escape(fsa.file_name)}</p>")
        try:
            frag = build_interactive_peak_plot_for_entry(e_combo)
            html_lines.append(frag if frag else "<p class='small'><em>Ingen data å vise.</em></p>")
        except Exception as ex:
            html_lines.append(f"<p class='small'><em>Kunne ikke lage plott: {escape(str(ex))}</em></p>")
        html_lines.append("</div>")
    html_lines.append("</div></div>")

def _render_sl_section(all_sl_entries: list[dict], html_lines: list[str]):
    """Renders the Size Ladder (DNA quality) section."""
    valid_entries = [e for e in all_sl_entries if e.get("sl_metrics")]
    if not valid_entries: return
    html_lines.append("<h2>Size Ladder (SL) – DNA-kvalitet</h2>")
    for e in sorted(valid_entries, key=lambda x: x.get("fsa").file_name if x.get("fsa") else ""):
        sl_metrics = e.get("sl_metrics")
        html_lines.append(f"<h3>SL-fil: {escape(e['fsa'].file_name)}</h3>")
        if not sl_metrics:
            html_lines.append("<p><em>Ingen SL-area-metrikker tilgjengelig.</em></p>")
            continue
        targets, areas, pcts = sl_metrics.get("targets_bp", []), sl_metrics.get("areas", []), sl_metrics.get("percents", [])
        total_area = sl_metrics.get("total_area", float("nan"))
        html_lines.append("<table><tr><th>Fragment (bp)</th><th>Area</th><th>% av total</th></tr>")
        for bp_val, area_val, pct_val in zip(targets, areas, pcts):
            area_str = f"{area_val:,.0f}".replace(",", " ") if not np.isnan(area_val) else "&mdash;"
            pct_str = f"{pct_val:.1f} %" if pct_val is not None and not np.isnan(pct_val) else "&mdash;"
            html_lines.append(f"<tr><td>{bp_val:.0f}</td><td>{area_str}</td><td>{pct_str}</td></tr>")
        tot_str = f"{total_area:,.0f}".replace(",", " ") if not np.isnan(total_area) else "&mdash;"
        html_lines.append(f"<tr><td><strong>Total</strong></td><td><strong>{tot_str}</strong></td><td></td></tr></table>")


def build_dit_html_reports(entries: list[dict], assay_outdir: Path):
    """Main entry for building per-patient DIT reports."""
    qc_sl_entries = []
    per_dit: dict[str, list[dict]] = defaultdict(list)
    from core.qc.qc_markers import control_id_from_filename
    for e in entries:
        ctrl_id = control_id_from_filename(e["fsa"].file_name)
        if ctrl_id in ("PK", "PK1", "PK2") and e.get("assay") == "SL":
            qc_sl_entries.append(e)
        if (dit := e.get("dit")): per_dit[dit].append(e)

    if not per_dit:
        print_warning("[DIT] Fant ingen DIT-nummer – ingen rapporter generert.")
        return

    assay_outdir.mkdir(exist_ok=True, parents=True)
    print_green(f"[DIT] Lager pasientrapporter i {assay_outdir}")
    display_name = _resolve_report_display_name(entries)

    for dit, dit_entries in sorted(per_dit.items()):
        year = dit_to_year(dit)
        assays: dict[str, list[dict]] = defaultdict(list)
        for e in dit_entries: assays[e["assay"]].append(e)

        html_lines: list[str] = []
        _create_html_header(dit, year, len(dit_entries), assay_outdir, html_lines, display_name=display_name)
        _render_file_summary_table(dit_entries, html_lines)

        html_lines.append("<h2>Assay-spesifikke oversikter</h2>")
        if "FLT3-ITD" in assays or "FLT3-D835" in assays or "NPM1" in assays:
            flt3_blocks = _flt3_report_blocks(assays)
            handled = {"FLT3-ITD", "FLT3-D835", "NPM1"}
            ordered = [a for a in ASSAY_DISPLAY_ORDER if a in assays and a not in handled] + [a for a in assays if a not in ASSAY_DISPLAY_ORDER and a not in handled]
            for assay_key, block_title, block_entries in flt3_blocks:
                _render_assay_block(block_title, block_entries, html_lines)

            for name in ordered:
                _render_assay_block(name, assays[name], html_lines)
                
                # Special Combination Sections
                if name == "TCRbC":
                    present = [a for a in ["TCRbA", "TCRbB", "TCRbC"] if a in assays]
                    sorted_rep = {a: sorted(assays[a], key=lambda x: x["fsa"].file_name) for a in present}
                    rep1 = [lst[0] for a, lst in sorted_rep.items() if len(lst) >= 1]
                    rep2 = [lst[1] for a, lst in sorted_rep.items() if len(lst) >= 2]
                    if rep1: html_lines.append("<h2>Kombinasjonsfigurer – TCRβ</h2>")
                    _render_tcrb_rep_block(rep1, "1", html_lines)
                    _render_tcrb_rep_block(rep2, "2", html_lines)
                
                if name == "TCRgB":
                    tcrg_all = []
                    for a in ["TCRgA", "TCRgB"]:
                        if a in assays: tcrg_all.extend(sorted(assays[a], key=lambda x: x["fsa"].file_name))
                    _render_tcrg_combo_block(tcrg_all, html_lines)
            _render_sl_section(dit_entries + qc_sl_entries, html_lines)
            
            html_lines.append("""
<div class="print-fab no-print">
  <button class="print-btn save-peaks-btn" onclick="PeakManager.downloadUpdatedHtml()">💾&nbsp; Save Peaks</button>
  <button class="print-btn" onclick="printReport()">🖨&nbsp; Print / PDF</button>
</div>
</body></html>""")
            
            out_html = assay_outdir / f"{dit}_{display_name}_Resultater.html"
            out_html.write_text("\n".join(html_lines), encoding="utf-8")
            print_green(f"[DIT] Lagret: {out_html}")
            continue

        ordered = [a for a in ASSAY_DISPLAY_ORDER if a in assays] + [a for a in assays if a not in ASSAY_DISPLAY_ORDER]
        
        for name in ordered:
            _render_assay_block(name, assays[name], html_lines)
            
            # Special Combination Sections
            if name == "TCRbC":
                present = [a for a in ["TCRbA", "TCRbB", "TCRbC"] if a in assays]
                sorted_rep = {a: sorted(assays[a], key=lambda x: x["fsa"].file_name) for a in present}
                rep1 = [lst[0] for a, lst in sorted_rep.items() if len(lst) >= 1]
                rep2 = [lst[1] for a, lst in sorted_rep.items() if len(lst) >= 2]
                if rep1: html_lines.append("<h2>Kombinasjonsfigurer – TCRβ</h2>")
                _render_tcrb_rep_block(rep1, "1", html_lines)
                _render_tcrb_rep_block(rep2, "2", html_lines)
            
            if name == "TCRgB":
                tcrg_all = []
                for a in ["TCRgA", "TCRgB"]:
                    if a in assays: tcrg_all.extend(sorted(assays[a], key=lambda x: x["fsa"].file_name))
                _render_tcrg_combo_block(tcrg_all, html_lines)

        _render_sl_section(dit_entries + qc_sl_entries, html_lines)
        
        html_lines.append("""
<div class="print-fab no-print">
  <button class="print-btn save-peaks-btn" onclick="PeakManager.downloadUpdatedHtml()">💾&nbsp; Save Peaks</button>
  <button class="print-btn" onclick="printReport()">🖨&nbsp; Print / PDF</button>
</div>
</body></html>""")
        
        out_html = assay_outdir / f"{dit}_{display_name}_Resultater.html"
        out_html.write_text("\n".join(html_lines), encoding="utf-8")
        print_green(f"[DIT] Lagret: {out_html}")


def interpret_sl_quality(percents, total_area):
    """Automatisk fortolkning av DNA-kvalitet basert på fragmentfordeling."""
    p100, p200, p300, p400, p600 = (percents[i] if i < len(percents) else float("nan") for i in range(5))
    if np.isnan(total_area) or total_area < 1e4: return "Materialet er uegnet (svært lite signal)."
    if np.isnan(p100) or p100 < 5: return "Materialet er uegnet (svært svak 100 bp-peak)."
    if p100 >= 85 and p200 <= 15 and p300 <= 5: return "Svært fragmentert materiale."
    sum_100_300, sum_100_200 = p100 + p200 + p300, p100 + p200
    if p100 >= 60 and sum_100_200 >= 80 and p300 <= 15: return "Mer enn 50 % fragmentert – redusert sensitivitet."
    if p100 >= 45 and sum_100_300 >= 70: return "Litt fragmentert – kan redusere sensitivitet."
    if p100 <= 50 and sum_100_200 <= 70 and p300 >= 10 and p400 >= 5: return "Bra kvalitet."
    return "Uvanlig fordeling – vurder manuelt."
