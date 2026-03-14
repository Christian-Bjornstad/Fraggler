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

def _create_html_header(dit: str, year: int | None, num_entries: int, dit_root: Path, html_lines: list[str]):
    """Appends the HTML head and page header to html_lines."""
    analysis_name = get_active_analysis_name()
    display_name = "Klonalitet" if analysis_name == "clonality" else analysis_name.capitalize()
    
    html_lines.extend(["<!DOCTYPE html>", "<html lang='no'>", "<head>", "<meta charset='utf-8'>"])
    html_lines.append(f"<title>{escape(dit)}_{display_name}_Resultater</title>")
    html_lines.append(REPORT_STYLE)
    html_lines.append('<script id="peak-data" type="application/json">{}</script>')
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
        var currentHtml = document.documentElement.outerHTML;
        var peakDataStr = JSON.stringify(allPeaks);
        var pattern = /<script id="peak-data" type="application\/json">[\\s\\S]*?<\/script>/;
        var newTag = '<script id="peak-data" type="application/json">\\n' + peakDataStr + '\\n<\/script>';
        var updatedHtml = currentHtml.replace(pattern, newTag);
        var blob = new Blob(['<!DOCTYPE html>\\n' + updatedHtml], {type: 'text/html'});
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a'); a.href = url; a.download = document.title + '.html'; a.click(); URL.revokeObjectURL(url);
    }
};
function printReport() { window.print(); }
</script>
""")
    html_lines.append(_local_plotly_tag(dit_root, version="2.35.2"))
    html_lines.extend(["</head>", "<body>"])

    gen_date = datetime.now().strftime("%Y-%m-%d %H:%M")
    meta = [f"År: {year}"] if year else []
    meta.extend([f"{num_entries} analyserte filer", f"Generert: {gen_date}"])
    meta_str = " &nbsp;&bull;&nbsp; ".join(meta)
    
    analysis_name = get_active_analysis_name()
    display_name = "Klonalitet" if analysis_name == "clonality" else analysis_name.capitalize()
    
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

def _render_flt3_results_table(dit_entries: list[dict], html_lines: list[str]):
    """Renders a specialized summary table for FLT3/NPM1 results."""
    html_lines.append("<h2>Analyseresultater – FLT3 / NPM1</h2>")
    html_lines.append(
        "<table><tr><th>Assay</th><th>Prøvetype</th><th>Parallell</th><th>Behandling</th>"
        "<th>Resultat</th></tr>"
    )
    
    for e in sorted(dit_entries, key=lambda x: (x["assay"], x["analysis_type"])):
        assay = e["assay"]
        g = e["group"]
        atype = e["analysis_type"]
        inj = e.get("injection_time", 0)
        ratio = e.get("ratio", 0.0)
        par = e.get("parallel") or "&mdash;"
        
        # Treatment description
        treatment = "Standard"
        if atype == "TKD_digested": treatment = "Digert (EcoRV)"
        elif atype == "10x_diluted": treatment = "Fortynnet 1:10"
        elif atype == "25x_diluted": treatment = "Fortynnet 1:25"
        elif atype == "ratio_quant": treatment = "Ratio-sett"
        elif atype == "undiluted": treatment = "Ufortynnet"

        # Extract peak areas
        wt_area, mut_area = 0.0, 0.0
        peaks = e["peaks_by_channel"].get(e["primary_peak_channel"], pd.DataFrame())
        if not peaks.empty:
            wt_row = peaks[peaks.label == "WT"]
            if not wt_row.empty: wt_area = wt_row.sort_values("peaks", ascending=False).iloc[0].area
            
            mut_row = peaks[peaks.label.isin(["MUT", "ITD"])]
            if not mut_row.empty: mut_area = mut_row.area.sum()

        ratio_str = f"{ratio:.4f}" if ratio > 0 else "&mdash;"
        
        # Simple interpretation
        label = "Villtype"
        if assay == "FLT3-ITD" and mut_area > 0: label = "MUTERT (ITD)"
        elif assay == "FLT3-D835" and ratio > 0.01: label = "MUTERT (D835)"
        elif assay == "NPM1" and mut_area > 0: label = "MUTERT (NPM1)"
        
        html_lines.append(
            f"<tr><td><strong>{escape(assay)}</strong></td>"
            f"<td>{escape(g)}</td><td>{escape(par)}</td>"
            f"<td>{escape(treatment)}<br><span class='small'>{inj}s inj.</span></td>"
            f"<td><strong>{label}</strong></td></tr>"
        )
    html_lines.append("</table>")

def _render_assay_block(assay_name: str, assay_entries: list[dict], html_lines: list[str]):
    """Renders a single assay block with plots for each file."""
    html_lines.append("<div class='assay-block'>")
    html_lines.append(f"<h3>{escape(assay_name)}</h3>")
    ref_ranges = ASSAY_REFERENCE_RANGES.get(assay_name)
    if ref_ranges:
        ranges_str = ", ".join(f"{int(a)}–{int(b)} bp" for (a, b) in ref_ranges)
        label_txt = ASSAY_REFERENCE_LABEL.get(assay_name, ranges_str)
        html_lines.append(f"<p class='small'><strong>Referanseområde:</strong> {escape(ranges_str)}<br>{escape(label_txt)}</p>")

    for e in sorted(assay_entries, key=lambda x: x["fsa"].file_name):
        fsa, primary_ch = e["fsa"], e["primary_peak_channel"]
        html_lines.append(f"<p class='sample-header'>{escape(fsa.file_name)} ({escape(primary_ch)})</p>")
        try:
            frag = build_interactive_peak_plot_for_entry(e)
            html_lines.append(frag if frag else "<p class='small'><em>Ingen data å vise.</em></p>")
        except Exception as ex:
            html_lines.append(f"<p class='small'><em>Kunne ikke lage plott: {escape(str(ex))}</em></p>")
            
        # Add detailed ratio tables directly below the plot for FLT3 targets
        if assay_name == "FLT3-ITD":
            html_lines.append(_build_itd_detailed_table(e))
            html_lines.append(_build_d835_detailed_table(e))

    # Add collapsible Comment Box for the overall assay
    html_lines.append(
        "<div class='comment-box-container'>"
        "<button class='comment-toggle-btn' onclick='toggleComment(this)'>"
        "💬 <span class='comment-label'>Legg til kommentar</span>"
        f" <em style='font-weight:400;opacity:0.7;'>({escape(assay_name)})</em>"
        "<i class='caret'>&#x25BC;</i>"
        "</button>"
        "<div class='comment-body'>"
        "<textarea class='report-comment' placeholder='Skriv inn eventuelle kommentarer her...'></textarea>"
        "</div>"
        "</div>"
    )

    html_lines.append("</div>")

def _build_itd_detailed_table(e: dict) -> str:
    """Builds the detailed Excel-like table for ITD (separated Green/Blue)."""
    peaks = e["peaks_by_channel"].get(e["primary_peak_channel"], pd.DataFrame())
    if peaks.empty: return ""
    
    wt_row = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False)
    mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])]
    
    # Extract channel 1 (Green) and channel 2 (Blue) if available
    # Default to 0 if missing
    wt_g, wt_b = 0.0, 0.0
    mut_g, mut_b = 0.0, 0.0
    
    if not wt_row.empty:
        r = wt_row.iloc[0]
        wt_g = r.get("area_DATA1", 0.0)
        wt_b = r.get("area_DATA2", 0.0)
        
    if not mut_rows.empty:
        mut_g = mut_rows.get("area_DATA1", pd.Series(0.0)).sum()
        mut_b = mut_rows.get("area_DATA2", pd.Series(0.0)).sum()
        
    ratio = e.get("ratio", 0.0)
    ratio_str = f"{ratio:.4f}" if ratio > 0 else "&mdash;"
    
    # Check if ratio suggests a positive result (e.g. > 0.01)
    # The summary already does this, but here it's clearly for the ratio
    label = "Negativ"
    if ratio > 0.01: label = "Positiv"
    
    t = []
    t.append("<div style='margin-top: 10px; margin-bottom: 24px;'>")
    t.append("<table style='width: 100%; border: 1px solid #e2e8f0;'>")
    t.append(
        "<tr><th colspan='2' style='text-align:center;'>Villtype (WT)</th>"
        "<th colspan='2' style='text-align:center;'>Mutert (ITD)</th>"
        "<th>Ratio</th><th>Tolkning</th></tr>"
    )
    t.append("<tr><th>Grønn (Area)</th><th>Blå (Area)</th><th>Grønn (Area)</th><th>Blå (Area)</th><th>ITD-ratio</th><th></th></tr>")
    
    t.append("<tr>")
    t.append(f"<td>{wt_g:,.0f}</td><td>{wt_b:,.0f}</td>")
    t.append(f"<td>{mut_g:,.0f}</td><td>{mut_b:,.0f}</td>")
    t.append(f"<td><strong>{ratio_str}</strong></td><td><strong>{label}</strong></td>")
    t.append("</tr>")
    t.append("</table></div>")
    return "".join(t)

def _build_d835_detailed_table(e: dict) -> str:
    """Builds the detailed Excel-like table for D835 (separate mutant peaks)."""
    peaks = e["peaks_by_channel"].get(e["primary_peak_channel"], pd.DataFrame())
    if peaks.empty: return ""
    
    wt_row = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False)
    mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].sort_values("peaks", ascending=False)
    
    wt_area = wt_row.iloc[0].area if not wt_row.empty else 0.0
    
    t = []
    t.append("<div style='margin-top: 10px; margin-bottom: 24px;'>")
    t.append("<table style='width: 100%; border: 1px solid #e2e8f0;'>")
    
    if mut_rows.empty:
        t.append("<tr><th>WT (Sort Area)</th><th>Mutert</th><th>Ratio</th><th>Tolkning</th></tr>")
        t.append(f"<tr><td>{wt_area:,.0f}</td><td>Ingen mutasjoner detektert</td><td>&mdash;</td><td><strong>Negativ</strong></td></tr>")
    else:
        # Build headers for up to N mutant peaks (usually 1 or 2)
        n_muts = len(mut_rows)
        headers = ["WT (Sort Area)"]
        for i in range(1, n_muts + 1):
            headers.append(f"Mutert Peak {i} (Area) [Størrelse]")
            headers.append(f"Ratio {i}")
        headers.append("Total Ratio")
        headers.append("Tolkning")
        
        t.append("<tr>")
        for h in headers: t.append(f"<th>{h}</th>")
        t.append("</tr>")
        
        t.append("<tr>")
        t.append(f"<td>{wt_area:,.0f}</td>")
        
        tot_mut = 0.0
        for _, r in mut_rows.iterrows():
            m_area = r.area
            m_bp = r.basepairs
            tot_mut += m_area
            r_val = m_area / wt_area if wt_area > 0 else 0.0
            t.append(f"<td>{m_area:,.0f} <span class='small'>[{m_bp:.1f} bp]</span></td>")
            t.append(f"<td>{r_val:.4f}</td>")
            
        tot_ratio = tot_mut / wt_area if wt_area > 0 else 0.0
        label = "Positiv" if tot_ratio > 0.01 else "Negativ"
        
        t.append(f"<td><strong>{tot_ratio:.4f}</strong></td>")
        t.append(f"<td><strong>{label}</strong></td>")
        t.append("</tr>")
        
    t.append("</table></div>")
    return "".join(t)

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

    for dit, dit_entries in sorted(per_dit.items()):
        year = dit_to_year(dit)
        assays: dict[str, list[dict]] = defaultdict(list)
        for e in dit_entries: assays[e["assay"]].append(e)

        html_lines: list[str] = []
        _create_html_header(dit, year, len(dit_entries), assay_outdir, html_lines)
        _render_file_summary_table(dit_entries, html_lines)

        if get_active_analysis_name() == "flt3":
            _render_flt3_results_table(dit_entries, html_lines)

        html_lines.append("<h2>Assay-spesifikke oversikter</h2>")
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
        
        analysis_name = get_active_analysis_name()
        display_name = "Klonalitet" if analysis_name == "clonality" else analysis_name.capitalize()
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
