"""Standalone HTML reporting for the general analysis."""
from __future__ import annotations

import re
from datetime import datetime
from html import escape
from pathlib import Path

import numpy as np

from core.html_reports import REPORT_STYLE, _build_plotly_reflow_script
from core.plotly_offline import local_plotly_tag as _local_plotly_tag
from core.plotting_plotly import build_interactive_peak_plot_for_entry


def _safe_report_name(label: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", label.strip())
    return cleaned.strip("_") or "general"


def _build_report_script() -> str:
    return """
<script>
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
    getAllPeaks: function() {
        var all = {};
        for (var id in this.plots) {
            if (Object.prototype.hasOwnProperty.call(this.plots, id)) {
                all[id] = this.plots[id].getPeaks();
            }
        }
        return all;
    },
    getInitialPeaksForPlot: function(id) {
        try {
            var data = JSON.parse(document.getElementById('peak-data').textContent || '{}');
            return data[id] || [];
        } catch (e) {
            return [];
        }
    },
    downloadUpdatedHtml: function() {
        var tas = document.querySelectorAll('textarea.report-comment');
        for (var i = 0; i < tas.length; i++) {
            var val = tas[i].value.trim();
            tas[i].innerHTML = val;
            var body = tas[i].closest('.comment-body');
            if (body) {
                if (val !== "") body.classList.add('open');
                else body.classList.remove('open');
            }
        }

        var allPeaks = this.getAllPeaks();
        var allPlotStates = (window.ReportPlotManager && window.ReportPlotManager.getAllStates)
            ? window.ReportPlotManager.getAllStates()
            : {};
        var currentHtml = document.documentElement.outerHTML;
        var peakDataStr = JSON.stringify(allPeaks);
        var plotStateStr = JSON.stringify(allPlotStates);
        var pattern = /<script id="peak-data" type="application\\/json">[\\s\\S]*?<\\/script>/;
        var newTag = '<script id="peak-data" type="application/json">\\n' + peakDataStr + '\\n<\\/script>';
        var plotPattern = /<script id="plot-state" type="application\\/json">[\\s\\S]*?<\\/script>/;
        var newPlotTag = '<script id="plot-state" type="application/json">\\n' + plotStateStr + '\\n<\\/script>';
        var updatedHtml = currentHtml.replace(pattern, newTag).replace(plotPattern, newPlotTag);
        var blob = new Blob(['<!DOCTYPE html>\\n' + updatedHtml], {type: 'text/html'});
        var url = URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = document.title + '.html';
        a.click();
        URL.revokeObjectURL(url);
    }
};

function printReport() { window.print(); }
</script>
"""


def _render_summary_table(entries: list[dict]) -> str:
    rows = [
        "<table><tr><th>Filnavn</th><th>Ladder</th><th>Trace-kanaler</th><th>bp-område</th><th>Ladder QC</th><th>R²</th></tr>"
    ]
    for entry in entries:
        r2 = entry.get("ladder_r2")
        r2_str = f"{r2:.4f}" if isinstance(r2, (int, float)) and np.isfinite(r2) else "&mdash;"
        status = entry.get("ladder_qc_status", "unknown")
        status_label = {
            "ok": "<span class='status-badge ok'>OK</span>",
            "manual_adjustment": "<span class='status-badge manual'>Manual</span>",
            "review_required": "<span class='status-badge warning'>Warning</span>",
            "ladder_qc_failed": "<span class='status-badge failed'>Failed</span>",
        }.get(status, "<span class='status-badge unknown'>Unknown</span>")
        channels = ", ".join(entry.get("trace_channels") or [])
        rows.append(
            "<tr>"
            f"<td>{escape(entry['fsa'].file_name)}</td>"
            f"<td>{escape(str(entry.get('ladder', '')))}</td>"
            f"<td>{escape(channels)}</td>"
            f"<td>{float(entry.get('bp_min', 0.0)):.0f}–{float(entry.get('bp_max', 0.0)):.0f} bp</td>"
            f"<td>{status_label}</td>"
            f"<td>{r2_str}</td>"
            "</tr>"
        )
    rows.append("</table>")
    return "\n".join(rows)


def _comment_block(label: str) -> str:
    return (
        "<div class='comment-box-container'>"
        "<button class='comment-toggle-btn' onclick='toggleComment(this)'>"
        "💬 <span class='comment-label'>Legg til kommentar</span>"
        f" <em style='font-weight:400;opacity:0.7;'>({escape(label)})</em>"
        "<i class='caret'>&#x25BC;</i>"
        "</button>"
        "<div class='comment-body'>"
        "<textarea class='report-comment' placeholder='Skriv inn eventuelle kommentarer her...'></textarea>"
        "</div>"
        "</div>"
    )


def build_general_html_report(entries: list[dict], assay_outdir: Path, run_label: str | None = None) -> Path | None:
    if not entries:
        return None

    assay_outdir.mkdir(parents=True, exist_ok=True)
    safe_label = _safe_report_name(run_label or assay_outdir.name or "general")
    out_html = assay_outdir / f"{safe_label}_General_Report.html"

    html_lines: list[str] = []
    html_lines.extend(["<!DOCTYPE html>", "<html lang='no'>", "<head>", "<meta charset='utf-8'>"])
    html_lines.append(f"<title>{escape(safe_label)}_General_Report</title>")
    html_lines.append(REPORT_STYLE)
    html_lines.append('<script id="peak-data" type="application/json">{}</script>')
    html_lines.append('<script id="plot-state" type="application/json">{}</script>')
    html_lines.append(_build_report_script())
    html_lines.append(_local_plotly_tag(assay_outdir, version="2.35.2"))
    html_lines.append(_build_plotly_reflow_script())
    html_lines.extend(["</head>", "<body>"])

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    html_lines.append(
        f"""
<div class='report-header no-print'>
  <h1>{escape(safe_label)}_General_Report</h1>
  <div class='meta'>{len(entries)} analyserte filer &nbsp;&bull;&nbsp; Generert: {generated}</div>
</div>
"""
    )
    html_lines.append("<h2>Oversikt</h2>")
    html_lines.append("<p class='small'>Generell analyse med valgte trace-kanaler, manuell peakredigering og lagret kommentarstøtte.</p>")
    html_lines.append(_render_summary_table(entries))
    html_lines.append(_comment_block("Run-level kommentar"))

    for entry in entries:
        fsa = entry["fsa"]
        html_lines.append("<div class='assay-block'>")
        html_lines.append(f"<h3>{escape(fsa.file_name)}</h3>")
        html_lines.append(
            "<p class='small'>"
            f"Ladder: {escape(str(entry.get('ladder', '')))} | "
            f"Trace-kanaler: {escape(', '.join(entry.get('trace_channels') or []))} | "
            f"Primærkanal: {escape(str(entry.get('primary_peak_channel', '')))}"
            "</p>"
        )
        fragment = build_interactive_peak_plot_for_entry(entry)
        html_lines.append(fragment if fragment else "<p class='small'><em>Ingen data å vise.</em></p>")
        html_lines.append(_comment_block(fsa.file_name))
        html_lines.append("</div>")

    html_lines.append(
        """
<div class="print-fab no-print">
  <button class="print-btn save-peaks-btn" onclick="PeakManager.downloadUpdatedHtml()">💾&nbsp; Save Peaks</button>
  <button class="print-btn" onclick="printReport()">🖨&nbsp; Print / PDF</button>
</div>
</body></html>"""
    )

    out_html.write_text("\n".join(html_lines), encoding="utf-8")
    return out_html
