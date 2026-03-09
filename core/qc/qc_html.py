"""
Fraggler QC — HTML QC report builder.
"""
from __future__ import annotations

from pathlib import Path
from html import escape as html_escape
from datetime import datetime
from typing import Any

import numpy as np

from core.plotly_offline import local_plotly_tag as _local_plotly_tag
from core.qc.qc_rules import QCRules, normalize_assay_qc
from core.qc.qc_markers import (
    control_id_from_filename,
    ladder_qc_grade,
    worst_grade,
)
from core.qc.qc_plots import build_interactive_peak_plot_for_entry_qc
from core.qc.qc_excel import update_excel_trends, apply_pk_excel_styling
from core.html_reports import interpret_sl_quality
import core.assay_config as master


def build_qc_html(entries: list[dict], out_html: Path, rules: QCRules, excel_path: Path) -> None:
    out_html.parent.mkdir(parents=True, exist_ok=True)

    by_assay: dict[str, list[dict]] = {}
    for e in entries:
        a = normalize_assay_qc(e.get("assay", "UNKNOWN"))
        by_assay.setdefault(a, []).append(e)

    ordered = []
    for a in getattr(master, "ASSAY_DISPLAY_ORDER", []):
        if a in by_assay:
            ordered.append(a)
    for a in sorted(by_assay.keys()):
        if a not in ordered:
            ordered.append(a)

    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    html: list[str] = []
    html.append("<!DOCTYPE html>")
    html.append("<html lang='no'>")
    html.append("<head>")
    html.append("<meta charset='utf-8'>")
    html.append("<title>QC – Fraggler</title>")
    html.append("""
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
tr:hover td { background: #f0fdfa; transition: background 0.15s ease; }

/* ── Cards ── */
.assay-block, .block {
    padding: 18px 22px;
    margin-bottom: 20px;
    background: #ffffff;
    border-radius: 12px;
    border: 1px solid rgba(226, 232, 240, 0.8);
    box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.05), 0 8px 10px -6px rgba(0, 0, 0, 0.01);
    transition: box-shadow 0.3s ease;
}
.assay-block:hover, .block:hover { box-shadow: 0 14px 28px -5px rgba(0, 0, 0, 0.08); }
.sample-header { font-size: 0.9rem; font-weight: 700; margin-top: 0.4rem; color: #0f172a; }
.small { font-size: 0.85rem; color: #64748b; font-weight: 500; }
.peak-editor-block { margin-top: 0.5rem; margin-bottom: 1.2rem; border-radius: 8px; overflow: hidden; }

/* ── Badges ── */
.badge { display: inline-block; padding: 0.1rem 0.45rem; border-radius: 999px; font-size: 0.8rem; font-weight: 600; }
.ok { background: #e6f4ea; color: #0f5132; border: 1px solid #b7e1c2; }
.warn { background: #fff4e5; color: #7a4b00; border: 1px solid #ffd59a; }
.fail { background: #fde7e9; color: #842029; border: 1px solid #f5b5bb; }
.na { background: #eeeeee; color: #444; border: 1px solid #d5d5d5; }

hr { border: none; border-top: 1px solid #e6e6e6; margin: 1.2rem 0; }
</style>
""")
    # Plotly CDN (samme versjon master bruker). [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    html.append(_local_plotly_tag(out_html.parent, version="2.35.2"))
    html.append("</head>")
    html.append("<body>")

    html.append("<h1>QC – Fraggler</h1>")
    html.append(f"<p class='small'>Generert: {now}. Antall QC-filer: <strong>{len(entries)}</strong>.</p>")
    html.append(f"<p class='small'>Excel trends: <code>{html_escape(excel_path.name)}</code> (oppdateres ved hver kjøring).</p>")

    # Oversiktstabell
    html.append("<h2>Oversikt</h2>")
    html.append("<table>")
    html.append(
        "<tr>"
        "<th>Control</th><th>Filnavn</th><th>Assay</th><th>Ladder</th>"
        "<th>Ladder QC</th><th>R²</th>"
        "</tr>"
    )

    all_entries = sorted(entries, key=lambda e: (e.get("assay", ""), e["fsa"].file_name))
    for e in all_entries:
        fsa_name = e["fsa"].file_name
        ctrl = control_id_from_filename(fsa_name)
        assay = normalize_assay_qc(e.get("assay", "UNKNOWN"))
        ladder = e.get("ladder", "—")
        r2 = e.get("ladder_r2", None)

        grade, note = ladder_qc_grade(r2, rules)
        badge_cls = grade.lower()
        r2_str = "—" if (r2 is None or not np.isfinite(r2)) else f"{float(r2):.4f}"

        sl_txt = "—"
        if assay == "SL":
            slm = e.get("sl_metrics")
            if slm and isinstance(slm, dict):
                pcts = slm.get("percents", [])
                total_area = slm.get("total_area", float("nan"))
                try:
                    from core.html_reports import interpret_sl_quality
                    sl_txt = interpret_sl_quality(pcts, total_area)
                except Exception:
                    sl_txt = "—"

        html.append(
            "<tr>"
            f"<td>{html_escape(ctrl)}</td>"
            f"<td>{html_escape(fsa_name)}</td>"
            f"<td>{html_escape(assay)}</td>"
            f"<td>{html_escape(str(ladder))}</td>"
            f"<td><span class='badge {badge_cls}' title='{html_escape(note)}'>{grade}</span> {html_escape(note)}</td>"
            f"<td>{html_escape(r2_str)}</td>"
            "</tr>"
        )
    html.append("</table>")

    # Per assay
    html.append("<h2>QC per assay</h2>")
    html.append("<p class='small'>Markører: lilla = sample-peak, gul = ladder-peak (vertikale linjer = expected bp).</p>")

    for assay in ordered:
        a_entries = sorted(by_assay.get(assay, []), key=lambda e: e["fsa"].file_name)
        if not a_entries:
            continue

        ref_label = getattr(master, "ASSAY_REFERENCE_LABEL", {}).get(assay, "")
        html.append("<div class='block'>")
        html.append(f"<h3>{html_escape(assay)}</h3>")
        if ref_label:
            html.append(f"<p class='small'><strong>Referanse:</strong> {html_escape(ref_label)}</p>")

        # Figurer per entry
        for e in a_entries:
            fsa_name = e["fsa"].file_name
            primary = e.get("primary_peak_channel", "")
            html.append(f"<p class='small'><strong>{html_escape(fsa_name)}</strong> ({html_escape(primary)})</p>")
            try:
                frag = build_interactive_peak_plot_for_entry_qc(e, rules)
                if frag is None:
                    html.append("<p class='small'><em>Ingen data å vise.</em></p>")
                else:
                    html.append(frag)
            except Exception as ex:
                html.append(f"<p class='small'><em>Kunne ikke lage plott: {html_escape(str(ex))}</em></p>")

        html.append("</div>")
        html.append("<hr>")

    # -----------------------------------------------------
    # 3) SIZE LADDER DATA (for QC files som er SL)
    # -----------------------------------------------------
    qc_sl_entries = [
        e for e in all_entries 
        if normalize_assay_qc(e.get("assay", "UNKNOWN")) == "SL"
        and control_id_from_filename(e["fsa"].file_name) in ("PK", "PK1", "PK2")
    ]
    if qc_sl_entries:
        html.append("<h2>Size Ladder (SL) – Fragmentfordeling</h2>")
        html.append("<p class='small'>Målte verdier per SL-fil i denne QC-kjøringen.</p>")

        for e in qc_sl_entries:
            sl_metrics = e.get("sl_metrics")
            fsa_name = html_escape(e["fsa"].file_name)
            html.append(f"<h3>SL-fil: {fsa_name}</h3>")

            if not sl_metrics:
                html.append("<p><em>Ingen SL-area-metrikker tilgjengelig for denne fila.</em></p>")
                continue

            targets_bp = sl_metrics.get("targets_bp", [])
            areas = sl_metrics.get("areas", [])
            percents = sl_metrics.get("percents", [])
            total_area = sl_metrics.get("total_area", float("nan"))

            html.append("<table>")
            html.append(
                "<tr>"
                "<th>Fragment (bp)</th>"
                "<th>Area</th>"
                "<th>% av total</th>"
                "</tr>"
            )

            for bp_val, area_val, pct_val in zip(targets_bp, areas, percents):
                area_str = "&mdash;" if np.isnan(area_val) else f"{area_val:,.0f}".replace(",", " ")
                pct_str = "&mdash;" if (pct_val is None or np.isnan(pct_val)) else f"{pct_val:.1f} %"
                html.append(
                    "<tr>"
                    f"<td>{bp_val} bp</td>"
                    f"<td>{area_str}</td>"
                    f"<td>{pct_str}</td>"
                    "</tr>"
                )

            tot_area_str = "&mdash;" if np.isnan(total_area) else f"{total_area:,.0f}".replace(",", " ")
            html.append(
                f"<tr style='background-color:#eef2f6; font-weight:bold'>"
                f"<td>TOTAL</td>"
                f"<td>{tot_area_str}</td>"
                f"<td>100 %</td>"
                "</tr>"
            )
            html.append("</table>")
            
    html.append("</body></html>")
    out_html.write_text("\n".join(html), encoding="utf-8")
