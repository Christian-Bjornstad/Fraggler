"""
Fraggler Diagnostics — DIT HTML Reports & SL Quality Interpretation.

Builds per-patient DIT HTML reports with embedded interactive Plotly figures.
"""
from __future__ import annotations

import re
import uuid
from pathlib import Path
from collections import defaultdict
from html import escape
from datetime import datetime

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
    OUTDIR_NAME,
)
from core.plotly_offline import local_plotly_tag as _local_plotly_tag
from core.plotting_plotly import (
    compute_group_ymax_for_entries,
    build_interactive_peak_plot_for_entry,
)
from config import APP_SETTINGS


DIT_PATTERN = re.compile(r"(\d{2}OUM\d{5})")

def extract_dit_from_name(name: str) -> str | None:
    """
    Finner første forekomst av 2-sifret år + 'OUM' + 5 siffer i filnavnet.
    Eksempler: '25OUM10166', '26OUM00042'.

    Returnerer DIT-strengen, eller None hvis ingen match.
    """
    m = DIT_PATTERN.search(name)
    if m:
        return m.group(1)
    return None

def dit_to_year(dit: str) -> int | None:
    """
    25OUM10166 -> 2025, 26OUMxxxxx -> 2026, etc.
    """
    if not dit or len(dit) < 2:
        return None
    try:
        yy = int(dit[:2])
    except ValueError:
        return None
    return 2000 + yy


def build_dit_html_reports(
    entries: list[dict],
    assay_outdir: Path,
):
    """
    Ny, forenklet DIT-rapport:

    - Én HTML per DIT.
    - Header + tabell over alle filer (QC, ladder-info osv.).
    - For hver assay og fil embeddes en interaktiv Plotly-figur direkte.
    - I tillegg: egne kombinasjonsblokker for TCRb og TCRg.
    """

    # -----------------------------

    # Gruppér alle entries per DIT
    qc_sl_entries = []
    per_dit: dict[str, list[dict]] = defaultdict(list)
    from core.qc.qc_markers import control_id_from_filename
    for e in entries:
        ctrl_id = control_id_from_filename(e["fsa"].file_name)
        if ctrl_id in ("PK", "PK1", "PK2") and e.get("assay") == "SL":
            qc_sl_entries.append(e)
            
        dit = e.get("dit")
        if not dit:
            continue
        per_dit[dit].append(e)

    if not per_dit:
        print_warning("[DIT] Fant ingen DIT-nummer i filnavnene – genererer ingen pasient-HTML.")
        return

    dit_root = assay_outdir
    dit_root.mkdir(exist_ok=True, parents=True)

    print_green(f"[DIT] Lager pasientrapporter i {dit_root}")

    for dit, dit_entries in sorted(per_dit.items()):
        year = dit_to_year(dit)

        # Gruppér innenfor DIT per assay
        assays: dict[str, list[dict]] = defaultdict(list)
        for e in dit_entries:
            assays[e["assay"]].append(e)

        html_lines: list[str] = []
        html_lines.append("<!DOCTYPE html>")
        html_lines.append("<html lang='no'>")
        html_lines.append("<head>")
        html_lines.append("<meta charset='utf-8'>")
        html_lines.append(f"<title>DIT {escape(dit)} – Fraggler-rapport</title>")
        html_lines.append(
            """
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
    /* Keep Plotly interactive modebar hidden in print */
    .modebar { display: none !important; }
    .js-plotly-plot .plotly .modebar { display: none !important; }
}
</style>
<script>
function printReport() {
    window.print();
}
</script>
            """
        )

        # Plotly inn her – ÉN gang
        
        html_lines.append(_local_plotly_tag(dit_root, version="2.35.2"))
        html_lines.append("</head>")
        html_lines.append("<body>")

        # HEADER
        gen_date = datetime.now().strftime("%Y-%m-%d %H:%M")
        meta_parts = []
        if year is not None:
            meta_parts.append(f"År: {year}")
        meta_parts.append(f"{len(dit_entries)} analyserte filer")
        meta_parts.append(f"Generert: {gen_date}")
        html_lines.append(f"""
<div class='report-header no-print'>
  <h1>Fraggler-rapport &mdash; DIT {escape(dit)}</h1>
  <div class='meta'>{" &nbsp;&bull;&nbsp; ".join(meta_parts)}</div>
</div>
<div style='display:none' class='print-only-header'>
  <h1>Fraggler-rapport &mdash; DIT {escape(dit)}</h1>
  <p>{" | ".join(meta_parts)}</p>
</div>
""")

        # 1) TABELL OVER ALLE FILER
        html_lines.append("<h2>Oversikt over analyserte filer</h2>")
        html_lines.append(
            "<p class='small'>Tabellen viser alle filer for dette DIT-nummeret, "
            "inkludert assay, ladder-type, bp-område og ladder-QC.</p>"
        )

        html_lines.append("<table>")
        html_lines.append(
            "<tr>"
            "<th>Filnavn</th>"
            "<th>Assay</th>"
            "<th>Ladder</th>"
            "<th>bp-område</th>"
            "<th>Ladder QC</th>"
            "<th>R²</th>"
            "</tr>"
        )
        for e in sorted(dit_entries, key=lambda x: (x["assay"], x["fsa"].file_name)):
            status = e.get("ladder_qc_status", "unknown")
            r2 = e.get("ladder_r2", None)
            if r2 is None or np.isnan(r2):
                r2_str = "&mdash;"
            else:
                r2_str = f"{r2:.4f}"

            html_lines.append(
                "<tr>"
                f"<td>{escape(e['fsa'].file_name)}</td>"
                f"<td>{escape(e['assay'])}</td>"
                f"<td>{escape(e['ladder'])}</td>"
                f"<td>{int(e['bp_min'])}–{int(e['bp_max'])} bp</td>"
                f"<td>{escape(status)}</td>"
                f"<td>{r2_str}</td>"
                "</tr>"
            )
        html_lines.append("</table>")

        # 2) ASSAY-SPESIFIKKE PLOTT (vanligt oppsett)
        html_lines.append("<h2>Assay-spesifikke oversikter</h2>")

        assays_present = list(assays.keys())
        ordered_assays = [a for a in ASSAY_DISPLAY_ORDER if a in assays_present]
        for a in assays_present:
            if a not in ASSAY_DISPLAY_ORDER:
                ordered_assays.append(a)

        # Vi bruker disse listene til å vite når vi skal putte inn kombinasjonsfigurer
        tcrb_assays = ["TCRbA", "TCRbB", "TCRbC"]
        tcrg_assays = ["TCRgA", "TCRgB"]

        for assay_name in ordered_assays:
            assay_entries = assays[assay_name]
            html_lines.append("<div class='assay-block'>")
            html_lines.append(f"<h3>{escape(assay_name)}</h3>")

            # Referanseområder / tekst
            ref_ranges = ASSAY_REFERENCE_RANGES.get(assay_name)
            if ref_ranges:
                ranges_str = ", ".join(
                    f"{int(a)}–{int(b)} bp" for (a, b) in ref_ranges
                )
                label_txt = ASSAY_REFERENCE_LABEL.get(assay_name, ranges_str)
                html_lines.append(
                    f"<p class='small'><strong>Referanseområde:</strong> "
                    f"{escape(ranges_str)}<br>{escape(label_txt)}</p>"
                )

            # Vanlige enkel-plot per fil (INGEN forced_ymax her)
            for e in sorted(assay_entries, key=lambda x: x["fsa"].file_name):
                fsa = e["fsa"]
                primary_ch = e["primary_peak_channel"]

                html_lines.append(
                    f"<p class='sample-header'>{escape(fsa.file_name)} "
                    f"({escape(primary_ch)})</p>"
                )

                try:
                    frag = build_interactive_peak_plot_for_entry(e)
                    if frag is None:
                        html_lines.append(
                            "<p class='small'><em>Ingen data å vise for denne fila.</em></p>"
                        )
                    else:
                        html_lines.append(frag)
                except Exception as ex:
                    html_lines.append(
                        f"<p class='small'><em>Kunne ikke lage interaktivt plott for "
                        f"{escape(fsa.file_name)}: {escape(str(ex))}</em></p>"
                    )

            html_lines.append("</div>")  # assay-block

            # --------------------------------------------------
            # Rett etter TCRbC-blokken: TCRβ-parallell 1 og 2
            # --------------------------------------------------
            if assay_name == "TCRbC":
                tcrb_present = [a for a in tcrb_assays if a in assays]
                if tcrb_present:
                    html_lines.append("<h2>Kombinasjonsfigurer – TCRβ</h2>")

                    # sorter per assay på filnavn (antatt 2 replikater)
                    tcrb_sorted = {
                        a: sorted(assays[a], key=lambda x: x["fsa"].file_name)
                        for a in tcrb_present
                    }

                    # repl 1 = første fil i hver assay, repl 2 = andre fil (hvis finnes)
                    rep1 = [lst[0] for a, lst in tcrb_sorted.items() if len(lst) >= 1]
                    rep2 = [lst[1] for a, lst in tcrb_sorted.items() if len(lst) >= 2]

                    def render_tcrb_rep_block(rep_entries: list[dict], rep_label: str):
                        if not rep_entries:
                            return
                        # felles ymax for denne replikat-blokken (bruk referansevindu + primærkanal)
                        group_y = compute_group_ymax_for_entries(rep_entries)
                        html_lines.append("<div class='assay-block'>")
                        html_lines.append(
                            f"<h3>TCRβ – parallell {escape(rep_label)} (mix A + B + C)</h3>"
                        )
                        html_lines.append("<div class='combo-grid'>")
                        for e in sorted(rep_entries, key=lambda x: x["assay"]):
                            fsa = e["fsa"]
                            primary_ch = e["primary_peak_channel"]
                            # bruk kopi slik at forced_ymax IKKE lekker til enkel-plottene
                            e_combo = dict(e)
                            e_combo["forced_ymax"] = group_y

                            html_lines.append("<div class='combo-item'>")
                            html_lines.append(
                                f"<p class='sample-header'>{escape(e_combo['assay'])} – "
                                f"{escape(fsa.file_name)} ({escape(primary_ch)})</p>"
                            )
                            try:
                                frag = build_interactive_peak_plot_for_entry(e_combo)
                                if frag is None:
                                    html_lines.append(
                                        "<p class='small'><em>Ingen data å vise.</em></p>"
                                    )
                                else:
                                    html_lines.append(frag)
                            except Exception as ex:
                                html_lines.append(
                                    f"<p class='small'><em>Kunne ikke lage plott: "
                                    f"{escape(str(ex))}</em></p>"
                                )
                            html_lines.append("</div>")  # combo-item

                        html_lines.append("</div>")  # combo-grid
                        html_lines.append("</div>")  # assay-block

                    # Parallell 1 og 2 (rep1/rep2)
                    render_tcrb_rep_block(rep1, "1")
                    render_tcrb_rep_block(rep2, "2")

            # --------------------------------------------------
            # Rett etter TCRgB-blokken: TCRγ mix A + B kombinasjon
            # --------------------------------------------------
            if assay_name == "TCRgB":
                # Samle alle TCRg A/B-entries for denne DIT
                tcrg_entries_all: list[dict] = []
                for a in tcrg_assays:
                    if a in assays:
                        tcrg_entries_all.extend(
                            sorted(assays[a], key=lambda x: x["fsa"].file_name)
                        )

                if tcrg_entries_all:
                    group_y = compute_group_ymax_for_entries(tcrg_entries_all)
                    html_lines.append("<h2>Kombinasjonsfigur – TCRγ</h2>")
                    html_lines.append("<div class='assay-block'>")
                    html_lines.append(
                        "<p class='small'>TCRγ mix A og mix B (begge paralleller) med felles y-akse for enkel sammenligning.</p>"
                    )
                    html_lines.append("<div class='combo-grid'>")

                    for e in tcrg_entries_all:
                        fsa = e["fsa"]
                        primary_ch = e["primary_peak_channel"]
                        e_combo = dict(e)
                        e_combo["forced_ymax"] = group_y

                        html_lines.append("<div class='combo-item'>")
                        html_lines.append(
                            f"<p class='sample-header'>{escape(e_combo['assay'])} – "
                            f"{escape(fsa.file_name)} ({escape(primary_ch)})</p>"
                        )
                        try:
                            frag = build_interactive_peak_plot_for_entry(e_combo)
                            if frag is None:
                                html_lines.append(
                                    "<p class='small'><em>Ingen data å vise.</em></p>"
                                )
                            else:
                                html_lines.append(frag)
                        except Exception as ex:
                            html_lines.append(
                                f"<p class='small'><em>Kunne ikke lage plott: "
                                f"{escape(str(ex))}</em></p>"
                            )
                        html_lines.append("</div>")  # combo-item

                    html_lines.append("</div>")  # combo-grid
                    html_lines.append("</div>")  # assay-block

        # 3) SIZE LADDER-SEKSJON (som før)
        sl_entries = [e for e in dit_entries if e.get("assay") == "SL"]
        all_sl_entries = sl_entries + qc_sl_entries
        if all_sl_entries:
            html_lines.append("<h2>Size Ladder (SL) – DNA-kvalitet</h2>")
            html_lines.append(
                "<p class='small'>Målte verdier per SL-fil.</p>"
            )

            for e in sorted(all_sl_entries, key=lambda x: x["fsa"].file_name):
                sl_metrics = e.get("sl_metrics")
                html_lines.append(
                    f"<h3>SL-fil: {escape(e['fsa'].file_name)}</h3>"
                )

                if not sl_metrics:
                    html_lines.append(
                        "<p><em>Ingen SL-area-metrikker tilgjengelig for denne fila.</em></p>"
                    )
                    continue

                pcts = sl_metrics.get("percents", [])
                total_area = sl_metrics.get("total_area", float("nan"))

                targets_bp = sl_metrics.get("targets_bp", [])
                areas = sl_metrics.get("areas", [])
                percents = sl_metrics.get("percents", [])
                total_area = sl_metrics.get("total_area", float("nan"))

                html_lines.append("<table>")
                html_lines.append(
                    "<tr>"
                    "<th>Fragment (bp)</th>"
                    "<th>Area</th>"
                    "<th>% av total</th>"
                    "</tr>"
                )

                for bp_val, area_val, pct_val in zip(targets_bp, areas, percents):
                    if np.isnan(area_val):
                        area_str = "&mdash;"
                    else:
                        area_str = f"{area_val:,.0f}".replace(",", " ")
                    if pct_val is None or np.isnan(pct_val):
                        pct_str = "&mdash;"
                    else:
                        pct_str = f"{pct_val:.1f} %"

                    html_lines.append(
                        "<tr>"
                        f"<td>{bp_val:.0f}</td>"
                        f"<td>{area_str}</td>"
                        f"<td>{pct_str}</td>"
                        "</tr>"
                    )

                if np.isnan(total_area):
                    tot_str = "&mdash;"
                else:
                    tot_str = f"{total_area:,.0f}".replace(",", " ")

                html_lines.append(
                    "<tr>"
                    "<td><strong>Total</strong></td>"
                    f"<td><strong>{tot_str}</strong></td>"
                    "<td></td>"
                    "</tr>"
                )
                html_lines.append("</table>")

        html_lines.append("""
<div class="print-fab no-print">
  <button class="print-btn" onclick="printReport()">
    🖨&nbsp; Print / Save PDF
  </button>
</div>
</body></html>""")

        out_html = dit_root / f"{dit}.html"
        out_html.write_text("\n".join(html_lines), encoding="utf-8")
        print_green(f"[DIT] Lagret pasientrapport for {dit}: {out_html}")


def interpret_sl_quality(percents, total_area):
    """
    Automatisk fortolkning av DNA-kvalitet basert på fragmentfordeling.

    percents = liste over %-andeler for [100, 200, 300, 400, 600] bp.
    total_area = summert area (kan brukes for å fange ekstremt svake prøver).
    """

    # Pakk ut med beskyttelse
    p100 = percents[0] if len(percents) > 0 else float("nan")
    p200 = percents[1] if len(percents) > 1 else float("nan")
    p300 = percents[2] if len(percents) > 2 else float("nan")
    p400 = percents[3] if len(percents) > 3 else float("nan")
    p600 = percents[4] if len(percents) > 4 else float("nan")

    # ------------------------------------------------------------------
    # 0) Tekniske "ugyldig/uegnet"-tilfeller
    # ------------------------------------------------------------------

    # Ekstremt lite signal (her kan du sette terskel etter behov)
    if np.isnan(total_area) or total_area < 1e4:
        return "Materialet er uegnet for PCR-analyser (svært lite signal / teknisk mislykket)."

    # Svært svak/negativ minstepeak
    if np.isnan(p100) or p100 < 5:
        return "Materialet er uegnet for PCR-analyser (svært svak/negativ 100 bp-peak)."

    # ------------------------------------------------------------------
    # 1) Svært fragmentert
    #    Tabell: 90 % på 100 bp, 10 % på 200 bp ⇒ nesten alt i 100–200 bp.
    # ------------------------------------------------------------------
    if p100 >= 85 and p200 <= 15 and p300 <= 5:
        return "Svært fragmentert materiale. Uegnet for de fleste analyser."

    # ------------------------------------------------------------------
    # 2) >50 % fragmentert
    #    Tabell: 60 % på 100 bp, 30 % på 200 bp, 10 % på 300 bp.
    #    Her tolker vi det som: veldig mye ligger ≤300 bp.
    # ------------------------------------------------------------------
    sum_100_300 = p100 + p200 + p300
    sum_100_200 = p100 + p200

    if p100 >= 60 and sum_100_200 >= 80 and p300 <= 15:
        return "Mer enn 50 % av materialet er fragmentert – redusert sensitivitet."

    # ------------------------------------------------------------------
    # 3) Litt fragmentert
    #    Tabell: 50 % / 40 % / 30 % + litt 400 bp.
    #    Vi tolker: mye i 100–300, men fortsatt noe "tyngre" materiale.
    # ------------------------------------------------------------------
    if p100 >= 45 and sum_100_300 >= 70:
        return "Litt fragmentert – kan redusere sensitiviteten på enkelte analyser."

    # ------------------------------------------------------------------
    # 4) Bra kvalitet
    #    Tabell: 40/30/20/10/5.
    #    Vi sier: moderat 100 bp, god andel i 200–400, noe i 600.
    # ------------------------------------------------------------------
    if (
        p100 <= 50 and
        sum_100_200 <= 70 and
        p300 >= 10 and
        p400 >= 5
    ):
        return "Bra kvalitet."

    # ------------------------------------------------------------------
    # 5) Fallback
    # ------------------------------------------------------------------
    return "Uvanlig fordeling – vurder manuelt."
