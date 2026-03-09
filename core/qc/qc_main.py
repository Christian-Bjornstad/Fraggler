"""
Fraggler QC — CLI entry point.
"""
from __future__ import annotations

from pathlib import Path

from core.qc.qc_rules import QCRules
from core.qc.qc_html import build_qc_html
from core.qc.qc_excel import update_excel_trends, apply_pk_excel_styling
from core.pipeline import run_pipeline
from core.assay_config import DEFAULT_FSA_DIR, OUTDIR_NAME
from datetime import datetime


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Fraggler QC report + marker tracking + Excel trends.")
    parser.add_argument("--fsa_dir", type=str, default=str(DEFAULT_FSA_DIR),
                        help="Mappe med .fsa-filer")
    parser.add_argument("--outdir", type=str, default="",
                        help="Output-mappe. Hvis tom: bruker fsa_dir/ASSAY_REPORTS")
    parser.add_argument("--outfile", type=str, default="QC_REPORT.html",
                        help="Navn på QC-HTML-filen")
    parser.add_argument("--excel", type=str, default="QC_TRENDS.xlsx",
                        help="Excel-filnavn (lagres i output-mappa)")

    # QC terskler
    parser.add_argument("--min_r2_ok", type=float, default=0.999)
    parser.add_argument("--min_r2_warn", type=float, default=0.995)
    parser.add_argument("--sample_peak_window_bp", type=float, default=2.0)
    # Peak søkevindu
    parser.add_argument("--w_sample", type=float, default=3.0, help="± bp-vindu for sample-markører")
    parser.add_argument("--w_ladder", type=float, default=3.0, help="± bp-vindu for ladder-markører")

    args = parser.parse_args()

    fsa_dir = Path(args.fsa_dir).expanduser()
    if not fsa_dir.exists():
        raise SystemExit(f"Finner ikke fsa_dir: {fsa_dir}")

    rules = QCRules(
        min_r2_ok=args.min_r2_ok,
        min_r2_warn=args.min_r2_warn,
        sample_peak_window_bp=args.sample_peak_window_bp,
        ladder_peak_window_bp=args.w_ladder,
    )

    # Output
    if args.outdir.strip():
        base_outdir = Path(args.outdir).expanduser()
    else:
        base_outdir = fsa_dir

    assay_folder = OUTDIR_NAME
    qc_outdir = base_outdir / assay_folder
    qc_outdir.mkdir(parents=True, exist_ok=True)

    out_html = qc_outdir / args.outfile
    excel_path = qc_outdir / args.excel

    # Kjør master i controls-modus, return entries (samme datastruktur som pasient-flyt). [1](https://hsorhf-my.sharepoint.com/personal/chrbj5_ous-hf_no/Documents/Microsoft%20Copilot%20Chat-filer/fraggler_master_assay_channels.py)
    entries = run_pipeline(
        fsa_dir=fsa_dir,
        base_outdir=base_outdir,
        return_entries=True,
        make_dit_reports=False,
        mode="controls",
    )

    if not entries:
        print("\n[QC] Ingen QC entries funnet. (Sjekk filnavn-prefix: PK/NK/RK/PK1/PK2_ ...)\n")
        return

    # Bygg HTML (den beregner også marker results per entry)
    build_qc_html(entries, out_html, rules, excel_path)

    # Oppdater Excel etter at marker results finnes i entries
    run_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    update_excel_trends(excel_path, entries, rules, run_ts)

    print(f"\n[QC] Lagret QC-rapport: {out_html}")
    print(f"[QC] Oppdatert Excel:   {excel_path}\n")


if __name__ == "__main__":
    main()