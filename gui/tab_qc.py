"""
Fraggler Diagnostics — QC Tab
"""
from __future__ import annotations

import panel as pn

from core.runner import executor, run_qc_job
from config import APP_SETTINGS, save_settings
from gui.components import make_card, VSpace, section_header


def make_qc_tab() -> pn.Column:
    s = APP_SETTINGS.setdefault("qc", {})

    input_dir = pn.widgets.TextInput(
        name="Input Folder (.fsa controls)",
        value=s.get("input_dir", ""),
        sizing_mode="stretch_width",
        placeholder="/path/to/controls/folder"
    )
    browse_in_btn = pn.widgets.Button(name="Browse...", width=90, height=32, align="end", margin=(0, 0, 4, 0))

    output_base = pn.widgets.TextInput(
        name="Output Folder (leave empty = same as input)",
        value=s.get("output_base", ""),
        sizing_mode="stretch_width",
    )
    browse_out_btn = pn.widgets.Button(name="Browse...", width=90, height=32, align="end", margin=(0, 0, 4, 0))
    smart_scan_btn = pn.widgets.Button(name="Smart Scan", width=110, height=32, align="end", button_type="default", margin=(0, 0, 4, 0))

    def _ask_dir(target_widget):
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder = filedialog.askdirectory(initialdir=target_widget.value or "~")
        root.destroy()
        if folder:
            target_widget.value = folder

    def _smart_scan(e):
        from pathlib import Path
        import re
        base = Path(input_dir.value or "~").expanduser()
        if not base.exists() or not base.is_dir():
            return
        
        # Look for folders containing .fsa files that might be controls
        found_dirs = []
        for d in base.rglob("*"):
            if d.is_dir() and not d.name.startswith("."):
                fsa_files = list(d.glob("*.fsa"))
                if not fsa_files:
                    continue
                
                # If folder name contains 'kontroll' or 'control'
                if re.search(r"kontroll|control", d.name, re.I):
                    found_dirs.append(d)
                # Or if any file starts with PK/NK/RK
                elif any(re.match(r"^(PK|NK|RK)", f.name, re.I) for f in fsa_files):
                    found_dirs.append(d)
        
        if found_dirs:
            # Pick the best match (shortest path usually)
            best = sorted(found_dirs, key=lambda x: len(str(x)))[0]
            input_dir.value = str(best)

    browse_in_btn.on_click(lambda e: _ask_dir(input_dir))
    browse_out_btn.on_click(lambda e: _ask_dir(output_base))
    smart_scan_btn.on_click(_smart_scan)
    outfile_html = pn.widgets.TextInput(name="HTML report filename", value=s.get("outfile_html", "QC_REPORT.html"), width=280)
    excel_name = pn.widgets.TextInput(name="Excel filename", value=s.get("excel_name", "QC_TRENDS.xlsx"), width=280)

    min_r2_ok    = pn.widgets.FloatInput(name="min R² (OK)",   value=float(s.get("min_r2_ok", 0.995)),  step=0.001, width=140)
    min_r2_warn  = pn.widgets.FloatInput(name="min R² (WARN)", value=float(s.get("min_r2_warn", 0.990)), step=0.001, width=140)
    max_mse_ok   = pn.widgets.FloatInput(name="max MSE (OK)",  value=float(s.get("max_mse_ok", 2.0)),   step=0.1,   width=140)
    max_mse_warn = pn.widgets.FloatInput(name="max MSE (WARN)",value=float(s.get("max_mse_warn", 5.0)), step=0.1,   width=140)
    nk_ymax      = pn.widgets.FloatInput(name="NK y-max floor",value=float(s.get("nk_ymax_floor", 250.0)), step=10.0, width=140)
    w_sample     = pn.widgets.FloatInput(name="±bp window (sample)", value=float(s.get("w_sample", 3.0)), step=0.5, width=140)
    w_ladder     = pn.widgets.FloatInput(name="±bp window (ladder)", value=float(s.get("w_ladder", 3.0)), step=0.5, width=140)

    run_btn  = pn.widgets.Button(name="Run QC", button_type="primary", width=140, height=40)
    open_btn = pn.widgets.Button(name="Open Output", button_type="default", width=140, height=40)
    spinner  = pn.indicators.LoadingSpinner(value=False, width=32, height=32, color="primary")
    
    status_md = pn.pane.HTML(
        '<div style="color:var(--muted); font-size:13px; margin:4px 0 0">Ready — set the controls folder and click Run QC.</div>',
        sizing_mode="stretch_width"
    )

    viewer = pn.pane.HTML(
        '<div style="display:flex; align-items:center; justify-content:center; height:500px; color:var(--muted); font-size:14px; font-weight:500; background:#ffffff">QC report will appear here after running.</div>',
        sizing_mode="stretch_both",
        min_height=600,
        styles={"border": "1px solid var(--border)", "border-radius": "8px", "overflow": "hidden"}
    )

    def on_run_clicked(event):
        APP_SETTINGS["qc"]["input_dir"] = input_dir.value
        APP_SETTINGS["qc"]["output_base"] = output_base.value
        APP_SETTINGS["qc"]["outfile_html"] = outfile_html.value
        APP_SETTINGS["qc"]["excel_name"] = excel_name.value
        APP_SETTINGS["qc"]["min_r2_ok"] = min_r2_ok.value
        APP_SETTINGS["qc"]["min_r2_warn"] = min_r2_warn.value
        APP_SETTINGS["qc"]["max_mse_ok"] = max_mse_ok.value
        APP_SETTINGS["qc"]["max_mse_warn"] = max_mse_warn.value
        APP_SETTINGS["qc"]["nk_ymax_floor"] = nk_ymax.value
        APP_SETTINGS["qc"]["w_sample"] = w_sample.value
        APP_SETTINGS["qc"]["w_ladder"] = w_ladder.value
        save_settings(APP_SETTINGS)

        from pathlib import Path
        in_path = Path(input_dir.value).expanduser() if input_dir.value else None
        if not in_path or not in_path.exists():
            status_md.object = '<div style="color:var(--red); font-size:13px">Invalid input directory.</div>'
            return

        out_path = Path(output_base.value).expanduser() if output_base.value else in_path

        from core.qc.qc_rules import QCRules
        rules = QCRules(
            min_r2_ok=min_r2_ok.value,
            min_r2_warn=min_r2_warn.value,
            max_mse_ok=max_mse_ok.value,
            max_mse_warn=max_mse_warn.value,
            nk_ymax_floor=nk_ymax.value,
            sample_peak_window_bp=w_sample.value,
            ladder_peak_window_bp=w_ladder.value,
        )

        status_md.object = '<div style="color:var(--amber); font-size:13px">Running QC analysis...</div>'
        spinner.value = True
        run_btn.disabled = True
        viewer.object = '<div style="display:flex;align-items:center;justify-content:center;height:300px;color:var(--muted)">Calculating...</div>'

        def job_wrapper():
            res = None
            try:
                res = run_qc_job(
                    fsa_dir=in_path,
                    base_outdir=out_path,
                    out_html_name=outfile_html.value,
                    excel_name=excel_name.value,
                    rules=rules
                )
            finally:
                pn.state.execute(lambda: _on_done(res))

        def _on_done(res_path):
            spinner.value = False
            run_btn.disabled = False
            if res_path and res_path.exists():
                status_md.object = f'<div style="color:var(--green); font-size:13px">QC complete. Report generated at {res_path}</div>'
                try:
                    content = res_path.read_text(encoding="utf-8").replace("'", "&#39;")
                    viewer.object = f"<iframe srcdoc='{content}' style='width:100%;height:700px;border:none;border-radius:8px;background:#fff'></iframe>"
                except Exception as e:
                    viewer.object = f'<div style="color:var(--red); padding:20px">Could not load preview: {e}</div>'
            else:
                status_md.object = '<div style="color:var(--red); font-size:13px">QC run failed. Check Log tab for details.</div>'
                viewer.object = '<div style="color:var(--red); padding:20px">Failed to generate report.</div>'

        started = executor.run_background(job_wrapper)
        if not started:
            status_md.object = '<div style="color:var(--amber); font-size:13px">A job is already running.</div>'
            spinner.value = False
            run_btn.disabled = False

    run_btn.on_click(on_run_clicked)

    def on_open_output(event):
        import subprocess, sys
        from pathlib import Path
        p = Path(output_base.value).expanduser() if output_base.value else Path(input_dir.value).expanduser()
        if p and p.exists():
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(p)])
            elif sys.platform == "win32":
                subprocess.Popen(["explorer", str(p)])
            else:
                subprocess.Popen(["xdg-open", str(p)])

    open_btn.on_click(on_open_output)

    return pn.Column(
        section_header("QC Analysis", "Run quality control analysis on PK/NK/RK control files."),
        
        pn.pane.HTML('<div class="fd-card-title">Folders</div>'),
        pn.Column(
            pn.Row(input_dir, browse_in_btn, smart_scan_btn, sizing_mode="stretch_width"),
            pn.Row(output_base, browse_out_btn, sizing_mode="stretch_width"),
            pn.Row(outfile_html, excel_name, styles={"gap": "12px", "flex-wrap": "wrap"}),
            sizing_mode="stretch_width",
            css_classes=["fd-card"]
        ),
        
        make_card(
            "QC Parameters",
            pn.Row(min_r2_ok, min_r2_warn, max_mse_ok, max_mse_warn, nk_ymax, styles={"gap": "10px", "flex-wrap": "wrap"}),
            VSpace(4),
            pn.Row(w_sample, w_ladder, styles={"gap": "10px"}),
            collapsed=True,
            css_classes=["fd-card"]
        ),
        
        pn.Row(run_btn, open_btn, spinner, styles={"gap": "10px", "align-items": "center", "margin": "10px 0 4px"}),
        status_md,
        VSpace(12),

        pn.pane.HTML('<div style="font-size:11px; font-weight:600; color:#94a3b8; text-transform:uppercase; letter-spacing:0.8px; margin-bottom:8px">Live Report Preview</div>'),
        viewer,

        sizing_mode="stretch_both",
        styles={"padding": "20px 24px", "gap": "0", "max-width": "1300px"},
    )
