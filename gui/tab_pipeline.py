"""
Fraggler Diagnostics — Pipeline Tab
"""
from __future__ import annotations

import panel as pn

from gui.components import make_card, VSpace, section_header, hero_row
from core.runner import executor, run_pipeline_job
from core.log import log
from config import APP_SETTINGS, save_settings


def make_pipeline_tab() -> pn.Column:
    s = APP_SETTINGS["pipeline"]
    active_analysis = APP_SETTINGS.get("active_analysis", "clonality")
    is_flt3 = active_analysis == "flt3"

    input_dir = pn.widgets.TextInput(
        name="Input Folder (.fsa files)",
        value=s.get("input_dir", ""),
        sizing_mode="stretch_width",
        placeholder="/path/to/fsa/folder"
    )
    output_base = pn.widgets.TextInput(
        name="Output Base Folder (leave empty = same as input)",
        value=s.get("output_base", ""),
        sizing_mode="stretch_width",
        placeholder="Leave empty to place results beside input"
    )
    out_folder_name = pn.widgets.TextInput(
        name="Reports Subfolder Name",
        value=s.get("out_folder_name", "ASSAY_REPORTS"),
        width=240
    )
    mode_select = pn.widgets.RadioButtonGroup(
        name="Scope",
        options=["all", "controls", "custom"],
        value=s.get("mode", "all"),
        button_type="primary",
    )
    assay_filter = pn.widgets.TextInput(
        name="Custom Assay Filter (only when scope=custom)",
        value=s.get("assay_filter_substring", ""),
        placeholder="e.g. ITD, D835, NPM1" if is_flt3 else "e.g. TCRgA, FR3",
        sizing_mode="stretch_width"
    )

    run_btn = pn.widgets.Button(
        name="▶  Run Pipeline", button_type="success", width=180, height=48,
        styles={"font-size": "15px", "font-weight": "600"}
    )
    open_btn = pn.widgets.Button(name="📂 Open Output", button_type="default", width=150, height=48)
    spinner = pn.indicators.LoadingSpinner(value=False, width=36, height=36, color="success")
    status_md = pn.pane.HTML(
        '<div style="color:#94a3b8; font-size:13px">Ready. Set the input folder and click Run.</div>',
        sizing_mode="stretch_width"
    )
    analysis_note = pn.pane.HTML(
        """
        <div style="padding:14px 16px; border-radius:12px; background:linear-gradient(135deg,#f8fafc 0%,#ecfeff 100%);
                    border:1px solid #cbd5e1; color:#0f172a;">
          <div style="font-size:12px; font-weight:700; text-transform:uppercase; letter-spacing:0.08em; color:#0f766e; margin-bottom:8px;">
            FLT3 / NPM1 workflow
          </div>
          <div style="font-size:14px; line-height:1.5;">
            Preferred injections from the current worksheet are <strong>1 s</strong> for FLT3-ITD and ratio runs,
            and <strong>3 s</strong> for FLT3-D835, TKD-digested, undiluted, and NPM1 runs.
            File names should still clearly contain <code>itd</code>, <code>d835</code>/<code>tkd</code>, or <code>npm1</code>.
          </div>
        </div>
        """ if is_flt3 else
        """
        <div style="padding:14px 16px; border-radius:12px; background:#f8fafc; border:1px solid #e2e8f0; color:#334155; font-size:14px; line-height:1.5;">
          The pipeline tab runs the active analysis on one folder at a time. Use <code>custom</code> scope if you want
          to stage a narrower subset by assay name.
        </div>
        """,
        sizing_mode="stretch_width",
    )

    def on_run_clicked(event):
        APP_SETTINGS["pipeline"]["input_dir"] = input_dir.value
        APP_SETTINGS["pipeline"]["output_base"] = output_base.value
        APP_SETTINGS["pipeline"]["out_folder_name"] = out_folder_name.value
        APP_SETTINGS["pipeline"]["mode"] = mode_select.value
        APP_SETTINGS["pipeline"]["assay_filter_substring"] = assay_filter.value
        save_settings(APP_SETTINGS)

        from pathlib import Path
        in_path = Path(input_dir.value).expanduser() if input_dir.value else None
        if not in_path or not in_path.exists():
            status_md.object = '<div style="color:#ef4444; font-size:13px">❌ Invalid input directory.</div>'
            return

        out_path = Path(output_base.value).expanduser() if output_base.value else in_path
        status_md.object = '<div style="color:#f59e0b; font-size:13px">⏳ Running pipeline... check Log tab for details.</div>'
        spinner.value = True
        run_btn.disabled = True

        def job_wrapper():
            try:
                run_pipeline_job(
                    fsa_dir=in_path,
                    base_outdir=out_path,
                    out_folder_name=out_folder_name.value,
                    scope=mode_select.value,
                    needle=assay_filter.value
                )
            finally:
                pn.state.execute(lambda: _on_done())

        def _on_done():
            spinner.value = False
            run_btn.disabled = False
            status_md.object = '<div style="color:#22c55e; font-size:13px">✅ Pipeline finished. Check Log tab for details.</div>'

        started = executor.run_background(job_wrapper)
        if not started:
            status_md.object = '<div style="color:#f59e0b; font-size:13px">⚠️ A job is already running — please wait.</div>'
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
        section_header(
            "Pipeline Run",
            "Run FLT3 / NPM1 fragment analysis on a single folder of .fsa files"
            if is_flt3 else
            "Run Fraggler assay analysis on a single folder of .fsa files"
        ),
        VSpace(8),
        analysis_note,
        VSpace(10),

        pn.Row(run_btn, open_btn, spinner, styles={"gap": "12px", "align-items": "center"}),
        status_md,
        VSpace(12),

        make_card(
            "Input / Output",
            input_dir,
            output_base,
            out_folder_name,
        ),
        VSpace(8),
        make_card(
            "Scope & Filtering",
            pn.Row(
                pn.pane.HTML('<div style="font-size:12px; color:#94a3b8; font-weight:500; margin-top:6px">Scope:</div>'),
                mode_select,
                styles={"align-items": "center", "gap": "12px"}
            ),
            assay_filter,
        ),

        sizing_mode="stretch_both",
        styles={"padding": "20px 24px", "gap": "0", "max-width": "900px"},
    )
