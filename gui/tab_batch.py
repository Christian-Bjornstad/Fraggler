"""
Fraggler Diagnostics — Batch Tab

Folder selection always visible. Job table and progress below. Settings collapsed.
"""
from __future__ import annotations

import panel as pn
import pandas as pd

from core.runner import executor
from core.batch import generate_jobs, run_batch_jobs
from core.log import log
from config import APP_SETTINGS, save_settings
from gui.components import badge, HSpace, VSpace, section_header


def make_batch_tab() -> pn.Column:
    s = APP_SETTINGS.get("batch", {})

    # ── Folder selection (always visible at top) ──────────────────────
    base_input_dir = pn.widgets.TextInput(
        name="Patient Folder (scan for subfolders with .fsa files)",
        value=s.get("base_input_dir", ""),
        placeholder="/path/to/patients",
        sizing_mode="stretch_width",
    )
    browse_in_btn = pn.widgets.Button(name="Browse...", width=90, height=32, align="end", margin=(0, 0, 4, 0))
    
    output_base = pn.widgets.TextInput(
        name="Output Folder (where results are saved)",
        value=s.get("output_base", ""),
        placeholder="/path/to/output  (leave empty = same as input)",
        sizing_mode="stretch_width",
    )
    browse_out_btn = pn.widgets.Button(name="Browse...", width=90, height=32, align="end", margin=(0, 0, 4, 0))

    def _asksystem_dir(target_widget):
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder = filedialog.askdirectory(initialdir=target_widget.value or "~")
        root.destroy()
        if folder:
            target_widget.value = folder

    browse_in_btn.on_click(lambda e: _asksystem_dir(base_input_dir))
    browse_out_btn.on_click(lambda e: _asksystem_dir(output_base))


    # ── Job options ───────────────────────────────────────────────────
    job_type = pn.widgets.Select(
        name="Job Type",
        options={"Pipeline — Run assays + DIT": "pipeline",
                 "QC — Controls only": "qc",
                 "DIT — HTML reports only": "dit"},
        value=s.get("job_type", "pipeline"),
        width=280,
    )
    aggregate_dit_reports = pn.widgets.Checkbox(
        name="Aggregate all patients into one DIT report",
        value=bool(s.get("aggregate_dit_reports", True)),
    )
    aggregate_by_patient = pn.widgets.Checkbox(
        name="Group files by Patient ID (regex)",
        value=bool(s.get("aggregate_by_patient", True)),
    )
    patient_id_regex = pn.widgets.TextInput(
        name="Patient ID Regex",
        value=s.get("patient_id_regex", r"\d{2}OUM\d{5}"),
        width=220,
        disabled=not s.get("aggregate_by_patient", True),
    )

    # ── Advanced (collapsed) ──────────────────────────────────────────
    out_folder_tmpl  = pn.widgets.TextInput(name="Folder template",   value=s.get("out_folder_tmpl", "ASSAY_REPORTS"),      width=220)
    mode_select      = pn.widgets.Select(name="Scope",                options=["all", "controls", "custom"],                 value=s.get("mode", "all"), width=140)
    assay_filter     = pn.widgets.TextInput(name="Assay filter",      value=s.get("assay_filter_substring", ""),             placeholder="e.g. TCRgA", width=180)
    outfile_html_tmpl= pn.widgets.TextInput(name="QC HTML template",  value=s.get("outfile_html_tmpl", "QC_REPORT_{name}.html"), width=240)
    excel_name_tmpl  = pn.widgets.TextInput(name="Excel template",    value=s.get("excel_name_tmpl", "QC_TRENDS_{name}.xlsx"),  width=240)
    min_r2_ok        = pn.widgets.FloatInput(name="min R² (OK)",      value=0.995, step=0.001, width=110)
    min_r2_warn      = pn.widgets.FloatInput(name="min R² (WARN)",    value=0.990, step=0.001, width=110)
    max_mse_ok       = pn.widgets.FloatInput(name="max MSE (OK)",     value=2.0,   step=0.1,   width=110)
    max_mse_warn     = pn.widgets.FloatInput(name="max MSE (WARN)",   value=5.0,   step=0.1,   width=110)
    nk_ymax          = pn.widgets.FloatInput(name="NK y-max",         value=250.0, step=10.0,  width=110)
    continue_on_error= pn.widgets.Checkbox(name="Continue batch on error", value=True)

    # ── Action buttons ────────────────────────────────────────────────
    scan_btn = pn.widgets.Button(name="Scan Jobs",      button_type="default",  width=140, height=40)
    run_btn  = pn.widgets.Button(name="Run Batch",      button_type="primary",  width=140, height=40, disabled=True)
    open_btn = pn.widgets.Button(name="Open Output",    button_type="default",  width=140, height=40)

    # ── Indicators ───────────────────────────────────────────────────
    spinner  = pn.indicators.LoadingSpinner(value=False, width=32, height=32, color="success")
    progress = pn.indicators.Progress(value=0, max=100, active=False,
                                      sizing_mode="stretch_width", bar_color="success")
    status_md = pn.pane.HTML(
        '<div style="color:#5a6a8a; font-size:13px; margin:4px 0 0">Ready — set folder and click Scan Jobs.</div>',
        sizing_mode="stretch_width",
    )

    # ── Stat counters (HTML panes so we can update them) ──────────────
    def _sc(val: str, label: str, color: str = "var(--text)") -> pn.pane.HTML:
        return pn.pane.HTML(
            f'<div class="stat-card"><div class="v" style="color:{color}">{val}</div>'
            f'<div class="l">{label}</div></div>',
            min_width=90,
        )

    stat_total   = _sc("—", "Total")
    stat_done    = _sc("—", "Done",    "var(--green)")
    stat_errors  = _sc("—", "Errors",  "var(--red)")
    stat_pending = _sc("—", "Pending", "var(--muted)")

    # ── Jobs table ────────────────────────────────────────────────────
    _EMPTY_DF = pd.DataFrame(columns=["Patient ID", "Source", "Files", "Status"])
    jobs_table = pn.widgets.Tabulator(
        _EMPTY_DF,
        sizing_mode="stretch_width",
        height=320,
        disabled=True,
        theme="bootstrap5",
        show_index=False,
        formatters={"Status": {"type": "html"}},
        styles={"border-radius": "10px", "overflow": "hidden"},
        selectable="checkbox",
    )

    select_all_btn = pn.widgets.Button(name="Select All", width=100, button_type="default")
    select_none_btn = pn.widgets.Button(name="Select None", width=100, button_type="default")

    def _select_all(e):
        if not jobs_table.value.empty:
            jobs_table.selection = list(range(len(jobs_table.value)))
    
    def _select_none(e):
        jobs_table.selection = []

    select_all_btn.on_click(_select_all)
    select_none_btn.on_click(_select_none)

    _detected_jobs: list[dict] = []
    _job_states: dict[str, str] = {}

    # ── Helpers ──────────────────────────────────────────────────────
    def _rebuild_table():
        rows = []
        for j in _detected_jobs:
            state = _job_states.get(j["name"], "pending")
            cls_map = {"done": "bd", "error": "be", "running": "br", "pending": "bp", "success": "bd"}
            cls = cls_map.get(state, "bp")
            badge_html = f'<span class="badge {cls}">{state.upper()}</span>'
            rows.append({
                "Patient ID": j["name"],
                "Source": str(j["path"]) if j.get("path") else "[Aggregated from multiple]",
                "Files": len(j["files"]) if j.get("files") else "auto",
                "Status": badge_html,
            })
        jobs_table.value = pd.DataFrame(rows) if rows else _EMPTY_DF

    def _update_stats():
        total = len(_detected_jobs)
        done  = sum(1 for s in _job_states.values() if s in ("done", "success"))
        errs  = sum(1 for s in _job_states.values() if s == "error")
        pend  = total - done - errs - sum(1 for s in _job_states.values() if s == "running")

        def _sc_html(val, label, color):
            return (f'<div class="stat-card"><div class="v" style="color:{color}">{val}</div>'
                    f'<div class="l">{label}</div></div>')

        stat_total.object   = _sc_html(total, "Total",   "var(--text)")
        stat_done.object    = _sc_html(done,  "Done",    "var(--green)")
        stat_errors.object  = _sc_html(errs,  "Errors",  "var(--red)")
        stat_pending.object = _sc_html(pend,  "Pending", "var(--muted)")

    # ── Watcher: patient regex toggle ────────────────────────────────
    aggregate_by_patient.param.watch(lambda ev: setattr(patient_id_regex, 'disabled', not ev.new), "value")
    job_type.param.watch(lambda ev: setattr(aggregate_dit_reports, 'disabled', ev.new != "pipeline"), "value")
    aggregate_dit_reports.disabled = (job_type.value != "pipeline")

    # ── Scan ─────────────────────────────────────────────────────────
    def on_scan(event):
        from pathlib import Path
        nonlocal _detected_jobs, _job_states
        spinner.value = True
        run_btn.disabled = True
        status_md.object = '<div style="color:#f59e0b;font-size:13px">Scanning for jobs...</div>'

        # Use base folder
        path_str = base_input_dir.value.strip()
        use_yaml = False
        path_obj = Path(path_str).expanduser() if path_str else None

        if not path_obj or not path_obj.exists():
            status_md.object = '<div style="color:#ef4444;font-size:13px">Invalid path — check Patient Folder or YAML path.</div>'
            spinner.value = False
            return

        try:
            _detected_jobs = generate_jobs(
                mode="yaml" if use_yaml else "subfolders",
                base_dir_or_yaml=path_obj,
                yaml_mode=use_yaml,
                aggregate_patients=aggregate_by_patient.value,
                patient_regex=patient_id_regex.value if aggregate_by_patient.value else "",
            )
            _job_states = {j["name"]: "pending" for j in _detected_jobs}

            if not _detected_jobs:
                status_md.object = '<div style="color:#f59e0b;font-size:13px">No jobs found — check that subfolders contain .fsa files.</div>'
            else:
                n = len(_detected_jobs)
                status_md.object = f'<div style="color:#22c55e;font-size:13px">Found <strong>{n}</strong> job{"s" if n!=1 else ""} — ready to run.</div>'
                run_btn.disabled = False
                progress.value = 0
                progress.max = n

            _rebuild_table()
            _update_stats()
        except Exception as e:
            status_md.object = f'<div style="color:#ef4444;font-size:13px">Scan error: {e}</div>'
        finally:
            spinner.value = False

    scan_btn.on_click(on_scan)

    # ── Run ──────────────────────────────────────────────────────────
    def on_run(event):
        if not _detected_jobs:
            return
        
        # Get selected jobs
        selected_indices = jobs_table.selection
        if not selected_indices:
            status_md.object = '<div style="color:#ef4444;font-size:13px">No jobs selected — check boxes in the list below.</div>'
            return
        
        actual_jobs_to_run = [_detected_jobs[i] for i in selected_indices]
        total_to_run = len(actual_jobs_to_run)

        from pathlib import Path

        # Persist settings
        _s = APP_SETTINGS.setdefault("batch", {})
        _s.update({
            "base_input_dir": base_input_dir.value,

            "job_type": job_type.value,
            "aggregate_dit_reports": aggregate_dit_reports.value,
            "aggregate_by_patient": aggregate_by_patient.value,
            "patient_id_regex": patient_id_regex.value,
            "output_base": output_base.value,
            "out_folder_tmpl": out_folder_tmpl.value,
            "outfile_html_tmpl": outfile_html_tmpl.value,
            "excel_name_tmpl": excel_name_tmpl.value,
            "mode": mode_select.value,
            "assay_filter_substring": assay_filter.value,
        })
        save_settings(APP_SETTINGS)

        # Resolve output path (default = same as input)
        out_path_str = output_base.value.strip() or base_input_dir.value.strip()
        out_path = Path(out_path_str).expanduser() if out_path_str else None
        if not out_path or not out_path.exists():
            status_md.object = '<div style="color:#ef4444;font-size:13px">Output folder does not exist — set it before running.</div>'
            return

        # Reset state for selected jobs
        for j in actual_jobs_to_run:
            _job_states[j["name"]] = "pending"
        _rebuild_table()
        _update_stats()

        progress.value = 0
        progress.max = total_to_run
        progress.active = True
        spinner.value = True
        run_btn.disabled = True
        scan_btn.disabled = True
        status_md.object = f'<div style="color:#f59e0b;font-size:13px">Running {total_to_run} selected jobs — see Log tab for details.</div>'

        from core.qc.qc_rules import QCRules
        rules = QCRules(
            min_r2_ok=min_r2_ok.value,
            min_r2_warn=min_r2_warn.value,
            max_mse_ok=max_mse_ok.value,
            max_mse_warn=max_mse_warn.value,
            nk_ymax_floor=nk_ymax.value,
        )

        def _update_progress(idx, total, name, state):
            _job_states[name] = state
            _rebuild_table()
            _update_stats()
            progress.value = idx
            if idx >= total:
                n_done = sum(1 for s in _job_states.values() if s in ("done", "success"))
                # Filter n_done to only those in the current run if needed, but the UI keeps track of all
                progress.active = False
                spinner.value = False
                run_btn.disabled = False
                scan_btn.disabled = False
                status_md.object = f'<div style="color:#22c55e;font-size:13px">Batch complete: {total} jobs processed.</div>'
            else:
                status_md.object = f'<div style="color:#f59e0b;font-size:13px">[{idx}/{total}] Running: <strong>{name}</strong></div>'

        def update_ui(idx, total, name, state):
            pn.state.execute(lambda: _update_progress(idx, total, name, state))

        def job_wrapper():
            try:
                run_batch_jobs(
                    jobs=actual_jobs_to_run,
                    job_type=job_type.value,
                    output_base=out_path,
                    out_folder_tmpl=out_folder_tmpl.value,
                    outfile_html_tmpl=outfile_html_tmpl.value,
                    excel_name_tmpl=excel_name_tmpl.value,
                    pipeline_scope=mode_select.value,
                    assay_filter=assay_filter.value,
                    qc_rules=rules,
                    continue_on_error=continue_on_error.value,
                    aggregate_dit_reports=aggregate_dit_reports.value,
                    update_callback=update_ui,
                )
            except Exception as e:
                pn.state.execute(lambda: _fail(e))

        def _fail(err):
            spinner.value = False
            progress.active = False
            run_btn.disabled = False
            scan_btn.disabled = False
            status_md.object = f'<div style="color:#ef4444;font-size:13px">Batch error: {err}</div>'

        if not executor.run_background(job_wrapper):
            status_md.object = '<div style="color:#f59e0b;font-size:13px">A job is already running — please wait.</div>'
            spinner.value = False
            run_btn.disabled = False
            scan_btn.disabled = False

    run_btn.on_click(on_run)

    def on_open_output(event):
        import subprocess, sys
        from pathlib import Path
        p = Path(output_base.value or base_input_dir.value).expanduser()
        if p.exists():
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(p)])
            elif sys.platform == "win32":
                subprocess.Popen(["explorer", str(p)])
            else:
                subprocess.Popen(["xdg-open", str(p)])

    open_btn.on_click(on_open_output)

    # ── Layout ────────────────────────────────────────────────────────
    return pn.Column(
        # Page heading
        section_header("Batch Processing", "Scan for patient folders and run the Fraggler pipeline."),

        # ── FOLDER CARD (always visible) ──────────────────────────────
        pn.pane.HTML('<div class="fd-card-title">Folders</div>'),
        pn.Column(
            pn.Row(base_input_dir, browse_in_btn, sizing_mode="stretch_width"),
            pn.Row(output_base, browse_out_btn, sizing_mode="stretch_width"),

            sizing_mode="stretch_width",
            css_classes=["fd-card"],
        ),

        # ── JOB OPTIONS (visible) ─────────────────────────────────────
        pn.pane.HTML('<div class="fd-card-title" style="margin-top:6px">Job Options</div>'),
        pn.Column(
            pn.Row(job_type, styles={"margin-bottom": "4px"}),
            pn.Row(aggregate_dit_reports, aggregate_by_patient, styles={"gap": "24px", "align-items": "center"}),
            pn.Row(pn.pane.HTML('<span style="font-size:11px;color:var(--muted)">Patient ID Regex:</span>'), patient_id_regex, styles={"gap": "10px", "align-items": "center"}),
            sizing_mode="stretch_width",
            css_classes=["fd-card"],
        ),

        # ── ACTIONS ───────────────────────────────────────────────────
        pn.Row(
            scan_btn, run_btn, open_btn, spinner,
            styles={"gap": "10px", "align-items": "center", "margin": "10px 0 4px"},
        ),
        progress,
        status_md,
        pn.Row(
            pn.pane.HTML('<div style="margin: 10px 0 6px; font-size:10px; font-weight:800; text-transform:uppercase; letter-spacing:0.8px; color:var(--muted)">Detected Jobs</div>'),
            HSpace(),
            select_all_btn, select_none_btn,
            styles={"align-items": "center", "margin-bottom": "4px"}
        ),

        # ── STAT ROW ─────────────────────────────────────────────────
        pn.Row(
            stat_total, stat_done, stat_errors, stat_pending,
            styles={"gap": "10px", "flex-wrap": "wrap", "margin-bottom": "8px"},
        ),

        # ── JOB TABLE ────────────────────────────────────────────────
        jobs_table,

        # ── ADVANCED (collapsed) ─────────────────────────────────────
        pn.pane.HTML("""
<details class="collapse-card" style="margin-top:14px">
  <summary>Advanced Settings</summary>
  <div class="collapse-body"></div>
</details>"""),
        pn.Column(
            pn.Row(out_folder_tmpl, mode_select, assay_filter, styles={"gap": "10px", "flex-wrap": "wrap"}),
            pn.Row(outfile_html_tmpl, excel_name_tmpl, styles={"gap": "10px", "flex-wrap": "wrap"}),
            pn.Row(min_r2_ok, min_r2_warn, max_mse_ok, max_mse_warn, nk_ymax, styles={"gap": "10px", "flex-wrap": "wrap"}),
            continue_on_error,
            sizing_mode="stretch_width",
            styles={"margin-top": "4px"},
            visible=False,
            name="_advanced",
        ),

        sizing_mode="stretch_width",
        styles={"padding": "24px 28px", "gap": "0", "background": "transparent"},
    )
