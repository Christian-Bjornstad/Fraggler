"""
HemaFrag Diagnostics — FLT3 Validation Tab
"""
from __future__ import annotations

import contextlib
import io
import subprocess
import sys
from pathlib import Path
from typing import Any

import panel as pn

from config import APP_SETTINGS, save_settings
from core.log import log, log_buffer
from core.runner import executor
from gui.components import VSpace, make_card, section_header


class _LogBufferStream(io.TextIOBase):
    def __init__(self, prefix: str = "[FLT3]") -> None:
        super().__init__()
        self._prefix = prefix
        self._buffer = ""

    def write(self, s: str) -> int:
        if not s:
            return 0
        self._buffer += s
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            if line.strip():
                log_buffer.write(f"{self._prefix} {line}")
        return len(s)

    def flush(self) -> None:
        if self._buffer.strip():
            log_buffer.write(f"{self._prefix} {self._buffer.strip()}")
        self._buffer = ""


def make_flt3_validation_tab() -> pn.Column:
    s = APP_SETTINGS.setdefault("analyses", {}).setdefault("flt3", {}).setdefault("validation", {})

    data_root = pn.widgets.TextInput(
        name="FLT3 Data Root",
        value=s.get("data_root", "/Volumes/T7 Shield/DATA/flt3"),
        sizing_mode="stretch_width",
        placeholder="/Volumes/T7 Shield/DATA/flt3",
    )
    output_root = pn.widgets.TextInput(
        name="Output Root",
        value=s.get("output_root", str(Path.home())),
        sizing_mode="stretch_width",
        placeholder="/path/to/output-root",
    )
    run_name = pn.widgets.TextInput(
        name="Run Name",
        value=s.get("run_name", ""),
        placeholder="Leave blank for auto timestamp",
        sizing_mode="stretch_width",
    )
    years = pn.widgets.MultiChoice(
        name="Years",
        options=["2024", "2025", "2026"],
        value=list(s.get("years", ["2025", "2026"])),
        sizing_mode="stretch_width",
    )
    require_run_name_contains = pn.widgets.TextInput(
        name="Run Name Must Contain",
        value=s.get("require_run_name_contains", "3730DNA"),
        placeholder="3730DNA",
        width=200,
    )
    workers = pn.widgets.IntInput(name="Workers", value=int(s.get("workers", 6)), step=1, start=1, width=120)
    limit = pn.widgets.IntInput(name="Limit", value=int(s.get("limit", 0)), step=100, start=0, width=120)
    timeout_seconds = pn.widgets.IntInput(
        name="Timeout / file (s)",
        value=int(s.get("timeout_seconds", 45)),
        step=5,
        start=1,
        width=140,
    )
    checkpoint_every = pn.widgets.IntInput(
        name="Checkpoint Every",
        value=int(s.get("checkpoint_every", 100)),
        step=25,
        start=0,
        width=140,
    )
    include_npm1 = pn.widgets.Checkbox(name="Include NPM1", value=bool(s.get("include_npm1", False)))
    dit_only = pn.widgets.Checkbox(name="DIT only", value=bool(s.get("dit_only", False)))

    browse_data_btn = pn.widgets.Button(name="Browse...", width=90, height=32, align="end", margin=(0, 0, 4, 0))
    browse_output_btn = pn.widgets.Button(name="Browse...", width=90, height=32, align="end", margin=(0, 0, 4, 0))
    run_btn = pn.widgets.Button(name="Run FLT3 Backfill Validation", button_type="primary", width=240, height=42)
    open_output_btn = pn.widgets.Button(name="Open Output", button_type="default", width=140, height=42)
    spinner = pn.indicators.LoadingSpinner(value=False, width=32, height=32, color="danger")

    status_md = pn.pane.HTML(
        '<div style="color:var(--muted); font-size:13px">Ready — configure the 3730DNA FLT3 backfill run and launch when you want.</div>',
        sizing_mode="stretch_width",
    )
    summary_md = pn.pane.Markdown(
        "Run summary will appear here after the first validation run.",
        sizing_mode="stretch_width",
    )
    command_preview = pn.pane.Markdown("", sizing_mode="stretch_width")

    def _ask_dir(target_widget: pn.widgets.TextInput) -> None:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        folder = filedialog.askdirectory(initialdir=target_widget.value or "~")
        root.destroy()
        if folder:
            target_widget.value = folder

    browse_data_btn.on_click(lambda e: _ask_dir(data_root))
    browse_output_btn.on_click(lambda e: _ask_dir(output_root))

    def _current_run_dir() -> Path | None:
        last_run_dir = str(APP_SETTINGS.setdefault("analyses", {}).setdefault("flt3", {}).setdefault("validation", {}).get("last_run_dir") or "").strip()
        if last_run_dir:
            p = Path(last_run_dir).expanduser()
            if p.exists():
                return p
        root = Path(output_root.value).expanduser() if output_root.value.strip() else None
        return root if root and root.exists() else None

    def _save_settings() -> None:
        target = APP_SETTINGS.setdefault("analyses", {}).setdefault("flt3", {}).setdefault("validation", {})
        target["data_root"] = data_root.value
        target["output_root"] = output_root.value
        target["run_name"] = run_name.value
        target["years"] = list(years.value)
        target["require_run_name_contains"] = require_run_name_contains.value
        target["workers"] = int(workers.value or 1)
        target["limit"] = int(limit.value or 0)
        target["timeout_seconds"] = int(timeout_seconds.value or 45)
        target["checkpoint_every"] = int(checkpoint_every.value or 100)
        target["include_npm1"] = bool(include_npm1.value)
        target["dit_only"] = bool(dit_only.value)
        save_settings(APP_SETTINGS)

    def _command_lines() -> list[str]:
        out_root = output_root.value.strip() or "<output-root>"
        lines = [
            "python3 scripts/run_flt3_backfill_validation.py",
            f"  --data-root '{data_root.value.strip() or '/Volumes/T7 Shield/DATA/flt3'}'",
            f"  --output-root '{out_root}'",
            f"  --require-run-name-contains '{require_run_name_contains.value.strip() or '3730DNA'}'",
            f"  --workers {int(workers.value or 1)}",
            f"  --timeout-seconds {int(timeout_seconds.value or 45)}",
            f"  --checkpoint-every {int(checkpoint_every.value or 100)}",
        ]
        if run_name.value.strip():
            lines.append(f"  --run-name '{run_name.value.strip()}'")
        for year in years.value:
            lines.append(f"  --year {year}")
        if int(limit.value or 0) > 0:
            lines.append(f"  --limit {int(limit.value)}")
        if include_npm1.value:
            lines.append("  --include-npm1")
        if dit_only.value:
            lines.append("  --dit-only")
        return lines

    def _refresh_command_preview(*_events: Any) -> None:
        _save_settings()
        command_preview.object = "```bash\n" + " \\\n".join(_command_lines()) + "\n```"

    for widget in (
        data_root,
        output_root,
        run_name,
        years,
        require_run_name_contains,
        workers,
        limit,
        timeout_seconds,
        checkpoint_every,
        include_npm1,
        dit_only,
    ):
        widget.param.watch(_refresh_command_preview, "value")

    def _load_summary(run_dir: Path) -> str:
        residual_md = run_dir / "residual_summary.md"
        summary_md_path = run_dir / "summary.md"
        chunks: list[str] = []
        if summary_md_path.exists():
            chunks.append(summary_md_path.read_text(encoding="utf-8", errors="replace"))
        if residual_md.exists():
            chunks.append(residual_md.read_text(encoding="utf-8", errors="replace"))
        return "\n\n---\n\n".join(chunk.strip() for chunk in chunks if chunk.strip()) or "Run completed, but no markdown summary was found."

    def _on_done(payload: dict[str, Any] | None, error_message: str | None = None) -> None:
        spinner.value = False
        run_btn.disabled = False
        if error_message:
            status_md.object = f'<div style="color:var(--red); font-size:13px">FLT3 validation failed: {error_message}</div>'
            return
        if not payload:
            status_md.object = '<div style="color:var(--red); font-size:13px">FLT3 validation finished without a payload.</div>'
            return

        run_dir = Path(str(payload.get("run_dir") or "")).expanduser()
        target = APP_SETTINGS.setdefault("analyses", {}).setdefault("flt3", {}).setdefault("validation", {})
        target["last_run_dir"] = str(run_dir)
        save_settings(APP_SETTINGS)

        status_counts = ((payload.get("validator_summary") or {}).get("status_counts") or {})
        status_md.object = (
            f'<div style="color:var(--green); font-size:13px">FLT3 validation complete. '
            f'Run folder: {run_dir}. Status counts: {status_counts}</div>'
        )
        summary_md.object = _load_summary(run_dir)

    def on_run(_event: Any) -> None:
        _save_settings()
        data_root_path = Path(data_root.value).expanduser()
        if not data_root_path.exists():
            status_md.object = '<div style="color:var(--red); font-size:13px">Data root does not exist.</div>'
            return

        output_root_path = Path(output_root.value).expanduser() if output_root.value.strip() else Path.home()
        output_root_path.mkdir(parents=True, exist_ok=True)

        status_md.object = '<div style="color:var(--amber); font-size:13px">Running FLT3 3730 backfill validation...</div>'
        summary_md.object = "Validation running. Summary will appear here when the job completes."
        spinner.value = True
        run_btn.disabled = True

        def job_wrapper() -> None:
            stream = _LogBufferStream("[FLT3-VALIDATE]")
            payload: dict[str, Any] | None = None
            error_message: str | None = None
            try:
                from scripts.run_flt3_backfill_validation import run_backfill_validation

                with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                    payload = run_backfill_validation(
                        data_root=data_root_path,
                        output_root=output_root_path,
                        run_name=run_name.value.strip() or None,
                        years=list(years.value),
                        workers=int(workers.value or 1),
                        limit=int(limit.value or 0),
                        input_manifest=None,
                        include_npm1=bool(include_npm1.value),
                        dit_only=bool(dit_only.value),
                        timeout_seconds=int(timeout_seconds.value or 45),
                        checkpoint_every=int(checkpoint_every.value or 100),
                        required_run_name_contains=require_run_name_contains.value.strip() or "3730DNA",
                    )
            except Exception as exc:  # pragma: no cover - defensive UI error handling
                error_message = str(exc)
                log(f"[ERROR] FLT3 validation failed: {exc}")
            finally:
                stream.flush()
                pn.state.execute(lambda payload=payload, error_message=error_message: _on_done(payload, error_message))

        started = executor.run_background(job_wrapper)
        if not started:
            spinner.value = False
            run_btn.disabled = False
            status_md.object = '<div style="color:var(--amber); font-size:13px">Another job is already running.</div>'

    def on_open_output(_event: Any) -> None:
        path = _current_run_dir()
        if path is None:
            status_md.object = '<div style="color:var(--amber); font-size:13px">No existing output directory to open yet.</div>'
            return
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(path)])
        elif sys.platform == "win32":
            subprocess.Popen(["explorer", str(path)])
        else:
            subprocess.Popen(["xdg-open", str(path)])

    run_btn.on_click(on_run)
    open_output_btn.on_click(on_open_output)

    _refresh_command_preview()

    return pn.Column(
        section_header("FLT3 Validation", "Run a 3730DNA-filtered FLT3 backfill validation with residual exports and review bundles."),
        pn.Column(
            pn.pane.HTML('<div class="fd-card-title">Data And Output</div>'),
            pn.Row(data_root, browse_data_btn, sizing_mode="stretch_width"),
            pn.Row(output_root, browse_output_btn, sizing_mode="stretch_width"),
            pn.Row(run_name, years, sizing_mode="stretch_width", styles={"gap": "12px"}),
            css_classes=["fd-card"],
            sizing_mode="stretch_width",
        ),
        make_card(
            "Validation Scope",
            pn.Row(
                require_run_name_contains,
                workers,
                limit,
                timeout_seconds,
                checkpoint_every,
                styles={"gap": "12px", "flex-wrap": "wrap"},
            ),
            pn.Row(include_npm1, dit_only, styles={"gap": "18px", "align-items": "center"}),
            pn.pane.Markdown(
                "This tab calls `run_flt3_backfill_validation.py`, which in turn runs the FLT3 ladder validator and writes residual summary files such as `residual_summary.json` and `residual_by_year.csv`."
            ),
            collapsed=False,
            css_classes=["fd-card"],
        ),
        pn.Row(run_btn, open_output_btn, spinner, styles={"gap": "10px", "align-items": "center", "margin": "10px 0 4px"}),
        status_md,
        VSpace(12),
        make_card(
            "Command Preview",
            command_preview,
            collapsed=False,
            css_classes=["fd-card"],
        ),
        make_card(
            "Last Summary",
            summary_md,
            collapsed=False,
            css_classes=["fd-card"],
        ),
        sizing_mode="stretch_both",
        styles={"padding": "20px 24px", "gap": "0", "max-width": "1300px"},
    )
