from __future__ import annotations

import argparse
import copy
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import panel as pn
import plotly.graph_objects as go
from scipy.signal import find_peaks


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from core.analysis import (  # noqa: E402
    apply_manual_ladder_mapping,
    compute_ladder_qc_metrics,
    estimate_running_baseline,
    get_ladder_candidates,
    load_ladder_adjustment,
    save_ladder_adjustment,
)
from core.analyses.flt3.pipeline import _analyse_fsa_candidate  # noqa: E402
from gui_qt.ladder_utils import detect_fsa_for_ladder, load_adjustable_fsa  # noqa: E402


pn.extension("plotly", sizing_mode="stretch_width")
logging.basicConfig(level=logging.INFO)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a browser-based FLT3 ladder review app."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Root directory containing FLT3 .fsa files.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5011,
        help="Port for the Panel server.",
    )
    parser.add_argument(
        "--case-manifest",
        type=Path,
        default=None,
        help="Optional CSV with FLT3 files to review. Must include a path column.",
    )
    return parser


def _markdown_escape(value: Any) -> str:
    text = str(value)
    return text.replace("|", "\\|")


def _preferred_injection_for_record(record: dict[str, Any]) -> int:
    assay = str(record.get("assay") or "")
    if assay == "FLT3-D835":
        return 3
    return 1


def _fit_grade(metrics: dict[str, float | int] | None) -> tuple[str, str]:
    if not metrics:
        return "check", "Map every ladder step to preview the fit."

    r2 = float(metrics.get("r2", float("nan")))
    max_abs = float(metrics.get("max_abs_error_bp", float("inf")))
    if not np.isfinite(r2):
        return "fail", "Preview fit failed."
    if r2 >= 0.9995 and max_abs <= 0.5:
        return "pass", "Stable fit."
    if r2 >= 0.9990 and max_abs <= 1.5:
        return "check", "Usable fit, but still worth reviewing."
    return "fail", "Fit still needs work."


def _bootstrap_candidates_from_trace(trace: np.ndarray) -> pd.DataFrame:
    if trace.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    try:
        baseline = estimate_running_baseline(trace, bin_size=200, quantile=0.10)
        corrected = np.clip(np.asarray(trace, dtype=float) - baseline, a_min=0, a_max=None)
    except Exception:
        corrected = np.asarray(trace, dtype=float)

    peak_idx, props = find_peaks(
        corrected,
        distance=12,
        prominence=35.0,
        height=40.0,
    )
    if peak_idx.size == 0:
        return pd.DataFrame(columns=["time", "intensity", "source"])

    rows = pd.DataFrame(
        {
            "time": peak_idx.astype(float),
            "intensity": corrected[peak_idx].astype(float),
            "prominence": np.asarray(props.get("prominences", []), dtype=float),
            "source": "trace_bootstrap",
        }
    )
    rows = rows[
        (rows["time"].astype(float) >= 1400.0)
        & (rows["time"].astype(float) <= 4700.0)
    ].copy()
    if rows.empty:
        return pd.DataFrame(columns=["time", "intensity", "source"])
    rows = rows.sort_values(["prominence", "intensity"], ascending=[False, False]).head(80)
    rows = rows.sort_values("time").reset_index(drop=True)
    return rows.loc[:, ["time", "intensity", "source"]]


class Flt3LadderReviewApp:
    def __init__(self, data_dir: Path | str, case_manifest: Path | str | None = None):
        self.data_dir = Path(data_dir).expanduser().resolve()
        self.case_manifest_path = Path(case_manifest).expanduser().resolve() if case_manifest else None
        self.case_manifest = self._load_case_manifest(self.case_manifest_path)
        self.records = self._scan_records()
        self.records_by_selection_key: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for record in self.records:
            self.records_by_selection_key[str(record.get("selection_key") or "")].append(record)

        self.filtered_records: list[dict[str, Any]] = []
        self.qc_cache: dict[str, dict[str, Any]] = {}
        self.current_record: dict[str, Any] | None = None
        self.current_fsa = None
        self.preview_fsa = None
        self.current_meta: dict[str, Any] | None = None
        self.current_candidates = pd.DataFrame(columns=["time", "intensity", "source"])
        self.ladder_steps = np.asarray([], dtype=float)
        self.mapping_times: dict[int, float] = {}
        self.manual_candidate_times: list[float] = []
        self.preview_metrics: dict[str, float | int] | None = None
        self.preview_reason = "No file selected."
        self.status_text = "Ready."
        self.last_clicked_time: float | None = None
        self.pending_record_path: str | None = None

        self.title = pn.pane.Markdown("## FLT3 Ladder Review")
        self.status = pn.pane.Markdown("")
        self.summary = pn.pane.Markdown("")
        self.instructions = pn.pane.Markdown(
            "Click the trace to inspect a position, then choose whether to keep the exact time or snap to a nearby peak. You can add manual candidates even when the pipeline did not propose them."
        )
        self.selection_help = pn.pane.Markdown("")
        self.mapping_view = pn.pane.Markdown("")
        self.plot = pn.pane.Plotly(height=580, config={"responsive": True, "scrollZoom": True})
        self.plot.param.watch(self._on_plot_click, "click_data")
        self.step_buttons: dict[int, pn.widgets.Button] = {}
        self.step_button_box = pn.FlexBox(sizing_mode="stretch_width")

        self.assay_filter = pn.widgets.RadioButtonGroup(
            name="Assay",
            options=["All", "FLT3-ITD", "FLT3-D835"],
            button_type="primary",
            value="All",
        )
        self.injection_filter = pn.widgets.RadioButtonGroup(
            name="Injection",
            options=["All", "1s", "3s"],
            button_type="default",
            value="All",
        )
        self.preferred_only = pn.widgets.Checkbox(
            name="Preferred injection only",
            value=False if self.case_manifest else True,
        )
        self.descending_steps = pn.widgets.Checkbox(name="Review high to low", value=True)
        self.fit_source = pn.widgets.RadioButtonGroup(
            name="Fit Source",
            options=["Pipeline rescue", "Raw adjustable"],
            button_type="default",
            value="Raw adjustable" if self.case_manifest else "Pipeline rescue",
        )
        self.search_input = pn.widgets.TextInput(name="Search", placeholder="Specimen, file, run, well")
        self.file_select = pn.widgets.Select(name="FLT3 File", options={})
        self.sort_mode = pn.widgets.RadioButtonGroup(
            name="Sort",
            options=["Default", "Worst residual", "Lowest r2"],
            button_type="default",
            value="Default",
        )
        self.rank_btn = pn.widgets.Button(name="Compute QC For Filter", button_type="primary")

        self.prev_btn = pn.widgets.Button(name="Previous", button_type="default")
        self.next_btn = pn.widgets.Button(name="Next", button_type="default")
        self.load_btn = pn.widgets.Button(name="Load Selected File", button_type="primary")
        self.reload_btn = pn.widgets.Button(name="Reload File", button_type="default")
        self.rescan_btn = pn.widgets.Button(name="Rescan Folder", button_type="default")

        self.step_select = pn.widgets.Select(name="Ladder Step", options={})
        self.candidate_select = pn.widgets.Select(name="Candidate Peak", options={})
        self.assign_btn = pn.widgets.Button(name="Assign Selected Peak", button_type="primary")
        self.manual_mode = pn.widgets.RadioButtonGroup(
            name="Manual Mode",
            options=["Snap to local peak", "Use exact time"],
            button_type="default",
            value="Snap to local peak",
        )
        self.manual_time_input = pn.widgets.FloatInput(name="Manual Peak Time", step=1.0, value=None)
        self.add_btn = pn.widgets.Button(name="Add Manual Candidate", button_type="default")
        self.snap_btn = pn.widgets.Button(name="Assign Manual Time", button_type="primary")
        self.clear_btn = pn.widgets.Button(name="Clear Step", button_type="warning")
        self.reset_btn = pn.widgets.Button(name="Reset Edits", button_type="default")
        self.save_btn = pn.widgets.Button(name="Save Adjustment", button_type="success")
        self.remove_btn = pn.widgets.Button(name="Delete Saved Adjustment", button_type="danger")
        self.step_table = pn.widgets.Tabulator(
            pd.DataFrame(columns=["step_idx", "step_bp", "assigned_time", "source", "residual", "status"]),
            selectable=1,
            pagination=None,
            disabled=True,
            height=300,
            show_index=False,
            sizing_mode="stretch_width",
        )
        self.candidate_table = pn.widgets.Tabulator(
            pd.DataFrame(columns=["candidate_idx", "time", "intensity", "source", "used_by"]),
            selectable=1,
            pagination=None,
            disabled=True,
            height=320,
            show_index=False,
            sizing_mode="stretch_width",
        )

        self.assay_filter.param.watch(self._on_filter_change, "value")
        self.injection_filter.param.watch(self._on_filter_change, "value")
        self.preferred_only.param.watch(self._on_filter_change, "value")
        self.descending_steps.param.watch(self._on_step_order_change, "value")
        self.fit_source.param.watch(self._on_fit_source_change, "value")
        self.search_input.param.watch(self._on_filter_change, "value")
        self.file_select.param.watch(self._on_file_change, "value")
        self.sort_mode.param.watch(self._on_filter_change, "value")
        self.step_select.param.watch(self._on_step_select_change, "value")
        self.candidate_select.param.watch(self._on_candidate_select_change, "value")
        self.manual_mode.param.watch(self._on_manual_mode_change, "value")
        self.step_table.param.watch(self._on_step_table_select, "selection")
        self.candidate_table.param.watch(self._on_candidate_table_select, "selection")

        self.rank_btn.on_click(lambda _event: self._compute_qc_for_filtered_records())
        self.prev_btn.on_click(lambda _event: self._move_file(-1))
        self.next_btn.on_click(lambda _event: self._move_file(1))
        self.load_btn.on_click(lambda _event: self._load_selected_record())
        self.reload_btn.on_click(lambda _event: self._reload_current_file())
        self.rescan_btn.on_click(lambda _event: self._rescan_files())
        self.assign_btn.on_click(lambda _event: self._assign_selected_candidate())
        self.add_btn.on_click(lambda _event: self._add_manual_candidate())
        self.snap_btn.on_click(lambda _event: self._assign_manual_time_input())
        self.clear_btn.on_click(lambda _event: self._clear_selected_step())
        self.reset_btn.on_click(lambda _event: self._load_record(self.current_record))
        self.save_btn.on_click(lambda _event: self._save_adjustment())
        self.remove_btn.on_click(lambda _event: self._delete_adjustment())

        sidebar = pn.Column(
            pn.Card(
                self.assay_filter,
                self.injection_filter,
                self.preferred_only,
                self.descending_steps,
                self.fit_source,
                self.search_input,
                self.sort_mode,
                self.rank_btn,
                title="Filter And Queue",
                collapsed=False,
            ),
            pn.Card(
                self.file_select,
                pn.Row(self.prev_btn, self.next_btn),
                self.load_btn,
                pn.Row(self.reload_btn, self.rescan_btn),
                title="Current File",
                collapsed=False,
            ),
            pn.Card(
                self.step_select,
                self.step_button_box,
                title="Quick Step Pad",
                collapsed=False,
            ),
            pn.Card(
                self.candidate_select,
                self.manual_mode,
                self.manual_time_input,
                pn.Row(self.assign_btn, self.add_btn),
                self.snap_btn,
                pn.Row(self.clear_btn, self.reset_btn),
                pn.Row(self.save_btn, self.remove_btn),
                title="Manual Peak Tools",
                collapsed=False,
            ),
            width=410,
        )

        review_tables = pn.Row(
            pn.Column("### Ladder Steps", self.step_table, sizing_mode="stretch_width"),
            pn.Column("### Candidates", self.candidate_table, sizing_mode="stretch_width"),
            sizing_mode="stretch_width",
        )

        main = pn.Column(
            self.title,
            self.status,
            pn.Row(
                pn.Card(self.summary, title="Case Summary", collapsed=False, sizing_mode="stretch_width"),
                pn.Card(self.selection_help, title="Selection State", collapsed=False, sizing_mode="stretch_width"),
                sizing_mode="stretch_width",
            ),
            pn.Card(self.instructions, title="Workflow", collapsed=False, sizing_mode="stretch_width"),
            pn.Card(self.plot, title="Ladder Trace", collapsed=False, sizing_mode="stretch_width"),
            pn.Card(review_tables, title="Steps And Candidates", collapsed=False, sizing_mode="stretch_width"),
            pn.Card(self.mapping_view, title="Current Mapping", collapsed=False, sizing_mode="stretch_width"),
            sizing_mode="stretch_width",
        )

        self.layout = pn.Row(sidebar, main, sizing_mode="stretch_width")

        self._update_filtered_records()
        if self.filtered_records:
            if self.case_manifest:
                self.file_select.param.update(value=None)
                self.current_record = None
                self.current_fsa = None
                self.preview_fsa = None
                self.current_meta = None
                self.current_candidates = pd.DataFrame(columns=["time", "intensity", "source"])
                self.ladder_steps = np.asarray([], dtype=float)
                self.mapping_times = {}
                self.manual_candidate_times = []
                self.preview_metrics = None
                self.preview_reason = "Select a file to load its ladder trace."
                self._set_status("Select one of the queued FLT3 files to load the trace.", level="info")
                self._refresh_views()
            else:
                first_path = str(self.filtered_records[0]["path"])
                self.file_select.param.update(value=first_path)
                selected = next(
                    (record for record in self.filtered_records if str(record["path"]) == first_path),
                    self.filtered_records[0],
                )
                self._load_record(selected)
        else:
            self._set_status("No FLT3 files found in this folder.", level="warning")
            self._refresh_views()

    def _load_case_manifest(self, manifest_path: Path | None) -> dict[str, dict[str, Any]]:
        if manifest_path is None:
            return {}
        df = pd.read_csv(manifest_path)
        path_col = next(
            (col for col in df.columns if str(col).lower() in {"path", "file_path", "fsa_path"}),
            None,
        )
        if path_col is None:
            raise ValueError(f"Case manifest {manifest_path} must include a path column.")

        manifest: dict[str, dict[str, Any]] = {}
        for idx, row in df.iterrows():
            path = Path(str(row[path_col])).expanduser()
            if not path.is_absolute():
                path = (manifest_path.parent / path).resolve()
            else:
                path = path.resolve()
            payload = row.to_dict()
            payload["queue_order"] = int(payload.get("queue_order", idx + 1) or idx + 1)
            manifest[str(path)] = payload
        return manifest

    def _scan_records(self) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for path in sorted(self.data_dir.rglob("*.fsa")):
            resolved_path = path.resolve()
            manifest_row = self.case_manifest.get(str(resolved_path)) if self.case_manifest else None
            if self.case_manifest and manifest_row is None:
                continue
            meta = detect_fsa_for_ladder(path, preferred_analysis="flt3")
            if not meta or meta.get("analysis") != "flt3":
                continue
            raw = meta.get("raw", {}) if isinstance(meta.get("raw"), dict) else {}
            record = {
                "path": resolved_path,
                "relative_path": str(resolved_path.relative_to(self.data_dir)),
                "assay": str(meta.get("assay") or raw.get("assay") or ""),
                "injection_time": int(raw.get("injection_time", 0) or 0),
                "specimen_id": str(raw.get("specimen_id") or ""),
                "selection_key": str(raw.get("selection_key") or ""),
                "source_run_dir": str(raw.get("source_run_dir") or ""),
                "well_id": str(raw.get("well_id") or ""),
                "analysis_type": str(raw.get("analysis_type") or ""),
                "ladder": str(meta.get("ladder") or raw.get("ladder") or ""),
                "meta": meta,
                "case_manifest": manifest_row or {},
            }
            records.append(record)
        return records

    def _set_status(self, text: str, *, level: str = "info") -> None:
        color = {"info": "#475569", "success": "#15803d", "warning": "#b45309", "error": "#b91c1c"}.get(
            level, "#475569"
        )
        self.status_text = text
        self.status.object = f"<div style='color:{color};font-weight:600'>{text}</div>"

    def _record_label(self, record: dict[str, Any]) -> str:
        specimen = record.get("specimen_id") or "Unknown"
        assay = record.get("assay") or "Unknown"
        inj = int(record.get("injection_time", 0) or 0)
        preferred = _preferred_injection_for_record(record)
        preferred_tag = "preferred" if inj == preferred else "non-preferred"
        label = (
            f"{specimen} | {assay} | {inj}s {preferred_tag} | "
            f"{record.get('well_id') or '?'} | {record['path'].name}"
        )
        manifest_row = record.get("case_manifest") or {}
        if manifest_row:
            label = (
                f"{int(manifest_row.get('queue_order', 0) or 0):03d} | "
                f"{manifest_row.get('queue_group', manifest_row.get('status', 'case'))} | {label}"
            )
        qc = self.qc_cache.get(str(record["path"]))
        if qc:
            label += f" | r2 {float(qc.get('r2', float('nan'))):.6f} | max {float(qc.get('max_abs_error_bp', float('nan'))):.2f} bp"
        return label

    def _step_button_action(self, step_idx: int) -> None:
        if not (0 <= int(step_idx) < len(self.ladder_steps)):
            return

        peak_time = self._selected_candidate_time()
        if peak_time is not None:
            candidate_idx = self._candidate_index_for_time(peak_time, tolerance=2.0)
            source = "manual"
            if candidate_idx is not None:
                source = str(self.current_candidates.iloc[candidate_idx].get("source", "auto"))
                peak_time = float(self.current_candidates.iloc[candidate_idx]["time"])
            self._assign_time_to_step(int(step_idx), float(peak_time), source=source)
            return

        if self.manual_time_input.value is not None:
            try:
                peak_time, source = self._resolve_manual_peak_time()
            except Exception:
                self.step_select.value = int(step_idx)
                return
            existing_idx = self._candidate_index_for_time(peak_time, tolerance=2.0)
            if existing_idx is not None:
                source = str(self.current_candidates.iloc[existing_idx].get("source", "auto"))
                peak_time = float(self.current_candidates.iloc[existing_idx]["time"])
            self._assign_time_to_step(int(step_idx), float(peak_time), source=source)
            return

        self.step_select.value = int(step_idx)

    def _rebuild_step_buttons(self) -> None:
        self.step_buttons = {}
        children: list[object] = []
        if self.ladder_steps.size == 0:
            self.step_button_box.objects = [
                pn.pane.Markdown("No ladder steps loaded.", styles={"color": "#64748b"})
            ]
            return

        selected_step = self._selected_step_index()
        ordered = self._ordered_step_indices()
        for idx in ordered:
            step_bp = float(self.ladder_steps[idx])
            assigned = self.mapping_times.get(idx)
            if assigned is None:
                name = f"{step_bp:.0f}"
                button_type = "default"
            else:
                name = f"{step_bp:.0f} *"
                button_type = "success"
            if selected_step == idx:
                button_type = "primary" if assigned is None else "warning"
            btn = pn.widgets.Button(
                name=name,
                button_type=button_type,
                width=66,
                height=40,
                margin=(0, 6, 6, 0),
            )
            btn.on_click(lambda _event, step_idx=idx: self._step_button_action(step_idx))
            self.step_buttons[int(idx)] = btn
            children.append(btn)
        self.step_button_box.objects = children

    def _ensure_record_qc(self, record: dict[str, Any]) -> dict[str, Any] | None:
        path_key = str(record["path"])
        cached = self.qc_cache.get(path_key)
        if cached is not None:
            return cached

        try:
            meta = record.get("meta") if isinstance(record.get("meta"), dict) else {}
            raw = meta.get("raw") if isinstance(meta.get("raw"), dict) else {}
            fsa = _analyse_fsa_candidate(
                record["path"],
                str(raw.get("primary_peak_channel") or meta.get("primary_peak_channel") or "DATA1"),
                str(record.get("assay") or ""),
                str(record.get("analysis_type") or ""),
            )
            metrics = compute_ladder_qc_metrics(fsa)
            cached = {
                "r2": float(metrics.get("r2", float("nan"))),
                "max_abs_error_bp": float(metrics.get("max_abs_error_bp", float("nan"))),
                "mean_abs_error_bp": float(metrics.get("mean_abs_error_bp", float("nan"))),
                "strategy": getattr(fsa, "ladder_fit_strategy", "auto_full"),
            }
        except Exception as exc:
            cached = {
                "r2": float("nan"),
                "max_abs_error_bp": float("inf"),
                "mean_abs_error_bp": float("inf"),
                "strategy": f"error: {exc}",
            }
        self.qc_cache[path_key] = cached
        return cached

    def _pipeline_context_for_record(self, record: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], str, str, str]:
        meta = record.get("meta") if isinstance(record.get("meta"), dict) else {}
        raw = meta.get("raw") if isinstance(meta.get("raw"), dict) else {}
        primary_channel = str(raw.get("primary_peak_channel") or meta.get("primary_peak_channel") or "DATA1")
        assay = str(record.get("assay") or raw.get("assay") or meta.get("assay") or "")
        analysis_type = str(record.get("analysis_type") or raw.get("analysis_type") or meta.get("analysis_type") or "")
        return meta, raw, primary_channel, assay, analysis_type

    def _load_pipeline_rescued_fsa(self, record: dict[str, Any]):
        meta, _raw, primary_channel, assay, analysis_type = self._pipeline_context_for_record(record)
        fsa = _analyse_fsa_candidate(
            record["path"],
            primary_channel,
            assay,
            analysis_type,
        )
        if fsa is None:
            raise ValueError("Pipeline analysis did not return a ladder fit.")
        return fsa, meta

    def _sort_filtered_records(self, records: list[dict[str, Any]]) -> list[dict[str, Any]]:
        if self.case_manifest and self.sort_mode.value == "Default":
            return sorted(
                records,
                key=lambda record: (
                    int((record.get("case_manifest") or {}).get("queue_order", 10**9) or 10**9),
                    record["path"].name,
                ),
            )
        if self.sort_mode.value == "Default":
            return records

        ranked: list[tuple[tuple[Any, ...], dict[str, Any]]] = []
        for record in records:
            qc = self.qc_cache.get(str(record["path"]))
            if qc is None:
                ranked.append(((1,), record))
                continue
            if self.sort_mode.value == "Worst residual":
                ranked.append(
                    (
                        (
                            0,
                            -float(qc.get("max_abs_error_bp", float("-inf"))),
                            float(qc.get("r2", float("inf"))),
                            record["path"].name,
                        ),
                        record,
                    )
                )
            else:
                ranked.append(
                    (
                        (
                            0,
                            float(qc.get("r2", float("inf"))),
                            -float(qc.get("max_abs_error_bp", float("-inf"))),
                            record["path"].name,
                        ),
                        record,
                    )
                )
        return [record for _, record in sorted(ranked, key=lambda item: item[0])]

    def _compute_qc_for_filtered_records(self) -> None:
        if not self.filtered_records:
            self._set_status("No FLT3 files match the current filters.", level="warning")
            return

        total = len(self.filtered_records)
        for idx, record in enumerate(self.filtered_records, start=1):
            self._set_status(f"Computing QC for filtered files: {idx}/{total} ...", level="info")
            self._ensure_record_qc(record)

        self._update_filtered_records()
        self._set_status(f"Computed QC for {total} filtered FLT3 files.", level="success")

    def _update_filtered_records(self) -> None:
        assay_filter = self.assay_filter.value
        injection_filter = self.injection_filter.value
        search_text = str(self.search_input.value or "").strip().lower()

        filtered: list[dict[str, Any]] = []
        for record in self.records:
            if assay_filter != "All" and record.get("assay") != assay_filter:
                continue
            record_injection = f"{int(record.get('injection_time', 0) or 0)}s"
            if injection_filter != "All" and record_injection != injection_filter:
                continue
            if self.preferred_only.value and int(record.get("injection_time", 0) or 0) != _preferred_injection_for_record(record):
                continue
            haystack = " ".join(
                [
                    str(record.get("specimen_id") or ""),
                    str(record.get("assay") or ""),
                    str(record.get("relative_path") or ""),
                    str(record.get("source_run_dir") or ""),
                    str(record.get("well_id") or ""),
                    str(record.get("analysis_type") or ""),
                    str(record["path"].name),
                ]
            ).lower()
            if search_text and search_text not in haystack:
                continue
            filtered.append(record)

        filtered = self._sort_filtered_records(filtered)
        self.filtered_records = filtered
        options = {self._record_label(record): str(record["path"]) for record in filtered}
        current_value = self.file_select.value
        self.file_select.options = options

        if current_value in options.values():
            self.file_select.value = current_value
        elif filtered:
            self.file_select.value = str(filtered[0]["path"])
        else:
            self.file_select.value = None

    def _on_filter_change(self, _event) -> None:
        current_value = self.file_select.value
        self._update_filtered_records()
        if current_value and any(str(record["path"]) == current_value for record in self.filtered_records):
            return
        if self.filtered_records:
            self._load_record(self.filtered_records[0])
        else:
            self.current_record = None
            self.current_fsa = None
            self.preview_fsa = None
            self.current_meta = None
            self.current_candidates = pd.DataFrame(columns=["time", "intensity", "source"])
            self.ladder_steps = np.asarray([], dtype=float)
            self.mapping_times = {}
            self.manual_candidate_times = []
            self.preview_metrics = None
            self.preview_reason = "No matching FLT3 files."
            self._set_status("No matching FLT3 files after filtering.", level="warning")
            self._refresh_views()

    def _on_file_change(self, event) -> None:
        value = event.new
        if not value:
            return
        self.pending_record_path = str(value)
        if self.case_manifest:
            self._set_status("Selection updated. Click 'Load Selected File' to open the trace.", level="info")
            return
        for record in self.filtered_records:
            if str(record["path"]) == str(value):
                self._load_record(record)
                return

    def _load_selected_record(self) -> None:
        value = self.pending_record_path or self.file_select.value
        if not value:
            self._set_status("Select a file first.", level="warning")
            return
        for record in self.filtered_records:
            if str(record["path"]) == str(value):
                self._load_record(record)
                return
        self._set_status("Could not find the selected file in the current queue.", level="error")

    def _move_file(self, delta: int) -> None:
        if not self.filtered_records or not self.current_record:
            return
        paths = [str(record["path"]) for record in self.filtered_records]
        try:
            idx = paths.index(str(self.current_record["path"]))
        except ValueError:
            idx = 0
        new_idx = max(0, min(len(paths) - 1, idx + delta))
        self.file_select.value = paths[new_idx]

    def _rescan_files(self) -> None:
        self.records = self._scan_records()
        self.records_by_selection_key = defaultdict(list)
        for record in self.records:
            self.records_by_selection_key[str(record.get("selection_key") or "")].append(record)
        self._update_filtered_records()
        if self.filtered_records:
            self._set_status(f"Rescanned folder and found {len(self.filtered_records)} filtered FLT3 files.", level="success")
            selected = next((r for r in self.filtered_records if str(r["path"]) == str(self.file_select.value)), self.filtered_records[0])
            self._load_record(selected)
        else:
            self._set_status("Rescanned folder, but no FLT3 files matched the current filters.", level="warning")
            self._refresh_views()

    def _reload_current_file(self) -> None:
        if not self.current_record:
            return
        self._load_record(self.current_record)

    def _on_fit_source_change(self, _event) -> None:
        self._reload_current_file()

    def _load_record(self, record: dict[str, Any] | None) -> None:
        if record is None:
            return
        self.current_record = record
        loaded_source_label = self.fit_source.value.lower()
        try:
            if self.fit_source.value == "Pipeline rescue":
                fsa, meta = self._load_pipeline_rescued_fsa(record)
            else:
                fsa, meta = load_adjustable_fsa(
                    record["path"],
                    preferred_analysis="flt3",
                    metadata=record["meta"],
                )
        except Exception as exc:
            if self.fit_source.value == "Pipeline rescue":
                try:
                    fsa, meta = load_adjustable_fsa(
                        record["path"],
                        preferred_analysis="flt3",
                        metadata=record["meta"],
                    )
                    loaded_source_label = "raw adjustable fallback"
                    self._set_status(
                        (
                            f"Pipeline fit was unavailable for {record['path'].name}. "
                            "Opened raw adjustable ladder state instead."
                        ),
                        level="warning",
                    )
                except Exception as raw_exc:
                    self.current_fsa = None
                    self.preview_fsa = None
                    self.current_meta = record.get("meta")
                    self.current_candidates = pd.DataFrame(columns=["time", "intensity", "source"])
                    self.ladder_steps = np.asarray([], dtype=float)
                    self.mapping_times = {}
                    self.manual_candidate_times = []
                    self.preview_metrics = None
                    self.preview_reason = str(raw_exc)
                    self._set_status(
                        f"Could not load ladder state for {record['path'].name}: {raw_exc}",
                        level="error",
                    )
                    self._refresh_views()
                    return
            else:
                self.current_fsa = None
                self.preview_fsa = None
                self.current_meta = record.get("meta")
                self.current_candidates = pd.DataFrame(columns=["time", "intensity", "source"])
                self.ladder_steps = np.asarray([], dtype=float)
                self.mapping_times = {}
                self.manual_candidate_times = []
                self.preview_metrics = None
                self.preview_reason = str(exc)
                self._set_status(f"Could not load ladder state for {record['path'].name}: {exc}", level="error")
                self._refresh_views()
                return

        self.current_fsa = fsa
        self.preview_fsa = None
        self.current_meta = meta
        self.ladder_steps = np.asarray(
            getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])),
            dtype=float,
        )

        payload = load_ladder_adjustment(fsa)
        base_mapping_times = self._infer_mapping_times_from_payload(payload, fsa)
        payload_manual_candidates = (payload or {}).get("manual_candidates", [])
        if not payload_manual_candidates:
            payload_manual_candidates = getattr(fsa, "ladder_review_manual_candidates", []) or []
        self.manual_candidate_times = sorted(
            {
                float(value)
                for value in payload_manual_candidates
                if value is not None and np.isfinite(float(value))
            }
        )
        self.mapping_times = {
            int(step_idx): float(time_value)
            for step_idx, time_value in dict(base_mapping_times).items()
            if time_value is not None and np.isfinite(float(time_value))
        }

        candidates = get_ladder_candidates(fsa).copy()
        if "source" not in candidates.columns:
            candidates["source"] = "auto"
        if candidates.empty:
            trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
            candidates = _bootstrap_candidates_from_trace(trace)
        self.current_candidates = candidates.reset_index(drop=True)

        for time_value in list(self.mapping_times.values()) + list(self.manual_candidate_times):
            self._ensure_candidate_time(time_value, source="manual")

        self._refresh_preview()
        if loaded_source_label != "raw adjustable fallback":
            self._set_status(f"Loaded {record['path'].name} using {loaded_source_label} fit.", level="success")
        self._refresh_views()

    def _infer_mapping_times_from_payload(self, payload: dict | None, fsa) -> dict[int, float]:
        if payload and payload.get("mapping_times"):
            return {
                int(k): float(v)
                for k, v in payload["mapping_times"].items()
                if v is not None and np.isfinite(float(v))
            }
        review_mapping = getattr(fsa, "ladder_review_mapping_times", None)
        if review_mapping:
            return {
                int(k): float(v)
                for k, v in dict(review_mapping).items()
                if v is not None and np.isfinite(float(v))
            }

        mapping: dict[int, float] = {}
        expected = np.atleast_1d(np.asarray(
            getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", [])),
            dtype=float,
        ))
        fitted_steps = np.atleast_1d(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float))
        fitted_peaks = np.atleast_1d(np.asarray(getattr(fsa, "best_size_standard", []), dtype=float))
        if expected.size == 0 or fitted_steps.size == 0 or fitted_peaks.size == 0:
            return mapping

        for step_bp, peak_time in zip(fitted_steps, fitted_peaks):
            if not np.isfinite(float(step_bp)) or not np.isfinite(float(peak_time)):
                continue
            matches = np.where(np.isclose(expected, float(step_bp), atol=1e-6))[0]
            if matches.size:
                mapping[int(matches[0])] = float(peak_time)
        return mapping

    def _candidate_index_for_time(self, peak_time: float, tolerance: float = 1e-6) -> int | None:
        if self.current_candidates.empty:
            return None
        diff = (self.current_candidates["time"].astype(float) - float(peak_time)).abs()
        matches = diff[diff <= tolerance]
        if matches.empty:
            return None
        return int(matches.index[0])

    def _ensure_candidate_time(self, peak_time: float, *, source: str) -> int:
        if not np.isfinite(float(peak_time)):
            return -1
        existing = self._candidate_index_for_time(peak_time)
        if existing is not None:
            return existing
        trace = np.asarray(getattr(self.current_fsa, "size_standard", []), dtype=float)
        peak_idx = int(round(float(peak_time)))
        intensity = float(trace[peak_idx]) if 0 <= peak_idx < trace.size else float("nan")
        row = pd.DataFrame(
            [
                {
                    "index": len(self.current_candidates),
                    "time": float(peak_time),
                    "intensity": intensity,
                    "source": source,
                }
            ]
        )
        self.current_candidates = pd.concat([self.current_candidates, row], ignore_index=True)
        self.current_candidates = self.current_candidates.sort_values(["time", "intensity"]).reset_index(drop=True)
        return int(self._candidate_index_for_time(peak_time) or 0)

    def _selected_step_index(self) -> int | None:
        value = self.step_select.value
        if value is None or value == "":
            return None
        return int(value)

    def _selected_candidate_time(self) -> float | None:
        value = self.candidate_select.value
        if value in {None, ""}:
            return None
        return float(value)

    def _find_local_peak_time(self, x_value: float, search_radius: int = 18) -> tuple[float, float]:
        trace = np.asarray(getattr(self.current_fsa, "size_standard", []), dtype=float)
        if trace.size == 0:
            raise ValueError("No ladder trace available.")
        center = int(round(float(x_value)))
        lo = max(center - search_radius, 0)
        hi = min(center + search_radius + 1, trace.size)
        if lo >= hi:
            raise ValueError("Selected ladder region could not be inspected.")
        window = trace[lo:hi]
        local_index = int(np.argmax(window))
        peak_index = lo + local_index
        return float(peak_index), float(trace[peak_index])

    def _resolve_manual_peak_time(self) -> tuple[float, str]:
        if self.manual_time_input.value is None:
            raise ValueError("Enter a peak time or click the plot first.")

        raw_time = float(self.manual_time_input.value)
        if self.manual_mode.value == "Use exact time":
            return float(raw_time), "manual"

        peak_time, _intensity = self._find_local_peak_time(raw_time)
        return float(peak_time), "manual"

    def _assign_time_to_step(self, step_idx: int, peak_time: float, *, source: str) -> None:
        if step_idx < 0 or step_idx >= len(self.ladder_steps):
            return
        self.mapping_times[int(step_idx)] = float(peak_time)
        if source == "manual":
            if not any(np.isclose(float(existing), float(peak_time), atol=1e-6) for existing in self.manual_candidate_times):
                self.manual_candidate_times.append(float(peak_time))
                self.manual_candidate_times.sort()
            self._ensure_candidate_time(peak_time, source="manual")
            self.candidate_select.value = float(peak_time)
        self.last_clicked_time = float(peak_time)
        self._refresh_preview()
        self._refresh_views()
        next_missing = self._next_missing_step()
        if next_missing is not None:
            self.step_select.value = next_missing
        self._set_status(
            f"Assigned {self.ladder_steps[step_idx]:.0f} bp to time {float(peak_time):.1f}.",
            level="success",
        )

    def _assign_selected_candidate(self) -> None:
        step_idx = self._selected_step_index()
        peak_time = self._selected_candidate_time()
        if step_idx is None or peak_time is None:
            self._set_status("Select both a ladder step and a candidate peak first.", level="warning")
            return
        candidate_idx = self._candidate_index_for_time(peak_time)
        source = "manual"
        if candidate_idx is not None:
            source = str(self.current_candidates.iloc[candidate_idx].get("source", "auto"))
        self._assign_time_to_step(step_idx, peak_time, source=source)

    def _assign_manual_time_input(self) -> None:
        step_idx = self._selected_step_index()
        if step_idx is None:
            self._set_status("Select a ladder step first.", level="warning")
            return
        try:
            peak_time, source = self._resolve_manual_peak_time()
        except Exception as exc:
            self._set_status(f"Could not resolve manual peak time: {exc}", level="error")
            return
        self.manual_time_input.value = float(peak_time)
        existing_idx = self._candidate_index_for_time(peak_time, tolerance=2.0)
        if existing_idx is not None:
            source = str(self.current_candidates.iloc[existing_idx].get("source", "auto"))
            peak_time = float(self.current_candidates.iloc[existing_idx]["time"])
        self._assign_time_to_step(step_idx, peak_time, source=source)

    def _add_manual_candidate(self) -> None:
        if self.current_fsa is None:
            self._set_status("No FLT3 file is loaded.", level="warning")
            return
        try:
            peak_time, _source = self._resolve_manual_peak_time()
        except Exception as exc:
            self._set_status(f"Could not add manual candidate: {exc}", level="error")
            return

        self.manual_time_input.value = float(peak_time)
        candidate_idx = self._ensure_candidate_time(float(peak_time), source="manual")
        if not any(np.isclose(float(existing), float(peak_time), atol=1e-6) for existing in self.manual_candidate_times):
            self.manual_candidate_times.append(float(peak_time))
            self.manual_candidate_times.sort()
        self.last_clicked_time = float(peak_time)
        self._refresh_preview()
        self._refresh_views()
        self.candidate_select.value = float(peak_time)
        candidate_df = self.candidate_table.value
        if candidate_df is not None and not candidate_df.empty:
            matches = candidate_df.index[(candidate_df["time"].astype(float) - float(peak_time)).abs() <= 1e-6]
            if len(matches):
                self.candidate_table.selection = [int(matches[0])]
        self._set_status(
            f"Added manual candidate at {float(peak_time):.1f} (candidate #{candidate_idx}).",
            level="success",
        )

    def _clear_selected_step(self) -> None:
        step_idx = self._selected_step_index()
        if step_idx is None:
            self._set_status("Select a ladder step first.", level="warning")
            return
        if step_idx in self.mapping_times:
            previous_time = float(self.mapping_times[step_idx])
            del self.mapping_times[step_idx]
            self.manual_time_input.value = previous_time
            self._refresh_preview()
            self._refresh_views()
            self._set_status(f"Cleared ladder step {self.ladder_steps[step_idx]:.0f} bp.", level="success")

    def _build_adjustment_payload(self) -> dict[str, Any]:
        return {
            "mapping": {},
            "mapping_times": {int(k): float(v) for k, v in self.mapping_times.items()},
            "manual_candidates": [float(v) for v in self.manual_candidate_times],
        }

    def _refresh_preview(self) -> None:
        self.preview_fsa = None
        self.preview_metrics = None
        if self.current_fsa is None or self.ladder_steps.size == 0:
            self.preview_reason = "No file selected."
            return
        if len(self.mapping_times) < len(self.ladder_steps):
            self.preview_reason = f"Mapped {len(self.mapping_times)}/{len(self.ladder_steps)} ladder steps."
            return
        try:
            preview = copy.deepcopy(self.current_fsa)
            preview.expected_ladder_steps = np.asarray(self.ladder_steps, dtype=float).copy()
            preview.ladder_steps = np.asarray(self.ladder_steps, dtype=float).copy()
            preview = apply_manual_ladder_mapping(preview, self._build_adjustment_payload())
            self.preview_fsa = preview
            self.preview_metrics = compute_ladder_qc_metrics(preview)
            _, grade_reason = _fit_grade(self.preview_metrics)
            self.preview_reason = grade_reason
        except Exception as exc:
            self.preview_fsa = None
            self.preview_metrics = None
            self.preview_reason = str(exc)

    def _save_adjustment(self) -> None:
        if self.current_fsa is None:
            self._set_status("No FLT3 file is loaded.", level="warning")
            return
        if len(self.mapping_times) < len(self.ladder_steps):
            self._set_status("Map every ladder step before saving.", level="warning")
            return
        try:
            save_ladder_adjustment(self.current_fsa, self._build_adjustment_payload())
            self._set_status(f"Saved ladder adjustment for {self.current_record['path'].name}.", level="success")
            self._load_record(self.current_record)
        except Exception as exc:
            self._set_status(f"Could not save ladder adjustment: {exc}", level="error")

    def _delete_adjustment(self) -> None:
        if not self.current_record:
            return
        adj_path = self.current_record["path"].with_suffix(".ladder_adj.json")
        if not adj_path.exists():
            self._set_status("No saved ladder adjustment exists for this file.", level="warning")
            return
        adj_path.unlink(missing_ok=True)
        self._set_status(f"Deleted saved ladder adjustment for {self.current_record['path'].name}.", level="success")
        self._load_record(self.current_record)

    def _on_plot_click(self, event) -> None:
        click_data = getattr(event, "new", None)
        if not click_data or not self.current_fsa:
            return
        points = click_data.get("points") or []
        if not points:
            return
        try:
            x_value = float(points[0].get("x"))
        except Exception as exc:
            self._set_status(f"Could not read clicked trace position: {exc}", level="error")
            return

        self.last_clicked_time = float(x_value)
        self.manual_time_input.value = float(x_value)
        existing_idx = self._candidate_index_for_time(x_value, tolerance=8.0)
        if existing_idx is not None:
            peak_time = float(self.current_candidates.iloc[existing_idx]["time"])
            self.candidate_select.value = peak_time
            candidate_df = self.candidate_table.value
            if candidate_df is not None and not candidate_df.empty:
                matches = candidate_df.index[(candidate_df["time"].astype(float) - peak_time).abs() <= 1e-6]
                if len(matches):
                    self.candidate_table.selection = [int(matches[0])]
            self._set_status(
                f"Clicked {x_value:.1f}. Nearest candidate at {peak_time:.1f} selected. Use exact mode if you want the clicked time itself.",
                level="info",
            )
        else:
            self._set_status(
                f"Clicked trace position {x_value:.1f}. Use 'Add Manual Candidate' or 'Assign Manual Time' to keep it.",
                level="info",
            )

    def _next_missing_step(self) -> int | None:
        for idx in self._ordered_step_indices():
            if idx not in self.mapping_times:
                return idx
        return None

    def _ordered_step_indices(self) -> list[int]:
        indices = list(range(len(self.ladder_steps)))
        if self.descending_steps.value:
            indices.reverse()
        return indices

    def _fit_residual_for_step(self, step_idx: int) -> float | None:
        if self.preview_fsa is None:
            return None
        peak_time = self.mapping_times.get(step_idx)
        if peak_time is None:
            return None
        ladder_model = getattr(self.preview_fsa, "ladder_model", None)
        if ladder_model is None:
            return None
        try:
            fitted_bp = float(ladder_model.predict(np.array([[float(peak_time)]], dtype=float))[0])
        except Exception:
            return None
        return float(fitted_bp - self.ladder_steps[step_idx])

    def _related_records_markdown(self) -> str:
        if not self.current_record:
            return ""
        selection_key = str(self.current_record.get("selection_key") or "")
        related = self.records_by_selection_key.get(selection_key, [])
        if len(related) <= 1:
            return ""
        lines = ["", "**Related Injection Twins**"]
        for record in sorted(related, key=lambda item: (int(item.get("injection_time", 0) or 0), item["path"].name)):
            inj = int(record.get("injection_time", 0) or 0)
            preferred = _preferred_injection_for_record(record)
            tag = "preferred" if inj == preferred else "other"
            current = "current" if record["path"] == self.current_record["path"] else ""
            lines.append(
                f"- `{inj}s` {tag} {current} | `{record['path'].name}`"
            )
        return "\n".join(lines)

    def _build_summary_markdown(self) -> str:
        if not self.current_record or self.current_fsa is None:
            return "No FLT3 file selected."

        raw = self.current_meta.get("raw", {}) if isinstance(self.current_meta, dict) else {}
        ladder_strategy = str(getattr(self.current_fsa, "ladder_fit_strategy", "auto_full")).replace("_", " ")
        review_required = bool(getattr(self.current_fsa, "ladder_review_required", False))
        missing_steps = list(map(float, getattr(self.current_fsa, "ladder_missing_expected_steps", [])))
        saved_adj = "yes" if load_ladder_adjustment(self.current_fsa) else "no"
        preview_grade, preview_label = _fit_grade(self.preview_metrics)
        preview_r2 = float(self.preview_metrics.get("r2", float("nan"))) if self.preview_metrics else float("nan")
        preview_max = (
            float(self.preview_metrics.get("max_abs_error_bp", float("nan")))
            if self.preview_metrics else float("nan")
        )

        selection_idx = 0
        if self.filtered_records:
            try:
                selection_idx = [str(record["path"]) for record in self.filtered_records].index(str(self.current_record["path"])) + 1
            except ValueError:
                selection_idx = 0

        lines = [
            f"### {self.current_record['path'].name}",
            "",
            f"- File `{selection_idx}/{len(self.filtered_records)}` in current filter set",
            f"- Display source `{self.fit_source.value}`",
            f"- Assay `{self.current_record.get('assay')}` | injection `{int(self.current_record.get('injection_time', 0) or 0)}s` | preferred `{_preferred_injection_for_record(self.current_record)}s`",
            f"- Specimen `{self.current_record.get('specimen_id') or 'unknown'}` | well `{self.current_record.get('well_id') or '?'}`",
            f"- Analysis type `{self.current_record.get('analysis_type') or 'unknown'}` | ladder `{self.current_record.get('ladder') or 'unknown'}`",
            f"- Selection key `{self.current_record.get('selection_key') or 'unknown'}`",
            f"- Current fit `{ladder_strategy}` | review required `{review_required}` | saved adjustment `{saved_adj}`",
            f"- Expected ladder steps `{len(self.ladder_steps)}` | mapped now `{len(self.mapping_times)}`",
            f"- Missing in loaded fit `{', '.join(f'{bp:.0f}' for bp in missing_steps) if missing_steps else 'none'}`",
        ]
        manifest_row = self.current_record.get("case_manifest") or {}
        if manifest_row:
            lines.extend(
                [
                    f"- Review queue `{manifest_row.get('queue_order')}` | group `{manifest_row.get('queue_group', manifest_row.get('status', 'case'))}`",
                    (
                        f"- Queue metrics R2 `{float(manifest_row.get('r2', float('nan'))):.6f}` | "
                        f"max residual `{float(manifest_row.get('max_bp_err', float('nan'))):.2f} bp` | "
                        f"mean residual `{float(manifest_row.get('mean_bp_err', float('nan'))):.2f} bp`"
                    ),
                ]
            )

        if self.preview_metrics:
            lines.append(
                f"- Preview `{preview_grade}` | R2 `{preview_r2:.6f}` | max residual `{preview_max:.2f} bp` | {preview_label}"
            )
        else:
            lines.append(f"- Preview pending | {self.preview_reason}")

        relative_path = self.current_record["relative_path"]
        lines.append(f"- Relative path `{relative_path}`")

        run_name = raw.get("run_name")
        if run_name:
            lines.append(f"- Run `{run_name}`")

        lines.append(self._related_records_markdown())
        return "\n".join(lines).strip()

    def _build_mapping_markdown(self) -> str:
        if self.ladder_steps.size == 0:
            return "No ladder mapping loaded."

        rows = [
            "| Ladder Step | Assigned Time | Source | Residual |",
            "| --- | ---: | --- | ---: |",
        ]
        for idx, step_bp in enumerate(self.ladder_steps):
            peak_time = self.mapping_times.get(idx)
            if peak_time is None:
                rows.append(f"| {step_bp:.0f} bp | — | missing | — |")
                continue
            candidate_idx = self._candidate_index_for_time(peak_time)
            source = "manual"
            if candidate_idx is not None:
                source = str(self.current_candidates.iloc[candidate_idx].get("source", "auto"))
            residual = self._fit_residual_for_step(idx)
            residual_text = "—" if residual is None else f"{residual:+.2f} bp"
            rows.append(
                f"| {step_bp:.0f} bp | {peak_time:.1f} | {_markdown_escape(source)} | {_markdown_escape(residual_text)} |"
            )
        rows.append("")
        rows.append(f"Preview note: {self.preview_reason}")
        return "\n".join(rows)

    def _step_table_df(self) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        for idx in self._ordered_step_indices():
            step_bp = self.ladder_steps[idx]
            peak_time = self.mapping_times.get(idx)
            if peak_time is None:
                rows.append(
                    {
                        "step_idx": idx,
                        "step_bp": float(step_bp),
                        "assigned_time": np.nan,
                        "source": "missing",
                        "residual": np.nan,
                        "status": "missing",
                    }
                )
                continue
            candidate_idx = self._candidate_index_for_time(peak_time, tolerance=2.0)
            source = "manual"
            if candidate_idx is not None:
                source = str(self.current_candidates.iloc[candidate_idx].get("source", "auto"))
            residual = self._fit_residual_for_step(idx)
            rows.append(
                {
                    "step_idx": idx,
                    "step_bp": float(step_bp),
                    "assigned_time": float(peak_time),
                    "source": source,
                    "residual": residual,
                    "status": "mapped",
                }
            )
        return pd.DataFrame(rows)

    def _candidate_table_df(self) -> pd.DataFrame:
        if self.current_candidates.empty:
            return pd.DataFrame(columns=["candidate_idx", "time", "intensity", "source", "used_by"])
        used_by: dict[float, float] = {}
        for step_idx, peak_time in self.mapping_times.items():
            used_by[float(peak_time)] = float(self.ladder_steps[step_idx])
        rows: list[dict[str, Any]] = []
        sorted_candidates = self.current_candidates.sort_values(["time", "intensity"], ascending=[True, False]).reset_index()
        for display_idx, row in sorted_candidates.iterrows():
            peak_time = float(row["time"])
            rows.append(
                {
                    "candidate_idx": int(display_idx),
                    "time": peak_time,
                    "intensity": float(row.get("intensity", float("nan"))),
                    "source": str(row.get("source", "auto")),
                    "used_by": used_by.get(peak_time, np.nan),
                }
            )
        return pd.DataFrame(rows)

    def _on_step_table_select(self, event) -> None:
        selection = list(event.new or [])
        if not selection:
            return
        df = self.step_table.value
        if df is None or df.empty:
            return
        row_idx = int(selection[0])
        if row_idx < 0 or row_idx >= len(df):
            return
        self.step_select.value = int(df.iloc[row_idx]["step_idx"])
        peak_time = df.iloc[row_idx].get("assigned_time")
        if np.isfinite(peak_time):
            self.manual_time_input.value = float(peak_time)
        self._refresh_selection_help()

    def _on_candidate_table_select(self, event) -> None:
        selection = list(event.new or [])
        if not selection:
            return
        df = self.candidate_table.value
        if df is None or df.empty:
            return
        row_idx = int(selection[0])
        if row_idx < 0 or row_idx >= len(df):
            return
        peak_time = float(df.iloc[row_idx]["time"])
        self.candidate_select.value = peak_time
        self.manual_time_input.value = peak_time
        self._refresh_selection_help()

    def _on_step_select_change(self, event) -> None:
        value = event.new
        if value is None:
            return
        step_idx = int(value)
        if step_idx in self.mapping_times:
            self.manual_time_input.value = float(self.mapping_times[step_idx])
        self._refresh_selection_help()

    def _on_candidate_select_change(self, event) -> None:
        value = event.new
        if value in {None, ""}:
            return
        self.manual_time_input.value = float(value)
        self._refresh_selection_help()

    def _on_manual_mode_change(self, _event) -> None:
        self._refresh_selection_help()

    def _on_step_order_change(self, _event) -> None:
        self._refresh_views()

    def _refresh_selection_help(self) -> None:
        step_idx = self._selected_step_index()
        peak_time = self._selected_candidate_time()
        step_text = "none"
        if step_idx is not None and 0 <= step_idx < len(self.ladder_steps):
            assigned = self.mapping_times.get(step_idx)
            step_text = f"{self.ladder_steps[step_idx]:.0f} bp"
            if assigned is not None:
                step_text += f" currently at {assigned:.1f}"
        cand_text = "none"
        if peak_time is not None:
            cand_text = f"{peak_time:.1f}"
        click_text = "none" if self.last_clicked_time is None else f"{self.last_clicked_time:.1f}"
        manual_mode = str(self.manual_mode.value or "Snap to local peak")
        self.selection_help.object = (
            f"- Selected step: `{step_text}`\n"
            f"- Selected candidate: `{cand_text}`\n"
            f"- Last clicked position: `{click_text}`\n"
            f"- Manual mode: `{manual_mode}`\n"
            "- `Add Manual Candidate` stores a peak without assigning it yet.\n"
            "- `Assign Manual Time` maps the current manual time directly to the selected step.\n"
            "- The step pad can assign the current candidate/manual time directly to `500`, `490`, `450`, etc."
        )

    def _step_options(self) -> dict[str, int]:
        options: dict[str, int] = {}
        for idx in self._ordered_step_indices():
            step_bp = self.ladder_steps[idx]
            peak_time = self.mapping_times.get(idx)
            if peak_time is None:
                label = f"{step_bp:.0f} bp -> missing"
            else:
                label = f"{step_bp:.0f} bp -> {peak_time:.1f}"
            options[label] = int(idx)
        return options

    def _candidate_options(self) -> dict[str, float]:
        if self.current_candidates.empty:
            return {}

        used_by: dict[float, float] = {}
        for step_idx, peak_time in self.mapping_times.items():
            used_by[float(peak_time)] = float(self.ladder_steps[step_idx])

        options: dict[str, float] = {}
        sorted_candidates = self.current_candidates.sort_values(["time", "intensity"], ascending=[True, False])
        for _, row in sorted_candidates.iterrows():
            peak_time = float(row["time"])
            intensity = float(row.get("intensity", float("nan")))
            source = str(row.get("source", "auto"))
            used_text = ""
            if peak_time in used_by:
                used_text = f" | used by {used_by[peak_time]:.0f} bp"
            label = f"{peak_time:.1f} | y={intensity:.0f} | {source}{used_text}"
            options[label] = peak_time
        return options

    def _update_select_widgets(self) -> None:
        step_options = self._step_options()
        candidate_options = self._candidate_options()
        current_step = self.step_select.value
        current_candidate = self.candidate_select.value

        self.step_select.options = step_options
        self.candidate_select.options = candidate_options

        if current_step in step_options.values():
            self.step_select.value = current_step
        elif self.mapping_times:
            ordered = self._ordered_step_indices()
            missing_steps = [idx for idx in ordered if idx not in self.mapping_times]
            self.step_select.value = (missing_steps[0] if missing_steps else ordered[0]) if ordered else None
        else:
            ordered = self._ordered_step_indices()
            self.step_select.value = ordered[0] if ordered else None

        if current_candidate in candidate_options.values():
            self.candidate_select.value = current_candidate
        elif candidate_options:
            self.candidate_select.value = next(iter(candidate_options.values()))
        else:
            self.candidate_select.value = None
        selected_step = self._selected_step_index()
        if selected_step is not None and selected_step in self.mapping_times:
            self.manual_time_input.value = float(self.mapping_times[selected_step])
        self._rebuild_step_buttons()
        self._refresh_selection_help()

    def _build_plot(self) -> go.Figure:
        fig = go.Figure()
        if self.current_fsa is None:
            fig.update_layout(title="No FLT3 file loaded.", height=580)
            return fig

        trace_raw = np.asarray(getattr(self.current_fsa, "size_standard", []), dtype=float)
        if trace_raw.size == 0:
            fig.update_layout(title="No ladder trace available.", height=580)
            return fig

        try:
            baseline = estimate_running_baseline(trace_raw, bin_size=200, quantile=0.10)
            y_arr = np.clip(trace_raw - baseline, a_min=0, a_max=None)
        except Exception:
            y_arr = trace_raw.astype(float)

        def corrected_intensity(time_value: float) -> float:
            idx = int(round(float(time_value)))
            if idx < 0 or idx >= len(y_arr):
                return 0.0
            return float(y_arr[idx])

        x_arr = np.arange(len(y_arr))
        fig.add_trace(
            go.Scatter(
                x=x_arr,
                y=y_arr,
                mode="lines",
                line=dict(color="#1d4ed8", width=2),
                fill="tozeroy",
                fillcolor="rgba(29,78,216,0.08)",
                name="Ladder Trace",
            )
        )

        sorted_candidates = self.current_candidates.sort_values(["time", "intensity"], ascending=[True, False])
        used_by_time = {float(time): int(step_idx) for step_idx, time in self.mapping_times.items()}
        auto_x: list[float] = []
        auto_y: list[float] = []
        manual_x: list[float] = []
        manual_y: list[float] = []
        assigned_x: list[float] = []
        assigned_y: list[float] = []
        assigned_text: list[str] = []

        for _, row in sorted_candidates.iterrows():
            peak_time = float(row["time"])
            intensity = corrected_intensity(peak_time)
            source = str(row.get("source", "auto"))
            step_idx = used_by_time.get(peak_time)
            if step_idx is not None:
                assigned_x.append(peak_time)
                assigned_y.append(intensity)
                assigned_text.append(f"{self.ladder_steps[step_idx]:.0f} bp")
                continue
            if source == "manual":
                manual_x.append(peak_time)
                manual_y.append(intensity)
            else:
                auto_x.append(peak_time)
                auto_y.append(intensity)

        if auto_x:
            fig.add_trace(
                go.Scatter(
                    x=auto_x,
                    y=auto_y,
                    mode="markers",
                    marker=dict(color="#dc2626", size=9, symbol="x"),
                    name="Auto Candidates",
                )
            )
        if manual_x:
            fig.add_trace(
                go.Scatter(
                    x=manual_x,
                    y=manual_y,
                    mode="markers",
                    marker=dict(color="#0f766e", size=10, symbol="diamond"),
                    name="Manual Candidates",
                )
            )
        if assigned_x:
            fig.add_trace(
                go.Scatter(
                    x=assigned_x,
                    y=assigned_y,
                    mode="markers+text",
                    marker=dict(color="#f59e0b", size=13, symbol="star"),
                    text=assigned_text,
                    textposition="top center",
                    name="Assigned Steps",
                )
            )

        selected_step = self._selected_step_index()
        if selected_step is not None and selected_step in self.mapping_times:
            peak_time = float(self.mapping_times[selected_step])
            fig.add_vline(
                x=peak_time,
                line_dash="dash",
                line_color="#059669",
                opacity=0.8,
            )

        x_candidates = list(sorted_candidates["time"].astype(float)) if not sorted_candidates.empty else []
        if x_candidates:
            x_min = max(0.0, min(x_candidates) - 120.0)
            x_max = min(float(len(y_arr) - 1), max(x_candidates) + 120.0)
            fig.update_xaxes(range=[x_min, x_max])
        else:
            fig.update_xaxes(range=[1400.0, min(float(len(y_arr) - 1), 4700.0)])

        fig.update_layout(
            title="FLT3 GS500ROX Ladder Trace",
            xaxis_title="Time (scan index)",
            yaxis_title="Corrected intensity",
            yaxis=dict(range=[0, 4000]),
            height=580,
            margin=dict(l=30, r=20, t=45, b=30),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0),
            template="simple_white",
        )
        return fig

    def _refresh_views(self) -> None:
        self._update_select_widgets()
        step_df = self._step_table_df()
        candidate_df = self._candidate_table_df()
        self.step_table.value = step_df
        self.candidate_table.value = candidate_df
        self.step_table.disabled = step_df.empty
        self.candidate_table.disabled = candidate_df.empty
        self.summary.object = self._build_summary_markdown()
        self.mapping_view.object = self._build_mapping_markdown()
        self.plot.object = self._build_plot()


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    def create_app():
        app = Flt3LadderReviewApp(args.data_dir, case_manifest=args.case_manifest)
        return app.layout

    pn.serve(create_app, port=args.port, show=True, title="FLT3 Ladder Review")


if __name__ == "__main__":
    main()
