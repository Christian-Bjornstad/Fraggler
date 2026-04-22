from __future__ import annotations

import argparse
import logging
import re
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import panel as pn
import plotly.graph_objects as go


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gui_qt.ladder_utils import detect_fsa_for_ladder, load_adjustable_fsa  # noqa: E402

pn.extension("plotly", sizing_mode="stretch_width")
logging.basicConfig(level=logging.INFO)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a browser-based ladder review annotator for the exported ladder review bundle."
    )
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        required=True,
        help="Directory containing ladder_review_cases.csv and ladder_review_candidates.csv.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Root directory where the original FSA run folders are stored.",
    )
    parser.add_argument("--port", type=int, default=5007, help="Port for the Panel server.")
    return parser


def _normalize_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


ASSAY_FILENAME_PATTERNS = {
    "IKZF1": ["ikzf1"],
    "Ktr-albumin": ["ktralbumin", "ktralbuminr", "ktralbuminr13"],
    "TCRbA": ["trbmixa", "tcrba", "tcrb_a"],
    "TCRbB": ["trbmixb", "tcrbb", "tcrb_b"],
    "TCRbC": ["trbmixc", "tcrbc", "tcrb_c"],
    "TCRgA": ["tcrga", "tcrg_a"],
    "TCRgB": ["tcrgb", "tcrg_b"],
    "FR1": ["fr1"],
    "FR2": ["fr2"],
    "FR3": ["fr3"],
    "DHJH_D": ["dhjhmixd", "dhjhd"],
    "DHJH_E": ["dhjhmixe", "dhjhe"],
    "IGK": ["igk"],
    "KDE": ["kde"],
}


def resolve_fsa_path(data_dir: Path, source_run_dir: str, assay: str, well: str, run_code: str) -> Path:
    run_dir = data_dir / source_run_dir
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Missing run directory: {run_dir}")

    assay_tokens = ASSAY_FILENAME_PATTERNS.get(assay, [_normalize_token(assay)])
    well_token = _normalize_token(well)
    run_code_token = _normalize_token(run_code)

    matches: list[Path] = []
    for candidate in sorted(run_dir.glob("*.fsa")):
        token = _normalize_token(candidate.name)
        if run_code_token and run_code_token not in token:
            continue
        if well_token and well_token not in token:
            continue
        if not any(pattern in token for pattern in assay_tokens):
            continue
        matches.append(candidate)

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        return matches[0]

    # Fallback: match by well + run code only, then prefer assay-ish names.
    fallback: list[Path] = []
    for candidate in sorted(run_dir.glob("*.fsa")):
        token = _normalize_token(candidate.name)
        if run_code_token and run_code_token not in token:
            continue
        if well_token and well_token not in token:
            continue
        fallback.append(candidate)
    if len(fallback) == 1:
        return fallback[0]
    if fallback:
        return fallback[0]

    raise FileNotFoundError(
        f"Could not resolve FSA for source_run_dir={source_run_dir}, assay={assay}, well={well}, run_code={run_code}"
    )


def load_fsa_full_analysis(fsa_path: Path, preferred_analysis: str = "clonality"):
    from core.analysis import compute_ladder_qc_metrics

    if not fsa_path.exists():
        raise FileNotFoundError(f"Missing FSA file: {fsa_path}")
    metadata = detect_fsa_for_ladder(fsa_path, preferred_analysis=preferred_analysis)
    if not metadata:
        raise ValueError(f"Could not classify {fsa_path.name} for ladder review.")
    fsa, metadata = load_adjustable_fsa(
        fsa_path,
        preferred_analysis=preferred_analysis,
        metadata=metadata,
    )
    metrics = compute_ladder_qc_metrics(fsa)
    return fsa, metadata, metrics


def build_candidate_rows_from_fsa(task: dict[str, Any], fsa) -> pd.DataFrame:
    from core.analysis import get_ladder_candidates

    candidate_table = get_ladder_candidates(fsa).copy()
    selected_map: dict[int, float] = {}
    best_size_standard = getattr(fsa, "best_size_standard", None)
    ladder_steps = getattr(fsa, "ladder_steps", None)
    best_size_standard = np.asarray([] if best_size_standard is None else best_size_standard, dtype=float)
    ladder_steps = np.asarray([] if ladder_steps is None else ladder_steps, dtype=float)
    for time_value, step_bp in zip(best_size_standard.tolist(), ladder_steps.tolist()):
        selected_map[int(round(float(time_value)))] = float(step_bp)

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    recovered_rows: list[dict[str, Any]] = []
    if best_size_standard.size > 0:
        existing_times = (
            candidate_table["time"].astype(float).to_numpy()
            if not candidate_table.empty and "time" in candidate_table.columns
            else np.asarray([], dtype=float)
        )
        next_index = (
            int(pd.to_numeric(candidate_table["index"], errors="coerce").max()) + 1
            if not candidate_table.empty and "index" in candidate_table.columns
            else 0
        )
        for time_value in best_size_standard.tolist():
            time_float = float(time_value)
            if existing_times.size and np.any(np.abs(existing_times - time_float) <= 1.5):
                continue
            idx = int(round(time_float))
            intensity = float(trace[idx]) if 0 <= idx < trace.size else float("nan")
            recovered_rows.append(
                {
                    "index": next_index,
                    "time": time_float,
                    "intensity": intensity,
                    "source": "fit_anchor",
                }
            )
            next_index += 1

    if recovered_rows:
        candidate_table = pd.concat([candidate_table, pd.DataFrame(recovered_rows)], ignore_index=True)
    if candidate_table.empty:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for row in candidate_table.to_dict("records"):
        candidate_index = int(row["index"])
        time_value = float(row["time"])
        rounded_time = int(round(time_value))
        rows.append(
            {
                "month": str(task.get("month", "")),
                "source_run_dir": str(task.get("source_run_dir", "")),
                "assay": str(task.get("assay", "")),
                "identity_key": str(task.get("identity_key", "")),
                "run_date": str(task.get("run_date", "")),
                "run_code": str(task.get("run_code", "")),
                "well": str(task.get("well", "")),
                "ladder": str(task.get("ladder", "")),
                "artifact_row_key": (
                    f"ladder:{task.get('identity_key', task.get('join_key', 'unknown'))}:candidate:{candidate_index}"
                ),
                "join_key": str(task.get("join_key", "")),
                "ladder_join_key": str(task.get("ladder_join_key", "")),
                "candidate_index": candidate_index,
                "candidate_time": time_value,
                "candidate_intensity": float(row["intensity"]),
                "candidate_source": "recovered",
                "selected_for_fit": rounded_time in selected_map,
                "selected_step_bp": selected_map.get(rounded_time, np.nan),
                "ladder_fit_strategy": str(task.get("ladder_fit_strategy", "")),
                "ladder_r2": task.get("ladder_r2", np.nan),
                "ladder_review_required": bool(task.get("ladder_qc") == "review_required"),
                "human_label": "",
                "human_note": "",
                "control": task.get("control", ""),
                "sample_kind": task.get("sample_kind", ""),
            }
        )
    return pd.DataFrame(rows)


def extract_ladder_trace(fsa, ladder: str) -> tuple[np.ndarray, str]:
    raw_keys = list(fsa.fsa.keys())

    def normalize(value: Any) -> str:
        if isinstance(value, bytes):
            value = value.decode("ascii", errors="ignore")
        return str(value).strip().upper()

    preferred = normalize(getattr(fsa, "size_standard_channel", ""))
    ladder_name = str(ladder or "").upper()
    aliases = {
        "LIZ": ["DATA105", "DATA5", "LIZ"],
        "LIZ500_250": ["DATA105", "DATA5", "LIZ"],
        "ROX": ["DATA4", "ROX"],
        "ROX400HD": ["DATA4", "ROX"],
    }
    search_targets: list[str] = []
    if preferred:
        search_targets.append(preferred)
    search_targets.extend(aliases.get(ladder_name, []))
    if ladder_name not in search_targets:
        search_targets.append(ladder_name)

    normalized_keys = [(key, normalize(key)) for key in raw_keys]
    for target in search_targets:
        for key, normalized in normalized_keys:
            if normalized == normalize(target):
                return np.asarray(fsa.fsa[key], dtype=float), normalized

    # Final fallback: if a LIZ case somehow lacks the explicit preferred channel, use the
    # high-number ABI ladder trace that the desktop tools expect.
    if ladder_name.startswith("LIZ"):
        for key, normalized in normalized_keys:
            if normalized in {"DATA105", "DATA5"}:
                return np.asarray(fsa.fsa[key], dtype=float), normalized
    else:
        for key, normalized in normalized_keys:
            if normalized == "DATA4":
                return np.asarray(fsa.fsa[key], dtype=float), normalized

    raise KeyError(f"Could not find ladder trace for {ladder}")


class LadderReviewAnnotator:
    def __init__(self, bundle_dir: Path, data_dir: Path):
        self.bundle_dir = bundle_dir
        self.data_dir = data_dir
        self.case_path = bundle_dir / "ladder_review_cases.csv"
        self.candidate_path = bundle_dir / "ladder_review_candidates.csv"
        self.summary_path = bundle_dir / "ladder_review_summary.json"

        self.case_df = pd.read_csv(self.case_path).fillna("")
        self.candidate_df = pd.read_csv(self.candidate_path).fillna("")
        if "label" not in self.case_df.columns:
            self.case_df["label"] = ""
        if "label_note" not in self.case_df.columns:
            self.case_df["label_note"] = ""
        if "reviewed_at_utc" not in self.case_df.columns:
            self.case_df["reviewed_at_utc"] = ""
        if "human_label" not in self.candidate_df.columns:
            self.candidate_df["human_label"] = ""
        if "human_note" not in self.candidate_df.columns:
            self.candidate_df["human_note"] = ""

        self.current_case: dict[str, Any] | None = None
        self.current_case_candidates = pd.DataFrame()
        self.current_candidate_selection: list[int] = []
        self.filtered_cases: list[dict[str, Any]] = []
        self.current_trace: np.ndarray = np.asarray([], dtype=float)
        self.current_trace_channel = ""
        self.current_plot_error = ""
        self.current_qc_metrics: dict[str, Any] | None = None
        self.current_preview_metrics: dict[str, Any] | None = None
        self.current_preview_error = ""
        self.current_trace_click_x: float | None = None
        self._updating_candidate_selection = False
        self._updating_step_selection = False
        self._dirty = False
        self._case_cache: dict[tuple[str, str, str, str], dict[str, Any]] = {}

        self.title = pn.pane.Markdown("## Clonality Ladder Review")
        self.status = pn.pane.Markdown("")
        self.review_state = pn.pane.Markdown("")
        self.progress = pn.pane.Markdown("")
        self.summary = pn.pane.Markdown("")
        self.meta = pn.pane.Markdown("")
        self.instructions = pn.pane.Markdown(
            "Use the sidebar like the FLT3 app: filter cases, open a file, click peaks in the plot or table, "
            "then save the case label and your kept ladder peaks."
        )
        self.selection_help = pn.pane.Markdown("")
        self.plot = pn.pane.Plotly(height=580, config={"responsive": True, "scrollZoom": True})
        self.plot.param.watch(self._on_plot_click, "click_data")

        assay_options = ["All"] + sorted(self.case_df["assay"].astype(str).dropna().unique().tolist())
        self.assay_filter = pn.widgets.Select(name="Assay", options=assay_options, value="All")
        self.ladder_filter = pn.widgets.RadioButtonGroup(
            name="Ladder",
            options=["All", "ROX", "LIZ"],
            button_type="default",
            value="All",
        )
        initial_queue = "pending"
        case_reviewed = (
            self.case_df["label"].astype(str).str.strip().ne("")
            | self.case_df["label_note"].astype(str).str.strip().ne("")
            | self.case_df["reviewed_at_utc"].astype(str).str.strip().ne("")
        )
        if len(self.case_df) and int((~case_reviewed).sum()) == 0:
            initial_queue = "all"

        self.queue_filter = pn.widgets.RadioButtonGroup(
            name="Queue",
            options={
                "Pending only": "pending",
                "All cases": "all",
                "Reviewed only": "reviewed",
            },
            button_type="primary",
            value=initial_queue,
        )
        self.sort_mode = pn.widgets.RadioButtonGroup(
            name="Sort",
            options=["Default", "Lowest r2", "Run/Well"],
            button_type="default",
            value="Default",
        )
        self.search_input = pn.widgets.TextInput(
            name="Search",
            placeholder="Run, assay, well, identity, run code",
        )
        self.file_select = pn.widgets.Select(name="Review Case", options={})

        self.prev_btn = pn.widgets.Button(name="Previous", button_type="default")
        self.next_btn = pn.widgets.Button(name="Next", button_type="default")
        self.reload_btn = pn.widgets.Button(name="Reload Case", button_type="default")
        self.use_fit_btn = pn.widgets.Button(name="Use Fit Peaks", button_type="primary")
        self.clear_peaks_btn = pn.widgets.Button(name="Clear Peak Picks", button_type="warning")
        self.add_peak_btn = pn.widgets.Button(name="Add Selected Peak", button_type="default")
        self.mark_keep_btn = pn.widgets.Button(name="Mark Selected Keep", button_type="success")
        self.mark_reject_btn = pn.widgets.Button(name="Mark Selected Reject", button_type="danger")
        self.clear_marks_btn = pn.widgets.Button(name="Clear Selected Marks", button_type="default")
        self.save_btn = pn.widgets.Button(name="Save", button_type="primary")
        self.save_next_btn = pn.widgets.Button(name="Save And Next", button_type="success")

        self.case_label = pn.widgets.RadioButtonGroup(
            name="Case Label",
            options={
                "Accept Current Fit": "accept_current_fit",
                "Needs Better Fit": "needs_better_fit",
                "Bad Signal": "bad_signal",
            },
            button_type="primary",
        )
        self.case_note = pn.widgets.TextAreaInput(name="Case Note", height=100)

        self.step_table = pn.widgets.Tabulator(
            pd.DataFrame(columns=["step_bp", "candidate_idx", "time", "intensity", "human_label"]),
            selectable=1,
            pagination=None,
            disabled=True,
            show_index=False,
            height=280,
            sizing_mode="stretch_width",
        )
        self.candidate_table = pn.widgets.Tabulator(
            pd.DataFrame(columns=["candidate_idx", "time", "intensity", "source", "fit_used", "step_bp", "human_label"]),
            selectable=True,
            pagination=None,
            disabled=False,
            show_index=False,
            height=320,
            sizing_mode="stretch_width",
        )

        self.assay_filter.param.watch(self._on_filter_change, "value")
        self.ladder_filter.param.watch(self._on_filter_change, "value")
        self.queue_filter.param.watch(self._on_filter_change, "value")
        self.sort_mode.param.watch(self._on_filter_change, "value")
        self.search_input.param.watch(self._on_filter_change, "value")
        self.file_select.param.watch(self._on_file_change, "value")
        self.candidate_table.param.watch(self._on_candidate_table_select, "selection")
        self.step_table.param.watch(self._on_step_table_select, "selection")
        self.case_label.param.watch(self._on_case_metadata_change, "value")
        self.case_note.param.watch(self._on_case_metadata_change, "value")

        self.prev_btn.on_click(lambda _event: self._move_case(-1))
        self.next_btn.on_click(lambda _event: self._move_case(1))
        self.reload_btn.on_click(lambda _event: self._reload_current_case())
        self.use_fit_btn.on_click(lambda _event: self._use_pipeline_fit_peaks())
        self.clear_peaks_btn.on_click(lambda _event: self._clear_peak_selection())
        self.add_peak_btn.on_click(lambda _event: self._add_selected_trace_peak())
        self.mark_keep_btn.on_click(lambda _event: self._mark_selected_peaks("keep_peak"))
        self.mark_reject_btn.on_click(lambda _event: self._mark_selected_peaks("reject_peak"))
        self.clear_marks_btn.on_click(lambda _event: self._clear_selected_peak_labels())
        self.save_btn.on_click(lambda _event: self._save_current_case(move_next=False))
        self.save_next_btn.on_click(lambda _event: self._save_current_case(move_next=True))

        sidebar = pn.Column(
            self.assay_filter,
            self.ladder_filter,
            self.queue_filter,
            self.sort_mode,
            self.search_input,
            self.file_select,
            pn.Row(self.prev_btn, self.next_btn),
            pn.Row(self.reload_btn, self.use_fit_btn),
            pn.Row(self.add_peak_btn, self.clear_peaks_btn),
            self.case_label,
            self.case_note,
            pn.Row(self.mark_keep_btn, self.mark_reject_btn),
            pn.Row(self.clear_marks_btn),
            self.save_btn,
            self.save_next_btn,
            width=400,
        )

        review_tables = pn.Row(
            pn.Column("### Fitted Ladder Steps", self.step_table, sizing_mode="stretch_width"),
            pn.Column("### Candidate Peaks", self.candidate_table, sizing_mode="stretch_width"),
            sizing_mode="stretch_width",
        )

        main = pn.Column(
            self.title,
            self.status,
            self.review_state,
            self.progress,
            self.summary,
            self.instructions,
            self.meta,
            self.selection_help,
            self.plot,
            review_tables,
            sizing_mode="stretch_width",
        )
        self.layout = pn.Row(sidebar, main, sizing_mode="stretch_width")

        self._update_filtered_cases()
        if self.filtered_cases:
            self.file_select.value = str(self.filtered_cases[0]["_case_index"])
            self._load_case(self.filtered_cases[0])
        else:
            self._set_status("No ladder review cases found for the current filters.", level="warning")
            self._clear_views()

    def _set_status(self, text: str, *, level: str = "info") -> None:
        color = {"info": "#475569", "success": "#15803d", "warning": "#b45309", "error": "#b91c1c"}.get(
            level, "#475569"
        )
        self.status.object = f"<div style='color:{color};font-weight:600'>{text}</div>"

    @staticmethod
    def _case_cache_key(task: dict[str, Any]) -> tuple[str, str, str, str]:
        return (
            str(task.get("source_run_dir", "")),
            str(task.get("assay", "")),
            str(task.get("well", "")),
            str(task.get("ladder", "")),
        )

    @staticmethod
    def _coerce_scalar_value(value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, (list, tuple, set)):
            if not value:
                return ""
            return str(next(iter(value)))
        return str(value)

    @staticmethod
    def _truthy(value: Any) -> bool:
        return str(value).strip().lower() in {"1", "true", "yes", "y"}

    @staticmethod
    def _safe_float(value: Any) -> float:
        try:
            return float(value)
        except Exception:
            return float("nan")

    def _set_dirty(self, value: bool = True) -> None:
        self._dirty = value
        self._update_review_state()

    def _update_review_state(self) -> None:
        if self.current_case is None:
            self.review_state.object = ""
            return
        reviewed_at = str(self.current_case.get("reviewed_at_utc", "")).strip()
        case_label = str(self.current_case.get("label", "")).strip()
        has_review = bool(reviewed_at or case_label)
        state_bits = []
        state_bits.append("**Unsaved changes:** yes" if self._dirty else "**Unsaved changes:** no")
        if has_review:
            stamp = reviewed_at or "saved"
            state_bits.append(f"**Saved review:** `{stamp}`")
        else:
            state_bits.append("**Saved review:** not yet")
        self.review_state.object = "  \n".join(state_bits)

    def _all_cases(self) -> pd.DataFrame:
        df = self.case_df.reset_index(drop=False).rename(columns={"index": "_case_index"})
        reviewed_case_keys = set()
        if not self.candidate_df.empty:
            reviewed_candidates = self.candidate_df[
                self.candidate_df["human_label"].astype(str).str.strip().ne("")
                | self.candidate_df["human_note"].astype(str).str.strip().ne("")
            ]
            reviewed_case_keys = set(
                zip(
                    reviewed_candidates["source_run_dir"].astype(str),
                    reviewed_candidates["assay"].astype(str),
                    reviewed_candidates["well"].astype(str),
                )
            )
        df["_reviewed"] = [
            bool(str(row.get("label", "")).strip())
            or bool(str(row.get("label_note", "")).strip())
            or bool(str(row.get("reviewed_at_utc", "")).strip())
            or (str(row.get("source_run_dir", "")), str(row.get("assay", "")), str(row.get("well", ""))) in reviewed_case_keys
            for row in df.to_dict("records")
        ]
        return df

    @staticmethod
    def _is_legacy_candidate_review_case(task: dict[str, Any]) -> bool:
        return (
            not bool(str(task.get("reviewed_at_utc", "")).strip())
            and not bool(str(task.get("label", "")).strip())
        )

    def _task_candidates(self, task: dict[str, Any]) -> pd.DataFrame:
        mask = (
            (self.candidate_df["source_run_dir"] == task["source_run_dir"])
            & (self.candidate_df["assay"] == task["assay"])
            & (self.candidate_df["well"] == task["well"])
        )
        rows = self.candidate_df.loc[mask].copy()
        if rows.empty:
            return pd.DataFrame(
                columns=[
                    "candidate_index",
                    "candidate_time",
                    "candidate_intensity",
                    "candidate_source",
                    "selected_for_fit",
                    "selected_step_bp",
                    "human_label",
                ]
            )
        rows["candidate_index"] = rows["candidate_index"].astype(int)
        rows["candidate_time"] = rows["candidate_time"].astype(float)
        rows["candidate_intensity"] = rows["candidate_intensity"].astype(float)
        rows["selected_step_bp"] = pd.to_numeric(rows["selected_step_bp"], errors="coerce")
        rows["human_label"] = rows["human_label"].astype(str)
        rows["human_note"] = rows["human_note"].astype(str)
        if self._is_legacy_candidate_review_case(task):
            rows.loc[rows["human_label"].str.strip().str.lower() == "reject_peak", "human_label"] = ""
        rows = rows.sort_values(["candidate_index"]).reset_index(drop=True)
        return rows

    @staticmethod
    def _empty_candidate_frame() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "month",
                "source_run_dir",
                "assay",
                "identity_key",
                "run_date",
                "run_code",
                "well",
                "ladder",
                "artifact_row_key",
                "join_key",
                "ladder_join_key",
                "candidate_index",
                "candidate_time",
                "candidate_intensity",
                "candidate_source",
                "selected_for_fit",
                "selected_step_bp",
                "ladder_fit_strategy",
                "ladder_r2",
                "ladder_review_required",
                "human_label",
                "human_note",
                "control",
                "sample_kind",
            ]
        )

    def _normalize_candidate_rows(self, rows: pd.DataFrame, task: dict[str, Any]) -> pd.DataFrame:
        if rows.empty:
            return self._empty_candidate_frame()
        normalized = rows.copy()
        defaults = {
            "month": str(task.get("month", "")),
            "source_run_dir": str(task.get("source_run_dir", "")),
            "assay": str(task.get("assay", "")),
            "identity_key": str(task.get("identity_key", "")),
            "run_date": str(task.get("run_date", "")),
            "run_code": str(task.get("run_code", "")),
            "well": str(task.get("well", "")),
            "ladder": str(task.get("ladder", "")),
            "artifact_row_key": "",
            "join_key": str(task.get("join_key", "")),
            "ladder_join_key": str(task.get("ladder_join_key", "")),
            "candidate_source": "auto",
            "selected_for_fit": False,
            "selected_step_bp": np.nan,
            "ladder_fit_strategy": str(task.get("ladder_fit_strategy", "")),
            "ladder_r2": task.get("ladder_r2", np.nan),
            "ladder_review_required": bool(task.get("ladder_qc") == "review_required"),
            "human_label": "",
            "human_note": "",
            "control": task.get("control", ""),
            "sample_kind": task.get("sample_kind", ""),
        }
        for column, default in defaults.items():
            if column not in normalized.columns:
                normalized[column] = default
        normalized["candidate_index"] = pd.to_numeric(normalized["candidate_index"], errors="coerce").astype(int)
        normalized["candidate_time"] = pd.to_numeric(normalized["candidate_time"], errors="coerce").astype(float)
        normalized["candidate_intensity"] = pd.to_numeric(normalized["candidate_intensity"], errors="coerce").astype(float)
        normalized["selected_step_bp"] = pd.to_numeric(normalized["selected_step_bp"], errors="coerce")
        normalized["human_label"] = normalized["human_label"].astype(str)
        normalized["human_note"] = normalized["human_note"].astype(str)
        normalized = normalized.sort_values(["candidate_index", "candidate_time"]).reset_index(drop=True)
        return normalized

    def _candidate_frame_from_fsa(self, task: dict[str, Any], fsa) -> pd.DataFrame:
        try:
            rows = build_candidate_rows_from_fsa(task, fsa)
        except Exception:
            logging.exception("Could not rebuild candidate rows from live FSA")
            return self._empty_candidate_frame()
        return self._normalize_candidate_rows(rows, task)

    def _merge_candidate_annotations(
        self,
        live_rows: pd.DataFrame,
        saved_rows: pd.DataFrame,
    ) -> pd.DataFrame:
        """Keep fresh live candidates, but preserve manual labels/notes when times align."""
        live = self._normalize_candidate_rows(live_rows, self.current_case or {}) if not live_rows.empty else self._empty_candidate_frame()
        if saved_rows.empty or live.empty:
            return live

        saved = self._normalize_candidate_rows(saved_rows, self.current_case or {})
        if saved.empty:
            return live

        for live_idx, live_row in live.iterrows():
            distances = (saved["candidate_time"] - float(live_row["candidate_time"])).abs()
            if distances.empty:
                continue
            best_idx = int(distances.idxmin())
            best_distance = float(distances.loc[best_idx])
            if best_distance > 3.0:
                continue
            saved_row = saved.loc[best_idx]
            label = str(saved_row.get("human_label", "")).strip()
            note = str(saved_row.get("human_note", "")).strip()
            if label:
                live.at[live_idx, "human_label"] = label
            if note:
                live.at[live_idx, "human_note"] = note

        return live

    def _add_manual_candidate_at_time(self, clicked_time: float) -> int | None:
        if self.current_case is None or self.current_trace.size == 0:
            return None
        center = int(round(clicked_time))
        if center < 0 or center >= len(self.current_trace):
            return None
        left = max(0, center - 12)
        right = min(len(self.current_trace), center + 13)
        window = self.current_trace[left:right]
        if window.size == 0:
            return None
        local_offset = int(np.argmax(window))
        peak_time = float(left + local_offset)
        peak_intensity = float(self.current_trace[int(peak_time)])

        if not self.current_case_candidates.empty:
            distances = (self.current_case_candidates["candidate_time"] - peak_time).abs()
            if not distances.empty and float(distances.min()) <= 3.0:
                return int(distances.idxmin())
            next_candidate_index = int(self.current_case_candidates["candidate_index"].max()) + 1
        else:
            next_candidate_index = 0

        task = self.current_case
        new_row = {
            "month": str(task.get("month", "")),
            "source_run_dir": str(task.get("source_run_dir", "")),
            "assay": str(task.get("assay", "")),
            "identity_key": str(task.get("identity_key", "")),
            "run_date": str(task.get("run_date", "")),
            "run_code": str(task.get("run_code", "")),
            "well": str(task.get("well", "")),
            "ladder": str(task.get("ladder", "")),
            "artifact_row_key": (
                f"ladder:{task.get('identity_key', task.get('join_key', 'unknown'))}:candidate:{next_candidate_index}"
            ),
            "join_key": str(task.get("join_key", "")),
            "ladder_join_key": str(task.get("ladder_join_key", "")),
            "candidate_index": next_candidate_index,
            "candidate_time": peak_time,
            "candidate_intensity": peak_intensity,
            "candidate_source": "manual",
            "selected_for_fit": False,
            "selected_step_bp": np.nan,
            "ladder_fit_strategy": str(task.get("ladder_fit_strategy", "")),
            "ladder_r2": task.get("ladder_r2", np.nan),
            "ladder_review_required": bool(task.get("ladder_qc") == "review_required"),
            "human_label": "",
            "human_note": "",
            "control": task.get("control", ""),
            "sample_kind": task.get("sample_kind", ""),
        }
        self.current_case_candidates = self._normalize_candidate_rows(
            pd.concat([self.current_case_candidates, pd.DataFrame([new_row])], ignore_index=True),
            task,
        )
        match = self.current_case_candidates.index[
            self.current_case_candidates["candidate_index"] == next_candidate_index
        ]
        if len(match) == 0:
            return None
        self._set_dirty(True)
        return int(match[0])

    def _recompute_preview_metrics(self) -> None:
        self.current_preview_metrics = None
        self.current_preview_error = ""
        if self.current_case is None:
            return
        if not self.current_candidate_selection:
            return
        cache_key = self._case_cache_key(self.current_case)
        cached = self._case_cache.get(cache_key)
        if cached is None or "fsa" not in cached:
            return

        from core.analysis import compute_ladder_qc_metrics, fit_size_standard_to_ladder

        selected_rows = self.current_case_candidates.iloc[self.current_candidate_selection].copy()
        selected_rows = selected_rows.sort_values("candidate_time")
        candidate_times = selected_rows["candidate_time"].astype(float).to_numpy()
        fsa = cached["fsa"]
        ladder_steps_value = getattr(fsa, "ladder_steps", None)
        ladder_steps = np.asarray([] if ladder_steps_value is None else ladder_steps_value, dtype=float)
        expected_count = int(len(ladder_steps))
        if expected_count == 0:
            self.current_preview_error = "No ladder steps available for preview."
            return
        if len(candidate_times) != expected_count:
            self.current_preview_error = f"Preview fit needs {expected_count} selected peaks, currently {len(candidate_times)}."
            return
        if np.any(np.diff(candidate_times) <= 0):
            self.current_preview_error = "Preview fit requires strictly increasing selected peak times."
            return

        try:
            trial = deepcopy(fsa)
            trial.best_size_standard = np.asarray(candidate_times, dtype=float)
            trial.expected_ladder_steps = ladder_steps.copy()
            trial.ladder_steps = ladder_steps.copy()
            trial = fit_size_standard_to_ladder(trial)
            self.current_preview_metrics = compute_ladder_qc_metrics(trial)
        except Exception as exc:
            self.current_preview_error = str(exc)

    def _case_label_text(self, task: dict[str, Any]) -> str:
        review_state = str(task.get("label") or ("reviewed" if bool(task.get("_reviewed")) else "pending"))
        run_date = str(task.get("run_date") or "")
        r2 = self._safe_float(task.get("ladder_r2"))
        r2_text = f"{r2:.6f}" if np.isfinite(r2) else "nan"
        return (
            f"{task.get('assay')} {task.get('well')} | {run_date} | {task.get('ladder')} | "
            f"{task.get('ladder_qc')} | r2 {r2_text} | {review_state} | {task.get('source_run_dir')}"
        )

    def _sort_cases(self, df: pd.DataFrame) -> pd.DataFrame:
        if self.sort_mode.value == "Lowest r2":
            sorted_df = df.copy()
            sorted_df["_ladder_r2_num"] = pd.to_numeric(sorted_df["ladder_r2"], errors="coerce")
            return sorted_df.sort_values(
                ["_ladder_r2_num", "source_run_dir", "assay", "well"],
                ascending=[True, True, True, True],
                na_position="last",
            ).drop(columns=["_ladder_r2_num"])
        if self.sort_mode.value == "Run/Well":
            return df.sort_values(["source_run_dir", "assay", "well", "_case_index"])
        return df.sort_values(["_case_index"])

    def _update_progress_text(self) -> None:
        all_cases = self._all_cases()
        total = len(all_cases)
        pending = int((~all_cases["_reviewed"].astype(bool)).sum())
        reviewed = total - pending
        filtered = len(self.filtered_cases)
        self.progress.object = (
            f"**Cases:** {filtered} filtered / {total} total  \n"
            f"**Pending:** {pending}  \n"
            f"**Reviewed:** {reviewed}"
        )

    def _update_filtered_cases(self) -> None:
        df = self._all_cases()
        queue_value = self.queue_filter.value
        if queue_value == "pending":
            df = df[~df["_reviewed"].astype(bool)]
        elif queue_value == "reviewed":
            df = df[df["_reviewed"].astype(bool)]

        if self.assay_filter.value != "All":
            df = df[df["assay"] == self.assay_filter.value]
        if self.ladder_filter.value != "All":
            df = df[df["ladder"] == self.ladder_filter.value]

        search_text = str(self.search_input.value or "").strip().lower()
        if search_text:
            haystacks = (
                df["source_run_dir"].astype(str)
                + " "
                + df["assay"].astype(str)
                + " "
                + df["well"].astype(str)
                + " "
                + df["run_code"].astype(str)
                + " "
                + df["identity_key"].astype(str)
                + " "
                + df["label"].astype(str)
            ).str.lower()
            df = df[haystacks.str.contains(search_text, na=False)]

        df = self._sort_cases(df)
        self.filtered_cases = df.to_dict("records")
        self._update_progress_text()

        current_value = self.file_select.value
        options = {self._case_label_text(task): str(task["_case_index"]) for task in self.filtered_cases}
        self.file_select.options = options

        if current_value in options.values():
            self.file_select.value = current_value
        elif self.filtered_cases:
            self.file_select.value = str(self.filtered_cases[0]["_case_index"])
        else:
            self.file_select.value = None

    def _find_case_by_index(self, case_index: int) -> dict[str, Any] | None:
        matches = self._all_cases()
        row = matches.loc[matches["_case_index"] == case_index]
        if row.empty:
            return None
        return row.iloc[0].to_dict()

    def _on_filter_change(self, _event) -> None:
        current_value = self.file_select.value
        self._update_filtered_cases()
        if current_value and any(str(task["_case_index"]) == str(current_value) for task in self.filtered_cases):
            return
        if self.filtered_cases:
            self._load_case(self.filtered_cases[0])
        else:
            self.current_case = None
            self._clear_views()
            self._set_status("No ladder review cases match the current filters.", level="warning")

    def _on_file_change(self, event) -> None:
        if not event.new:
            return
        try:
            case_index = int(str(event.new))
        except ValueError:
            return
        task = next((task for task in self.filtered_cases if int(task["_case_index"]) == case_index), None)
        if task is not None:
            self._load_case(task)

    def _move_case(self, delta: int) -> None:
        if not self.filtered_cases or self.current_case is None:
            return
        indices = [int(task["_case_index"]) for task in self.filtered_cases]
        try:
            current_pos = indices.index(int(self.current_case["_case_index"]))
        except ValueError:
            current_pos = 0
        new_pos = max(0, min(len(indices) - 1, current_pos + delta))
        self.file_select.value = str(indices[new_pos])

    def _reload_current_case(self) -> None:
        if self.current_case is None:
            return
        fresh = self._find_case_by_index(int(self.current_case["_case_index"]))
        if fresh is not None:
            self._load_case(fresh)

    def _clear_views(self) -> None:
        self.summary.object = ""
        self.meta.object = ""
        self.selection_help.object = ""
        self.review_state.object = ""
        self.plot.object = go.Figure()
        self.step_table.value = pd.DataFrame(columns=["step_bp", "candidate_idx", "time", "intensity", "human_label"])
        self.candidate_table.value = pd.DataFrame(
            columns=["candidate_idx", "time", "intensity", "source", "fit_used", "step_bp", "human_label"]
        )
        self.current_case_candidates = pd.DataFrame()
        self.current_candidate_selection = []
        self.current_qc_metrics = None
        self.current_preview_metrics = None
        self.current_preview_error = ""
        self.current_trace_click_x = None
        self._dirty = False

    def _load_case(self, task: dict[str, Any]) -> None:
        self.current_case = task
        self.case_label.value = task.get("label") or None
        saved_case_candidates = self._task_candidates(task)
        self.current_case_candidates = saved_case_candidates.copy()
        note_value = str(task.get("label_note") or "")
        if not note_value and not self.current_case_candidates.empty:
            candidate_notes = [
                str(value).strip()
                for value in self.current_case_candidates["human_note"].tolist()
                if str(value).strip()
            ]
            if candidate_notes:
                note_value = max(set(candidate_notes), key=candidate_notes.count)
        self.case_note.value = note_value
        selected_rows: list[int] = []
        if not self.current_case_candidates.empty:
            for row_idx, row in self.current_case_candidates.iterrows():
                if str(row.get("human_label", "")).strip().lower() in {"keep_peak", "relevant_peak", "reject_peak"}:
                    selected_rows.append(int(row_idx))
        self.current_candidate_selection = sorted(selected_rows)
        self.current_preview_metrics = None
        self.current_preview_error = ""

        try:
            cache_key = self._case_cache_key(task)
            cached = self._case_cache.get(cache_key)
            if cached is None:
                fsa_path = resolve_fsa_path(
                    self.data_dir,
                    str(task["source_run_dir"]),
                    str(task["assay"]),
                    str(task["well"]),
                    str(task["run_code"]),
                )
                fsa, _metadata, metrics = load_fsa_full_analysis(fsa_path, preferred_analysis="clonality")
                trace_raw, trace_channel = extract_ladder_trace(fsa, str(task["ladder"]))
                from core.analysis import estimate_running_baseline

                try:
                    baseline = estimate_running_baseline(trace_raw, bin_size=200, quantile=0.10)
                    trace = np.clip(trace_raw - baseline, a_min=0, a_max=None)
                except Exception:
                    trace = trace_raw.astype(float)
                cached = {
                    "trace": np.asarray(trace, dtype=float),
                    "trace_channel": trace_channel,
                    "plot_error": "",
                    "metrics": metrics,
                    "fsa": fsa,
                }
                self._case_cache[cache_key] = cached
            self.current_trace = np.asarray(cached["trace"], dtype=float)
            self.current_trace_channel = str(cached["trace_channel"])
            self.current_plot_error = str(cached["plot_error"])
            self.current_qc_metrics = cached["metrics"]
            self.current_trace_click_x = None
            live_candidates = self._candidate_frame_from_fsa(task, cached["fsa"])
            self.current_case_candidates = self._merge_candidate_annotations(live_candidates, saved_case_candidates)
            self._recompute_preview_metrics()
            self._set_status(f"Loaded {task['assay']} {task['well']} for review.", level="success")
        except Exception as exc:
            logging.exception("Failed to load ladder review task")
            self.current_trace = np.asarray([], dtype=float)
            self.current_trace_channel = ""
            self.current_plot_error = str(exc)
            self.current_qc_metrics = None
            self._set_status(f"Could not load case: {exc}", level="error")

        self._set_dirty(False)
        self._refresh_views()

    def _refresh_views(self) -> None:
        if self.current_case is None:
            self._clear_views()
            return
        task = self.current_case
        selected_count = len(self.current_candidate_selection)
        fit_selected = int(
            self.current_case_candidates["selected_for_fit"].astype(str).str.strip().str.lower().eq("true").sum()
        ) if not self.current_case_candidates.empty else 0
        metric_r2 = self._safe_float((self.current_qc_metrics or {}).get("r2"))
        workbook_r2 = self._safe_float(task.get("ladder_r2"))
        mean_abs_error = self._safe_float((self.current_qc_metrics or {}).get("mean_abs_error_bp"))
        max_abs_error = self._safe_float((self.current_qc_metrics or {}).get("max_abs_error_bp"))
        preview_r2 = self._safe_float((self.current_preview_metrics or {}).get("r2"))
        preview_mean_abs_error = self._safe_float((self.current_preview_metrics or {}).get("mean_abs_error_bp"))
        preview_max_abs_error = self._safe_float((self.current_preview_metrics or {}).get("max_abs_error_bp"))
        r2_text = f"{metric_r2:.6f}" if np.isfinite(metric_r2) else (f"{workbook_r2:.6f}" if np.isfinite(workbook_r2) else "nan")
        mean_text = f"{mean_abs_error:.2f} bp" if np.isfinite(mean_abs_error) else "n/a"
        max_text = f"{max_abs_error:.2f} bp" if np.isfinite(max_abs_error) else "n/a"
        if self.current_preview_metrics is not None and np.isfinite(preview_r2):
            preview_text = (
                f"`preview R2 {preview_r2:.6f}` | "
                f"`preview mean residual {preview_mean_abs_error:.2f} bp` | "
                f"`preview max residual {preview_max_abs_error:.2f} bp`"
            )
        elif self.current_preview_error:
            preview_text = self.current_preview_error
        else:
            preview_text = "Select peaks to preview updated residuals."
        self.summary.object = (
            f"**Current case:** `{task['source_run_dir']}`  \n"
            f"**Assay / Well:** `{task['assay']}` / `{task['well']}`  \n"
            f"**Pipeline fit:** `{task['ladder_fit_strategy']}` | `{task['ladder_fitted_step_count']}` / "
            f"`{task['ladder_expected_step_count']}` steps | `R2 {r2_text}` | "
            f"`mean residual {mean_text}` | `max residual {max_text}`  \n"
            f"**Peak picks:** {selected_count} selected now | {fit_selected} used by current fit  \n"
            f"**Live preview:** {preview_text}"
        )
        self.meta.object = (
            f"**Run date:** `{task['run_date']}`  \n"
            f"**Run code:** `{task['run_code']}`  \n"
            f"**Ladder:** `{task['ladder']}`  \n"
            f"**QC:** `{task['ladder_qc']}`  \n"
            f"**Identity:** `{task['identity_key']}`  \n"
            f"**Residuals:** mean `{mean_text}` | max `{max_text}`  \n"
            f"**Preview status:** {preview_text}"
        )
        self._refresh_step_table()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()

    def _refresh_step_table(self) -> None:
        if self.current_case_candidates.empty:
            df = pd.DataFrame(columns=["step_bp", "candidate_idx", "time", "intensity", "human_label"])
        else:
            step_rows = self.current_case_candidates[
                self.current_case_candidates["selected_for_fit"].astype(str).str.strip().str.lower() == "true"
            ].copy()
            step_rows = step_rows.sort_values(["selected_step_bp", "candidate_time"], na_position="last")
            df = pd.DataFrame(
                {
                    "step_bp": step_rows["selected_step_bp"].round(1),
                    "candidate_idx": step_rows["candidate_index"].astype(int),
                    "time": step_rows["candidate_time"].round(1),
                    "intensity": step_rows["candidate_intensity"].round(0).astype(int),
                    "human_label": step_rows["human_label"].replace("", "-"),
                }
            )
        self.step_table.value = df
        self.step_table.disabled = df.empty

    def _refresh_candidate_table(self) -> None:
        if self.current_case_candidates.empty:
            df = pd.DataFrame(columns=["candidate_idx", "time", "intensity", "source", "fit_used", "step_bp", "human_label"])
        else:
            df = pd.DataFrame(
                {
                    "candidate_idx": self.current_case_candidates["candidate_index"].astype(int),
                    "time": self.current_case_candidates["candidate_time"].round(1),
                    "intensity": self.current_case_candidates["candidate_intensity"].round(0).astype(int),
                    "source": self.current_case_candidates["candidate_source"],
                    "fit_used": self.current_case_candidates["selected_for_fit"],
                    "step_bp": self.current_case_candidates["selected_step_bp"].round(1),
                    "human_label": self.current_case_candidates["human_label"].replace("", "-"),
                }
            )
        self.candidate_table.value = df
        self._updating_candidate_selection = True
        self.candidate_table.selection = list(self.current_candidate_selection)
        self._updating_candidate_selection = False

    def _plot_peak_intensity(self, time_value: float, radius: int = 4) -> float:
        if self.current_trace.size == 0:
            return float("nan")
        center = int(round(time_value))
        if center < 0 or center >= self.current_trace.size:
            return float("nan")
        left = max(0, center - radius)
        right = min(self.current_trace.size, center + radius + 1)
        window = self.current_trace[left:right]
        if window.size == 0:
            return float(self.current_trace[center])
        return float(np.max(window))

    def _refresh_plot(self) -> None:
        if self.current_case is None:
            self.plot.object = go.Figure()
            return
        if self.current_plot_error:
            self.plot.object = go.Figure()
            return

        task = self.current_case
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=np.arange(len(self.current_trace)),
                y=self.current_trace,
                mode="lines",
                line=dict(color="#1f77b4", width=2),
                fill="tozeroy",
                fillcolor="rgba(31,119,180,0.08)",
                name=f"{task['ladder']} trace",
            )
        )

        for row_idx, row in self.current_case_candidates.iterrows():
            candidate_idx = int(row["candidate_index"])
            time_value = float(row["candidate_time"])
            source_intensity = float(row["candidate_intensity"])
            intensity_value = self._plot_peak_intensity(time_value)
            if not np.isfinite(intensity_value):
                intensity_value = source_intensity
            fit_used = self._truthy(row.get("selected_for_fit"))
            human_label = str(row.get("human_label", "")).strip().lower()
            is_selected = row_idx in self.current_candidate_selection

            color = "#ef476f"
            symbol = "circle"
            if fit_used:
                color = "#ffd166"
                symbol = "diamond"
            if human_label in {"keep_peak", "relevant_peak"}:
                color = "#06d6a0"
                symbol = "diamond"
            elif human_label == "reject_peak":
                color = "#ef476f"
                symbol = "x"
            size = 11
            line_width = 0
            if is_selected:
                size = 18
                line_width = 2

            fig.add_trace(
                go.Scatter(
                    x=[time_value],
                    y=[intensity_value],
                    mode="markers+text",
                    marker=dict(
                        color=color,
                        size=size,
                        symbol=symbol,
                        line=dict(color="#0f172a", width=line_width),
                    ),
                    text=[str(candidate_idx)],
                    textposition="top center",
                    name=f"candidate {candidate_idx}",
                    customdata=[[int(row_idx), candidate_idx]],
                    hovertemplate=(
                        f"candidate #{candidate_idx}<br>time={time_value:.1f}<br>plot_intensity={intensity_value:.0f}"
                        f"<br>candidate_intensity={source_intensity:.0f}"
                        f"<br>fit_used={fit_used}<br>human_label={human_label or '-'}<extra></extra>"
                    ),
                )
            )

        fig.update_layout(
            title=(
                f"{task['source_run_dir']} | {task['assay']} {task['well']} | "
                f"strategy={task['ladder_fit_strategy']} | fitted={task['ladder_fitted_step_count']}/"
                f"{task['ladder_expected_step_count']}"
            ),
            xaxis_title=f"Ladder Time / Scan Index ({self.current_trace_channel})",
            yaxis_title="Signal",
            template="plotly_white",
            height=580,
            showlegend=False,
            margin=dict(l=40, r=20, t=70, b=40),
        )
        self.plot.object = fig

    def _refresh_selection_help(self) -> None:
        if self.current_case is None:
            self.selection_help.object = ""
            return
        if self.current_plot_error:
            self.selection_help.object = f"**Plot unavailable:** {self.current_plot_error}"
            return
        selected_bits: list[str] = []
        for row_idx in self.current_candidate_selection:
            if row_idx < 0 or row_idx >= len(self.current_case_candidates):
                continue
            row = self.current_case_candidates.iloc[row_idx]
            selected_bits.append(
                f"`#{int(row['candidate_index'])}` @ `{float(row['candidate_time']):.1f}`"
            )
        selected_text = ", ".join(selected_bits) if selected_bits else "none"
        click_text = (
            f"**Last clicked trace position:** `{self.current_trace_click_x:.1f}`  \n"
            if self.current_trace_click_x is not None
            else ""
        )
        self.selection_help.object = (
            f"**Selected peaks:** {selected_text}  \n"
            f"{click_text}"
            "Click a plotted candidate or select rows in the candidate table. "
            "Click bare trace to set an add-peak position. Only explicit keep/reject labels are saved."
        )

    def _on_case_metadata_change(self, _event) -> None:
        if self.current_case is None:
            return
        self._set_dirty(True)

    def _on_candidate_table_select(self, event) -> None:
        if self._updating_candidate_selection:
            return
        self.current_candidate_selection = sorted(int(idx) for idx in (event.new or []))
        self._recompute_preview_metrics()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()

    def _on_step_table_select(self, event) -> None:
        if self._updating_step_selection or self.current_case_candidates.empty or not event.new:
            return
        step_row_idx = int(event.new[0])
        step_df = self.step_table.value
        if step_row_idx >= len(step_df):
            return
        candidate_idx = int(step_df.iloc[step_row_idx]["candidate_idx"])
        match = self.current_case_candidates.index[self.current_case_candidates["candidate_index"] == candidate_idx]
        if len(match) == 0:
            return
        self.current_candidate_selection = [int(match[0])]
        self._recompute_preview_metrics()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()

    def _nearest_candidate_row(self, clicked_time: float) -> int | None:
        if self.current_case_candidates.empty:
            return None
        distances = (self.current_case_candidates["candidate_time"] - clicked_time).abs()
        if distances.empty:
            return None
        nearest_idx = int(distances.idxmin())
        if float(distances.loc[nearest_idx]) > 25.0:
            return None
        return nearest_idx

    def _on_plot_click(self, event) -> None:
        click_data = getattr(event, "new", None)
        if not click_data:
            return
        point = (click_data.get("points") or [{}])[0]
        if point.get("x") is not None:
            self.current_trace_click_x = float(point["x"])
        if point.get("customdata"):
            nearest_row = int(point["customdata"][0])
        elif point.get("x") is not None:
            nearest_row = self._nearest_candidate_row(float(point["x"]))
            if nearest_row is None:
                self._refresh_selection_help()
                self._set_status(
                    f"Selected trace position at {float(point['x']):.1f}. Use 'Add Selected Peak' to create a manual candidate.",
                    level="info",
                )
                return
        else:
            return

        selection = set(self.current_candidate_selection)
        if nearest_row in selection:
            selection.remove(nearest_row)
        else:
            selection.add(nearest_row)
        self.current_candidate_selection = sorted(selection)
        self._recompute_preview_metrics()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()

    def _add_selected_trace_peak(self) -> None:
        if self.current_trace_click_x is None:
            self._set_status("Click the trace where you want to add a peak first.", level="warning")
            return
        added_row = self._add_manual_candidate_at_time(self.current_trace_click_x)
        if added_row is None:
            self._set_status("Could not add a manual peak at that trace position.", level="error")
            return
        self.current_candidate_selection = [added_row]
        self._recompute_preview_metrics()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()
        self._set_status("Added a manual candidate peak from the selected trace position.", level="success")

    def _use_pipeline_fit_peaks(self) -> None:
        if self.current_case_candidates.empty:
            return
        selection = self.current_case_candidates.index[
            self.current_case_candidates["selected_for_fit"].astype(str).str.strip().str.lower() == "true"
        ]
        self.current_candidate_selection = sorted(int(idx) for idx in selection.tolist())
        self._recompute_preview_metrics()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()

    def _clear_peak_selection(self) -> None:
        self.current_candidate_selection = []
        self._recompute_preview_metrics()
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()

    def _mark_selected_peaks(self, label: str) -> None:
        if self.current_case_candidates.empty or not self.current_candidate_selection:
            self._set_status("Select one or more candidate peaks first.", level="warning")
            return
        note_value = self._coerce_scalar_value(self.case_note.value)
        for row_idx in self.current_candidate_selection:
            if row_idx < 0 or row_idx >= len(self.current_case_candidates):
                continue
            self.current_case_candidates.at[row_idx, "human_label"] = label
            self.current_case_candidates.at[row_idx, "human_note"] = note_value
        self._set_dirty(True)
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()
        action = "keep" if label == "keep_peak" else "reject"
        self._set_status(f"Marked {len(self.current_candidate_selection)} selected peaks as {action}.", level="success")

    def _clear_selected_peak_labels(self) -> None:
        if self.current_case_candidates.empty or not self.current_candidate_selection:
            self._set_status("Select one or more candidate peaks first.", level="warning")
            return
        for row_idx in self.current_candidate_selection:
            if row_idx < 0 or row_idx >= len(self.current_case_candidates):
                continue
            self.current_case_candidates.at[row_idx, "human_label"] = ""
            self.current_case_candidates.at[row_idx, "human_note"] = ""
        self._set_dirty(True)
        self._refresh_candidate_table()
        self._refresh_plot()
        self._refresh_selection_help()
        self._set_status(f"Cleared explicit labels for {len(self.current_candidate_selection)} selected peaks.", level="success")

    def _persist_current_task(self) -> str:
        if self.current_case is None:
            return ""
        task = self.current_case
        case_index = int(task["_case_index"])
        timestamp = datetime.now(timezone.utc).isoformat()

        case_label_value = self._coerce_scalar_value(self.case_label.value)
        case_note_value = self._coerce_scalar_value(self.case_note.value)
        has_peak_feedback = False
        if not self.current_case_candidates.empty:
            has_peak_feedback = bool(
                self.current_case_candidates["human_label"].astype(str).str.strip().ne("").any()
                or self.current_case_candidates["human_note"].astype(str).str.strip().ne("").any()
            )
        self.case_df.at[case_index, "label"] = case_label_value
        self.case_df.at[case_index, "label_note"] = case_note_value
        self.case_df.at[case_index, "reviewed_at_utc"] = (
            timestamp if (case_label_value or case_note_value or has_peak_feedback) else ""
        )
        self.case_df.to_csv(self.case_path, index=False)

        mask = (
            (self.candidate_df["source_run_dir"] == task["source_run_dir"])
            & (self.candidate_df["assay"] == task["assay"])
            & (self.candidate_df["well"] == task["well"])
        )
        current_by_candidate = {
            int(row["candidate_index"]): row
            for _, row in self.current_case_candidates.iterrows()
        }
        existing_candidate_indices = (
            {int(value) for value in self.candidate_df.loc[mask, "candidate_index"].tolist()}
            if not self.candidate_df.empty
            else set()
        )
        for idx in self.candidate_df.loc[mask].index:
            candidate_index = int(self.candidate_df.at[idx, "candidate_index"])
            current_row = current_by_candidate.get(candidate_index)
            if current_row is None:
                continue
            self.candidate_df.at[idx, "human_label"] = str(current_row.get("human_label", "")).strip()
            self.candidate_df.at[idx, "human_note"] = str(current_row.get("human_note", "")).strip()
        new_rows = []
        for candidate_index, current_row in current_by_candidate.items():
            if candidate_index in existing_candidate_indices:
                continue
            new_rows.append({column: current_row.get(column, "") for column in self.candidate_df.columns})
        if new_rows:
            self.candidate_df = pd.concat([self.candidate_df, pd.DataFrame(new_rows)], ignore_index=True)
        self.candidate_df.to_csv(self.candidate_path, index=False)

        if self.summary_path.exists():
            summary = json_load(self.summary_path)
        else:
            summary = {}
        summary["last_reviewed_at_utc"] = timestamp
        self.summary_path.write_text(json_dumps(summary), encoding="utf-8")
        refreshed = self._find_case_by_index(case_index)
        if refreshed is not None:
            self.current_case = refreshed
        self._set_dirty(False)
        return timestamp

    def _save_current_case(self, *, move_next: bool) -> None:
        if self.current_case is None:
            return
        current_index = int(self.current_case["_case_index"])
        prior_indices = [int(task["_case_index"]) for task in self.filtered_cases]
        prior_pos = prior_indices.index(current_index) if current_index in prior_indices else 0
        saved_at = self._persist_current_task()
        self._update_filtered_cases()
        current_indices = [int(task["_case_index"]) for task in self.filtered_cases]

        target_index: int | None = None
        if move_next:
            if current_index in current_indices:
                current_pos = current_indices.index(current_index)
                if current_pos < len(current_indices) - 1:
                    target_index = current_indices[current_pos + 1]
            elif current_indices:
                target_index = current_indices[min(prior_pos, len(current_indices) - 1)]
        else:
            if current_index in current_indices:
                target_index = current_index
            elif current_indices and self.queue_filter.value == "pending":
                target_index = current_indices[min(prior_pos, len(current_indices) - 1)]

        if target_index is not None:
            self.file_select.value = str(target_index)
            if str(target_index) == str(current_index):
                refreshed = self._find_case_by_index(current_index)
                if refreshed is not None:
                    self._load_case(refreshed)
        elif self.filtered_cases:
            self.file_select.value = str(self.filtered_cases[0]["_case_index"])
        else:
            self.current_case = None
            self._clear_views()

        if self.queue_filter.value == "pending" and current_index not in current_indices:
            self._set_status(
                f"Saved ladder review case at {saved_at}. It is now hidden from the pending queue.",
                level="success",
            )
        elif move_next and target_index is None:
            self._set_status("Saved current case. No more cases in the current filter.", level="success")
        else:
            self._set_status(f"Saved ladder review case at {saved_at}.", level="success")


def json_load(path: Path) -> dict[str, Any]:
    import json

    return json.loads(path.read_text(encoding="utf-8"))


def json_dumps(payload: dict[str, Any]) -> str:
    import json

    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    def create_app():
        return LadderReviewAnnotator(
            args.bundle_dir.expanduser().resolve(),
            args.data_dir.expanduser().resolve(),
        ).layout

    pn.serve(create_app, port=args.port, show=True, title="Ladder Review Annotator")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
