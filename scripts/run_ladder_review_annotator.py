from __future__ import annotations

import argparse
import logging
import re
import sys
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


def load_fsa_full_analysis(fsa_path: Path, ladder: str):
    from core.analysis import analyse_fsa_liz, analyse_fsa_rox

    if not fsa_path.exists():
        raise FileNotFoundError(f"Missing FSA file: {fsa_path}")

    ladder_name = str(ladder or "").upper()
    if ladder_name == "LIZ":
        fsa = analyse_fsa_liz(fsa_path, sample_channel="DATA1")
    else:
        fsa = analyse_fsa_rox(fsa_path, sample_channel="DATA1")

    if fsa is None:
        raise ValueError(f"Could not fit ladder for {fsa_path.name}")
    return fsa


def extract_ladder_trace(fsa, ladder: str) -> tuple[np.ndarray, str]:
    target = "DATA5" if str(ladder).upper() == "LIZ" else "DATA4"
    possible = [target, str(ladder).upper()]
    for key in fsa.fsa.keys():
        key_str = key.decode("ascii", errors="ignore") if isinstance(key, bytes) else str(key)
        if key_str.upper() in possible:
            return np.asarray(fsa.fsa[key], dtype=float), key_str
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
        self.current_idx = 0
        self.current_task_candidates: pd.DataFrame | None = None
        self.candidate_option_map: dict[str, int] = {}

        self.pending_cases = self.case_df[
            self.case_df["label"].astype(str).str.strip() == ""
        ].reset_index(drop=False).rename(columns={"index": "_case_index"})
        self.tasks = self.pending_cases.to_dict("records")

        self.title = pn.pane.Markdown("### Ladder Review Annotator")
        self.progress = pn.pane.Markdown("")
        self.meta = pn.pane.Markdown("")
        self.status = pn.pane.Markdown("")
        self.plot = pn.pane.Plotly(height=500)

        self.case_label = pn.widgets.RadioButtonGroup(
            name="Case Label",
            options={
                "Accept Current Fit": "accept_current_fit",
                "Needs Better Fit": "needs_better_fit",
                "Bad Signal": "bad_signal",
            },
            button_type="primary",
        )
        self.case_note = pn.widgets.TextAreaInput(name="Case Note", height=90)
        self.candidate_selector = pn.widgets.CheckBoxGroup(name="Relevant Ladder Peaks", options=[])
        self.prev_btn = pn.widgets.Button(name="Previous", button_type="default")
        self.save_btn = pn.widgets.Button(name="Save And Next", button_type="success")
        self.skip_btn = pn.widgets.Button(name="Skip", button_type="warning")

        self.prev_btn.on_click(lambda _event: self.prev_task())
        self.save_btn.on_click(lambda _event: self.save_and_next())
        self.skip_btn.on_click(lambda _event: self.next_task())

        self.layout = pn.Column(
            pn.Row(self.title, self.progress),
            self.meta,
            self.status,
            self.plot,
            self.case_label,
            self.case_note,
            self.candidate_selector,
            pn.Row(self.prev_btn, self.skip_btn, self.save_btn),
            sizing_mode="stretch_width",
        )

        if not self.tasks:
            self.progress.object = "**No pending ladder cases.**"
            self.save_btn.disabled = True
            self.skip_btn.disabled = True
            self.prev_btn.disabled = True
        else:
            self.load_task(0)

    def _task_candidates(self, task: dict[str, Any]) -> pd.DataFrame:
        mask = (
            (self.candidate_df["source_run_dir"] == task["source_run_dir"])
            & (self.candidate_df["assay"] == task["assay"])
            & (self.candidate_df["well"] == task["well"])
        )
        rows = self.candidate_df.loc[mask].copy()
        rows["candidate_index"] = rows["candidate_index"].astype(int)
        rows["candidate_time"] = rows["candidate_time"].astype(float)
        rows["candidate_intensity"] = rows["candidate_intensity"].astype(float)
        rows = rows.sort_values(["candidate_index"]).reset_index(drop=True)
        return rows

    def load_task(self, idx: int) -> None:
        if idx < 0 or idx >= len(self.tasks):
            return
        self.current_idx = idx
        task = self.tasks[idx]
        candidates = self._task_candidates(task)
        self.current_task_candidates = candidates
        self.progress.object = f"**{idx + 1} / {len(self.tasks)} pending cases**"
        self.case_label.value = task.get("label") or None
        self.case_note.value = task.get("label_note") or ""

        try:
            fsa_path = resolve_fsa_path(
                self.data_dir,
                str(task["source_run_dir"]),
                str(task["assay"]),
                str(task["well"]),
                str(task["run_code"]),
            )
            fsa = load_fsa_full_analysis(fsa_path, str(task["ladder"]))
            trace_raw, trace_channel = extract_ladder_trace(fsa, str(task["ladder"]))
            from core.analysis import estimate_running_baseline

            try:
                baseline = estimate_running_baseline(trace_raw, bin_size=200, quantile=0.10)
                y_arr = np.clip(trace_raw - baseline, a_min=0, a_max=None)
            except Exception:
                y_arr = trace_raw.astype(float)

            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=np.arange(len(y_arr)),
                    y=y_arr,
                    mode="lines",
                    line=dict(color="#1f77b4", width=2),
                    fill="tozeroy",
                    fillcolor="rgba(31,119,180,0.08)",
                    name=f"{task['ladder']} trace",
                )
            )

            self.candidate_option_map = {}
            options: list[tuple[str, str]] = []
            selected_values: list[str] = []
            for _, row in candidates.iterrows():
                idx_value = int(row["candidate_index"])
                time_value = float(row["candidate_time"])
                intensity_value = float(row["candidate_intensity"])
                option_value = str(idx_value)
                label = (
                    f"#{idx_value} | t={time_value:.1f} | y={intensity_value:.0f}"
                    f" | selected={row['selected_for_fit']}"
                )
                options.append((label, option_value))
                self.candidate_option_map[option_value] = idx_value
                if str(row.get("human_label", "")).strip().lower() in {"keep_peak", "relevant_peak"}:
                    selected_values.append(option_value)

                is_selected = str(row["selected_for_fit"]).strip().lower() == "true"
                color = "#ffd166" if is_selected else "#ef476f"
                size = 16 if is_selected else 10
                fig.add_trace(
                    go.Scatter(
                        x=[time_value],
                        y=[intensity_value],
                        mode="markers+text",
                        marker=dict(color=color, size=size, symbol="diamond" if is_selected else "circle"),
                        text=[str(idx_value)],
                        textposition="top center",
                        name=f"candidate {idx_value}",
                    )
                )

            self.candidate_selector.options = options
            self.candidate_selector.value = selected_values

            fig.update_layout(
                title=(
                    f"{task['source_run_dir']} | {task['assay']} {task['well']} | "
                    f"strategy={task['ladder_fit_strategy']} | fitted={task['ladder_fitted_step_count']}/{task['ladder_expected_step_count']}"
                ),
                xaxis_title=f"Ladder Time / Scan Index ({trace_channel})",
                yaxis_title="Signal",
                template="plotly_white",
                height=520,
                showlegend=False,
                margin=dict(l=40, r=20, t=70, b=40),
            )
            self.plot.object = fig
            self.meta.object = (
                f"**Run:** `{task['source_run_dir']}`  \n"
                f"**Assay:** `{task['assay']}`  \n"
                f"**Well:** `{task['well']}`  \n"
                f"**Ladder:** `{task['ladder']}`  \n"
                f"**QC:** `{task['ladder_qc']}`  \n"
                f"**R2:** `{task['ladder_r2']}`"
            )
            self.status.object = ""
        except Exception as exc:
            logging.exception("Failed to load ladder review task")
            self.plot.object = go.Figure()
            self.meta.object = ""
            self.status.object = f"**Error loading task:** {exc}"

    def _persist_current_task(self) -> None:
        if not self.tasks:
            return
        task = self.tasks[self.current_idx]
        case_index = int(task["_case_index"])
        timestamp = datetime.now(timezone.utc).isoformat()

        case_label_value = self._coerce_scalar_value(self.case_label.value)
        case_note_value = self._coerce_scalar_value(self.case_note.value)
        self.case_df.at[case_index, "label"] = case_label_value
        self.case_df.at[case_index, "label_note"] = case_note_value
        self.case_df.to_csv(self.case_path, index=False)

        if self.current_task_candidates is not None:
            selected_values = self.candidate_selector.value or []
            selected: set[int] = set()
            for value in selected_values:
                normalized_value = self._coerce_option_value(value)
                if normalized_value in self.candidate_option_map:
                    selected.add(self.candidate_option_map[normalized_value])
            mask = (
                (self.candidate_df["source_run_dir"] == task["source_run_dir"])
                & (self.candidate_df["assay"] == task["assay"])
                & (self.candidate_df["well"] == task["well"])
            )
            for idx in self.candidate_df.loc[mask].index:
                candidate_index = int(self.candidate_df.at[idx, "candidate_index"])
                self.candidate_df.at[idx, "human_label"] = "keep_peak" if candidate_index in selected else "reject_peak"
                self.candidate_df.at[idx, "human_note"] = case_note_value
            self.candidate_df.to_csv(self.candidate_path, index=False)

        if self.summary_path.exists():
            summary = json_load(self.summary_path)
        else:
            summary = {}
        summary["last_reviewed_at_utc"] = timestamp
        self.summary_path.write_text(json_dumps(summary), encoding="utf-8")

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
    def _coerce_option_value(value: Any) -> str:
        if isinstance(value, tuple):
            if len(value) >= 2:
                return str(value[1])
            if len(value) == 1:
                return str(value[0])
            return ""
        return str(value)

    def save_and_next(self) -> None:
        self._persist_current_task()
        self.next_task()

    def next_task(self) -> None:
        if self.current_idx < len(self.tasks) - 1:
            self.load_task(self.current_idx + 1)
        else:
            self._persist_current_task()
            self.progress.object = "**All ladder cases reviewed.**"
            self.status.object = "You’re done with the current ladder review bundle."

    def prev_task(self) -> None:
        if self.current_idx > 0:
            self.load_task(self.current_idx - 1)


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
