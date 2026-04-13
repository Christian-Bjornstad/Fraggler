from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


from core.analysis import apply_manual_ladder_mapping, compute_ladder_qc_metrics  # noqa: E402
from core.analyses.flt3.pipeline import FLT3_TEMPLATE_STEPS  # noqa: E402
from gui_qt.ladder_utils import detect_fsa_for_ladder, load_adjustable_fsa  # noqa: E402


DEFAULT_DATA_ROOT = Path("/Volumes/T7 Shield/DATA/flt3")
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "validation_outputs"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze saved FLT3 manual ladder adjustments and summarize timing/intensity patterns."
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help=f"FLT3 data root. Default: {DEFAULT_DATA_ROOT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to validation_outputs/flt3_manual_adjustment_analysis_<timestamp>.",
    )
    return parser


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def _csv_value(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    if isinstance(value, (list, tuple)):
        return ";".join(str(_csv_value(v)) for v in value)
    return value


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key in seen:
                continue
            seen.add(key)
            fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def _discover_adjustments(data_root: Path) -> list[Path]:
    return sorted(
        p
        for p in data_root.rglob("*.ladder_adj.json")
        if not p.name.startswith("._")
    )


def _fsa_path_for_adjustment(adj_path: Path) -> Path:
    suffix = ".ladder_adj.json"
    name = adj_path.name
    if name.endswith(suffix):
        return adj_path.with_name(name[: -len(suffix)] + ".fsa")
    return adj_path.with_suffix(".fsa")


def _median(values: list[float]) -> float | None:
    finite = [float(v) for v in values if np.isfinite(v)]
    if not finite:
        return None
    return float(statistics.median(finite))


def _percentile(values: list[float], q: float) -> float | None:
    finite = sorted(float(v) for v in values if np.isfinite(v))
    if not finite:
        return None
    idx = int(round((len(finite) - 1) * q))
    return float(finite[max(0, min(len(finite) - 1, idx))])


def _fwhm_and_asymmetry(trace: np.ndarray, peak_time: float) -> tuple[float | None, float | None]:
    if trace.size == 0:
        return None, None
    idx = int(round(float(peak_time)))
    if idx < 0 or idx >= trace.size:
        return None, None
    lo = max(0, idx - 30)
    hi = min(trace.size - 1, idx + 30)
    peak_height = float(trace[idx])
    baseline = float(np.min(trace[lo : hi + 1]))
    half_level = baseline + (peak_height - baseline) * 0.5

    left = idx
    while left > lo and float(trace[left]) > half_level:
        left -= 1
    right = idx
    while right < hi and float(trace[right]) > half_level:
        right += 1

    width = float(right - left) if right > left else None
    left_span = max(1.0, float(idx - left))
    right_span = max(1.0, float(right - idx))
    asymmetry = float(right_span / left_span) if width is not None else None
    return width, asymmetry


def _row_common_fields(fsa_path: Path, metadata: dict[str, Any], qc: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    rel = ""
    try:
        rel = str(fsa_path.relative_to(DEFAULT_DATA_ROOT))
    except Exception:
        rel = str(fsa_path)
    return {
        "path": str(fsa_path),
        "relative_path": rel,
        "year": rel.split("/")[0] if rel else "",
        "file": fsa_path.name,
        "run_dir": fsa_path.parent.name,
        "assay": str(metadata.get("assay") or ""),
        "analysis": str(metadata.get("analysis") or ""),
        "ladder": str(metadata.get("ladder") or ""),
        "manual_candidate_count": len(payload.get("manual_candidates", []) or []),
        "manual_mapping_time_count": len(payload.get("mapping_times", {}) or {}),
        "r2": float(qc.get("r2", float("nan"))),
        "mean_abs_error_bp": float(qc.get("mean_abs_error_bp", float("nan"))),
        "max_abs_error_bp": float(qc.get("max_abs_error_bp", float("nan"))),
    }


def analyze_adjustments(data_root: Path, output_dir: Path) -> dict[str, Any]:
    adjustment_paths = _discover_adjustments(data_root)
    step_rows: list[dict[str, Any]] = []
    gap_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    manual_candidate_rows: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []

    expected_steps = list(map(float, FLT3_TEMPLATE_STEPS))

    for adj_path in adjustment_paths:
        fsa_path = _fsa_path_for_adjustment(adj_path)
        try:
            payload = json.loads(adj_path.read_text(encoding="utf-8"))
            metadata = detect_fsa_for_ladder(fsa_path, preferred_analysis="flt3")
            if not metadata or metadata.get("analysis") != "flt3":
                raise ValueError("File could not be classified as FLT3.")

            fsa, metadata = load_adjustable_fsa(fsa_path, preferred_analysis="flt3", metadata=metadata)
            manual_fsa = apply_manual_ladder_mapping(fsa, payload)
            qc = compute_ladder_qc_metrics(manual_fsa)
            trace = np.asarray(getattr(manual_fsa, "size_standard", []), dtype=float)
            best_size_standard = np.asarray(getattr(manual_fsa, "best_size_standard", []), dtype=float)
            ladder_steps = np.asarray(getattr(manual_fsa, "ladder_steps", expected_steps), dtype=float)

            common = _row_common_fields(fsa_path, metadata, qc, payload)
            file_rows.append(common)

            time_by_step: dict[int, float] = {}
            intensity_by_step: dict[int, float] = {}
            norm_by_step: dict[int, float] = {}
            widths_by_step: dict[int, float | None] = {}
            asym_by_step: dict[int, float | None] = {}

            intensities = []
            for step_bp, peak_time in zip(ladder_steps.tolist(), best_size_standard.tolist()):
                idx = int(round(float(step_bp)))
                peak_index = int(round(float(peak_time)))
                intensity = float(trace[peak_index]) if 0 <= peak_index < trace.size else float("nan")
                intensities.append(intensity)
                time_by_step[idx] = float(peak_time)
                intensity_by_step[idx] = intensity

            valid_anchor_intensities = [v for v in intensities if np.isfinite(v) and v > 0]
            median_intensity = float(statistics.median(valid_anchor_intensities)) if valid_anchor_intensities else float("nan")

            for step_bp in expected_steps:
                idx = int(round(step_bp))
                peak_time = time_by_step.get(idx)
                if peak_time is None:
                    continue
                intensity = intensity_by_step.get(idx, float("nan"))
                width, asymmetry = _fwhm_and_asymmetry(trace, peak_time)
                norm_intensity = (
                    float(intensity / median_intensity)
                    if np.isfinite(intensity) and np.isfinite(median_intensity) and median_intensity > 0
                    else float("nan")
                )
                widths_by_step[idx] = width
                asym_by_step[idx] = asymmetry
                norm_by_step[idx] = norm_intensity
                step_rows.append(
                    {
                        **common,
                        "step_bp": idx,
                        "peak_time": float(peak_time),
                        "peak_intensity": intensity,
                        "peak_intensity_vs_median": norm_intensity,
                        "peak_fwhm_points": width,
                        "peak_asymmetry_ratio": asymmetry,
                    }
                )

            for left_step, right_step in zip(expected_steps[:-1], expected_steps[1:]):
                left_idx = int(round(left_step))
                right_idx = int(round(right_step))
                left_time = time_by_step.get(left_idx)
                right_time = time_by_step.get(right_idx)
                if left_time is None or right_time is None:
                    continue
                gap_rows.append(
                    {
                        **common,
                        "left_step_bp": left_idx,
                        "right_step_bp": right_idx,
                        "gap_time": float(right_time - left_time),
                        "left_intensity_vs_median": norm_by_step.get(left_idx),
                        "right_intensity_vs_median": norm_by_step.get(right_idx),
                        "left_width": widths_by_step.get(left_idx),
                        "right_width": widths_by_step.get(right_idx),
                    }
                )

            for candidate_time in payload.get("manual_candidates", []) or []:
                candidate_time = float(candidate_time)
                peak_index = int(round(candidate_time))
                intensity = float(trace[peak_index]) if 0 <= peak_index < trace.size else float("nan")
                nearest_step = ""
                if time_by_step:
                    nearest_step = min(time_by_step.keys(), key=lambda step: abs(time_by_step[step] - candidate_time))
                manual_candidate_rows.append(
                    {
                        **common,
                        "candidate_time": candidate_time,
                        "candidate_intensity": intensity,
                        "nearest_step_bp": nearest_step,
                        "delta_to_nearest_step_time": (
                            float(candidate_time - time_by_step[nearest_step])
                            if nearest_step != ""
                            else None
                        ),
                    }
                )
        except Exception as exc:
            errors.append({"path": str(fsa_path), "error": str(exc)})

    step_summary_rows: list[dict[str, Any]] = []
    by_step: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in step_rows:
        by_step[int(row["step_bp"])].append(row)
    for step_bp in sorted(by_step):
        rows = by_step[step_bp]
        step_summary_rows.append(
            {
                "step_bp": step_bp,
                "n": len(rows),
                "median_peak_time": _median([float(r["peak_time"]) for r in rows]),
                "median_peak_intensity": _median([float(r["peak_intensity"]) for r in rows]),
                "median_peak_intensity_vs_median": _median([float(r["peak_intensity_vs_median"]) for r in rows]),
                "p10_peak_intensity_vs_median": _percentile([float(r["peak_intensity_vs_median"]) for r in rows], 0.10),
                "p90_peak_intensity_vs_median": _percentile([float(r["peak_intensity_vs_median"]) for r in rows], 0.90),
                "median_peak_fwhm_points": _median([float(r["peak_fwhm_points"]) for r in rows if r["peak_fwhm_points"] is not None]),
                "median_peak_asymmetry_ratio": _median([float(r["peak_asymmetry_ratio"]) for r in rows if r["peak_asymmetry_ratio"] is not None]),
            }
        )

    gap_summary_rows: list[dict[str, Any]] = []
    by_gap: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in gap_rows:
        by_gap[(int(row["left_step_bp"]), int(row["right_step_bp"]))].append(row)
    for left_step, right_step in sorted(by_gap):
        rows = by_gap[(left_step, right_step)]
        gap_summary_rows.append(
            {
                "left_step_bp": left_step,
                "right_step_bp": right_step,
                "n": len(rows),
                "median_gap_time": _median([float(r["gap_time"]) for r in rows]),
                "p10_gap_time": _percentile([float(r["gap_time"]) for r in rows], 0.10),
                "p90_gap_time": _percentile([float(r["gap_time"]) for r in rows], 0.90),
                "median_left_intensity_vs_median": _median([float(r["left_intensity_vs_median"]) for r in rows if r["left_intensity_vs_median"] is not None]),
                "median_right_intensity_vs_median": _median([float(r["right_intensity_vs_median"]) for r in rows if r["right_intensity_vs_median"] is not None]),
            }
        )

    focus_pairs = {(139, 150), (150, 160), (340, 350), (450, 490), (490, 500)}
    focus_gap_rows = [row for row in gap_summary_rows if (row["left_step_bp"], row["right_step_bp"]) in focus_pairs]

    file_rows_sorted = sorted(file_rows, key=lambda row: (str(row["year"]), str(row["file"])))
    status = {
        "adjustment_file_count": len(adjustment_paths),
        "analyzed_file_count": len(file_rows),
        "error_count": len(errors),
        "year_counts": dict(sorted(Counter(str(row["year"]) for row in file_rows).items())),
        "assay_counts": dict(sorted(Counter(str(row["assay"]) for row in file_rows).items())),
        "median_r2": _median([float(row["r2"]) for row in file_rows]),
        "median_max_abs_error_bp": _median([float(row["max_abs_error_bp"]) for row in file_rows]),
        "median_mean_abs_error_bp": _median([float(row["mean_abs_error_bp"]) for row in file_rows]),
        "manual_candidate_entry_count": len(manual_candidate_rows),
    }

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_root": str(data_root),
        "output_dir": str(output_dir),
        "summary": status,
        "focus_gap_rows": focus_gap_rows,
        "errors": errors,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "manual_adjustment_files.csv", file_rows_sorted)
    _write_csv(output_dir / "manual_adjustment_step_metrics.csv", step_rows)
    _write_csv(output_dir / "manual_adjustment_step_summary.csv", step_summary_rows)
    _write_csv(output_dir / "manual_adjustment_gap_metrics.csv", gap_rows)
    _write_csv(output_dir / "manual_adjustment_gap_summary.csv", gap_summary_rows)
    _write_csv(output_dir / "manual_adjustment_manual_candidates.csv", manual_candidate_rows)
    (output_dir / "summary.json").write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# FLT3 Manual Ladder Adjustment Analysis",
        "",
        f"Generated UTC: {payload['generated_at_utc']}",
        "",
        f"Manual adjustment files: **{status['adjustment_file_count']}**",
        f"Analyzed files: **{status['analyzed_file_count']}**",
        f"Errors: **{status['error_count']}**",
        f"Median R2: **{status['median_r2']}**",
        f"Median max residual bp: **{status['median_max_abs_error_bp']}**",
        "",
        "## Focus Gap Pairs",
        "",
        "| Pair | Median gap | P10 | P90 |",
        "| --- | --- | --- | --- |",
    ]
    for row in focus_gap_rows:
        lines.append(
            f"| {int(row['left_step_bp'])}-{int(row['right_step_bp'])} | "
            f"{row['median_gap_time']} | {row['p10_gap_time']} | {row['p90_gap_time']} |"
        )
    lines.extend(
        [
            "",
            "## Step Intensity Summary",
            "",
            "| Step | Median intensity vs file median | P10 | P90 | Median FWHM | Median asymmetry |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    for row in step_summary_rows:
        if int(row["step_bp"]) not in {139, 150, 160, 340, 350, 400, 450, 490, 500}:
            continue
        lines.append(
            f"| {int(row['step_bp'])} | {row['median_peak_intensity_vs_median']} | "
            f"{row['p10_peak_intensity_vs_median']} | {row['p90_peak_intensity_vs_median']} | "
            f"{row['median_peak_fwhm_points']} | {row['median_peak_asymmetry_ratio']} |"
        )
    (output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return payload


def main() -> None:
    args = build_arg_parser().parse_args()
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else (DEFAULT_OUTPUT_ROOT / f"flt3_manual_adjustment_analysis_{_timestamp()}").resolve()
    )
    payload = analyze_adjustments(args.data_root.expanduser().resolve(), output_dir)
    print(json.dumps(_json_safe(payload), indent=2, sort_keys=True))
    print(f"Summary JSON: {output_dir / 'summary.json'}")
    print(f"Summary Markdown: {output_dir / 'summary.md'}")


if __name__ == "__main__":
    main()
