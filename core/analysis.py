"""
Fraggler Diagnostics — Analysis Functions.

Ladder fitting (LIZ / ROX), SL peak detection, ladder QC metrics,
SL area metrics, local-maxima helpers, and running-baseline estimation.
"""
from __future__ import annotations

from pathlib import Path
import copy
import time

import numpy as np
import pandas as pd
from scipy import signal
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import mean_squared_error, r2_score

import json
from fraggler.fraggler import (
    FsaFile,
    estimate_combination_count,
    find_size_standard_peaks,
    return_maxium_allowed_distance_between_size_standard_peaks,
    generate_combinations,
    calculate_best_combination_of_size_standard_peaks,
    fit_size_standard_to_ladder,
    baseline_arPLS,
    print_green,
    print_warning,
)

from core.assay_config import (
    LIZ_LADDER,
    ROX_LADDER,
    MIN_DISTANCE_BETWEEN_PEAKS_LIZ,
    MIN_SIZE_STANDARD_HEIGHT_LIZ,
    MIN_DISTANCE_BETWEEN_PEAKS_ROX,
    MIN_SIZE_STANDARD_HEIGHT_ROX,
    SL_WINDOW_BP,
)

# --------------------------------------------------------------
# Analysis Constants (extracted magic numbers)
# --------------------------------------------------------------
LADDER_MAX_ITERATIONS = 15
BASELINE_BIN_SIZE = 200
BASELINE_QUANTILE = 0.10
PEAK_MIN_HEIGHT = 800.0
YMAX_PADDING_FACTOR = 1.15
MIN_R2_QUALITY = 0.998  # Post-fit quality gate
LADDER_CANDIDATE_COUNT = 5
HIGH_END_RESCUE_R2 = 0.9999
LOW_INTENSITY_RATIO_FLOOR = 0.55
MEDIAN_INTENSITY_TARGET_RATIO = 0.80
EARLY_PEAK_INTENSITY_WEIGHT = 1.2
GLOBAL_PEAK_INTENSITY_WEIGHT = 0.45
SEVERE_WEAK_PEAK_PENALTY = 0.30
DESCENDING_RECOVERY_R2_FLOOR = 0.9985
DESCENDING_RECOVERY_MAX_ABS_ERROR = 3.0
DESCENDING_RECOVERY_MEAN_ABS_ERROR = 1.6
ASCENDING_RECOVERY_R2_FLOOR = 0.9985
ASCENDING_RECOVERY_MAX_ABS_ERROR = 3.0
ASCENDING_RECOVERY_MEAN_ABS_ERROR = 1.6
ASCENDING_RECOVERY_MIN_INTENSITY = 300.0
GENERAL_COMPLETION_R2_FLOOR = 0.9985
GENERAL_COMPLETION_MAX_ABS_ERROR = 3.0
GENERAL_COMPLETION_MEAN_ABS_ERROR = 1.8
GENERAL_COMPLETION_MIN_INTENSITY = 300.0
CORE_COMPLETION_MIN_ASSIGNED = 12


def _log_ladder_timing(label: str, phase: str, fsa_path: Path, elapsed_seconds: float, **details: object) -> None:
    detail_text = ""
    if details:
        detail_text = " | " + ", ".join(f"{key}={value}" for key, value in details.items())
    print_green(
        f"[{label}][TIMING] {phase} for {fsa_path.name} took {elapsed_seconds:.3f}s{detail_text}"
    )
DESCENDING_RECOVERY_MIN_INTENSITY = 350.0
PARTIAL_RESCUE_MISSING_STEP_PENALTY = 0.50
ROX_DYEBLOB_HEIGHT_MULTIPLIER = 2.0
ROX_DYEBLOB_EARLY_INDEX = 2000
ROX_DYEBLOB_TIGHT_GAP = 40
ROX_DYEBLOB_CLUSTER_GAP = 70
ROX_BASELINE_FALLBACK_MIN_HEIGHT = 50.0
ROX_BASELINE_FALLBACK_MIN_PEAKS = 3
LADDER_RESCORING_MAX_COMBINATIONS = 250
ROX_COMBINATION_ESTIMATE_LIMIT = 10_000
ROX_ALLOWED_EXTRA_SIZE_STANDARD_PEAKS = 2
ROX_BEAM_WIDTH = 64
ROX_BEAM_KEEP_FINISHED = 5
ROX_BEAM_MIN_COMPLETION_RATIO = 0.60
EARLY_ACCEPT_R2 = 0.99995
EARLY_ACCEPT_MEAN_ABS_ERROR = 0.35
EARLY_ACCEPT_MAX_ABS_ERROR = 0.90
ROX_PREFERRED_TIME_MIN = 1500.0
ROX_PREFERRED_TIME_MAX = 4000.0
ROX_PREFERRED_TIME_MARGIN = 75.0
ROX_HARD_FILTER_TIME_MIN = 1450.0
ROX_HARD_FILTER_TIME_MAX = 4050.0
ROX_APEX_SNAP_RADIUS = 6
ROX_PREFERRED_SUPPLEMENT_MIN_HEIGHT = 250.0
ROX_PREFERRED_SUPPLEMENT_DISTANCE = 15
ROX_PROFILE_TIME_WEIGHT = 0.10
ROX_PROFILE_LOW_INTENSITY_WEIGHT = 0.04
ROX_PROFILE_SEVERE_WEAK_PENALTY = 0.08
ROX_PROFILE_SEVERE_WEAK_INTENSITY = 80.0
ROX_BEAM_EXPECTED_GAP_WEIGHT = 0.60
ROX_EDGE_MISSING_STEP_PENALTY = 0.85
ROX_NEAR_EDGE_MISSING_STEP_PENALTY = 0.60
ROX_MIDDLE_MISSING_STEP_PENALTY = 0.35
ROX_TAIL_EXPANSION_STEPS = 5
ROX_TAIL_GAP_MULTIPLIER = 3.10
ROX_PARTIAL_ALIGNMENT_VARIANTS = 8


def _project_root_for(path: Path) -> Path | None:
    parts = path.resolve().parts
    try:
        desktop_idx = parts.index("Desktop")
    except ValueError:
        return None
    if desktop_idx + 1 >= len(parts):
        return None
    return Path(*parts[: desktop_idx + 2])


def _sibling_fsa_paths(fsa_path: Path) -> list[Path]:
    root = _project_root_for(fsa_path)
    if root is None:
        return []

    desktop_root = root.parent
    try:
        rel = fsa_path.resolve().relative_to(root.resolve())
    except Exception:
        return []

    siblings: list[Path] = []
    for candidate_root in desktop_root.iterdir():
        if candidate_root == root or not candidate_root.is_dir():
            continue
        name = candidate_root.name.lower()
        if name.startswith(".") or "backup" in name:
            continue
        candidate = candidate_root / rel
        if candidate.exists():
            siblings.append(candidate)
    return siblings


def _get_expected_ladder_steps(fsa: FsaFile) -> np.ndarray:
    expected = getattr(fsa, "expected_ladder_steps", None)
    if expected is None:
        return np.asarray(fsa.ladder_steps, dtype=float)
    return np.asarray(expected, dtype=float)


def _missing_expected_ladder_steps(fsa: FsaFile) -> list[float]:
    expected = _get_expected_ladder_steps(fsa)
    current = np.asarray(getattr(fsa, "ladder_steps", expected), dtype=float)
    missing = [float(bp) for bp in expected if not np.any(np.isclose(current, bp, atol=1e-6))]
    return missing


def _count_missing_low_end_steps(fsa: FsaFile) -> int:
    expected = _get_expected_ladder_steps(fsa)
    current = np.asarray(getattr(fsa, "ladder_steps", expected), dtype=float)
    missing = 0
    for bp in expected:
        if np.any(np.isclose(current, bp, atol=1e-6)):
            break
        missing += 1
    return missing


def _set_ladder_fit_metadata(fsa: FsaFile, strategy: str, note: str | None = None) -> FsaFile:
    fsa.ladder_fit_strategy = strategy
    fsa.ladder_missing_expected_steps = _missing_expected_ladder_steps(fsa)
    fsa.ladder_review_required = bool(fsa.ladder_missing_expected_steps)
    fsa.ladder_expected_step_count = int(len(_get_expected_ladder_steps(fsa)))
    fsa.ladder_fitted_step_count = int(len(getattr(fsa, "ladder_steps", [])))
    if note is None:
        if fsa.ladder_missing_expected_steps:
            missing_txt = ", ".join(f"{bp:.0f}" for bp in fsa.ladder_missing_expected_steps)
            note = f"Missing expected ladder steps: {missing_txt} bp"
        else:
            note = "All expected ladder steps were fitted."
    fsa.ladder_fit_note = note
    return fsa


def _finalize_auto_fit_metadata(fsa: FsaFile) -> FsaFile:
    existing = getattr(fsa, "ladder_fit_strategy", None)
    if existing:
        if not hasattr(fsa, "ladder_missing_expected_steps"):
            fsa.ladder_missing_expected_steps = _missing_expected_ladder_steps(fsa)
        fsa.ladder_review_required = bool(getattr(fsa, "ladder_missing_expected_steps", []))
        fsa.ladder_expected_step_count = int(len(_get_expected_ladder_steps(fsa)))
        fsa.ladder_fitted_step_count = int(len(getattr(fsa, "ladder_steps", [])))
        if not getattr(fsa, "ladder_fit_note", None):
            _set_ladder_fit_metadata(fsa, existing)
        return fsa
    strategy = "auto_full" if not _missing_expected_ladder_steps(fsa) else "auto_partial"
    return _set_ladder_fit_metadata(fsa, strategy)


def _map_step_indices(source_steps: np.ndarray, target_steps: np.ndarray) -> dict[int, int]:
    mapping: dict[int, int] = {}
    used: set[int] = set()
    for source_idx, source_bp in enumerate(np.asarray(source_steps, dtype=float)):
        matches = np.where(np.isclose(target_steps, source_bp, atol=1e-6))[0]
        for target_idx in matches:
            if int(target_idx) in used:
                continue
            mapping[int(source_idx)] = int(target_idx)
            used.add(int(target_idx))
            break
    return mapping


def _clone_fsa_for_ladder_trial(
    fsa: FsaFile,
    *,
    strip_candidate_table: bool = True,
) -> FsaFile:
    """Clone only ladder-mutable state while reusing heavy trace payloads."""
    trial = copy.copy(fsa)

    for attr_name in (
        "ladder_steps",
        "expected_ladder_steps",
        "size_standard_peaks",
        "best_size_standard",
        "ladder_missing_expected_steps",
    ):
        value = getattr(fsa, attr_name, None)
        if isinstance(value, np.ndarray):
            setattr(trial, attr_name, value.copy())
        elif isinstance(value, list):
            setattr(trial, attr_name, list(value))

    if hasattr(fsa, "best_size_standard_combinations"):
        combinations = getattr(fsa, "best_size_standard_combinations")
        if strip_candidate_table:
            setattr(trial, "best_size_standard_combinations", None)
        elif isinstance(combinations, pd.DataFrame):
            setattr(trial, "best_size_standard_combinations", combinations.copy(deep=False))

    return trial




def _candidate_combination_arrays(
    fsa: FsaFile,
    *,
    max_combinations: int = LADDER_RESCORING_MAX_COMBINATIONS,
) -> list[np.ndarray]:
    """Return a bounded candidate subset from the generated ladder combinations."""
    combinations = getattr(fsa, "best_size_standard_combinations", None)
    if combinations is None or getattr(combinations, "empty", True):
        return []

    total = int(getattr(combinations, "shape", [0])[0])
    if total <= max_combinations:
        source = combinations["combinations"].tolist()
    else:
        sampled_idx = np.linspace(0, total - 1, num=max_combinations, dtype=int)
        sampled_idx = np.unique(sampled_idx)
        print_warning(
            f"[LADDER] Sampling {len(sampled_idx)}/{total} ladder combinations for "
            f"{Path(getattr(fsa, 'file', 'unknown')).name} to avoid a long stall."
        )
        source = combinations.iloc[sampled_idx]["combinations"].tolist()

    return [np.asarray(combo, dtype=float) for combo in source]


def _rank_size_standard_combinations(fsa: FsaFile) -> list[np.ndarray]:
    """Return the smoothest ladder candidates using a bounded candidate subset."""
    ranked: list[tuple[float, np.ndarray]] = []
    ladder_steps = np.asarray(fsa.ladder_steps, dtype=float)

    for combo in _candidate_combination_arrays(fsa):
        if combo.size != ladder_steps.size:
            continue
        try:
            derivative = UnivariateSpline(ladder_steps, combo, s=0).derivative(n=2)
            score = float(max(abs(derivative(ladder_steps))))
        except Exception:
            continue
        ranked.append((score, combo))

    ranked.sort(key=lambda item: item[0])
    return [combo for _, combo in ranked[:LADDER_CANDIDATE_COUNT]]


def _fit_score_tuple(
    metrics: dict[str, float | int],
    intensity_penalty: float,
    *,
    missing_penalty: float = 0.0,
) -> tuple[float, float, float, float]:
    return (
        float(metrics.get("mean_abs_error_bp", float("inf")))
        + float(intensity_penalty)
        + float(missing_penalty),
        float(metrics.get("max_abs_error_bp", float("inf"))),
        -float(metrics.get("r2", float("-inf"))),
        float(intensity_penalty),
    )


def _is_early_accept_candidate(
    metrics: dict[str, float | int],
    *,
    missing_count: int = 0,
) -> bool:
    return (
        missing_count == 0
        and float(metrics.get("r2", float("-inf"))) >= EARLY_ACCEPT_R2
        and float(metrics.get("mean_abs_error_bp", float("inf"))) <= EARLY_ACCEPT_MEAN_ABS_ERROR
        and float(metrics.get("max_abs_error_bp", float("inf"))) <= EARLY_ACCEPT_MAX_ABS_ERROR
    )


def _missing_step_penalty(fsa: FsaFile) -> float:
    missing_steps = _missing_expected_ladder_steps(fsa)
    if not missing_steps:
        return 0.0

    ladder_name = str(getattr(fsa, "ladder", "") or "").upper()
    expected = _get_expected_ladder_steps(fsa)
    if "ROX" not in ladder_name or expected.size == 0:
        return len(missing_steps) * PARTIAL_RESCUE_MISSING_STEP_PENALTY

    penalty = 0.0
    for missing_bp in missing_steps:
        idx_matches = np.where(np.isclose(expected, missing_bp, atol=1e-6))[0]
        if idx_matches.size == 0:
            penalty += PARTIAL_RESCUE_MISSING_STEP_PENALTY
            continue
        idx = int(idx_matches[0])
        if idx <= 1 or idx >= (len(expected) - 3):
            penalty += ROX_EDGE_MISSING_STEP_PENALTY
        elif idx <= 3 or idx >= (len(expected) - 5):
            penalty += ROX_NEAR_EDGE_MISSING_STEP_PENALTY
        else:
            penalty += ROX_MIDDLE_MISSING_STEP_PENALTY
    return penalty


def _estimate_size_standard_combination_count(
    fsa: FsaFile,
    *,
    cap: int = ROX_COMBINATION_ESTIMATE_LIMIT + 1,
) -> int:
    peaks = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    if peaks.size == 0:
        return 0
    length = int(getattr(fsa, "n_ladder_peaks", len(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float))))
    distance = float(getattr(fsa, "maxium_allowed_distance_between_size_standard_peaks", 0.0) or 0.0)
    return estimate_combination_count(peaks, length, distance, cap=cap)


def _should_use_bounded_rox_search(
    fsa: FsaFile,
    *,
    combination_estimate: int | None = None,
) -> tuple[bool, int]:
    peak_count = int(len(np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)))
    expected_count = int(len(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)))
    combination_estimate = (
        _estimate_size_standard_combination_count(fsa)
        if combination_estimate is None
        else int(combination_estimate)
    )
    use_bounded = (
        combination_estimate > ROX_COMBINATION_ESTIMATE_LIMIT
        or peak_count > (expected_count + ROX_ALLOWED_EXTRA_SIZE_STANDARD_PEAKS)
    )
    return use_bounded, combination_estimate


def _build_bounded_rox_candidate_specs(
    fsa: FsaFile,
    *,
    beam_width: int = ROX_BEAM_WIDTH,
    keep_finished: int = ROX_BEAM_KEEP_FINISHED,
    allow_partial: bool = True,
) -> list[dict[str, object]]:
    peaks = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    expected_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    if peaks.size == 0 or expected_steps.size == 0:
        return []

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    peak_idx = np.rint(peaks).astype(int)
    valid_idx = np.clip(peak_idx, 0, max(len(trace) - 1, 0))
    intensities = trace[valid_idx] if trace.size else np.ones_like(peaks, dtype=float)
    positive_intensities = intensities[intensities > 0]
    global_intensity = float(np.median(positive_intensities)) if positive_intensities.size else 1.0
    global_gap = float(np.median(np.diff(peaks))) if peaks.size > 1 else 1.0
    max_gap = float(getattr(fsa, "maxium_allowed_distance_between_size_standard_peaks", 0.0) or 0.0)
    expected_bp_gaps = np.diff(expected_steps) if expected_steps.size > 1 else np.array([], dtype=float)
    target_len = int(expected_steps.size)
    min_partial_len = max(10, int(np.ceil(target_len * ROX_BEAM_MIN_COMPLETION_RATIO)))

    states: list[tuple[float, list[int], float, float]] = []
    for i in range(len(peaks)):
        start_ratio = float(intensities[i]) / max(global_intensity, 1.0)
        start_penalty = max(0.0, 0.90 - start_ratio) * 0.60
        start_penalty += _rox_peak_time_penalty(float(peaks[i])) * 0.90
        states.append((start_penalty, [i], global_gap, global_gap))
    states.sort(key=lambda item: item[0])
    states = states[:beam_width]

    best_states = list(states)
    best_depth = 1

    for _depth in range(1, target_len):
        next_states: list[tuple[float, list[int], float, float]] = []
        for score, path_indices, mean_gap, last_gap in states:
            last_idx = path_indices[-1]
            remaining_needed = target_len - len(path_indices) - 1
            adaptive_gap_limit = max_gap
            is_tail_expansion = len(path_indices) >= max(1, target_len - ROX_TAIL_EXPANSION_STEPS - 1)
            if is_tail_expansion:
                adaptive_gap_limit = max(
                    adaptive_gap_limit,
                    float(last_gap) * ROX_TAIL_GAP_MULTIPLIER,
                    global_gap * ROX_TAIL_GAP_MULTIPLIER,
                )
            for next_idx in range(last_idx + 1, len(peaks)):
                if adaptive_gap_limit > 0 and (peaks[next_idx] - peaks[last_idx]) > adaptive_gap_limit:
                    break
                if (len(peaks) - (next_idx + 1)) < max(0, remaining_needed):
                    continue
                gap = float(peaks[next_idx] - peaks[last_idx])
                local_target = last_gap if len(path_indices) > 1 else global_gap
                late_relaxation = 0.35 if is_tail_expansion else 1.0
                smooth_penalty = (abs(gap - local_target) / max(local_target, 1.0)) * late_relaxation
                drift_penalty = (abs(gap - mean_gap) / max(mean_gap, 1.0)) * late_relaxation
                intensity_ratio = float(intensities[next_idx]) / max(global_intensity, 1.0)
                intensity_penalty = max(0.0, 0.90 - intensity_ratio)
                edge_penalty = (0.05 if is_tail_expansion else 0.15) if gap > (max_gap * 0.90) else 0.0
                time_penalty = _rox_peak_time_penalty(float(peaks[next_idx])) * 0.90
                expected_gap_penalty = 0.0
                gap_position = len(path_indices) - 1
                if 0 <= gap_position < len(expected_bp_gaps):
                    expected_gap = float(expected_bp_gaps[gap_position])
                    if gap_position > 0:
                        previous_expected_gap = float(expected_bp_gaps[gap_position - 1])
                        expected_ratio = expected_gap / max(previous_expected_gap, 1.0)
                        observed_ratio = gap / max(last_gap, 1.0)
                        expected_gap_penalty = abs(observed_ratio - expected_ratio) * late_relaxation
                next_score = (
                    float(score)
                    + (smooth_penalty * 1.10)
                    + (drift_penalty * 0.45)
                    + (intensity_penalty * 0.80)
                    + (expected_gap_penalty * ROX_BEAM_EXPECTED_GAP_WEIGHT)
                    + time_penalty
                    + edge_penalty
                )
                next_mean_gap = gap if len(path_indices) == 1 else ((mean_gap * (len(path_indices) - 1)) + gap) / len(path_indices)
                next_states.append((next_score, path_indices + [next_idx], next_mean_gap, gap))

        if not next_states:
            break

        next_states.sort(
            key=lambda item: (
                item[0],
                -len(item[1]),
                -float(np.sum(intensities[np.asarray(item[1], dtype=int)])),
            )
        )
        states = next_states[:beam_width]
        current_depth = len(states[0][1])
        if current_depth >= best_depth:
            best_depth = current_depth
            best_states = list(states)
        if current_depth >= target_len:
            break

    candidate_states = [
        state for state in best_states
        if len(state[1]) >= target_len
    ]
    if candidate_states:
        specs: list[dict[str, object]] = []
        for score, path_indices, _mean_gap, _last_gap in candidate_states[:keep_finished]:
            specs.append(
                {
                    "times": peaks[np.asarray(path_indices, dtype=int)],
                    "ladder_steps": expected_steps.copy(),
                    "beam_score": float(score),
                    "complete": True,
                    "bounded": True,
                }
            )
        return specs

    if not allow_partial:
        return []

    partial_depth = max((len(state[1]) for state in best_states), default=0)
    if partial_depth < min_partial_len:
        return []

    partial_states = [state for state in best_states if len(state[1]) == partial_depth]
    specs = []
    for score, path_indices, _mean_gap, _last_gap in partial_states[:keep_finished]:
        times = peaks[np.asarray(path_indices, dtype=int)]
        for ladder_steps in _build_partial_rox_step_assignments(
            expected_steps,
            times,
            max_variants=ROX_PARTIAL_ALIGNMENT_VARIANTS,
        ):
            specs.append(
                {
                    "times": times,
                    "ladder_steps": ladder_steps,
                    "beam_score": float(score),
                    "complete": False,
                    "bounded": True,
                }
            )
    return specs


def _round_to_monotonic_indices(position_values: np.ndarray, *, size: int) -> np.ndarray:
    positions = np.asarray(position_values, dtype=float)
    if positions.size == 0:
        return np.array([], dtype=int)

    rounded = np.rint(positions).astype(int)
    min_allowed = np.arange(positions.size, dtype=int)
    max_allowed = size - (positions.size - np.arange(positions.size, dtype=int))
    rounded = np.clip(rounded, min_allowed, max_allowed)

    for idx in range(1, rounded.size):
        rounded[idx] = max(rounded[idx], rounded[idx - 1] + 1)
    for idx in range(rounded.size - 2, -1, -1):
        rounded[idx] = min(rounded[idx], rounded[idx + 1] - 1)
        rounded[idx] = max(rounded[idx], idx)
    return rounded


def _build_partial_rox_step_assignments(
    expected_steps: np.ndarray,
    observed_times: np.ndarray,
    *,
    max_variants: int = ROX_PARTIAL_ALIGNMENT_VARIANTS,
) -> list[np.ndarray]:
    expected = np.asarray(expected_steps, dtype=float)
    times = np.asarray(observed_times, dtype=float)
    target_len = expected.size
    observed_len = times.size
    if observed_len == 0 or target_len == 0 or observed_len > target_len:
        return []
    if observed_len == target_len:
        return [expected.copy()]

    assignments: list[np.ndarray] = []
    seen: set[tuple[int, ...]] = set()

    def add_indices(indices: np.ndarray) -> None:
        key = tuple(int(value) for value in np.asarray(indices, dtype=int))
        if len(key) != observed_len or any(b <= a for a, b in zip(key, key[1:])):
            return
        if key in seen:
            return
        seen.add(key)
        assignments.append(expected[np.asarray(key, dtype=int)].copy())

    # Baseline contiguous windows remain useful for truly truncated ladders.
    max_start = max(0, target_len - observed_len)
    for start_idx in range(max_start + 1):
        add_indices(np.arange(start_idx, start_idx + observed_len, dtype=int))

    # Add sparse assignments that span wider ROX ranges so partial paths can
    # keep valid low-end and high-end peaks without forcing a contiguous bp window.
    span_pairs = [
        (0, target_len - 1),
        (0, max(observed_len - 1, target_len - 2)),
        (1, target_len - 1),
        (1, max(observed_len, target_len - 2)),
    ]
    if observed_len > 1 and times[-1] > times[0]:
        obs_norm = (times - times[0]) / max(times[-1] - times[0], 1.0)
    else:
        obs_norm = np.linspace(0.0, 1.0, observed_len)

    for start_idx, end_idx in span_pairs:
        start_idx = max(0, int(start_idx))
        end_idx = min(target_len - 1, int(end_idx))
        if end_idx - start_idx + 1 < observed_len:
            continue

        span_steps = expected[start_idx : end_idx + 1]
        if span_steps.size == observed_len:
            add_indices(np.arange(start_idx, end_idx + 1, dtype=int))
            continue

        if span_steps[-1] > span_steps[0]:
            step_norm = (span_steps - span_steps[0]) / max(span_steps[-1] - span_steps[0], 1.0)
        else:
            step_norm = np.linspace(0.0, 1.0, span_steps.size)
        approx_positions = np.interp(obs_norm, step_norm, np.arange(start_idx, end_idx + 1, dtype=float))
        add_indices(_round_to_monotonic_indices(approx_positions, size=target_len))

    return assignments[:max_variants]


def _select_best_bounded_ladder_fit(
    fsa: FsaFile,
    candidate_specs: list[dict[str, object]],
    *,
    rescue_mode: bool = False,
) -> FsaFile | None:
    if not candidate_specs:
        return None

    expected_steps = _get_expected_ladder_steps(fsa)
    best_complete_fit = None
    best_complete_score = None
    best_partial_fit = None
    best_partial_score = None
    evaluated_specs = 0
    matched_specs = 0
    search_start = time.perf_counter()

    for spec in candidate_specs:
        evaluated_specs += 1
        times = np.asarray(spec.get("times", []), dtype=float)
        ladder_steps = np.asarray(spec.get("ladder_steps", []), dtype=float)
        if times.size == 0 or ladder_steps.size == 0 or times.size != ladder_steps.size:
            continue

        trial = _clone_fsa_for_ladder_trial(fsa)
        trial.expected_ladder_steps = expected_steps.copy()
        trial.ladder_steps = ladder_steps
        trial.n_ladder_peaks = int(ladder_steps.size)
        trial.best_size_standard = times

        try:
            trial = fit_size_standard_to_ladder(trial)
        except Exception:
            continue
        if not getattr(trial, "fitted_to_model", False):
            continue
        matched_specs += 1

        metrics = compute_ladder_qc_metrics(trial)
        intensity_penalty = _candidate_intensity_penalty(trial)
        profile_penalty = _candidate_rox_profile_penalty(trial)
        missing_count = len(_missing_expected_ladder_steps(trial))
        score = _fit_score_tuple(
            metrics,
            intensity_penalty + profile_penalty,
            missing_penalty=_missing_step_penalty(trial),
        )
        used_bounded = bool(spec.get("bounded", False))
        strategy = "auto_full" if missing_count == 0 else "auto_partial"
        note = (
            f"Bounded ROX beam search selected a {'full' if missing_count == 0 else 'partial'} ladder fit "
            f"from explosive candidate space ({spec.get('beam_score', 0.0):.3f})."
            if used_bounded
            else None
        )
        trial = _set_ladder_fit_metadata(trial, strategy, note)

        if _is_early_accept_candidate(metrics, missing_count=missing_count):
            _log_ladder_timing(
                "ROX" if not rescue_mode else "ROX-RESCUE",
                "bounded candidate selection",
                Path(str(getattr(fsa, "file", "unknown.fsa"))),
                time.perf_counter() - search_start,
                candidates=evaluated_specs,
                fitted=matched_specs,
                complete=missing_count == 0,
                rescue=rescue_mode,
            )
            return trial
        if missing_count == 0:
            if best_complete_score is None or score < best_complete_score:
                best_complete_fit = trial
                best_complete_score = score
            continue

        if best_partial_score is None or score < best_partial_score:
            best_partial_fit = trial
            best_partial_score = score

    selected = best_complete_fit or best_partial_fit
    if selected is not None:
        _log_ladder_timing(
            "ROX" if not rescue_mode else "ROX-RESCUE",
            "bounded candidate selection",
            Path(str(getattr(fsa, "file", "unknown.fsa"))),
            time.perf_counter() - search_start,
            candidates=evaluated_specs,
            fitted=matched_specs,
            complete=best_complete_fit is not None,
            rescue=rescue_mode,
        )
    return selected


def _candidate_fit_score(fsa: FsaFile) -> tuple[float, float, float, float]:
    metrics = compute_ladder_qc_metrics(fsa)
    intensity_penalty = _candidate_intensity_penalty(fsa) + _candidate_rox_profile_penalty(fsa)
    return _fit_score_tuple(metrics, intensity_penalty, missing_penalty=_missing_step_penalty(fsa))


def _rescue_fit_score(fsa: FsaFile) -> tuple[float, float, float, float]:
    metrics = compute_ladder_qc_metrics(fsa)
    intensity_penalty = _candidate_intensity_penalty(fsa) + _candidate_rox_profile_penalty(fsa)
    return _fit_score_tuple(metrics, intensity_penalty, missing_penalty=_missing_step_penalty(fsa))


def _candidate_intensity_penalty(fsa: FsaFile) -> float:
    best = getattr(fsa, "best_size_standard", None)
    if best is None:
        return float("inf")

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return 0.0

    peak_idx = np.rint(np.asarray(best, dtype=float)).astype(int)
    valid = (peak_idx >= 0) & (peak_idx < trace.size)
    if not np.any(valid):
        return float("inf")

    intensities = trace[peak_idx[valid]]
    if intensities.size == 0:
        return float("inf")

    median_intensity = float(np.median(intensities))
    if median_intensity <= 0:
        return 0.0

    target_floor = median_intensity * MEDIAN_INTENSITY_TARGET_RATIO
    global_deficit = np.clip((target_floor - intensities) / median_intensity, a_min=0.0, a_max=None)

    early_count = max(1, int(np.ceil(len(intensities) * 0.25)))
    early_intensities = intensities[:early_count]
    early_deficit = np.clip((target_floor - early_intensities) / median_intensity, a_min=0.0, a_max=None)

    severe_weak_count = int(np.sum(intensities < (median_intensity * LOW_INTENSITY_RATIO_FLOOR)))

    return (
        float(np.sum(global_deficit)) * GLOBAL_PEAK_INTENSITY_WEIGHT
        + float(np.sum(early_deficit)) * EARLY_PEAK_INTENSITY_WEIGHT
        + (severe_weak_count * SEVERE_WEAK_PEAK_PENALTY)
    )


def _clean_rox_size_standard_peaks(all_found: np.ndarray, rox_data: np.ndarray) -> np.ndarray:
    if all_found is None or len(all_found) == 0:
        return np.array([], dtype=int)

    heights = np.array([rox_data[p] for p in all_found], dtype=float)
    median_h = float(np.median(heights)) if heights.size else 0.0
    cleaned: list[int] = []
    for idx, peak in enumerate(np.asarray(all_found, dtype=int)):
        height = float(rox_data[peak])
        if height > 31000 or peak < 1000:
            continue

        if median_h > 0 and height > (median_h * ROX_DYEBLOB_HEIGHT_MULTIPLIER):
            prev_gap = float("inf") if idx == 0 else peak - int(all_found[idx - 1])
            next_gap = float("inf") if idx == len(all_found) - 1 else int(all_found[idx + 1]) - peak
            crowded = min(prev_gap, next_gap) < ROX_DYEBLOB_TIGHT_GAP or (
                prev_gap < ROX_DYEBLOB_CLUSTER_GAP and next_gap < ROX_DYEBLOB_CLUSTER_GAP
            )
            if peak < ROX_DYEBLOB_EARLY_INDEX or crowded:
                continue

        cleaned.append(int(peak))

    return np.asarray(cleaned, dtype=int)


def _snap_peak_times_to_local_apexes(
    peak_times: np.ndarray,
    trace: np.ndarray,
    *,
    radius: int = ROX_APEX_SNAP_RADIUS,
) -> np.ndarray:
    """Snap candidate times to the nearest local apex in a small neighborhood."""
    peaks = np.asarray(peak_times, dtype=float)
    signal_trace = np.asarray(trace, dtype=float)
    if peaks.size == 0 or signal_trace.size == 0:
        return peaks

    snapped: list[float] = []
    max_index = signal_trace.size - 1
    for time_value in peaks:
        idx = int(np.clip(np.rint(time_value), 0, max_index))
        lo = max(0, idx - radius)
        hi = min(signal_trace.size, idx + radius + 1)
        window = signal_trace[lo:hi]
        if window.size == 0:
            snapped.append(float(idx))
            continue
        apex_value = float(np.max(window))
        apex_indices = np.flatnonzero(window == apex_value) + lo
        if apex_indices.size == 0:
            snapped.append(float(idx))
            continue
        best_idx = int(apex_indices[np.argmin(np.abs(apex_indices - idx))])
        snapped.append(float(best_idx))

    if not snapped:
        return np.array([], dtype=float)

    snapped_array = np.asarray(sorted(set(snapped)), dtype=float)
    return snapped_array


def _prepare_rox_size_standard_peaks(
    peak_times: np.ndarray,
    trace: np.ndarray,
    *,
    expected_count: int,
) -> np.ndarray:
    """Snap ROX peaks to apices and hard-trim to the human-good region when safe."""
    snapped = _snap_peak_times_to_local_apexes(peak_times, trace)
    if snapped.size == 0:
        return snapped

    preferred_mask = (
        (snapped >= ROX_HARD_FILTER_TIME_MIN)
        & (snapped <= ROX_HARD_FILTER_TIME_MAX)
    )
    preferred = snapped[preferred_mask]
    if preferred.size < max(1, int(expected_count)):
        lo = int(max(0, np.floor(ROX_HARD_FILTER_TIME_MIN)))
        hi = int(min(np.asarray(trace, dtype=float).size, np.ceil(ROX_HARD_FILTER_TIME_MAX) + 1))
        if hi > lo:
            supplemental_peaks, _ = signal.find_peaks(
                np.asarray(trace, dtype=float)[lo:hi],
                height=ROX_PREFERRED_SUPPLEMENT_MIN_HEIGHT,
                distance=ROX_PREFERRED_SUPPLEMENT_DISTANCE,
            )
            if supplemental_peaks.size > 0:
                merged = np.concatenate([snapped, supplemental_peaks.astype(float) + float(lo)])
                snapped = _snap_peak_times_to_local_apexes(merged, trace)
                preferred_mask = (
                    (snapped >= ROX_HARD_FILTER_TIME_MIN)
                    & (snapped <= ROX_HARD_FILTER_TIME_MAX)
                )
                preferred = snapped[preferred_mask]
    if preferred.size >= max(1, int(expected_count)):
        return preferred
    return snapped


def _supplement_rox_preferred_region_peaks(
    peak_times: np.ndarray,
    trace: np.ndarray,
    *,
    expected_count: int,
    min_distance: float,
) -> np.ndarray:
    """Recover moderate ROX peaks in the reviewed time window before later artifacts crowd them out."""
    peaks = np.asarray(peak_times, dtype=float)
    signal_trace = np.asarray(trace, dtype=float)
    if signal_trace.size == 0:
        return peaks

    preferred_mask = (
        (peaks >= ROX_HARD_FILTER_TIME_MIN)
        & (peaks <= ROX_HARD_FILTER_TIME_MAX)
    )
    if np.count_nonzero(preferred_mask) >= max(1, int(expected_count)):
        return peaks

    lo = int(max(0, np.floor(ROX_HARD_FILTER_TIME_MIN)))
    hi = int(min(signal_trace.size, np.ceil(ROX_HARD_FILTER_TIME_MAX) + 1))
    if hi <= lo:
        return peaks

    supplemental_peaks, props = signal.find_peaks(
        signal_trace[lo:hi],
        height=ROX_PREFERRED_SUPPLEMENT_MIN_HEIGHT,
        distance=max(1, int(round(float(min_distance)))),
    )
    if supplemental_peaks.size == 0:
        return peaks

    supplemental = supplemental_peaks.astype(float) + float(lo)
    if peaks.size == 0:
        return supplemental

    merged = list(peaks.tolist())
    for candidate in supplemental.tolist():
        if any(abs(float(existing) - float(candidate)) <= max(2.0, float(min_distance) * 0.35) for existing in merged):
            continue
        merged.append(float(candidate))
    if not merged:
        return np.array([], dtype=float)
    return np.asarray(sorted(merged), dtype=float)


def _rox_peak_time_penalty(time_value: float) -> float:
    if time_value < ROX_PREFERRED_TIME_MIN:
        return max(0.0, (ROX_PREFERRED_TIME_MIN - time_value) / 300.0)
    if time_value > ROX_PREFERRED_TIME_MAX:
        return max(0.0, (time_value - ROX_PREFERRED_TIME_MAX) / 300.0)
    if time_value < (ROX_PREFERRED_TIME_MIN + ROX_PREFERRED_TIME_MARGIN):
        return ((ROX_PREFERRED_TIME_MIN + ROX_PREFERRED_TIME_MARGIN) - time_value) / 900.0
    if time_value > (ROX_PREFERRED_TIME_MAX - ROX_PREFERRED_TIME_MARGIN):
        return (time_value - (ROX_PREFERRED_TIME_MAX - ROX_PREFERRED_TIME_MARGIN)) / 900.0
    return 0.0


def _candidate_rox_profile_penalty(fsa: FsaFile) -> float:
    best = getattr(fsa, "best_size_standard", None)
    if best is None:
        return 0.0
    ladder_name = str(getattr(fsa, "ladder", "") or "").upper()
    if "ROX" not in ladder_name:
        return 0.0

    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    if trace.size == 0:
        return 0.0

    peak_times = np.asarray(best, dtype=float)
    peak_idx = np.rint(peak_times).astype(int)
    valid = (peak_idx >= 0) & (peak_idx < trace.size)
    if not np.any(valid):
        return 0.0

    peak_times = peak_times[valid]
    intensities = trace[peak_idx[valid]]
    if intensities.size == 0:
        return 0.0

    time_penalty = sum(_rox_peak_time_penalty(float(time_value)) for time_value in peak_times)
    low_intensity = np.clip((250.0 - intensities) / 250.0, a_min=0.0, a_max=None)
    severe_weak = int(np.sum(intensities < ROX_PROFILE_SEVERE_WEAK_INTENSITY))

    return (
        float(time_penalty) * ROX_PROFILE_TIME_WEIGHT
        + float(np.sum(low_intensity)) * ROX_PROFILE_LOW_INTENSITY_WEIGHT
        + (severe_weak * ROX_PROFILE_SEVERE_WEAK_PENALTY)
    )


def _recover_rox_size_standard_peaks_from_baseline(fsa: FsaFile, raw_trace: np.ndarray) -> bool:
    """Retry ROX peak detection on a baseline-corrected trace when the raw pass fails."""
    corrected = np.maximum(np.asarray(raw_trace, dtype=float) - baseline_arPLS(raw_trace), 0.0)
    fallback_height = max(20.0, min(float(fsa.min_size_standard_height), ROX_BASELINE_FALLBACK_MIN_HEIGHT))
    found_peaks, _ = signal.find_peaks(
        corrected,
        height=fallback_height,
        distance=fsa.min_distance_between_peaks,
    )
    supplemented = _supplement_rox_preferred_region_peaks(
        np.asarray(found_peaks, dtype=float),
        corrected,
        expected_count=int(len(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float))),
        min_distance=float(getattr(fsa, "min_distance_between_peaks", 1.0) or 1.0),
    )
    cleaned = _clean_rox_size_standard_peaks(np.asarray(supplemented, dtype=int), corrected)
    if len(cleaned) < ROX_BASELINE_FALLBACK_MIN_PEAKS:
        return False

    expected_count = int(len(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)))
    prepared = _prepare_rox_size_standard_peaks(
        np.asarray(cleaned, dtype=float),
        corrected,
        expected_count=expected_count,
    )
    fsa.size_standard = corrected
    fsa.size_standard_peaks = np.asarray(prepared, dtype=float)
    fsa.size_standard_baseline_corrected = True
    return True

def _select_best_ladder_candidate(fsa: FsaFile, ranked_combinations: list[np.ndarray] | None = None) -> FsaFile | None:
    """Fit the top smooth candidates and keep the best actual ladder fit."""
    if ranked_combinations is None:
        ranked_combinations = _rank_size_standard_combinations(fsa)
    if not ranked_combinations:
        return None

    best_fit = None
    best_score = None

    for combo in ranked_combinations:
        trial = _clone_fsa_for_ladder_trial(fsa)
        trial.best_size_standard = combo
        try:
            trial = fit_size_standard_to_ladder(trial)
        except Exception:
            continue
        if not getattr(trial, "fitted_to_model", False):
            continue

        metrics = compute_ladder_qc_metrics(trial)
        intensity_penalty = _candidate_intensity_penalty(trial)
        profile_penalty = _candidate_rox_profile_penalty(trial)
        missing_count = len(_missing_expected_ladder_steps(trial))
        score = _fit_score_tuple(
            metrics,
            intensity_penalty + profile_penalty,
            missing_penalty=_missing_step_penalty(trial),
        )
        if _is_early_accept_candidate(metrics, missing_count=missing_count):
            return trial
        if best_score is None or score < best_score:
            best_score = score
            best_fit = trial

    return best_fit


def _build_rox_candidate_specs(
    fsa: FsaFile,
    *,
    label: str,
    fsa_path: Path,
    allow_partial: bool = True,
) -> tuple[list[dict[str, object]], bool, int]:
    combination_estimate = 0
    warned_bounded = False
    for _ in range(LADDER_MAX_ITERATIONS):
        combination_estimate = _estimate_size_standard_combination_count(fsa)
        use_bounded, combination_estimate = _should_use_bounded_rox_search(
            fsa,
            combination_estimate=combination_estimate,
        )
        if use_bounded:
            bounded_start = time.perf_counter()
            if not warned_bounded:
                peak_count = int(len(np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)))
                print_warning(
                    f"[{label}] Using bounded ladder beam search for {fsa_path.name} "
                    f"({peak_count} detected peaks, estimated {combination_estimate} combinations)."
                )
                warned_bounded = True
            specs = _build_bounded_rox_candidate_specs(fsa, allow_partial=allow_partial)
            _log_ladder_timing(
                label,
                "bounded ladder search",
                fsa_path,
                time.perf_counter() - bounded_start,
                peaks=int(len(np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float))),
                estimate=combination_estimate,
                specs=len(specs),
                partial=allow_partial,
            )
            return specs, True, combination_estimate

        fsa = generate_combinations(fsa)
        best = getattr(fsa, "best_size_standard_combinations", None)
        if best is not None and best.shape[0] > 0:
            break
        fsa.maxium_allowed_distance_between_size_standard_peaks += 10

    best = getattr(fsa, "best_size_standard_combinations", None)
    if best is None or best.shape[0] == 0:
        return [], False, combination_estimate

    specs = [
        {
            "times": np.asarray(combo, dtype=float),
            "ladder_steps": np.asarray(getattr(fsa, "ladder_steps", []), dtype=float).copy(),
            "beam_score": 0.0,
            "complete": True,
            "bounded": False,
        }
        for combo in _rank_size_standard_combinations(fsa)
    ]
    return specs, False, combination_estimate


def _should_attempt_high_end_rox_rescue(fsa: FsaFile, qc: dict[str, float | int]) -> bool:
    """Only attempt the expensive high-end rescue when the current fit is incomplete."""
    if float(qc.get("r2", float("-inf"))) >= HIGH_END_RESCUE_R2:
        return False
    return bool(_missing_expected_ladder_steps(fsa))


def _try_high_end_ladder_rescue(fsa: FsaFile, label: str, fsa_path: Path) -> FsaFile | None:
    rescue_start = time.perf_counter()
    full_steps = _get_expected_ladder_steps(fsa)
    if full_steps.size < 12:
        return None

    best_fit = None
    best_score = None
    max_skip = min(6, max(0, int(full_steps.size) - 8))
    if max_skip < 1:
        return None
    low_end_missing = _count_missing_low_end_steps(fsa)
    max_skip = min(max_skip, max(1, low_end_missing + 1))
    attempted_skips = 0
    bounded_attempts = 0

    for skip_low in range(1, max_skip + 1):
        attempted_skips += 1
        trial = _clone_fsa_for_ladder_trial(fsa)
        trial.expected_ladder_steps = full_steps.copy()
        trial.ladder_steps = np.asarray(full_steps[skip_low:], dtype=float)
        trial.n_ladder_peaks = trial.ladder_steps.size
        trial.max_peaks_allow_in_size_standard = trial.n_ladder_peaks + 15

        ss_peaks = getattr(trial, "size_standard_peaks", None)
        if ss_peaks is None or len(ss_peaks) < trial.n_ladder_peaks:
            continue

        try:
            trial = return_maxium_allowed_distance_between_size_standard_peaks(trial, multiplier=2)
            candidate_specs, used_bounded, _estimate = _build_rox_candidate_specs(
                trial,
                label=label,
                fsa_path=fsa_path,
                allow_partial=True,
            )
            if used_bounded:
                bounded_attempts += 1
            if not candidate_specs:
                continue

            selected_fit = _select_best_bounded_ladder_fit(trial, candidate_specs, rescue_mode=True)
            if selected_fit is None:
                if used_bounded:
                    continue
                trial = calculate_best_combination_of_size_standard_peaks(trial)
                trial = fit_size_standard_to_ladder(trial)
            else:
                trial = selected_fit

            if not getattr(trial, "fitted_to_model", False):
                continue

            metrics = compute_ladder_qc_metrics(trial)
            score = _rescue_fit_score(trial)
            if best_score is None or score < best_score:
                best_fit = trial
                best_score = score
            if _is_early_accept_candidate(metrics, missing_count=len(_missing_expected_ladder_steps(trial))):
                break
        except Exception:
            continue

    if best_fit is not None:
        kept = len(getattr(best_fit, "ladder_steps", []))
        total = len(getattr(best_fit, "expected_ladder_steps", getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", []))))
        best_fit = _set_ladder_fit_metadata(
            best_fit,
            "high_end_rescue",
            f"High-end rescue used the stable top {kept}/{total} ladder steps because the lower ROX region was unreliable.",
        )
    _log_ladder_timing(
        label,
        "high-end rescue",
        fsa_path,
        time.perf_counter() - rescue_start,
        skip_trials=attempted_skips,
        low_end_missing=low_end_missing,
        bounded_trials=bounded_attempts,
        rescued=best_fit is not None,
    )
    return best_fit


def _try_descending_low_end_completion(fsa: FsaFile, label: str, fsa_path: Path) -> FsaFile | None:
    expected = _get_expected_ladder_steps(fsa)
    current_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    current_times = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    candidate_times = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    ladder_model = getattr(fsa, "ladder_model", None)

    if (
        expected.size == 0
        or current_steps.size == 0
        or current_times.size != current_steps.size
        or candidate_times.size == 0
        or ladder_model is None
    ):
        return None

    full_times = np.full(expected.size, np.nan, dtype=float)
    step_map = _map_step_indices(current_steps, expected)
    for current_idx, full_idx in step_map.items():
        full_times[full_idx] = current_times[current_idx]

    missing_indices = [idx for idx, value in enumerate(full_times) if np.isnan(value)]
    if not missing_indices:
        return None

    xs = np.arange(trace.size, dtype=float)
    predicted_bp = np.asarray(ladder_model.predict(xs.reshape(-1, 1)), dtype=float)
    anchor_intensities = trace[np.rint(current_times).astype(int)]
    median_anchor_intensity = float(np.median(anchor_intensities)) if anchor_intensities.size else 0.0
    used_times = {round(float(t), 6) for t in current_times}
    added_steps: list[float] = []

    for step_idx in reversed(missing_indices):
        higher_indices = [idx for idx in range(step_idx + 1, expected.size) if not np.isnan(full_times[idx])]
        if not higher_indices:
            continue

        next_higher_idx = higher_indices[0]
        next_higher_time = float(full_times[next_higher_idx])
        target_bp = float(expected[step_idx])
        target_time = int(np.argmin(np.abs(predicted_bp - target_bp)))
        gap_to_next = max(18.0, abs(next_higher_time - target_time))
        search_radius = min(120.0, max(30.0, gap_to_next * 0.8))
        lo = max(0.0, target_time - search_radius)
        hi = min(next_higher_time - 1.0, target_time + search_radius)
        if hi <= lo:
            continue

        candidates_in_window: list[tuple[float, float]] = []
        for candidate_time in candidate_times:
            candidate_time = float(candidate_time)
            if round(candidate_time, 6) in used_times:
                continue
            if not (lo <= candidate_time <= hi):
                continue
            intensity = float(trace[int(round(candidate_time))])
            candidates_in_window.append((candidate_time, intensity))

        if not candidates_in_window:
            continue

        def candidate_score(item: tuple[float, float]) -> tuple[float, float, float]:
            candidate_time, intensity = item
            distance_penalty = abs(candidate_time - target_time)
            if median_anchor_intensity > 0:
                relative_intensity = intensity / median_anchor_intensity
            else:
                relative_intensity = 1.0
            weak_penalty = max(0.0, 0.22 - relative_intensity)
            return (
                distance_penalty,
                weak_penalty,
                -intensity,
            )

        chosen_time, chosen_intensity = min(candidates_in_window, key=candidate_score)
        if chosen_intensity < DESCENDING_RECOVERY_MIN_INTENSITY:
            continue

        full_times[step_idx] = chosen_time
        used_times.add(round(chosen_time, 6))
        added_steps.append(target_bp)

    if not added_steps:
        return None

    assigned_mask = ~np.isnan(full_times)
    assigned_times = full_times[assigned_mask]
    assigned_steps = expected[assigned_mask]
    if assigned_times.size < current_times.size or np.any(np.diff(assigned_times) <= 0):
        return None

    trial = _clone_fsa_for_ladder_trial(fsa)
    trial.expected_ladder_steps = expected.copy()
    trial.ladder_steps = np.asarray(assigned_steps, dtype=float)
    trial.best_size_standard = np.asarray(assigned_times, dtype=float)
    trial.n_ladder_peaks = trial.ladder_steps.size

    try:
        trial = fit_size_standard_to_ladder(trial)
        if not getattr(trial, "fitted_to_model", False):
            return None
        qc = compute_ladder_qc_metrics(trial)
    except Exception:
        return None

    if (
        qc["r2"] < DESCENDING_RECOVERY_R2_FLOOR
        or qc["max_abs_error_bp"] > DESCENDING_RECOVERY_MAX_ABS_ERROR
        or qc["mean_abs_error_bp"] > DESCENDING_RECOVERY_MEAN_ABS_ERROR
    ):
        return None

    missing_after = [
        float(bp) for bp in expected if not np.any(np.isclose(trial.ladder_steps, bp, atol=1e-6))
    ]
    added_text = ", ".join(f"{bp:.0f}" for bp in added_steps)
    if missing_after:
        note = (
            f"High-end rescue recovered lower ladder steps {added_text} bp using a descending search. "
            f"Remaining missing steps: {', '.join(f'{bp:.0f}' for bp in missing_after)} bp."
        )
    else:
        note = (
            f"High-end rescue recovered all lower ladder steps using a descending search "
            f"({added_text} bp)."
        )

    return _set_ladder_fit_metadata(trial, "high_end_rescue", note)


def _try_ascending_high_end_completion(fsa: FsaFile, label: str, fsa_path: Path) -> FsaFile | None:
    expected = _get_expected_ladder_steps(fsa)
    current_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    current_times = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    candidate_times = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    ladder_model = getattr(fsa, "ladder_model", None)

    if (
        expected.size == 0
        or current_steps.size == 0
        or current_times.size != current_steps.size
        or candidate_times.size == 0
        or ladder_model is None
    ):
        return None

    full_times = np.full(expected.size, np.nan, dtype=float)
    step_map = _map_step_indices(current_steps, expected)
    for current_idx, full_idx in step_map.items():
        full_times[full_idx] = current_times[current_idx]

    missing_indices = [idx for idx, value in enumerate(full_times) if np.isnan(value)]
    if not missing_indices:
        return None

    highest_present = max((idx for idx, value in enumerate(full_times) if not np.isnan(value)), default=-1)
    high_end_missing = [idx for idx in missing_indices if idx > highest_present]
    if not high_end_missing:
        return None

    xs = np.arange(trace.size, dtype=float)
    predicted_bp = np.asarray(ladder_model.predict(xs.reshape(-1, 1)), dtype=float)
    anchor_intensities = trace[np.rint(current_times).astype(int)]
    median_anchor_intensity = float(np.median(anchor_intensities)) if anchor_intensities.size else 0.0
    used_times = {round(float(t), 6) for t in current_times}
    added_steps: list[float] = []

    for step_idx in high_end_missing:
        lower_indices = [idx for idx in range(step_idx - 1, -1, -1) if not np.isnan(full_times[idx])]
        if not lower_indices:
            continue

        prev_idx = lower_indices[0]
        prev_time = float(full_times[prev_idx])
        target_bp = float(expected[step_idx])
        target_time = int(np.argmin(np.abs(predicted_bp - target_bp)))
        gap_from_prev = max(18.0, abs(target_time - prev_time))
        search_radius = min(140.0, max(35.0, gap_from_prev * 0.8))
        lo = max(prev_time + 1.0, target_time - search_radius)
        hi = min(float(trace.size - 1), target_time + search_radius)
        if hi <= lo:
            continue

        candidates_in_window: list[tuple[float, float]] = []
        for candidate_time in candidate_times:
            candidate_time = float(candidate_time)
            if round(candidate_time, 6) in used_times:
                continue
            if not (lo <= candidate_time <= hi):
                continue
            intensity = float(trace[int(round(candidate_time))])
            candidates_in_window.append((candidate_time, intensity))

        if not candidates_in_window:
            continue

        def candidate_score(item: tuple[float, float]) -> tuple[float, float, float]:
            candidate_time, intensity = item
            distance_penalty = abs(candidate_time - target_time)
            if median_anchor_intensity > 0:
                relative_intensity = intensity / median_anchor_intensity
            else:
                relative_intensity = 1.0
            weak_penalty = max(0.0, 0.30 - relative_intensity)
            return (
                distance_penalty,
                weak_penalty,
                -intensity,
            )

        chosen_time, chosen_intensity = min(candidates_in_window, key=candidate_score)
        if chosen_intensity < ASCENDING_RECOVERY_MIN_INTENSITY:
            continue

        full_times[step_idx] = chosen_time
        used_times.add(round(chosen_time, 6))
        added_steps.append(target_bp)

    if not added_steps:
        return None

    assigned_mask = ~np.isnan(full_times)
    assigned_times = full_times[assigned_mask]
    assigned_steps = expected[assigned_mask]
    if assigned_times.size < current_times.size or np.any(np.diff(assigned_times) <= 0):
        return None

    trial = _clone_fsa_for_ladder_trial(fsa)
    trial.expected_ladder_steps = expected.copy()
    trial.ladder_steps = np.asarray(assigned_steps, dtype=float)
    trial.best_size_standard = np.asarray(assigned_times, dtype=float)
    trial.n_ladder_peaks = trial.ladder_steps.size

    try:
        trial = fit_size_standard_to_ladder(trial)
        if not getattr(trial, "fitted_to_model", False):
            return None
        qc = compute_ladder_qc_metrics(trial)
    except Exception:
        return None

    if (
        qc["r2"] < ASCENDING_RECOVERY_R2_FLOOR
        or qc["max_abs_error_bp"] > ASCENDING_RECOVERY_MAX_ABS_ERROR
        or qc["mean_abs_error_bp"] > ASCENDING_RECOVERY_MEAN_ABS_ERROR
    ):
        return None

    note = (
        f"Ascending high-end completion recovered upper ladder steps "
        f"({', '.join(f'{bp:.0f}' for bp in added_steps)} bp)."
    )
    return _set_ladder_fit_metadata(trial, "high_end_rescue", note)


def _candidate_time_window_for_missing_step(
    full_times: np.ndarray,
    step_idx: int,
    target_time: float,
    trace_size: int,
) -> tuple[float, float]:
    lower_time = None
    upper_time = None
    for idx in range(step_idx - 1, -1, -1):
        if not np.isnan(full_times[idx]):
            lower_time = float(full_times[idx])
            break
    for idx in range(step_idx + 1, len(full_times)):
        if not np.isnan(full_times[idx]):
            upper_time = float(full_times[idx])
            break

    if lower_time is not None and upper_time is not None:
        left_gap = max(24.0, target_time - lower_time)
        right_gap = max(24.0, upper_time - target_time)
        lo = max(lower_time + 1.0, target_time - (left_gap * 0.85))
        hi = min(upper_time - 1.0, target_time + (right_gap * 0.85))
    elif lower_time is not None:
        gap = max(35.0, target_time - lower_time)
        lo = max(lower_time + 1.0, target_time - (gap * 0.60))
        hi = min(float(trace_size - 1), target_time + min(150.0, gap * 0.90))
    elif upper_time is not None:
        gap = max(35.0, upper_time - target_time)
        lo = max(0.0, target_time - min(150.0, gap * 0.90))
        hi = min(upper_time - 1.0, target_time + (gap * 0.60))
    else:
        lo = max(0.0, target_time - 120.0)
        hi = min(float(trace_size - 1), target_time + 120.0)
    return lo, hi


def _estimate_missing_step_time_from_assigned(
    expected_steps: np.ndarray,
    full_times: np.ndarray,
    step_idx: int,
    fallback_time: float,
) -> float:
    assigned = [
        (idx, float(expected_steps[idx]), float(full_times[idx]))
        for idx in range(len(expected_steps))
        if not np.isnan(full_times[idx])
    ]
    if len(assigned) < 2:
        return float(fallback_time)

    lower = [item for item in assigned if item[0] < step_idx]
    upper = [item for item in assigned if item[0] > step_idx]
    target_bp = float(expected_steps[step_idx])

    if lower and upper:
        left_idx, left_bp, left_time = lower[-1]
        right_idx, right_bp, right_time = upper[0]
        if right_bp > left_bp and right_time > left_time:
            ratio = (target_bp - left_bp) / max(right_bp - left_bp, 1.0)
            return float(left_time + (ratio * (right_time - left_time)))

    if len(lower) >= 2:
        left0 = lower[-2]
        left1 = lower[-1]
        bp_delta = left1[1] - left0[1]
        time_delta = left1[2] - left0[2]
        if bp_delta > 0 and time_delta > 0:
            slope = time_delta / bp_delta
            return float(left1[2] + ((target_bp - left1[1]) * slope))

    if len(upper) >= 2:
        right0 = upper[0]
        right1 = upper[1]
        bp_delta = right1[1] - right0[1]
        time_delta = right1[2] - right0[2]
        if bp_delta > 0 and time_delta > 0:
            slope = time_delta / bp_delta
            return float(right0[2] - ((right0[1] - target_bp) * slope))

    return float(fallback_time)


def _try_complete_missing_steps_by_prediction(fsa: FsaFile, label: str, fsa_path: Path) -> FsaFile | None:
    expected = _get_expected_ladder_steps(fsa)
    current_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    current_times = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    candidate_times = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    ladder_model = getattr(fsa, "ladder_model", None)

    if (
        expected.size == 0
        or current_steps.size == 0
        or current_times.size != current_steps.size
        or candidate_times.size == 0
        or ladder_model is None
        or trace.size == 0
    ):
        return None

    full_times = np.full(expected.size, np.nan, dtype=float)
    step_map = _map_step_indices(current_steps, expected)
    for current_idx, full_idx in step_map.items():
        full_times[full_idx] = current_times[current_idx]

    missing_indices = [idx for idx, value in enumerate(full_times) if np.isnan(value)]
    if not missing_indices:
        return None

    xs = np.arange(trace.size, dtype=float)
    predicted_bp = np.asarray(ladder_model.predict(xs.reshape(-1, 1)), dtype=float)
    anchor_intensities = trace[np.rint(current_times).astype(int)]
    median_anchor_intensity = float(np.median(anchor_intensities)) if anchor_intensities.size else 0.0
    used_times = {round(float(t), 6) for t in current_times}
    added_steps: list[float] = []

    for step_idx in missing_indices:
        target_bp = float(expected[step_idx])
        fallback_time = float(int(np.argmin(np.abs(predicted_bp - target_bp))))
        target_time = _estimate_missing_step_time_from_assigned(
            expected,
            full_times,
            step_idx,
            fallback_time,
        )
        lo, hi = _candidate_time_window_for_missing_step(full_times, step_idx, target_time, trace.size)
        if hi <= lo:
            continue

        candidates_in_window: list[tuple[float, float]] = []
        for candidate_time in candidate_times:
            candidate_time = float(candidate_time)
            if round(candidate_time, 6) in used_times:
                continue
            if not (lo <= candidate_time <= hi):
                continue
            intensity = float(trace[int(round(candidate_time))])
            candidates_in_window.append((candidate_time, intensity))

        if not candidates_in_window:
            continue

        def candidate_score(item: tuple[float, float]) -> tuple[float, float, float]:
            candidate_time, intensity = item
            distance_penalty = abs(candidate_time - target_time)
            if median_anchor_intensity > 0:
                relative_intensity = intensity / median_anchor_intensity
            else:
                relative_intensity = 1.0
            weak_penalty = max(0.0, 0.28 - relative_intensity)
            return (
                distance_penalty,
                weak_penalty,
                -intensity,
            )

        chosen_time, chosen_intensity = min(candidates_in_window, key=candidate_score)
        if chosen_intensity < GENERAL_COMPLETION_MIN_INTENSITY:
            continue

        full_times[step_idx] = chosen_time
        used_times.add(round(chosen_time, 6))
        added_steps.append(target_bp)

    if not added_steps:
        return None

    assigned_mask = ~np.isnan(full_times)
    assigned_times = full_times[assigned_mask]
    assigned_steps = expected[assigned_mask]
    if assigned_times.size < current_times.size or np.any(np.diff(assigned_times) <= 0):
        return None

    trial = _clone_fsa_for_ladder_trial(fsa)
    trial.expected_ladder_steps = expected.copy()
    trial.ladder_steps = np.asarray(assigned_steps, dtype=float)
    trial.best_size_standard = np.asarray(assigned_times, dtype=float)
    trial.n_ladder_peaks = trial.ladder_steps.size

    try:
        trial = fit_size_standard_to_ladder(trial)
        if not getattr(trial, "fitted_to_model", False):
            return None
        qc = compute_ladder_qc_metrics(trial)
    except Exception:
        return None

    if (
        qc["r2"] < GENERAL_COMPLETION_R2_FLOOR
        or qc["max_abs_error_bp"] > GENERAL_COMPLETION_MAX_ABS_ERROR
        or qc["mean_abs_error_bp"] > GENERAL_COMPLETION_MEAN_ABS_ERROR
    ):
        return None

    note = (
        f"Predicted-step completion recovered missing ladder steps "
        f"({', '.join(f'{bp:.0f}' for bp in added_steps)} bp)."
    )
    return _set_ladder_fit_metadata(trial, "high_end_rescue", note)


def _try_core_anchored_step_completion(fsa: FsaFile, label: str, fsa_path: Path) -> FsaFile | None:
    expected = _get_expected_ladder_steps(fsa)
    current_steps = np.asarray(getattr(fsa, "ladder_steps", []), dtype=float)
    current_times = np.asarray(getattr(fsa, "best_size_standard", []), dtype=float)
    if current_steps.size < CORE_COMPLETION_MIN_ASSIGNED or current_times.size != current_steps.size:
        return None

    full_times = np.full(expected.size, np.nan, dtype=float)
    step_map = _map_step_indices(current_steps, expected)
    for current_idx, full_idx in step_map.items():
        full_times[full_idx] = current_times[current_idx]

    assigned_indices = [idx for idx, value in enumerate(full_times) if not np.isnan(value)]
    if len(assigned_indices) < CORE_COMPLETION_MIN_ASSIGNED:
        return None

    # Anchor around the longest assigned contiguous run and then grow outward.
    runs: list[list[int]] = []
    current_run: list[int] = []
    for idx in assigned_indices:
        if not current_run or idx == current_run[-1] + 1:
            current_run.append(idx)
        else:
            runs.append(current_run)
            current_run = [idx]
    if current_run:
        runs.append(current_run)
    if not runs:
        return None

    core_run = max(runs, key=len)
    core_start = core_run[0]
    core_end = core_run[-1]
    if len(core_run) < 8:
        return None

    candidate_times = np.asarray(getattr(fsa, "size_standard_peaks", []), dtype=float)
    trace = np.asarray(getattr(fsa, "size_standard", []), dtype=float)
    ladder_model = getattr(fsa, "ladder_model", None)
    if candidate_times.size == 0 or trace.size == 0 or ladder_model is None:
        return None

    xs = np.arange(trace.size, dtype=float)
    predicted_bp = np.asarray(ladder_model.predict(xs.reshape(-1, 1)), dtype=float)
    anchor_intensities = trace[np.rint(current_times).astype(int)]
    median_anchor_intensity = float(np.median(anchor_intensities)) if anchor_intensities.size else 0.0
    used_times = {round(float(t), 6) for t in current_times}
    added_steps: list[float] = []

    expansion_order = list(range(core_start - 1, -1, -1)) + list(range(core_end + 1, len(expected)))
    for step_idx in expansion_order:
        if not np.isnan(full_times[step_idx]):
            continue
        target_bp = float(expected[step_idx])
        fallback_time = float(int(np.argmin(np.abs(predicted_bp - target_bp))))
        target_time = _estimate_missing_step_time_from_assigned(expected, full_times, step_idx, fallback_time)
        lo, hi = _candidate_time_window_for_missing_step(full_times, step_idx, target_time, trace.size)
        if hi <= lo:
            continue

        candidates_in_window: list[tuple[float, float]] = []
        for candidate_time in candidate_times:
            candidate_time = float(candidate_time)
            if round(candidate_time, 6) in used_times:
                continue
            if not (lo <= candidate_time <= hi):
                continue
            intensity = float(trace[int(round(candidate_time))])
            candidates_in_window.append((candidate_time, intensity))

        if not candidates_in_window:
            continue

        def candidate_score(item: tuple[float, float]) -> tuple[float, float, float]:
            candidate_time, intensity = item
            distance_penalty = abs(candidate_time - target_time)
            if median_anchor_intensity > 0:
                relative_intensity = intensity / median_anchor_intensity
            else:
                relative_intensity = 1.0
            weak_penalty = max(0.0, 0.30 - relative_intensity)
            return (distance_penalty, weak_penalty, -intensity)

        chosen_time, chosen_intensity = min(candidates_in_window, key=candidate_score)
        if chosen_intensity < GENERAL_COMPLETION_MIN_INTENSITY:
            continue

        full_times[step_idx] = chosen_time
        used_times.add(round(chosen_time, 6))
        added_steps.append(target_bp)

    if not added_steps:
        return None

    assigned_mask = ~np.isnan(full_times)
    assigned_times = full_times[assigned_mask]
    assigned_steps = expected[assigned_mask]
    if assigned_times.size < current_times.size or np.any(np.diff(assigned_times) <= 0):
        return None

    trial = _clone_fsa_for_ladder_trial(fsa)
    trial.expected_ladder_steps = expected.copy()
    trial.ladder_steps = np.asarray(assigned_steps, dtype=float)
    trial.best_size_standard = np.asarray(assigned_times, dtype=float)
    trial.n_ladder_peaks = trial.ladder_steps.size

    try:
        trial = fit_size_standard_to_ladder(trial)
        if not getattr(trial, "fitted_to_model", False):
            return None
        qc = compute_ladder_qc_metrics(trial)
    except Exception:
        return None

    if (
        qc["r2"] < GENERAL_COMPLETION_R2_FLOOR
        or qc["max_abs_error_bp"] > GENERAL_COMPLETION_MAX_ABS_ERROR
        or qc["mean_abs_error_bp"] > GENERAL_COMPLETION_MEAN_ABS_ERROR
    ):
        return None

    note = (
        f"Core-anchored completion recovered missing ladder steps "
        f"({', '.join(f'{bp:.0f}' for bp in added_steps)} bp)."
    )
    return _set_ladder_fit_metadata(trial, "high_end_rescue", note)


# ==================================================================
# ==================== ANALYSEFUNKSJONER ===========================
# ==================================================================

def _normalize_ladder_adjustment_payload(adjustment: dict | None) -> dict | None:
    """Normalizes legacy and enriched ladder adjustment payloads."""
    if not adjustment:
        return None

    if "mapping" in adjustment or "mapping_times" in adjustment or "manual_candidates" in adjustment:
        mapping_raw = adjustment.get("mapping", {})
        mapping_times_raw = adjustment.get("mapping_times", {})
        manual_candidates_raw = adjustment.get("manual_candidates", [])
        return {
            "mapping": {int(k): int(v) for k, v in mapping_raw.items()},
            "mapping_times": {int(k): float(v) for k, v in mapping_times_raw.items()},
            "manual_candidates": [float(v) for v in manual_candidates_raw],
        }

    return {
        "mapping": {int(k): int(v) for k, v in adjustment.items()},
        "mapping_times": {},
        "manual_candidates": [],
    }


def save_ladder_adjustment(
    fsa: FsaFile,
    adjustment: dict[int, int] | dict,
    *,
    manual_candidates: list[float] | None = None,
    mapping_times: dict[int, float] | None = None,
) -> None:
    """Saves a manual mapping payload to a .json file alongside the .fsa file."""
    adj_path = Path(fsa.file).resolve().with_suffix(".ladder_adj.json")
    try:
        if manual_candidates is not None or mapping_times is not None:
            payload = {
                "mapping": {int(k): int(v) for k, v in adjustment.items()},
                "mapping_times": {int(k): float(v) for k, v in (mapping_times or {}).items()},
                "manual_candidates": [float(v) for v in (manual_candidates or [])],
            }
        else:
            payload = _normalize_ladder_adjustment_payload(adjustment) or {
                "mapping": {},
                "mapping_times": {},
                "manual_candidates": [],
            }
        with open(adj_path, "w") as f:
            json.dump(payload, f)
        print_green(f"Saved ladder adjustment to {adj_path.name}")
    except Exception as e:
        print_warning(f"Could not save ladder adjustment: {e}")


def load_ladder_adjustment(fsa: FsaFile) -> dict | None:
    """Loads a manual mapping payload from a .json file if it exists."""
    candidate_files: list[Path] = [Path(fsa.file)]
    try:
        resolved = Path(fsa.file).resolve()
    except Exception:
        resolved = None
    if resolved is not None and resolved not in candidate_files:
        candidate_files.append(resolved)

    for candidate_file in candidate_files:
        adj_path = candidate_file.with_suffix(".ladder_adj.json")
        if not adj_path.exists():
            continue
        try:
            with open(adj_path, "r", encoding="utf-8", errors="replace") as f:
                payload = json.load(f)
                return _normalize_ladder_adjustment_payload(payload)
        except Exception as e:
            print_warning(f"Could not load ladder adjustment {adj_path.name}: {e}")
    return None


def _try_apply_saved_ladder_adjustment(fsa: FsaFile, adjustment: dict | None, label: str) -> FsaFile | None:
    """Applies a saved ladder adjustment if valid, otherwise warns and falls back to auto-fit."""
    if not adjustment:
        return None
    try:
        print_green(f"[{label}] Applying manual ladder adjustment for {fsa.file_name}")
        return _set_ladder_fit_metadata(
            apply_manual_ladder_mapping(fsa, adjustment),
            "manual_adjustment",
            "Manual ladder adjustment applied from saved sidecar.",
        )
    except Exception as exc:
        print_warning(
            f"[{label}] Ignoring invalid saved ladder adjustment for {fsa.file_name}: {exc}. Falling back to auto-fit."
        )
        return None


def analyse_fsa_liz(
    fsa_path: Path,
    sample_channel: str,
    *,
    ladder_name: str | None = None,
    min_distance_between_peaks: float | None = None,
    min_size_standard_height: float | None = None,
) -> FsaFile | None:
    """Ladder-fit for LIZ (TCRg/IGK/KDE).
    
    Uses multi-config search, dye-blob filtering, polynomial sizing,
    and iterative outlier removal for robust ladder fitting.
    """
    ladder_name = ladder_name or LIZ_LADDER
    base_min_distance = float(
        MIN_DISTANCE_BETWEEN_PEAKS_LIZ if min_distance_between_peaks is None else min_distance_between_peaks
    )
    base_min_height = float(
        MIN_SIZE_STANDARD_HEIGHT_LIZ if min_size_standard_height is None else min_size_standard_height
    )

    print_green(
        f"=== Analysing {fsa_path} ({ladder_name}, sample {sample_channel}, Python API) ==="
    )

    configs = [
        {"min_h": base_min_height, "min_d": base_min_distance},
        {"min_h": 200, "min_d": 20},
        {"min_h": 100, "min_d": 15},
        {"min_h": 50, "min_d": 10},
    ]

    base_fsa = FsaFile(
        file=str(fsa_path),
        ladder=ladder_name,
        sample_channel=sample_channel,
        min_distance_between_peaks=configs[0]["min_d"],
        min_size_standard_height=configs[0]["min_h"],
        size_standard_channel="DATA105",
    )
    base_fsa = find_size_standard_peaks(base_fsa)
    
    applied = _try_apply_saved_ladder_adjustment(base_fsa, load_ladder_adjustment(base_fsa), "LIZ")
    if applied is not None:
        return applied

    best_fallback_fsa = None

    for cfg in configs:
        fsa = FsaFile(
            file=str(fsa_path),
            ladder=ladder_name,
            sample_channel=sample_channel,
            min_distance_between_peaks=cfg["min_d"],
            min_size_standard_height=cfg["min_h"],
            size_standard_channel="DATA105",
        )
        liz_data = np.asarray(fsa.fsa["DATA105"]).astype(float)
        fsa = find_size_standard_peaks(fsa)

        all_found = getattr(fsa, "size_standard_peaks", None)
        if all_found is not None:
            # Dye-blob detection via median height
            heights = np.array([liz_data[p] for p in all_found])
            median_h = np.median(heights)
            cleaned = []
            for p in all_found:
                h = liz_data[p]
                if h > 31000 or p < 1000 or h > 3.0 * median_h:
                    continue
                cleaned.append(p)
            
            # Fall-through Logic: Only apply cleaned if we have a reasonable amount of peaks left
            expected_steps = len(getattr(fsa, "expected_ladder_steps", []))
            min_required = max(10, int(expected_steps * 0.6)) if expected_steps > 0 else 10
            if len(cleaned) >= min_required:
                fsa.size_standard_peaks = np.array(cleaned)

        ss_peaks = getattr(fsa, "size_standard_peaks", None)
        if ss_peaks is None or getattr(ss_peaks, "shape", [0])[0] < 2:
            continue

        try:
            fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=2)
            for _ in range(LADDER_MAX_ITERATIONS):
                fsa = generate_combinations(fsa)
                best = getattr(fsa, "best_size_standard_combinations", None)
                if best is not None and best.shape[0] > 0:
                    break
                fsa.maxium_allowed_distance_between_size_standard_peaks += 10
                
            best = getattr(fsa, "best_size_standard_combinations", None)
            if best is None or best.shape[0] == 0:
                continue
                
            selected_fit = _select_best_ladder_candidate(fsa)
            if selected_fit is not None:
                fsa = selected_fit
            else:
                fsa = calculate_best_combination_of_size_standard_peaks(fsa)
            if best_fallback_fsa is None:
                best_fallback_fsa = _clone_fsa_for_ladder_trial(fsa)
            
            try:
                if not getattr(fsa, "fitted_to_model", False):
                    fsa = fit_size_standard_to_ladder(fsa)
                if getattr(fsa, "fitted_to_model", False):
                    qc = compute_ladder_qc_metrics(fsa)
                    if qc["r2"] >= 0.9995 and not _missing_expected_ladder_steps(fsa):
                        return _finalize_auto_fit_metadata(fsa)
                    if _should_attempt_high_end_rox_rescue(fsa, qc):
                        rescued = _try_high_end_ladder_rescue(fsa, "LIZ", fsa_path)
                        if rescued is not None and _rescue_fit_score(rescued) < _rescue_fit_score(fsa):
                            kept = len(getattr(rescued, "ladder_steps", []))
                            total = len(getattr(rescued, "expected_ladder_steps", getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", []))))
                            print_green(
                                f"[LIZ] High-end ladder rescue selected for {fsa_path.name} using the top {kept}/{total} ladder steps."
                            )
                            fsa = rescued
                            qc = compute_ladder_qc_metrics(fsa)
                    return _finalize_auto_fit_metadata(fsa)
            except ValueError:
                pass
        except ValueError:
            continue

    if best_fallback_fsa is not None:
        return _finalize_auto_fit_metadata(best_fallback_fsa)

    print_warning(f"[LIZ] Fant ingen gyldige size-standard kombinasjoner for {fsa_path.name}")
    return None

def analyse_fsa_rox(
    fsa_path: Path,
    sample_channel: str,
    *,
    ladder_name: str | None = None,
    min_distance_between_peaks: float | None = None,
    min_size_standard_height: float | None = None,
) -> FsaFile | None:
    """Ladder-fit for ROX (FR1–3, TCRbA/B/C, SL, DHJH_D/E).
    
    Uses multi-config search, dye-blob filtering, polynomial sizing,
    and iterative outlier removal for robust ladder fitting.
    """
    ladder_name = ladder_name or ROX_LADDER
    base_min_distance = float(
        MIN_DISTANCE_BETWEEN_PEAKS_ROX if min_distance_between_peaks is None else min_distance_between_peaks
    )
    base_min_height = float(
        MIN_SIZE_STANDARD_HEIGHT_ROX if min_size_standard_height is None else min_size_standard_height
    )

    print_green(
        f"=== Analysing {fsa_path} ({ladder_name}, sample {sample_channel}, Python API) ==="
    )

    configs = [
        {"min_h": base_min_height, "min_d": base_min_distance},
        {"min_h": 100, "min_d": 15},
        {"min_h": 50, "min_d": 10},
        {"min_h": 20, "min_d": 8},
    ]

    base_fsa = FsaFile(
        file=str(fsa_path),
        ladder=ladder_name,
        sample_channel=sample_channel,
        min_distance_between_peaks=configs[0]["min_d"],
        min_size_standard_height=configs[0]["min_h"],
        size_standard_channel="DATA4",
    )
    base_fsa = find_size_standard_peaks(base_fsa)
    base_raw_rox = np.asarray(base_fsa.fsa["DATA4"], dtype=float)
    base_found = np.asarray(getattr(base_fsa, "size_standard_peaks", []), dtype=float)
    base_supplemented = _supplement_rox_preferred_region_peaks(
        base_found,
        base_raw_rox,
        expected_count=int(len(np.asarray(getattr(base_fsa, "ladder_steps", []), dtype=float))),
        min_distance=float(getattr(base_fsa, "min_distance_between_peaks", 1.0) or 1.0),
    )
    base_cleaned = _clean_rox_size_standard_peaks(
        np.asarray(base_supplemented, dtype=int),
        base_raw_rox,
    )
    if len(base_cleaned) >= ROX_BASELINE_FALLBACK_MIN_PEAKS:
        base_fsa.size_standard_peaks = _prepare_rox_size_standard_peaks(
            np.asarray(base_cleaned, dtype=float),
            base_raw_rox,
            expected_count=int(len(np.asarray(getattr(base_fsa, "ladder_steps", []), dtype=float))),
        )
    else:
        _recover_rox_size_standard_peaks_from_baseline(base_fsa, base_raw_rox)
    
    applied = _try_apply_saved_ladder_adjustment(base_fsa, load_ladder_adjustment(base_fsa), "ROX")
    if applied is not None:
        return applied

    best_fallback_fsa = None

    for cfg in configs:
        fsa = FsaFile(
            file=str(fsa_path),
            ladder=ladder_name,
            sample_channel=sample_channel,
            min_distance_between_peaks=cfg["min_d"],
            min_size_standard_height=cfg["min_h"],
            size_standard_channel="DATA4",
        )
        rox_data = np.asarray(fsa.fsa["DATA4"]).astype(float)
        fsa = find_size_standard_peaks(fsa)

        all_found = getattr(fsa, "size_standard_peaks", None)
        if all_found is not None:
            supplemented = _supplement_rox_preferred_region_peaks(
                np.asarray(all_found, dtype=float),
                rox_data,
                expected_count=int(len(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float))),
                min_distance=float(getattr(fsa, "min_distance_between_peaks", 1.0) or 1.0),
            )
            cleaned = _clean_rox_size_standard_peaks(np.asarray(supplemented, dtype=int), rox_data)
            if len(cleaned) >= ROX_BASELINE_FALLBACK_MIN_PEAKS:
                fsa.size_standard_peaks = _prepare_rox_size_standard_peaks(
                    np.asarray(cleaned, dtype=float),
                    rox_data,
                    expected_count=int(len(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float))),
                )
            elif _recover_rox_size_standard_peaks_from_baseline(fsa, rox_data):
                print_green(f"[ROX] Baseline-corrected ladder detection used for {fsa_path.name}")
        elif _recover_rox_size_standard_peaks_from_baseline(fsa, rox_data):
            print_green(f"[ROX] Baseline-corrected ladder detection used for {fsa_path.name}")

        ss_peaks = getattr(fsa, "size_standard_peaks", None)
        if ss_peaks is None or getattr(ss_peaks, "shape", [0])[0] < 2:
            continue

        try:
            fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=2)
            candidate_specs, used_bounded, _estimate = _build_rox_candidate_specs(
                fsa,
                label="ROX",
                fsa_path=fsa_path,
                allow_partial=True,
            )
            if not candidate_specs:
                if not used_bounded and getattr(getattr(fsa, "best_size_standard_combinations", None), "shape", [0])[0] > 0:
                    fsa = calculate_best_combination_of_size_standard_peaks(fsa)
                else:
                    continue

            selected_fit = _select_best_bounded_ladder_fit(fsa, candidate_specs, rescue_mode=False) if candidate_specs else None
            if selected_fit is not None:
                fsa = selected_fit
            elif candidate_specs:
                if used_bounded:
                    continue
                fsa = calculate_best_combination_of_size_standard_peaks(fsa)
            if best_fallback_fsa is None:
                best_fallback_fsa = _clone_fsa_for_ladder_trial(fsa)
            
            try:
                if not getattr(fsa, "fitted_to_model", False):
                    fsa = fit_size_standard_to_ladder(fsa)
                if getattr(fsa, "fitted_to_model", False):
                    qc = compute_ladder_qc_metrics(fsa)
                    if qc["r2"] >= 0.9995:
                        return _finalize_auto_fit_metadata(fsa)
                    if _should_attempt_high_end_rox_rescue(fsa, qc):
                        rescued = _try_high_end_ladder_rescue(fsa, "ROX", fsa_path)
                        if rescued is not None and _rescue_fit_score(rescued) < _rescue_fit_score(fsa):
                            kept = len(getattr(rescued, "ladder_steps", []))
                            total = len(getattr(rescued, "expected_ladder_steps", getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", []))))
                            print_green(
                                f"[ROX] High-end ladder rescue selected for {fsa_path.name} using the top {kept}/{total} ladder steps."
                            )
                            fsa = rescued
                            qc = compute_ladder_qc_metrics(fsa)
                            completed = _try_descending_low_end_completion(fsa, "ROX", fsa_path)
                            if completed is not None:
                                rescued_score = _rescue_fit_score(fsa)
                                completed_score = _rescue_fit_score(completed)
                                rescued_steps = len(getattr(fsa, "ladder_steps", []))
                                completed_steps = len(getattr(completed, "ladder_steps", []))
                                if completed_steps > rescued_steps or (
                                    completed_steps == rescued_steps and completed_score < rescued_score
                                ):
                                    added_steps = [
                                        float(bp)
                                        for bp in np.asarray(getattr(completed, "ladder_steps", []), dtype=float)
                                        if not np.any(np.isclose(np.asarray(getattr(fsa, "ladder_steps", []), dtype=float), bp, atol=1e-6))
                                    ]
                                    if added_steps:
                                        print_green(
                                            f"[ROX] Descending low-end recovery accepted for {fsa_path.name}: "
                                            f"{', '.join(f'{bp:.0f}' for bp in added_steps)} bp"
                                        )
                                    fsa = completed
                                    qc = compute_ladder_qc_metrics(fsa)
                    high_completed = _try_ascending_high_end_completion(fsa, "ROX", fsa_path)
                    if high_completed is not None:
                        current_score = _rescue_fit_score(fsa)
                        completed_score = _rescue_fit_score(high_completed)
                        current_steps = len(getattr(fsa, "ladder_steps", []))
                        completed_steps = len(getattr(high_completed, "ladder_steps", []))
                        if completed_steps > current_steps or (
                            completed_steps == current_steps and completed_score < current_score
                        ):
                            fsa = high_completed
                    general_completed = _try_complete_missing_steps_by_prediction(fsa, "ROX", fsa_path)
                    if general_completed is not None:
                        current_score = _rescue_fit_score(fsa)
                        completed_score = _rescue_fit_score(general_completed)
                        current_steps = len(getattr(fsa, "ladder_steps", []))
                        completed_steps = len(getattr(general_completed, "ladder_steps", []))
                        if completed_steps > current_steps or (
                            completed_steps == current_steps and completed_score < current_score
                        ):
                            fsa = general_completed
                    core_completed = _try_core_anchored_step_completion(fsa, "ROX", fsa_path)
                    if core_completed is not None:
                        current_score = _rescue_fit_score(fsa)
                        completed_score = _rescue_fit_score(core_completed)
                        current_steps = len(getattr(fsa, "ladder_steps", []))
                        completed_steps = len(getattr(core_completed, "ladder_steps", []))
                        if completed_steps > current_steps or (
                            completed_steps == current_steps and completed_score < current_score
                        ):
                            fsa = core_completed
                    return _finalize_auto_fit_metadata(fsa)
            except ValueError:
                pass
        except ValueError:
            continue

    if best_fallback_fsa is not None:
        return _finalize_auto_fit_metadata(best_fallback_fsa)

    print_warning(f"[ROX] Fant ingen gyldige size-standard kombinasjoner for {fsa_path.name}")
    return None




# ==================================================================
# ===================== EGEN PEAK-DETEKTOR =========================
# ==================================================================

def _find_local_maxima(y: np.ndarray) -> np.ndarray:
    """Enkel lokal maks-deteksjon."""
    if y.size < 3:
        return np.array([], dtype=int)
    idx = np.arange(1, y.size - 1)
    left = y[idx - 1]
    mid = y[idx]
    right = y[idx + 1]
    mask = (mid > left) & (mid >= right)
    return idx[mask]


def estimate_running_baseline(
    trace: np.ndarray,
    bin_size: int = BASELINE_BIN_SIZE,
    quantile: float = BASELINE_QUANTILE,
) -> np.ndarray:
    """Rask tilnærming til rullende baseline (O(N), for plotting)."""
    n = trace.size
    if n == 0:
        return np.zeros_like(trace, dtype=float)
    if bin_size < 20:
        bin_size = 20

    n_bins = int(np.ceil(n / bin_size))
    centers = []
    base_vals = []

    for b in range(n_bins):
        start = b * bin_size
        end = min((b + 1) * bin_size, n)
        seg = trace[start:end]
        if seg.size == 0:
            continue
        centers.append(0.5 * (start + end - 1))
        base_vals.append(np.quantile(seg, quantile))

    centers = np.asarray(centers, dtype=float)
    base_vals = np.asarray(base_vals, dtype=float)

    if centers.size == 0:
        return np.zeros_like(trace, dtype=float)

    idx = np.arange(n, dtype=float)
    baseline = np.interp(idx, centers, base_vals,
                         left=base_vals[0], right=base_vals[-1])
    return baseline


# ==================================================================
# ================= LADDER-QC: METRIKKER ===========================
# ==================================================================

def compute_ladder_qc_metrics(fsa: FsaFile) -> dict[str, float | int]:
    """Beregner QC-metrikker for ladder-fit using actual basepair mapping."""
    ladder_size = np.array(fsa.ladder_steps, dtype=float)
    best_combination = np.array(fsa.best_size_standard, dtype=float)

    ladder_model = getattr(fsa, "ladder_model", None)
    if ladder_model is not None:
        predicted = np.asarray(ladder_model.predict(best_combination.reshape(-1, 1)), dtype=float).reshape(-1)
    else:
        df = getattr(fsa, "sample_data_with_basepairs", None)
        if df is not None and "basepairs" in df.columns and "time" in df.columns:
            lookup = (
                df.loc[:, ["time", "basepairs"]]
                .drop_duplicates(subset=["time"], keep="last")
                .set_index("time")["basepairs"]
                .to_dict()
            )
            predicted = np.array(
                [float(lookup.get(int(idx), np.nan)) for idx in best_combination],
                dtype=float,
            )
        else:
            predicted = np.array([], dtype=float)

    if predicted is None or len(predicted) == 0:
        return {
            "r2": float("nan"),
            "mean_abs_error_bp": float("inf"),
            "max_abs_error_bp": float("inf"),
            "n_ladder_steps": 0,
            "n_size_standard_peaks": 0,
        }

    if np.any(np.isnan(predicted)):
        return {
            "r2": float("nan"),
            "mean_abs_error_bp": float("inf"),
            "max_abs_error_bp": float("inf"),
            "n_ladder_steps": int(ladder_size.size),
            "n_size_standard_peaks": int(best_combination.size),
        }

    r2 = float(r2_score(ladder_size, predicted))
    abs_errors = np.abs(ladder_size - predicted)
    mean_abs_error = float(np.mean(abs_errors)) if abs_errors.size else float("inf")
    max_abs_error = float(np.max(abs_errors)) if abs_errors.size else float("inf")

    return {
        "r2": r2,
        "mean_abs_error_bp": mean_abs_error,
        "max_abs_error_bp": max_abs_error,
        "n_ladder_steps": int(ladder_size.size),
        "n_size_standard_peaks": int(best_combination.size),
    }


# ==================================================================
# ================= MANUAL LADDER ADJUSTMENT =======================
# ==================================================================

def get_ladder_candidates(fsa: FsaFile) -> pd.DataFrame:
    """
    Returns all detected peaks in the size standard channel as a DataFrame.
    Useful for manual selection.
    """
    ss_peaks = getattr(fsa, "size_standard_peaks", None)
    if ss_peaks is None:
        return pd.DataFrame(columns=["index", "time", "intensity"])

    peak_times = np.asarray(ss_peaks, dtype=float)
    trace = np.asarray(fsa.size_standard, dtype=float)
    peak_indices = np.rint(peak_times).astype(int)
    valid_mask = (peak_indices >= 0) & (peak_indices < len(trace))
    peak_times = peak_times[valid_mask]
    peak_indices = peak_indices[valid_mask]

    manual_candidates = set(
        float(v) for v in getattr(fsa, "manual_ladder_candidates", []) or []
    )
    sources = [
        "manual" if any(abs(float(time_value) - manual) <= 1e-6 for manual in manual_candidates) else "auto"
        for time_value in peak_times
    ]

    return pd.DataFrame({
        "index": np.arange(len(peak_times)),
        "time": peak_times,
        "intensity": trace[peak_indices],
        "source": sources,
    })


def apply_manual_ladder_mapping(fsa: FsaFile, adjustment: dict[int, int] | dict) -> FsaFile:
    """
    Applies a manual mapping of ladder steps to candidate peak indices.
    
    mapping: {ladder_step_index: candidate_peak_index}
    """
    payload = _normalize_ladder_adjustment_payload(adjustment)
    if payload is None:
        raise ValueError("No ladder adjustment payload provided.")

    mapping = payload["mapping"]
    mapping_times = payload["mapping_times"]
    manual_candidates = payload["manual_candidates"]
    ladder_steps = _get_expected_ladder_steps(fsa)
    current_ladder_steps = np.asarray(getattr(fsa, "ladder_steps", ladder_steps), dtype=float)
    ss_peaks = fsa.size_standard_peaks

    if ss_peaks is None:
        raise ValueError("No size standard peaks found in FsaFile.")

    if manual_candidates:
        merged = list(np.asarray(ss_peaks, dtype=float))
        for time_value in manual_candidates:
            if not any(abs(float(existing) - float(time_value)) <= 1e-6 for existing in merged):
                merged.append(float(time_value))
        merged.sort()
        fsa.size_standard_peaks = np.asarray(merged, dtype=float)
        ss_peaks = fsa.size_standard_peaks
    fsa.manual_ladder_candidates = [float(v) for v in manual_candidates]
    
    selected_peaks = np.full(len(ladder_steps), np.nan, dtype=float)
    current = getattr(fsa, "best_size_standard", None)
    if current is not None and len(current) == len(current_ladder_steps):
        current = np.asarray(current, dtype=float)
        step_map = _map_step_indices(current_ladder_steps, ladder_steps)
        for current_idx, full_idx in step_map.items():
            selected_peaks[full_idx] = current[current_idx]

    for step_idx, peak_time in mapping_times.items():
        if step_idx < 0 or step_idx >= len(ladder_steps):
            continue
        selected_peaks[step_idx] = float(peak_time)

    for step_idx, peak_idx in mapping.items():
        if step_idx < 0 or step_idx >= len(ladder_steps):
            continue
        if step_idx in mapping_times:
            continue
        if peak_idx < 0 or peak_idx >= len(ss_peaks):
            continue
        selected_peaks[step_idx] = ss_peaks[peak_idx]

    missing = np.isnan(selected_peaks)
    if np.any(missing):
        raise ValueError("Manual ladder mapping is incomplete. Start from an auto-fit or map every missing ladder step.")

    if np.any(np.diff(selected_peaks) <= 0):
        raise ValueError("Selected ladder peaks must be strictly increasing in time.")

    fsa.expected_ladder_steps = ladder_steps.copy()
    fsa.ladder_steps = ladder_steps.copy()
    fsa.best_size_standard = selected_peaks

    # Re-run fitting
    fsa = fit_size_standard_to_ladder(fsa)
    if not getattr(fsa, "fitted_to_model", False):
        raise ValueError("Manual ladder mapping did not produce a valid fit.")

    return fsa


# ==================================================================
# =========== SL AREA METRICS =====================================
# ==================================================================

def compute_sl_area_metrics(
    fsa: FsaFile,
    trace_channel: str,
    targets_bp: list[float],
    window_bp: float = SL_WINDOW_BP,
) -> dict[str, list[float] | float]:
    """Beregner area for SL-fragmenter ved å integrere råtrace i et bp-vindu."""
    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        raise ValueError("sample_data_with_basepairs er tom/None – kan ikke beregne SL-area.")

    if trace_channel not in fsa.fsa:
        raise ValueError(f"Fant ikke kanal {trace_channel} i FSA-filen.")

    trace = np.asarray(fsa.fsa[trace_channel])

    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        raise ValueError("sample_data_with_basepairs mangler 'time' og/eller 'basepairs'.")

    time_arr = raw_df["time"].astype(int).to_numpy()
    bp_arr = raw_df["basepairs"].to_numpy()
    
    # Calculate baseline to avoid integrating background/negative noise
    from core.analysis import estimate_running_baseline
    baseline = estimate_running_baseline(trace, bin_size=5000, quantile=0.01)
    trace_corr = np.maximum(trace - baseline, 0.0)

    results = []
    for target_bp in targets_bp:
        mask = (bp_arr >= (target_bp - window_bp)) & (bp_arr <= (target_bp + window_bp))
        if not np.any(mask):
            area_val = 0.0
        else:
            time_idx = time_arr[mask]
            time_idx = time_idx[(time_idx >= 0) & (time_idx < len(trace))]
            if time_idx.size == 0:
                area_val = 0.0
            else:
                area_val = float(trace_corr[time_idx].sum())

        results.append({"bp": float(target_bp), "area": area_val})

    total_area = float(sum(r["area"] for r in results))
    for r in results:
        if total_area > 0:
            r["percent"] = (r["area"] / total_area) * 100.0
        else:
            r["percent"] = float("nan")

    return {
        "targets_bp": [r["bp"] for r in results],
        "areas": [r["area"] for r in results],
        "percents": [r["percent"] for r in results],
        "total_area": total_area,
    }


# ==================================================================
# =========== SL AUTO PEAK DETECTION ===============================
# ==================================================================

def auto_detect_sl_peaks(
    fsa: FsaFile,
    peak_channels: list[str],
    targets_bp: list[float],
    window_bp: float,
    min_height: float = PEAK_MIN_HEIGHT,
) -> dict[str, pd.DataFrame]:
    """Automatisk peak-detection for SL."""
    peaks_by_channel: dict[str, pd.DataFrame] = {}

    raw_df = getattr(fsa, "sample_data_with_basepairs", None)
    if raw_df is None or raw_df.empty:
        print_warning(
            f"[SL_PEAKS] sample_data_with_basepairs er tom/None for {fsa.file_name}"
        )
        for ch in peak_channels:
            peaks_by_channel[ch] = pd.DataFrame(columns=["basepairs", "peaks", "keep"])
        return peaks_by_channel

    if "time" not in raw_df.columns or "basepairs" not in raw_df.columns:
        print_warning(
            f"[SL_PEAKS] sample_data_with_basepairs mangler 'time'/'basepairs' for {fsa.file_name}"
        )
        for ch in peak_channels:
            peaks_by_channel[ch] = pd.DataFrame(columns=["basepairs", "peaks", "keep"])
        return peaks_by_channel

    time_all = raw_df["time"].astype(int).to_numpy()
    bp_all = raw_df["basepairs"].to_numpy()

    for ch in peak_channels:
        if ch not in fsa.fsa:
            print_warning(f"[SL_PEAKS] Kanal {ch} finnes ikke i {fsa.file_name}")
            peaks_by_channel[ch] = pd.DataFrame(columns=["basepairs", "peaks", "keep"])
            continue

        trace = np.asarray(fsa.fsa[ch])
        bp_list: list[float] = []
        height_list: list[float] = []

        # 1) Hoved-fragmenter
        for target_bp in targets_bp:
            local_window = float(window_bp)
            if abs(target_bp - 600.0) <= 1.0:
                local_window = 40.0

            win_min = float(target_bp) - local_window
            win_max = float(target_bp) + local_window

            win_mask = (bp_all >= win_min) & (bp_all <= win_max)
            if not np.any(win_mask):
                continue

            bp_win = bp_all[win_mask]
            time_win = time_all[win_mask]

            valid_mask = (time_win >= 0) & (time_win < len(trace))
            if not np.any(valid_mask):
                continue

            bp_win = bp_win[valid_mask]
            time_win = time_win[valid_mask]
            y_win = trace[time_win]

            if y_win.size == 0 or not np.any(np.isfinite(y_win)):
                continue

            j = int(np.nanargmax(y_win))
            bp_peak = float(bp_win[j])
            height_peak = float(y_win[j])

            if height_peak >= min_height:
                bp_list.append(bp_peak)
                height_list.append(height_peak)

        # 2) Ekstra skulder-peak rundt ~90 bp
        extra_center = 90.0
        extra_halfwidth = 5.0
        win_min = extra_center - extra_halfwidth
        win_max = extra_center + extra_halfwidth

        win_mask = (bp_all >= win_min) & (bp_all <= win_max)
        if np.any(win_mask):
            bp_win = bp_all[win_mask]
            time_win = time_all[win_mask]
            valid_mask = (time_win >= 0) & (time_win < len(trace))
            if np.any(valid_mask):
                bp_win = bp_win[valid_mask]
                time_win = time_win[valid_mask]
                y_win = trace[time_win]
                if y_win.size > 0 and np.any(np.isfinite(y_win)):
                    j = int(np.nanargmax(y_win))
                    bp_peak = float(bp_win[j])
                    height_peak = float(y_win[j])
                    if height_peak >= min_height:
                        bp_list.append(bp_peak)
                        height_list.append(height_peak)

        # 3) Bygg DataFrame
        if bp_list:
            df = pd.DataFrame({
                "basepairs": bp_list,
                "peaks": height_list,
                "keep": [True] * len(bp_list),
            })
        else:
            df = pd.DataFrame(columns=["basepairs", "peaks", "keep"])

        peaks_by_channel[ch] = df

    return peaks_by_channel
