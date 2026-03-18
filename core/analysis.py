"""
Fraggler Diagnostics — Analysis Functions.

Ladder fitting (LIZ / ROX), SL peak detection, ladder QC metrics,
SL area metrics, local-maxima helpers, and running-baseline estimation.
"""
from __future__ import annotations

from pathlib import Path
import copy

import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from sklearn.metrics import mean_squared_error, r2_score

import json
from fraggler.fraggler import (
    FsaFile,
    find_size_standard_peaks,
    return_maxium_allowed_distance_between_size_standard_peaks,
    generate_combinations,
    calculate_best_combination_of_size_standard_peaks,
    fit_size_standard_to_ladder,
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
LADDER_MAX_ITERATIONS = 40
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
DESCENDING_RECOVERY_MIN_INTENSITY = 350.0
PARTIAL_RESCUE_MISSING_STEP_PENALTY = 0.50
ROX_DYEBLOB_HEIGHT_MULTIPLIER = 2.0
ROX_DYEBLOB_EARLY_INDEX = 2000
ROX_DYEBLOB_TIGHT_GAP = 40
ROX_DYEBLOB_CLUSTER_GAP = 70


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


def _refine_polynomial(fsa: FsaFile, label: str, fsa_path) -> FsaFile | None:
    """Apply polynomial degree-3 fit with iterative outlier removal.

    Used as a fallback when standard fitting fails, and also as a
    post-fit refinement when R² < MIN_R2_QUALITY.
    """
    try:
        from sklearn.linear_model import LinearRegression
        if not hasattr(fsa, "expected_ladder_steps") or getattr(fsa, "expected_ladder_steps", None) is None:
            fsa.expected_ladder_steps = np.asarray(fsa.ladder_steps, dtype=float).copy()
        X_vals = np.array(fsa.best_size_standard, dtype=float)
        Y_vals = np.array(fsa.ladder_steps, dtype=float)

        X_iter = X_vals.copy()
        Y_iter = Y_vals.copy()
        for _it in range(5):
            if len(X_iter) < 8:
                break
            coeffs = np.polyfit(X_iter, Y_iter, 3)
            Y_pred = np.polyval(coeffs, X_iter)
            residuals = np.abs(Y_iter - Y_pred)
            max_idx = np.argmax(residuals)
            if residuals[max_idx] < 2.0:
                break
            X_iter = np.delete(X_iter, max_idx)
            Y_iter = np.delete(Y_iter, max_idx)

        coeffs = np.polyfit(X_iter, Y_iter, 3)
        all_time = np.arange(len(fsa.sample_data))
        bp_all = np.polyval(coeffs, all_time)

        df = (pd.DataFrame({"peaks": fsa.sample_data})
              .reset_index().rename(columns={"index": "time"}))
        df["basepairs"] = np.round(bp_all, 2)
        df = df.loc[df.basepairs >= 0]

        fsa.sample_data_with_basepairs = df
        fsa.fitted_to_model = True
        fsa.best_size_standard = X_iter
        fsa.ladder_steps = Y_iter

        dummy_model = LinearRegression()
        dummy_model.coef_ = np.array([coeffs[-2]])
        dummy_model.intercept_ = coeffs[-1]
        fsa.ladder_model = dummy_model
        strategy = "polynomial_refine_partial" if _missing_expected_ladder_steps(fsa) else "polynomial_refine_full"
        fsa = _set_ladder_fit_metadata(fsa, strategy)

        print_green(f"[{label}] Brukte polynomial refinement (degree-3, outlier removal) for {fsa_path.name}")
        return fsa
    except Exception as e:
        print_warning(f"[{label}] Polynomial refinement feilet for {fsa_path.name}: {e}")
        return None


def _rank_size_standard_combinations(fsa: FsaFile) -> list[np.ndarray]:
    """Return the smoothest ladder candidates, matching the original spline heuristic."""
    combinations = getattr(fsa, "best_size_standard_combinations", None)
    if combinations is None or getattr(combinations, "empty", True):
        return []

    ranked = (
        combinations.assign(
            der=lambda x: [
                UnivariateSpline(fsa.ladder_steps, y, s=0).derivative(n=2)
                for y in x.combinations
            ]
        )
        .assign(max_value=lambda x: [max(abs(y(fsa.ladder_steps))) for y in x.der])
        .sort_values("max_value", ascending=True)
    )
    return [
        np.asarray(row.combinations, dtype=float)
        for row in ranked.head(LADDER_CANDIDATE_COUNT).itertuples()
    ]


def _candidate_fit_score(fsa: FsaFile) -> tuple[float, float, float, float]:
    metrics = compute_ladder_qc_metrics(fsa)
    intensity_penalty = _candidate_intensity_penalty(fsa)
    return (
        float(metrics.get("mean_abs_error_bp", float("inf"))) + float(intensity_penalty),
        float(metrics.get("max_abs_error_bp", float("inf"))),
        -float(metrics.get("r2", float("-inf"))),
        float(intensity_penalty),
    )


def _rescue_fit_score(fsa: FsaFile) -> tuple[float, float, float, float]:
    metrics = compute_ladder_qc_metrics(fsa)
    missing_count = len(_missing_expected_ladder_steps(fsa))
    intensity_penalty = _candidate_intensity_penalty(fsa)
    return (
        float(metrics.get("mean_abs_error_bp", float("inf")))
        + (missing_count * PARTIAL_RESCUE_MISSING_STEP_PENALTY)
        + float(intensity_penalty),
        float(metrics.get("max_abs_error_bp", float("inf"))),
        -float(metrics.get("r2", float("-inf"))),
        float(intensity_penalty),
    )


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

def _select_best_ladder_candidate(fsa: FsaFile, ranked_combinations: list[np.ndarray] | None = None) -> FsaFile | None:
    """Fit the top smooth candidates and keep the best actual ladder fit."""
    if ranked_combinations is None:
        ranked_combinations = _rank_size_standard_combinations(fsa)
    if not ranked_combinations:
        return None

    best_fit = None
    best_score = None

    for combo in ranked_combinations:
        trial = copy.deepcopy(fsa)
        trial.best_size_standard = combo
        try:
            trial = fit_size_standard_to_ladder(trial)
        except Exception:
            continue
        if not getattr(trial, "fitted_to_model", False):
            continue

        score = _candidate_fit_score(trial)
        if best_score is None or score < best_score:
            best_score = score
            best_fit = trial

    return best_fit


def _try_high_end_ladder_rescue(fsa: FsaFile, label: str, fsa_path: Path) -> FsaFile | None:
    full_steps = _get_expected_ladder_steps(fsa)
    if full_steps.size < 12:
        return None

    best_fit = None
    best_score = None
    max_skip = min(6, max(0, int(full_steps.size) - 8))
    if max_skip < 1:
        return None

    for skip_low in range(1, max_skip + 1):
        trial = copy.deepcopy(fsa)
        trial.expected_ladder_steps = full_steps.copy()
        trial.ladder_steps = np.asarray(full_steps[skip_low:], dtype=float)
        trial.n_ladder_peaks = trial.ladder_steps.size
        trial.max_peaks_allow_in_size_standard = trial.n_ladder_peaks + 15

        ss_peaks = getattr(trial, "size_standard_peaks", None)
        if ss_peaks is None or len(ss_peaks) < trial.n_ladder_peaks:
            continue

        try:
            trial = return_maxium_allowed_distance_between_size_standard_peaks(trial, multiplier=2)
            for _ in range(LADDER_MAX_ITERATIONS):
                trial = generate_combinations(trial)
                best = getattr(trial, "best_size_standard_combinations", None)
                if best is not None and best.shape[0] > 0:
                    break
                trial.maxium_allowed_distance_between_size_standard_peaks += 10

            best = getattr(trial, "best_size_standard_combinations", None)
            if best is None or best.shape[0] == 0:
                continue

            selected_fit = _select_best_ladder_candidate(trial)
            if selected_fit is not None:
                trial = selected_fit
            else:
                trial = calculate_best_combination_of_size_standard_peaks(trial)
                trial = fit_size_standard_to_ladder(trial)

            if not getattr(trial, "fitted_to_model", False):
                continue

            score = _rescue_fit_score(trial)
            if best_score is None or score < best_score:
                best_fit = trial
                best_score = score
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

    trial = copy.deepcopy(fsa)
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
    adj_path = fsa.file.with_suffix(".ladder_adj.json")
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
        mirrored = 0
        for sibling in _sibling_fsa_paths(fsa.file):
            sibling_adj = sibling.with_suffix(".ladder_adj.json")
            try:
                with open(sibling_adj, "w") as f:
                    json.dump(payload, f)
                mirrored += 1
            except Exception:
                continue
        print_green(f"Saved ladder adjustment to {adj_path.name}")
        if mirrored:
            print_green(f"Mirrored ladder adjustment to {mirrored} sibling project copie(s)")
    except Exception as e:
        print_warning(f"Could not save ladder adjustment: {e}")


def load_ladder_adjustment(fsa: FsaFile) -> dict | None:
    """Loads a manual mapping payload from a .json file if it exists."""
    adj_path = fsa.file.with_suffix(".ladder_adj.json")
    if adj_path.exists():
        try:
            with open(adj_path, "r") as f:
                payload = json.load(f)
                return _normalize_ladder_adjustment_payload(payload)
        except Exception as e:
            print_warning(f"Could not load ladder adjustment {adj_path.name}: {e}")

    for sibling in _sibling_fsa_paths(fsa.file):
        sibling_adj = sibling.with_suffix(".ladder_adj.json")
        if not sibling_adj.exists():
            continue
        try:
            with open(sibling_adj, "r") as f:
                payload = json.load(f)
            print_green(f"Using sibling ladder adjustment from {sibling_adj}")
            return _normalize_ladder_adjustment_payload(payload)
        except Exception as e:
            print_warning(f"Could not load sibling ladder adjustment {sibling_adj.name}: {e}")
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
                if h > 31000 or p < 1000 or h > 2.0 * median_h:
                    continue
                cleaned.append(p)
            if len(cleaned) >= 3:
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
                import copy
                best_fallback_fsa = copy.deepcopy(fsa) 
            
            try:
                if not getattr(fsa, "fitted_to_model", False):
                    fsa = fit_size_standard_to_ladder(fsa)
                if getattr(fsa, "fitted_to_model", False):
                    qc = compute_ladder_qc_metrics(fsa)
                    if qc["r2"] < HIGH_END_RESCUE_R2:
                        rescued = _try_high_end_ladder_rescue(fsa, "LIZ", fsa_path)
                        if rescued is not None and _rescue_fit_score(rescued) < _rescue_fit_score(fsa):
                            kept = len(getattr(rescued, "ladder_steps", []))
                            total = len(getattr(rescued, "expected_ladder_steps", getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", []))))
                            print_green(
                                f"[LIZ] High-end ladder rescue selected for {fsa_path.name} using the top {kept}/{total} ladder steps."
                            )
                            fsa = rescued
                            qc = compute_ladder_qc_metrics(fsa)
                    if qc["r2"] >= MIN_R2_QUALITY:
                        return _finalize_auto_fit_metadata(fsa)
                    else:
                        refined = _refine_polynomial(fsa, "LIZ", fsa_path)
                        if refined:
                            return refined
                        return _finalize_auto_fit_metadata(fsa)
            except ValueError:
                pass
        except ValueError:
            continue

    if best_fallback_fsa is not None:
        refined = _refine_polynomial(best_fallback_fsa, "LIZ", fsa_path)
        if refined:
            return refined

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
            cleaned = _clean_rox_size_standard_peaks(np.asarray(all_found), rox_data)
            if len(cleaned) >= 3:
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
                import copy
                best_fallback_fsa = copy.deepcopy(fsa)
            
            try:
                if not getattr(fsa, "fitted_to_model", False):
                    fsa = fit_size_standard_to_ladder(fsa)
                if getattr(fsa, "fitted_to_model", False):
                    qc = compute_ladder_qc_metrics(fsa)
                    if qc["r2"] < HIGH_END_RESCUE_R2:
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
                    if qc["r2"] >= MIN_R2_QUALITY:
                        return _finalize_auto_fit_metadata(fsa)
                    else:
                        refined = _refine_polynomial(fsa, "ROX", fsa_path)
                        if refined:
                            return refined
                        return _finalize_auto_fit_metadata(fsa)
            except ValueError:
                pass
        except ValueError:
            continue

    if best_fallback_fsa is not None:
        refined = _refine_polynomial(best_fallback_fsa, "ROX", fsa_path)
        if refined:
            return refined

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
    
    # Compute predicted basepairs from the actual basepair mapping
    df = getattr(fsa, "sample_data_with_basepairs", None)
    if df is not None and "basepairs" in df.columns and "time" in df.columns:
        # Look up predicted bp for each size standard peak index
        predicted = []
        for idx in best_combination:
            row = df.loc[df["time"] == int(idx)]
            if len(row) > 0:
                predicted.append(float(row["basepairs"].iloc[0]))
            else:
                # Fallback: use model prediction
                predicted.append(float(fsa.ladder_model.predict(np.array([[idx]]))[0]))
        predicted = np.array(predicted)
    else:
        predicted = fsa.ladder_model.predict(best_combination.reshape(-1, 1))

    if predicted is None or len(predicted) == 0:
        return {
            "r2": float("nan"),
            "mean_abs_error_bp": float("inf"),
            "max_abs_error_bp": float("inf"),
            "n_ladder_steps": 0,
            "n_size_standard_peaks": 0,
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
