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


def _refine_polynomial(fsa: FsaFile, label: str, fsa_path) -> FsaFile | None:
    """Apply polynomial degree-3 fit with iterative outlier removal.

    Used as a fallback when standard fitting fails, and also as a
    post-fit refinement when R² < MIN_R2_QUALITY.
    """
    try:
        from sklearn.linear_model import LinearRegression
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


def _candidate_fit_score(fsa: FsaFile) -> tuple[float, float, float]:
    metrics = compute_ladder_qc_metrics(fsa)
    return (
        float(metrics.get("mean_abs_error_bp", float("inf"))),
        float(metrics.get("max_abs_error_bp", float("inf"))),
        -float(metrics.get("r2", float("-inf"))),
    )

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


# ==================================================================
# ==================== ANALYSEFUNKSJONER ===========================
# ==================================================================

def save_ladder_adjustment(fsa: FsaFile, mapping: dict[int, int]) -> None:
    """Saves a manual mapping to a .json file alongside the .fsa file."""
    adj_path = fsa.file.with_suffix(".ladder_adj.json")
    try:
        with open(adj_path, "w") as f:
            json.dump(mapping, f)
        print_green(f"Saved ladder adjustment to {adj_path.name}")
    except Exception as e:
        print_warning(f"Could not save ladder adjustment: {e}")


def load_ladder_adjustment(fsa: FsaFile) -> dict[int, int] | None:
    """Loads a manual mapping from a .json file if it exists."""
    adj_path = fsa.file.with_suffix(".ladder_adj.json")
    if adj_path.exists():
        try:
            with open(adj_path, "r") as f:
                mapping = json.load(f)
                return {int(k): int(v) for k, v in mapping.items()}
        except Exception as e:
            print_warning(f"Could not load ladder adjustment {adj_path.name}: {e}")
    return None


def _try_apply_saved_ladder_adjustment(fsa: FsaFile, mapping: dict[int, int] | None, label: str) -> FsaFile | None:
    """Applies a saved ladder adjustment if valid, otherwise warns and falls back to auto-fit."""
    if not mapping:
        return None
    try:
        print_green(f"[{label}] Applying manual ladder adjustment for {fsa.file_name}")
        return apply_manual_ladder_mapping(fsa, mapping)
    except Exception as exc:
        print_warning(
            f"[{label}] Ignoring invalid saved ladder adjustment for {fsa.file_name}: {exc}. Falling back to auto-fit."
        )
        return None


def analyse_fsa_liz(fsa_path: Path, sample_channel: str) -> FsaFile | None:
    """Ladder-fit for LIZ (TCRg/IGK/KDE).
    
    Uses multi-config search, dye-blob filtering, polynomial sizing,
    and iterative outlier removal for robust ladder fitting.
    """
    print_green(
        f"=== Analysing {fsa_path} (LIZ, sample {sample_channel}, Python API) ==="
    )

    configs = [
        {"min_h": MIN_SIZE_STANDARD_HEIGHT_LIZ, "min_d": MIN_DISTANCE_BETWEEN_PEAKS_LIZ},
        {"min_h": 200, "min_d": 20},
        {"min_h": 100, "min_d": 15},
        {"min_h": 50, "min_d": 10},
    ]

    base_fsa = FsaFile(
        file=str(fsa_path),
        ladder=LIZ_LADDER,
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
            ladder=LIZ_LADDER,
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
                    if qc["r2"] >= MIN_R2_QUALITY:
                        return fsa
                    else:
                        refined = _refine_polynomial(fsa, "LIZ", fsa_path)
                        if refined:
                            return refined
                        return fsa
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

def analyse_fsa_rox(fsa_path: Path, sample_channel: str) -> FsaFile | None:
    """Ladder-fit for ROX (FR1–3, TCRbA/B/C, SL, DHJH_D/E).
    
    Uses multi-config search, dye-blob filtering, polynomial sizing,
    and iterative outlier removal for robust ladder fitting.
    """
    print_green(
        f"=== Analysing {fsa_path} (ROX, sample {sample_channel}, Python API) ==="
    )

    configs = [
        {"min_h": MIN_SIZE_STANDARD_HEIGHT_ROX, "min_d": MIN_DISTANCE_BETWEEN_PEAKS_ROX},
        {"min_h": 100, "min_d": 15},
        {"min_h": 50, "min_d": 10},
        {"min_h": 20, "min_d": 8},
    ]

    base_fsa = FsaFile(
        file=str(fsa_path),
        ladder=ROX_LADDER,
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
            ladder=ROX_LADDER,
            sample_channel=sample_channel,
            min_distance_between_peaks=cfg["min_d"],
            min_size_standard_height=cfg["min_h"],
            size_standard_channel="DATA4",
        )
        rox_data = np.asarray(fsa.fsa["DATA4"]).astype(float)
        fsa = find_size_standard_peaks(fsa)

        all_found = getattr(fsa, "size_standard_peaks", None)
        if all_found is not None:
            # Calculate median height for dye-blob detection
            heights = np.array([rox_data[p] for p in all_found])
            median_h = np.median(heights)
            cleaned = []
            for p in all_found:
                h = rox_data[p]
                # Filter: saturated (>31000), early noise (<1000 idx),
                # or dye blobs (>3x median height)
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
                    if qc["r2"] >= MIN_R2_QUALITY:
                        return fsa
                    else:
                        refined = _refine_polynomial(fsa, "ROX", fsa_path)
                        if refined:
                            return refined
                        return fsa
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
    
    trace = fsa.size_standard
    return pd.DataFrame({
        "index": np.arange(len(ss_peaks)),
        "time": ss_peaks,
        "intensity": trace[ss_peaks]
    })


def apply_manual_ladder_mapping(fsa: FsaFile, mapping: dict[int, int]) -> FsaFile:
    """
    Applies a manual mapping of ladder steps to candidate peak indices.
    
    mapping: {ladder_step_index: candidate_peak_index}
    """
    ladder_steps = fsa.ladder_steps
    ss_peaks = fsa.size_standard_peaks
    
    if ss_peaks is None:
        raise ValueError("No size standard peaks found in FsaFile.")
    
    current = getattr(fsa, "best_size_standard", None)
    if current is not None and len(current) == len(ladder_steps):
        selected_peaks = np.asarray(current, dtype=float).copy()
    else:
        selected_peaks = np.full(len(ladder_steps), np.nan, dtype=float)
    
    for step_idx, peak_idx in mapping.items():
        if step_idx < 0 or step_idx >= len(ladder_steps):
            continue
        if peak_idx < 0 or peak_idx >= len(ss_peaks):
            continue
        selected_peaks[step_idx] = ss_peaks[peak_idx]

    missing = np.isnan(selected_peaks)
    if np.any(missing):
        raise ValueError("Manual ladder mapping is incomplete. Start from an auto-fit or map every missing ladder step.")

    if np.any(np.diff(selected_peaks) <= 0):
        raise ValueError("Selected ladder peaks must be strictly increasing in time.")

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
