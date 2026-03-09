"""
Fraggler Diagnostics — Analysis Functions.

Ladder fitting (LIZ / ROX), SL peak detection, ladder QC metrics,
SL area metrics, local-maxima helpers, and running-baseline estimation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score

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


# ==================================================================
# ==================== ANALYSEFUNKSJONER ===========================
# ==================================================================

def analyse_fsa_liz(fsa_path: Path, sample_channel: str):
    """Ladder-fit for LIZ (TCRg/IGK/KDE)."""
    print_green(
        f"=== Analysing {fsa_path} (LIZ, sample {sample_channel}, Python API) ==="
    )

    fsa = FsaFile(
        file=str(fsa_path),
        ladder=LIZ_LADDER,
        sample_channel=sample_channel,
        min_distance_between_peaks=MIN_DISTANCE_BETWEEN_PEAKS_LIZ,
        min_size_standard_height=MIN_SIZE_STANDARD_HEIGHT_LIZ,
        size_standard_channel="DATA105",
    )

    fsa = find_size_standard_peaks(fsa)

    ss_peaks = getattr(fsa, "size_standard_peaks", None)
    n_ss = 0
    if ss_peaks is not None and hasattr(ss_peaks, "shape"):
        n_ss = ss_peaks.shape[0]

    if n_ss < 2:
        print_warning(
            f"[LIZ] Fant {n_ss} ladder-peaks for {fsa_path.name} – for lite til å "
            f"beregne avstander. Skipper fila."
        )
        return None

    try:
        fsa = return_maxium_allowed_distance_between_size_standard_peaks(
            fsa, multiplier=2
        )
    except ValueError as e:
        print_warning(
            f"[LIZ] return_maxium_allowed_distance_between_size_standard_peaks feilet "
            f"for {fsa_path.name}: {e} – skipper fila."
        )
        return None

    for _ in range(20):
        fsa = generate_combinations(fsa)
        best = getattr(fsa, "best_size_standard_combinations", None)
        if best is not None and best.shape[0] > 0:
            break
        fsa.maxium_allowed_distance_between_size_standard_peaks += 10

    best = getattr(fsa, "best_size_standard_combinations", None)
    if best is None or best.shape[0] == 0:
        print_warning(
            f"[LIZ] Ingen gyldige size-standard kombinasjoner for {fsa_path.name} – skipper."
        )
        return None

    try:
        fsa = calculate_best_combination_of_size_standard_peaks(fsa)
        fsa = fit_size_standard_to_ladder(fsa)
    except ValueError as e:
        print_warning(
            f"[LIZ] Ladder-fit feilet for {fsa_path.name}: {e} – skipper fila."
        )
        return None

    if not getattr(fsa, "fitted_to_model", False):
        print_warning(
            f"[LIZ] Modell kunne ikke fit'es for {fsa_path.name} – hopper over denne."
        )
        return None

    return fsa


def analyse_fsa_rox(fsa_path: Path, sample_channel: str):
    """Ladder-fit for ROX (FR1–3, TCRbA/B/C, SL, DHJH_D/E)."""
    print_green(
        f"=== Analysing {fsa_path} (ROX, sample {sample_channel}, Python API) ==="
    )

    fsa = FsaFile(
        file=str(fsa_path),
        ladder=ROX_LADDER,
        sample_channel=sample_channel,
        min_distance_between_peaks=MIN_DISTANCE_BETWEEN_PEAKS_ROX,
        min_size_standard_height=MIN_SIZE_STANDARD_HEIGHT_ROX,
        size_standard_channel="DATA4",
    )

    fsa = find_size_standard_peaks(fsa)

    ss_peaks = getattr(fsa, "size_standard_peaks", None)
    n_ss = 0
    if ss_peaks is not None and hasattr(ss_peaks, "shape"):
        n_ss = ss_peaks.shape[0]

    if n_ss < 2:
        print_warning(
            f"[ROX] Fant {n_ss} ladder-peaks for {fsa_path.name} – for lite til å "
            f"beregne avstander. Skipper fila."
        )
        return None

    try:
        fsa = return_maxium_allowed_distance_between_size_standard_peaks(
            fsa, multiplier=2
        )
    except ValueError as e:
        print_warning(
            f"[ROX] return_maxium_allowed_distance_between_size_standard_peaks feilet "
            f"for {fsa_path.name}: {e} – skipper fila."
        )
        return None

    for _ in range(20):
        fsa = generate_combinations(fsa)
        best = getattr(fsa, "best_size_standard_combinations", None)
        if best is not None and best.shape[0] > 0:
            break
        fsa.maxium_allowed_distance_between_size_standard_peaks += 10

    best = getattr(fsa, "best_size_standard_combinations", None)
    if best is None or best.shape[0] == 0:
        print_warning(
            f"[ROX] Fant ingen gyldige size-standard kombinasjoner for {fsa_path.name} "
            f"– mulig dårlig eller feil ladder. Skipper fila."
        )
        return None

    try:
        fsa = calculate_best_combination_of_size_standard_peaks(fsa)
        fsa = fit_size_standard_to_ladder(fsa)
    except ValueError as e:
        print_warning(
            f"[ROX] Ladder-fit feilet for {fsa_path.name}: {e} – skipper fila."
        )
        return None

    if not getattr(fsa, "fitted_to_model", False):
        print_warning(
            f"[ROX] Modell ble ikke fit'et for {fsa_path.name} – skipper fila."
        )
        return None

    return fsa


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
    bin_size: int = 200,
    quantile: float = 0.10,
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

def compute_ladder_qc_metrics(fsa) -> dict:
    """Beregner MSE og R² for ladder-fit."""
    ladder_size = fsa.ladder_steps
    best_combination = fsa.best_size_standard

    predicted = fsa.ladder_model.predict(best_combination.reshape(-1, 1))

    if predicted is None:
        return {"r2": float("nan"), "n_ladder_steps": 0, "n_size_standard_peaks": 0}

    r2 = float(r2_score(ladder_size, predicted))

    return {
        "r2": r2,
        "n_ladder_steps": int(ladder_size.size),
        "n_size_standard_peaks": int(best_combination.size),
    }


# ==================================================================
# =========== SL AREA METRICS =====================================
# ==================================================================

def compute_sl_area_metrics(
    fsa,
    trace_channel: str,
    targets_bp: list[float],
    window_bp: float = SL_WINDOW_BP,
) -> dict:
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
                area_val = float(trace[time_idx].sum())

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
    fsa,
    peak_channels: list[str],
    targets_bp: list[float],
    window_bp: float,
    min_height: float = 200.0,
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
