import numpy as np
import pandas as pd
from lmfit.models import GaussianModel
from fraggler.fraggler import baseline_arPLS

def compute_peak_area_gaussian(
    trace: np.ndarray,
    time_all: np.ndarray,
    bp_all: np.ndarray,
    center_bp: float,
    half_width_bp: float,
) -> float:
    """
    Computes the area under a peak using robust baseline correction (arPLS)
    followed by a Gaussian fit using lmfit. If the fit fails or R-squared is poor,
    falls back to a simple sum of the baseline-corrected trace.
    
    Args:
        trace: Full trace array.
        time_all: Corresponding time indices.
        bp_all: Corresponding basepairs.
        center_bp: Expected center of the peak in bp.
        half_width_bp: Window half-width in bp around center_bp.
        
    Returns:
        float: Calculated area of the peak.
    """
    if trace.size == 0 or time_all.size == 0 or bp_all.size == 0:
        return 0.0

    # 1) Isolate the region of interest based on basepairs
    mask = (bp_all >= center_bp - half_width_bp) & (bp_all <= center_bp + half_width_bp)
    
    if not np.any(mask):
        return 0.0
        
    roi_time = time_all[mask].astype(int)
    
    # Ensure valid indices
    valid_mask = (roi_time >= 0) & (roi_time < trace.size)
    roi_time = roi_time[valid_mask]
    
    if roi_time.size < 3: # Need at least a few points for any fit/sum
        return 0.0
        
    roi_trace = trace[roi_time]
    
    # 2) Extract a slightly larger window for robust baseline estimation
    # We expand by a factor (e.g., 3x) to give arPLS enough background to work with
    extended_half_width = max(half_width_bp * 3.0, 10.0) # at least 10 bp window
    ext_mask = (bp_all >= center_bp - extended_half_width) & (bp_all <= center_bp + extended_half_width)
    ext_time = time_all[ext_mask].astype(int)
    ext_valid = (ext_time >= 0) & (ext_time < trace.size)
    ext_time = ext_time[ext_valid]
    
    if len(ext_time) < 10:
        # If the extended window is too small, fallback to roi_time
        ext_time = roi_time
        ext_trace = roi_trace
    else:
        ext_trace = trace[ext_time]
        
    # Apply arPLS for baseline calculation
    # We pass it the extended trace subset rather than the entire 10,000+ point trace
    # to avoid extreme performance hits and handle local variations effectively.
    # Note: baseline_arPLS defaults are ratio=0.99, lam=100 in fraggler, we use lam=1e4 for smoother baseline
    try:
        baseline = baseline_arPLS(ext_trace, ratio=0.01, lam=1e4)
    except Exception:
        # Fallback to simple min if arPLS fails
        baseline = np.min(ext_trace)

    ext_corr = np.maximum(ext_trace - baseline, 0.0)
    
    # Map back the corrected trace to our exact ROI
    # ext_time is sorted, roi_time is a subset
    roi_indices_in_ext = np.searchsorted(ext_time, roi_time)
    
    # Safety check
    roi_indices_in_ext = roi_indices_in_ext[roi_indices_in_ext < ext_corr.size]
    if roi_indices_in_ext.size == 0:
        return 0.0
        
    y = ext_corr[roi_indices_in_ext]
    
    if len(y) < 3:
        return float(np.sum(y))
        
    raw_sum_area = float(np.sum(y))

    # 3) Gaussian Fit using lmfit
    model = GaussianModel()
    try:
        # We fit Gaussian over the array index.
        # This keeps the integral (amplitude) in the same scale as np.sum() over time_points
        x_idx = np.arange(len(y))
        params_idx = model.guess(y, x=x_idx)
        out_idx = model.fit(y, params_idx, x=x_idx)
        
        # Check fit quality
        r_sq_idx = out_idx.rsquared if hasattr(out_idx, 'rsquared') else 1 - out_idx.residual.var() / np.var(y)
        
        # If fit is extremely poor, or the amplitude is negative, fallback to simple sum
        if r_sq_idx < 0.5 or out_idx.values["amplitude"] <= 0:
            return raw_sum_area
            
        return float(out_idx.values["amplitude"])
        
    except Exception:
        # Fallback to simple sum
        return raw_sum_area
