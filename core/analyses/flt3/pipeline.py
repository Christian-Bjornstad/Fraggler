from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

from fraggler.fraggler import (
    print_green, 
    print_warning, 
    FsaFile,
    find_size_standard_peaks,
    return_maxium_allowed_distance_between_size_standard_peaks,
    generate_combinations,
    calculate_best_combination_of_size_standard_peaks,
    fit_size_standard_to_ladder
)

from core.analyses.flt3.config import (
    ASSAY_CONFIG,
    PREFERRED_INJECTION_TIME,
)
from core.analyses.flt3.classification import classify_fsa
from core.analysis import (
    compute_ladder_qc_metrics,
    estimate_running_baseline,
)
from core.plotting_mpl import compute_zoom_ymax
from core.html_reports import (
    extract_dit_from_name,
    build_dit_html_reports,
)
from core.utils import is_water_file

def _scan_files(fsa_dir: Path) -> list[Path]:
    """Scans for .fsa files, filtering out water files."""
    if not fsa_dir.exists():
        print_warning(f"FSA-katalog finnes ikke: {fsa_dir}")
        return []

    fsa_files = [
        p for p in sorted(fsa_dir.glob("*.fsa"))
        if not is_water_file(p.name)
    ]
    return fsa_files

def _filter_best_injections(classified_files: list[tuple[Path, dict]]) -> list[tuple[Path, dict]]:
    """
    Keep the best injection for each (DIT, assay, parallel) group based on PREFERRED_INJECTION_TIME.
    Falls back to any available injection if preferred is missing or fails.
    """
    groups = defaultdict(list)
    for p, c in classified_files:
        dit = extract_dit_from_name(p.name) or "unknown"
        par = c.get("parallel") or "none"
        atype = c.get("analysis_type") or "standard"
        groups[(dit, c["assay"], par, atype)].append((p, c))
    
    final_files = []
    for (dit, assay, par, atype), files in groups.items():
        # 1) Try to find preferred injection time
        # Check both assay and analysis_type for preference
        pref_time = PREFERRED_INJECTION_TIME.get(assay)
        if pref_time is None:
            # Check secondary classification if assay is generic (e.g. D8365_kutting)
            for f_p, f_c in files:
                if f_c["analysis_type"] in PREFERRED_INJECTION_TIME:
                    pref_time = PREFERRED_INJECTION_TIME[f_c["analysis_type"]]
                    break
        
        if pref_time is not None:
            best_match = [f for f in files if f[1]["injection_time"] == pref_time]
            if best_match:
                # If multiple with same time, just pick first for now
                final_files.append(best_match[0])
                continue
        
        # 2) Fallback: just pick the first one available
        final_files.append(files[0])
        
    return final_files

def _calculate_auc(trace: np.ndarray, time_idx: np.ndarray) -> float:
    """Calculates the area under curve (AUC) for given time indices."""
    if time_idx.size == 0:
        return 0.0
    return float(trace[time_idx].sum())

def _detect_peaks(fsa: FsaFile, assay: str, wt_bp: float, trace: np.ndarray, mut_bp: float | None = None, analysis_type: str | None = None) -> pd.DataFrame:
    """Detects WT and Mutant peaks and calculates their AUC."""
    sample_data = getattr(fsa, "sample_data_with_basepairs", None)
    if sample_data is None or sample_data.empty:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])
    
    time_all = sample_data["time"].astype(int).to_numpy()
    bp_all = sample_data["basepairs"].to_numpy()
    
    # Broad search window
    bp_min, bp_max = 50.0, 1000.0
    auc_window = 5.0 
    peaks = []
    
    from scipy.signal import find_peaks
    
    # 1) Detect ALL peaks in range 50-1000
    mask = (bp_all >= bp_min) & (bp_all <= bp_max)
    if not mask.any():
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])
        
    y_win = trace[time_all[mask]]
    bp_win = bp_all[mask]
    
    if y_win.size < 3:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])
        
    p_idx, _ = find_peaks(y_win, height=200, distance=20)
    
    for idx in p_idx:
        p_bp = bp_win[idx]
        p_h = y_win[idx]
        
        # AUC
        a_mask = (bp_all >= p_bp - auc_window) & (bp_all <= p_bp + auc_window)
        p_a = _calculate_auc(trace, time_all[a_mask])
        
        # Determine label based on assay and bp
        label = "unspecific"
        # Simple labeling logic
        # For ITD: be stricter with WT window (usually +-4bp)
        wt_tol = 4.0 if assay == "FLT3-ITD" else 5.0
        
        if wt_bp and abs(p_bp - wt_bp) < wt_tol:
            label = "WT"
        elif assay == "FLT3-ITD" and p_bp >= wt_bp + 4.9: # Stricter ITD start
            label = "ITD"
        elif mut_bp and abs(p_bp - mut_bp) < (wt_tol + 2):
            label = "MUT"
            
        peak_info = {
            "basepairs": p_bp,
            "peaks": p_h,
            "area": p_a,
            "label": label,
            "keep": True,
        }
        
        # Calculate AUC for individual channels if they exist
        for ch in ["DATA1", "DATA2", "DATA3"]:
            if ch in fsa.fsa:
                ch_y = np.asarray(fsa.fsa[ch]).astype(float)
                # Apply a generic baseline for the individual channel AUC to be somewhat clean
                base = estimate_running_baseline(ch_y, bin_size=5000, quantile=0.01)
                ch_y_corr = np.maximum(ch_y - base, 0.0)
                peak_info[f"area_{ch}"] = _calculate_auc(ch_y_corr, time_all[a_mask])
        
        peaks.append(peak_info)

    if not peaks:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])

    return pd.DataFrame(peaks)

def _robust_analyse_fsa_rox(fsa_path: Path, sample_channel: str, assay: str) -> FsaFile | None:
    # Try GS500ROX only for now
    ladders = ["GS500ROX"]
    
    # Try multiple parameter sets
    configs = [
        {"min_h": 200, "min_d": 20},
        {"min_h": 100, "min_d": 15},
        {"min_h": 50, "min_d": 10},
        {"min_h": 20, "min_d": 8},
    ]
    
    for ladder_name in ladders:
        for cfg in configs:
            try:
                fsa = FsaFile(
                    file=str(fsa_path),
                    ladder=ladder_name,
                    sample_channel=sample_channel,
                    min_distance_between_peaks=cfg["min_d"],
                    min_size_standard_height=cfg["min_h"],
                    size_standard_channel="DATA4",
                )
                # --- Improved ROX Signal Cleaning ---
                rox_data = np.asarray(fsa.fsa["DATA4"]).astype(float)
                fsa = find_size_standard_peaks(fsa)
                all_found_peaks = fsa.size_standard_peaks
                
                # Filter: Exclude clearly saturated artifacts or very early primer peaks
                cleaned_peaks = []
                for p in all_found_peaks:
                    h = rox_data[p]
                    if h < 31000 and p > 1480: # Filter artifacts before 35bp fragment
                        cleaned_peaks.append(p)
                
                if len(cleaned_peaks) >= 3:
                    fsa.size_standard_peaks = np.array(cleaned_peaks)
                # ------------------------------------
                
                ss_peaks = getattr(fsa, "size_standard_peaks", None)
                if ss_peaks is None or ss_peaks.shape[0] < 3:
                    continue
                    
                fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=1.5)
                for _ in range(30):
                    fsa = generate_combinations(fsa)
                    if getattr(fsa, "best_size_standard_combinations", None) is not None:
                        if fsa.best_size_standard_combinations.shape[0] > 0:
                            break
                    fsa.maxium_allowed_distance_between_size_standard_peaks += 10
                
                fsa = calculate_best_combination_of_size_standard_peaks(fsa)
                
                # Manual linear sizing for absolute robustness
                X_vals = fsa.best_size_standard
                y_vals = fsa.ladder_steps
                
                x1, x2 = float(X_vals[0]), float(X_vals[-1])
                y1, y2 = float(y_vals[0]), float(y_vals[-1])
                
                m = (y2 - y1) / (x2 - x1)
                c = y1 - m * x1
                
                df = (
                    pd.DataFrame({"peaks": fsa.sample_data})
                    .reset_index()
                    .rename(columns={"index": "time"})
                    .assign(basepairs=lambda x: (m * x.time + c).round(2))
                    .loc[lambda x: x.basepairs >= 0]
                )
                fsa.sample_data_with_basepairs = df
                fsa.fitted_to_model = True
                
                # Create a dummy model for compatibility if needed
                from sklearn.linear_model import LinearRegression
                dummy_model = LinearRegression()
                dummy_model.coef_ = np.array([m])
                dummy_model.intercept_ = c
                fsa.ladder_model = dummy_model
                
                return fsa
            except Exception:
                continue
            
    return None

def _calculate_ratios(entries: list[dict]):
    """Calculates all mutant ratios (ITD, D835, NPM1)."""
    for e in entries:
        peaks = e["peaks_by_channel"][e["primary_peak_channel"]]
        if peaks.empty:
            e["ratio"] = 0.0
            continue
        
        wt_peak = peaks[peaks.label == "WT"]
        mut_peaks = peaks[peaks.label.isin(["MUT", "ITD"])]
        
        if not wt_peak.empty:
            # Pick highest WT
            max_wt = wt_peak.sort_values("peaks", ascending=False).iloc[0]
            
            if not max_wt.get("keep", True):
                e["ratio"] = 0.0
                continue
                
            wt_area = max_wt.area
            # Sum all ITD/MUT areas
            mut_area = mut_peaks.area.sum()
            
            if wt_area > 0:
                e["ratio"] = mut_area / wt_area
            else:
                e["ratio"] = 0.0
        else:
            e["ratio"] = 0.0

def generate_flt3_qc_report(entries: list[dict], outdir: Path):
    """Generates a CSV report summarizing FLT3 QC controls."""
    qc_data = []
    
    for e in entries:
        group = e.get("group")
        if group not in ["negative_control", "positive_control", "reactive_control"]:
            continue
            
        peaks = e["peaks_by_channel"][e["primary_peak_channel"]]
        wt_peaks = peaks[peaks.label == "WT"]
        mut_peaks = peaks[peaks.label.isin(["MUT", "ITD"])]
        
        status = "FAIL"
        details = ""
        
        if group == "negative_control":
            # Expect nominal signal (WT might be present but very low, or no mutants)
            if mut_peaks.empty:
                status = "PASS"
            else:
                status = "FAIL"
                details = f"Unexpected mutant peaks found: {mut_peaks.basepairs.tolist()}"
        elif group == "reactive_control":
            # ivs-0000: Expect WT peaks
            if not wt_peaks.empty:
                status = "PASS"
            else:
                status = "FAIL"
                details = "No WT peak detected"
        elif group == "positive_control":
            # ivs-p001: Expect mutant peaks
            if not mut_peaks.empty:
                status = "PASS"
            else:
                status = "FAIL"
                details = "No mutant/ITD peak detected"
                
        qc_data.append({
            "File": e["fsa"].file_name,
            "ControlGroup": group,
            "Assay": e["assay"],
            "InjectionTime": e["injection_time"],
            "Status": status,
            "Details": details
        })
        
    if qc_data:
        df = pd.DataFrame(qc_data)
        outdir.mkdir(parents=True, exist_ok=True)
        csv_path = outdir / "QC_FLT3_Injections.csv"
        df.to_csv(csv_path, index=False)
        print_green(f"FLT3 QC report saved to {csv_path}")

def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
) -> list[dict] | None:
    """
    Kjør FLT3-pipeline på alle .fsa-filer i fsa_dir.
    """
    fsa_dir = Path(fsa_dir).expanduser()
    base_outdir = Path(base_outdir or fsa_dir).expanduser()
    assay_dir = base_outdir / (assay_folder_name or "REPORTS")

    raw_files = _scan_files(fsa_dir)
    classified = []
    for p in raw_files:
        c = classify_fsa(p)
        if c: classified.append((p, c))
            
    if not classified:
        return [] if return_entries else None

    filtered = _filter_best_injections(classified)
    
    entries = []
    for fsa_path, c in filtered:
        fsa = _robust_analyse_fsa_rox(fsa_path, c["primary_peak_channel"], c["assay"])
        if fsa is None:
            print_warning(f"Ladder-fit failed for {fsa_path.name}")
            continue
            
        # 1) Combine traces if multiple peak_channels (e.g. DATA1 + DATA2)
        peak_channels = c.get("peak_channels", [c["primary_peak_channel"]])
        combined_trace = None
        for ch in peak_channels:
            if ch in fsa.fsa:
                raw_y = np.asarray(fsa.fsa[ch]).astype(float)
                # Use less aggressive baseline for FLT3 to avoid masking WT in high-mutant samples.
                baseline = estimate_running_baseline(raw_y, bin_size=5000, quantile=0.01)
                corr_y = np.maximum(raw_y - baseline, 0.0)
                if combined_trace is None:
                    combined_trace = corr_y
                else:
                    combined_trace += corr_y
        
        if combined_trace is None:
            combined_trace = np.zeros(len(next(iter(fsa.fsa.values()))))

        # 2) Detect peaks on combined/corrected trace
        peaks = _detect_peaks(fsa, c["assay"], c["wt_bp"], combined_trace, c.get("mut_bp"), c.get("analysis_type"))
        
        # 3) Use combined trace for zoom max
        ymax = float(np.max(combined_trace)) * 1.1 if combined_trace.any() else 1000.0
        metrics = compute_ladder_qc_metrics(fsa)
        
        entries.append({
            "fsa": fsa,
            "peaks_by_channel": {c["primary_peak_channel"]: peaks},
            "trace_channels": c["trace_channels"],
            "primary_peak_channel": c["primary_peak_channel"],
            "ymax": ymax,
            "assay": c["assay"],
            "analysis_type": c["analysis_type"],
            "parallel": c.get("parallel"),
            "group": c["group"],
            "ladder": fsa.ladder,
            "bp_min": c["bp_min"],
            "bp_max": c["bp_max"],
            "dit": extract_dit_from_name(fsa.file_name),
            "ladder_qc_status": "ok" if metrics["r2"] > 0.99 else "ladder_qc_failed",
            "ladder_r2": metrics["r2"],
            "n_ladder_steps": metrics["n_ladder_steps"],
            "n_size_standard_peaks": metrics["n_size_standard_peaks"],
            "injection_time": c["injection_time"],
        })

    if not entries:
        return [] if return_entries else None

    # Calculate ratios for all entries that have mutant/wt peaks
    _calculate_ratios(entries)

    # Generate QC Report
    generate_flt3_qc_report(entries, assay_dir)

    if make_dit_reports and mode != "controls":
        build_dit_html_reports(entries, assay_dir)

    return entries if return_entries else None
