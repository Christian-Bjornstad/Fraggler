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
)
from core.analyses.shared_pipeline import (
    finalize_pipeline_run,
    normalize_pipeline_paths,
    scan_fsa_files,
)

def _scan_files(fsa_dir: Path, mode: str = "all") -> list[Path]:
    """Scans for .fsa files, filtering out water files."""
    return scan_fsa_files(fsa_dir, mode=mode)

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

def _assay_positive_ratio(assay: str) -> float:
    return float(ASSAY_CONFIG.get(assay, {}).get("positive_ratio", 0.01))

def _bp_in_ranges(bp: float, ranges: list[tuple[float, float]] | None) -> bool:
    if not ranges:
        return False
    return any(start <= bp <= end for start, end in ranges)

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
    assay_cfg = ASSAY_CONFIG.get(assay, {})
    
    from scipy.signal import find_peaks
    
    # 1) Detect ALL peaks in range 50-1000
    mask = (bp_all >= bp_min) & (bp_all <= bp_max)
    if not mask.any():
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])
        
    y_win = trace[time_all[mask]]
    bp_win = bp_all[mask]
    
    if y_win.size < 3:
        return pd.DataFrame(columns=["basepairs", "peaks", "area", "keep", "label"])
        
    peak_height_min = float(assay_cfg.get("peak_height_min", 200))
    peak_distance = int(assay_cfg.get("peak_distance", 20))
    p_idx, _ = find_peaks(y_win, height=peak_height_min, distance=peak_distance)
    
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
        wt_tol = 4.0 if assay == "FLT3-ITD" else 2.0 if assay in {"FLT3-D835", "NPM1"} else 5.0
        itd_min_bp = float(ASSAY_CONFIG.get("FLT3-ITD", {}).get("itd_min_bp", wt_bp + 4.9)) if wt_bp else None
        wt_range = assay_cfg.get("wt_range")
        mut_ranges = assay_cfg.get("mut_ranges")
        
        if _bp_in_ranges(p_bp, [wt_range] if wt_range else None) or (wt_bp and abs(p_bp - wt_bp) < wt_tol):
            label = "WT"
        elif assay == "FLT3-ITD" and itd_min_bp is not None and p_bp >= itd_min_bp:
            label = "ITD"
        elif _bp_in_ranges(p_bp, mut_ranges) or (mut_bp and abs(p_bp - mut_bp) < (wt_tol + 2)):
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
            # D835/NPM1 behave better with the dominant mutant peak, while ITD often has several true mutant peaks.
            if e["assay"] == "FLT3-ITD":
                mut_area = mut_peaks.area.sum()
            elif not mut_peaks.empty:
                mut_area = float(mut_peaks.sort_values("area", ascending=False).iloc[0].area)
            else:
                mut_area = 0.0
            
            if wt_area > 0:
                e["ratio"] = mut_area / wt_area
            else:
                e["ratio"] = 0.0
        else:
            e["ratio"] = 0.0

def _summarize_peak_areas(entry: dict) -> tuple[float, float]:
    peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
    if peaks.empty:
        return 0.0, 0.0

    wt_peak = peaks[peaks.label == "WT"]
    mut_peaks = peaks[peaks.label.isin(["MUT", "ITD"])]

    wt_area = 0.0
    if not wt_peak.empty:
        wt_area = float(wt_peak.sort_values("peaks", ascending=False).iloc[0].area)
    if entry["assay"] == "FLT3-ITD":
        mut_area = float(mut_peaks.area.sum()) if not mut_peaks.empty else 0.0
    elif not mut_peaks.empty:
        mut_area = float(mut_peaks.sort_values("area", ascending=False).iloc[0].area)
    else:
        mut_area = 0.0
    return wt_area, mut_area

def _summarize_detected_peaks(entry: dict) -> dict:
    peaks = entry["peaks_by_channel"][entry["primary_peak_channel"]]
    wt_rows = peaks[peaks.label == "WT"].sort_values("peaks", ascending=False) if not peaks.empty else pd.DataFrame()
    mut_rows = peaks[peaks.label.isin(["MUT", "ITD"])].sort_values("basepairs") if not peaks.empty else pd.DataFrame()

    wt_bp = float(wt_rows.iloc[0].basepairs) if not wt_rows.empty else np.nan
    wt_area = float(wt_rows.iloc[0].area) if not wt_rows.empty else 0.0
    mut_bps = [round(float(v), 2) for v in mut_rows.basepairs.tolist()] if not mut_rows.empty else []
    mut_areas = [round(float(v), 2) for v in mut_rows.area.tolist()] if not mut_rows.empty else []

    return {
        "wt_bp": wt_bp,
        "wt_area": wt_area,
        "mut_bps": mut_bps,
        "mut_areas": mut_areas,
        "mut_area_total": float(sum(mut_areas)),
    }

def _interpret_entry(entry: dict) -> str:
    assay = entry["assay"]
    ratio = float(entry.get("ratio", 0.0))
    peak_summary = _summarize_detected_peaks(entry)
    positive_ratio = _assay_positive_ratio(assay)

    if assay == "FLT3-ITD":
        if ratio >= positive_ratio:
            return "Positiv FLT3-ITD"
        if peak_summary["mut_bps"]:
            return "Negativ FLT3-ITD - lavniva dokumentert"
        return "Ingen FLT3-ITD pavist"
    if assay == "FLT3-D835":
        if ratio >= positive_ratio:
            return "Positiv FLT3-D835"
        if peak_summary["mut_bps"]:
            return "FLT3-D835 under positiv grense - dokumentert"
        return "Ingen FLT3-D835 pavist"
    if assay == "NPM1":
        if ratio >= positive_ratio:
            return "Positiv NPM1"
        if peak_summary["mut_bps"]:
            return "Mulig NPM1 - vurder manuelt"
        return "Ingen NPM1-mutasjon pavist"
    return "Ingen tolkning"

def generate_flt3_peak_report(entries: list[dict], outdir: Path):
    """Write a flat CSV that mirrors the worksheet-style assay summary."""
    rows = []
    for entry in entries:
        peak_summary = _summarize_detected_peaks(entry)
        rows.append({
            "DIT": entry.get("dit") or "",
            "File": entry["fsa"].file_name,
            "Assay": entry["assay"],
            "Group": entry.get("group") or "",
            "Parallel": entry.get("parallel") or "",
            "Treatment": entry.get("analysis_type") or "",
            "InjectionTime": entry.get("injection_time"),
            "WT_bp": round(float(peak_summary["wt_bp"]), 2) if not np.isnan(peak_summary["wt_bp"]) else "",
            "WT_Area": round(peak_summary["wt_area"], 2),
            "Mutant_bp": ", ".join(f"{bp:.2f}" for bp in peak_summary["mut_bps"]),
            "Mutant_Area": ", ".join(f"{area:.2f}" for area in peak_summary["mut_areas"]),
            "Mutant_Area_Total": round(peak_summary["mut_area_total"], 2),
            "Ratio": round(float(entry.get("ratio", 0.0)), 4),
            "Interpretation": _interpret_entry(entry),
            "LadderQC": entry.get("ladder_qc_status", ""),
            "LadderR2": round(float(entry.get("ladder_r2", np.nan)), 4) if not np.isnan(entry.get("ladder_r2", np.nan)) else "",
        })

    if rows:
        df = pd.DataFrame(rows)
        outdir.mkdir(parents=True, exist_ok=True)
        csv_path = outdir / "Final_Detailed_Peak_Report.csv"
        df.to_csv(csv_path, index=False)
        print_green(f"FLT3 detailed peak report saved to {csv_path}")

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
        wt_area, mut_area = _summarize_peak_areas(e)
        ratio = float(e.get("ratio", 0.0))
        assay_cfg = ASSAY_CONFIG.get(e["assay"], {})
        min_wt_area = float(assay_cfg.get("control_wt_min_area", 0.0))
        min_ratio = _assay_positive_ratio(e["assay"])
        
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
            if not wt_peaks.empty and wt_area >= min_wt_area:
                status = "PASS"
            else:
                status = "FAIL"
                if wt_peaks.empty:
                    details = "No WT peak detected"
                else:
                    details = f"WT area below threshold ({wt_area:.0f} < {min_wt_area:.0f})"
        elif group == "positive_control":
            # ivs-p001: Expect mutant peaks
            if not mut_peaks.empty and (ratio >= min_ratio or wt_peaks.empty):
                status = "PASS"
            else:
                status = "FAIL"
                if mut_peaks.empty:
                    details = "No mutant/ITD peak detected"
                else:
                    details = f"Mutant ratio below threshold ({ratio:.4f} < {min_ratio:.4f})"
                
        qc_data.append({
            "File": e["fsa"].file_name,
            "ControlGroup": group,
            "Assay": e["assay"],
            "InjectionTime": e["injection_time"],
            "WT_Area": round(wt_area, 2),
            "Mutant_Area": round(mut_area, 2),
            "Ratio": round(ratio, 4),
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
    fsa_dir, assay_dir = normalize_pipeline_paths(fsa_dir, base_outdir, assay_folder_name)

    raw_files = _scan_files(fsa_dir, mode=mode)
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
            "protocol_injection_time": c.get("protocol_injection_time", c["injection_time"]),
        })

    if not entries:
        return [] if return_entries else None

    # Calculate ratios for all entries that have mutant/wt peaks
    _calculate_ratios(entries)

    # Generate QC Report
    generate_flt3_qc_report(entries, assay_dir)
    generate_flt3_peak_report(entries, assay_dir)

    return finalize_pipeline_run(
        entries,
        assay_dir,
        return_entries=return_entries,
        make_dit_reports=make_dit_reports,
        mode=mode,
    )
