import os
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from fraggler.fraggler import FsaFile, find_size_standard_peaks, return_maxium_allowed_distance_between_size_standard_peaks, generate_combinations, calculate_best_combination_of_size_standard_peaks, fit_size_standard_to_ladder

def _robust_analyse_fsa_rox(fsa_path: Path, sample_channel: str) -> FsaFile | None:
    ladder_name = "GS500ROX"
    configs = [
        {"min_h": 200, "min_d": 20},
        {"min_h": 100, "min_d": 15},
    ]
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
            fsa = find_size_standard_peaks(fsa)
            ss_peaks = getattr(fsa, "size_standard_peaks", None)
            if ss_peaks is None or ss_peaks.shape[0] < 3:
                continue
            fsa = return_maxium_allowed_distance_between_size_standard_peaks(fsa, multiplier=3)
            for _ in range(30):
                fsa = generate_combinations(fsa)
                if getattr(fsa, "best_size_standard_combinations", None) is not None:
                    if fsa.best_size_standard_combinations.shape[0] > 0:
                        break
                fsa.maxium_allowed_distance_between_size_standard_peaks += 15
            best = getattr(fsa, "best_size_standard_combinations", None)
            if best is None or best.shape[0] == 0:
                continue
            fsa = calculate_best_combination_of_size_standard_peaks(fsa)
            fsa = fit_size_standard_to_ladder(fsa)
            if getattr(fsa, "fitted_to_model", False):
                return fsa
        except:
            continue
    return None

def analyze_peaks(fsa: FsaFile, channel: str):
    trace = np.asarray(fsa.fsa[channel]).astype(float)
    sample_data = getattr(fsa, "sample_data_with_basepairs", None)
    if sample_data is None:
        return None
    
    time_all = sample_data["time"].astype(int).to_numpy()
    bp_all = sample_data["basepairs"].to_numpy()
    
    # D835 window
    mask = (bp_all >= 50.0) & (bp_all <= 200.0)
    y_win = trace[time_all[mask]]
    bp_win = bp_all[mask]
    
    from scipy.signal import find_peaks
    p_idx, _ = find_peaks(y_win, height=200, distance=20)
    
    peaks = []
    for idx in p_idx:
        peaks.append({"bp": bp_win[idx], "height": y_win[idx]})
    return peaks

dirs = [
    "/Users/christian/Desktop/OUS/data/flt3/data flt3",
    "/Users/christian/Desktop/OUS/data/flt3/data2 flt3"
]

results = []
for d in dirs:
    d_path = Path(d)
    if not d_path.exists(): 
        print(f"Directory {d} does not exist")
        continue
    print(f"Scanning directory: {d_path.name}")
    for f in sorted(d_path.glob("*D835*.fsa")):
        print(f"  Processing file: {f.name}")
        fsa = _robust_analyse_fsa_rox(f, "DATA3")
        if fsa:
            peaks = analyze_peaks(fsa, "DATA3")
            results.append({
                "file": f.name,
                "dir": d_path.name,
                "peaks": peaks
            })
        else:
            print(f"    Failed to analyze {f.name}")

for res in results:
    print(f"File: {res['file']} ({res['dir']})")
    for p in res['peaks']:
        print(f"  Peak: {p['bp']:.2f} bp, Height: {p['height']:.0f}")
