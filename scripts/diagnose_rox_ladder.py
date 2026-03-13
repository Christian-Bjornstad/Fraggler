import os
from pathlib import Path
import numpy as np
from Bio import SeqIO
from fraggler.fraggler import FsaFile, find_size_standard_peaks, return_maxium_allowed_distance_between_size_standard_peaks

failing_files = [
    "/Users/christian/Desktop/OUS/data/flt3/data flt3/25OUM00002_D835_050326_E03_H920GFLJ.fsa",
    "/Users/christian/Desktop/OUS/data/flt3/data2 flt3/IVS-0000_D835_050326_C03_H920GFLJ.fsa",
    "/Users/christian/Desktop/OUS/data/flt3/data2 flt3/IVS-P001_D835_050326_E03_H920GFLJ.fsa"
]

print(f"{'File':<50} | {'Peaks Found':<12} | {'Max Peak H':<12} | {'Avg Peak H':<12}")
print("-" * 100)

for f_str in failing_files:
    f_path = Path(f_str)
    if not f_path.exists():
        print(f"File not found: {f_path.name}")
        continue
    
    try:
        # Try multiple thresholds to see what's happening
        results = []
        for thresh in [200, 100, 50, 20]:
            fsa = FsaFile(
                file=str(f_path),
                ladder="GS500ROX",
                sample_channel="DATA3",
                min_distance_between_peaks=15,
                min_size_standard_height=thresh,
                size_standard_channel="DATA4",
            )
            fsa = find_size_standard_peaks(fsa)
            ss_peaks = getattr(fsa, "size_standard_peaks", None)
            
            n_peaks = len(ss_peaks) if ss_peaks is not None else 0
            
            # Get peak heights
            rox_data = np.asarray(fsa.fsa["DATA4"]).astype(float)
            peak_heights = rox_data[ss_peaks] if n_peaks > 0 else [0]
            max_h = np.max(peak_heights) if n_peaks > 0 else 0
            avg_h = np.mean(peak_heights) if n_peaks > 0 else 0
            
            results.append((thresh, n_peaks, max_h, avg_h))
        
        # Print info for the highest threshold that found anything, or just the 200 one
        t, n, mh, ah = results[0] # Default to 200
        print(f"{f_path.name:<50} | {n:<12} | {mh:<12.0f} | {ah:<12.0f} (Thresh: 200)")
        for t, n, mh, ah in results[1:]:
            print(f"{'':<50} | {n:<12} | {mh:<12.0f} | {ah:<12.0f} (Thresh: {t})")

    except Exception as e:
        print(f"{f_path.name:<50} | Error: {e}")
