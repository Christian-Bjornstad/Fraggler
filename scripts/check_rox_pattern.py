import os
from pathlib import Path
import numpy as np
from Bio import SeqIO
from scipy.signal import find_peaks

f = "/Users/christian/Desktop/OUS/data/flt3/data2 flt3/IVS-0000_ITD_050326_C01_H920GFLJ.fsa"

record = SeqIO.read(f, "abi")
tags = record.annotations.get("abif_raw", {})
rox = np.asarray(tags.get("DATA4", [])).astype(float)

peaks, props = find_peaks(rox, height=200, distance=20)
peak_heights = props["peak_heights"]

print(f"File: {Path(f).name}")
print(f"{'Index':<10} | {'Height':<10} | {'Distance to next':<10}")
print("-" * 40)
for i in range(len(peaks)-1):
    print(f"{peaks[i]:<10} | {peak_heights[i]:<10.0f} | {peaks[i+1]-peaks[i]:<10}")
print(f"{peaks[-1]:<10} | {peak_heights[-1]:<10.0f} | -")
