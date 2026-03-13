import os
from pathlib import Path
import numpy as np
from Bio import SeqIO

dirs = [
    "/Users/christian/Desktop/OUS/data/flt3/data flt3",
    "/Users/christian/Desktop/OUS/data/flt3/data2 flt3"
]

print(f"{'Directory':<15} | {'File':<50} | {'Max DATA4':<10}")
print("-" * 80)

for d in dirs:
    d_path = Path(d)
    if not d_path.exists(): continue
    for f in sorted(d_path.glob("*D835*.fsa")):
        try:
            record = SeqIO.read(str(f), "abi")
            tags = record.annotations.get("abif_raw", {})
            rox = np.asarray(tags.get("DATA4", [])).astype(float)
            if len(rox) > 0:
                print(f"{d_path.name:<15} | {f.name:<50} | {np.max(rox):<10.1f}")
            else:
                print(f"{d_path.name:<15} | {f.name:<50} | No DATA4")
        except Exception as e:
            print(f"{d_path.name:<15} | {f.name:<50} | Error: {e}")
