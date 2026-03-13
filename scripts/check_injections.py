import os
from pathlib import Path
from Bio import SeqIO

def get_injection_metadata(fsa_path: Path) -> dict:
    try:
        record = SeqIO.read(str(fsa_path), "abi")
        tags = record.annotations.get("abif_raw", {})
        return {
            "injection_time": tags.get("InSc1", 0),
            "injection_voltage": tags.get("InVt1", 0),
        }
    except Exception as e:
        return {"injection_time": 0, "injection_voltage": 0, "error": str(e)}

dirs = [
    "/Users/christian/Desktop/OUS/data/flt3/data flt3",
    "/Users/christian/Desktop/OUS/data/flt3/data2 flt3"
]

print(f"{'Directory':<40} | {'File':<50} | {'Inj Time (s)':<12}")
print("-" * 110)

for d in dirs:
    d_path = Path(d)
    if not d_path.exists():
        continue
    for f in sorted(d_path.glob("*.fsa")):
        meta = get_injection_metadata(f)
        print(f"{d_path.name:<40} | {f.name:<50} | {meta['injection_time']:<12}")
