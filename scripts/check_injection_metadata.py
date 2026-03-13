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
    "/Users/christian/Desktop/OUS/data/flt3/data/rerun_all",
    "/Users/christian/Desktop/OUS/data/flt3/data/rerun_all_2"
]

for d in dirs:
    d_path = Path(d)
    if not d_path.exists():
        print(f"Directory {d} does not exist")
        continue
    f = next(d_path.glob("*.fsa"))
    meta = get_injection_metadata(f)
    print(f"Directory: {d_path.name} | File: {f.name} | Inj Time: {meta['injection_time']}s")
