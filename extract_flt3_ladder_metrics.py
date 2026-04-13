import os
import sys
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, "/Users/christian/Desktop/OUS")

from config import APP_SETTINGS
APP_SETTINGS["active_analysis"] = "flt3"

from fraggler.fraggler import FsaFile
from core.analysis import analyse_fsa_rox

def process_file(fsa_path_str):
    fsa_path = Path(fsa_path_str)
    try:
        fsa = FsaFile(fsa_path_str, ladder="GS500ROX")
        if not hasattr(fsa, "assay") or not fsa.assay:
            fsa.assay = "FLT3-ITD"
            
        fsa = analyse_fsa_rox(
            fsa,
            ladder_fit_profile="flt3_gs500rox",
            ladder_candidate_count=1,
            ladder_fit_rescue=False
        )
        
        metrics = getattr(fsa, "ladder_fit_metrics", {})
        
        return {
            "file_name": fsa_path.name,
            "folder_path": str(fsa_path.parent),
            "r2": metrics.get("r2"),
            "max_residual_bp": metrics.get("max_abs_error_bp"),
            "mean_residual_bp": metrics.get("mean_abs_error_bp"),
            "fitted_steps": getattr(fsa, "ladder_fitted_step_count", 0),
            "status": "success"
        }
    except Exception as e:
        return {
            "file_name": fsa_path.name,
            "folder_path": str(fsa_path.parent),
            "r2": None,
            "max_residual_bp": None,
            "mean_residual_bp": None,
            "fitted_steps": 0,
            "status": f"error: {str(e)}"
        }

def main():
    root = "/Volumes/T7 Shield/DATA/flt3"
    files = []
    print(f"Scanning {root} for FSA files...")
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.lower().endswith(".fsa") and "V_" not in f and "Vann" not in f:
                files.append(os.path.join(dirpath, f))
    
    print(f"Found {len(files)} FSA files. Processing...")
    
    results = []
    # Avoid creating too many processes if there are tons of files
    max_w = min(12, os.cpu_count() or 4)
    with ProcessPoolExecutor(max_workers=max_w) as executor:
        futures = {executor.submit(process_file, f): f for f in files}
        for i, future in enumerate(as_completed(futures), 1):
            results.append(future.result())
            if i % 100 == 0:
                print(f"Processed {i}/{len(files)} files...")
                
    df = pd.DataFrame(results)
    out_path = "/Users/christian/Desktop/flt3_ladder_metrics_t7.xlsx"
    df.to_excel(out_path, index=False)
    print(f"Done! Saved results to {out_path}")

if __name__ == "__main__":
    main()
