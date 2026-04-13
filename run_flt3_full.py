import os
import sys
from pathlib import Path
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# Add source code to path
sys.path.insert(0, "/Users/christian/Desktop/OUS")
from config import APP_SETTINGS

# Force FLT3 analysis active
APP_SETTINGS["active_analysis"] = "flt3"
from core.runner import run_pipeline_job_collect

def _process_chunk(chunk_files):
    # This runs in a worker process
    try:
        entries = run_pipeline_job_collect(
            fsa_dir=None,  # Or parent if needed, but we pass files
            base_outdir=Path("/tmp"),
            out_folder_name="tmp",
            scope="all",
            needle="",
            files=chunk_files,
            chunk_files=False,
            update_tracking_workbook=False
        )
        
        results = []
        for e in entries:
            fsa = e.get("fsa")
            file_name = getattr(fsa, "file_name", "Unknown") if fsa else "Unknown"
            
            # Start flattening the data.
            row = {
                "file_name": file_name,
                "assay": e.get("assay"),
                "status": "success",
            }
            
            # Ladder metrics
            if fsa:
                ladder_metrics = getattr(fsa, "ladder_fit_metrics", {})
                row["r2"] = ladder_metrics.get("r2")
                row["max_residual_bp"] = ladder_metrics.get("max_abs_error_bp")
                row["ladder_fitted_step_count"] = getattr(fsa, "ladder_fitted_step_count", 0)
                row["ladder_review_required"] = getattr(fsa, "ladder_review_required", False)

            # Ratio / Mutant data
            if "ratio_resolution" in e:
                res = e["ratio_resolution"]
                row["ratio_mode"] = res.get("ratio_mode")
                row["ratio"] = res.get("ratio")
                row["mutant_fraction"] = res.get("mutant_fraction")
                row["wt_bp"] = res.get("selected_wt_bp")
                row["wt_area"] = res.get("selected_wt_area")
                
                # Combine mutant bps if multiple
                mut_bps = res.get("selected_mutant_bps", [])
                row["mutant_bps"] = ", ".join(map(str, mut_bps)) if mut_bps else None
                row["mutant_area"] = res.get("selected_mutant_area")
                
            # Generic pipeline output
            for k, v in e.items():
                # skip heavy stuff
                if k in ["fsa", "trace_data", "peaks", "peaks_by_channel", "ladder_trace_data", "ratio_resolution", "manual_ratio_selection"]:
                    continue
                if k not in row:
                    row[k] = str(v) if isinstance(v, (dict, list)) else v
                    
            results.append(row)
        return results
    except Exception as exc:
        print(f"Error processing chunk: {exc}")
        return [{"file_name": p.name, "status": f"failed: {exc}"} for p in chunk_files]

def main():
    root = "/Volumes/T7 Shield/DATA/flt3"
    files = []
    print(f"Scanning {root} for FSA files...")
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.lower().endswith(".fsa") and "V_" not in f and "Vann" not in f:
                path_obj = Path(os.path.join(dirpath, f))
                # Quick skip of water/controls to focus on actual samples if needed,
                # but "V_" and "Vann" filter achieves most of it.
                files.append(path_obj)
    
    # Process in chunks
    CHUNK_SIZE = 100
    chunks = [files[i:i + CHUNK_SIZE] for i in range(0, len(files), CHUNK_SIZE)]
    print(f"Found {len(files)} FSA files. Divided into {len(chunks)} chunks.")
    
    results = []
    # Using ProcessPoolExecutor to distribute the workload
    # We use fewer workers to avoid memory overload with heavy FSA processing
    max_w = min(8, os.cpu_count() or 4)
    with ProcessPoolExecutor(max_workers=max_w) as executor:
        futures = {executor.submit(_process_chunk, chunk): chunk for chunk in chunks}
        for i, future in enumerate(as_completed(futures), 1):
            chunk_res = future.result()
            results.extend(chunk_res)
            print(f"Processed chunk {i}/{len(chunks)}...")
            
    df = pd.DataFrame(results)
    out_path = "/Users/christian/Desktop/flt3_full_pipeline_t7.xlsx"
    df.to_excel(out_path, index=False)
    print(f"Done! Pipeline extracted {len(results)} rows. Saved to {out_path}")

if __name__ == "__main__":
    main()
