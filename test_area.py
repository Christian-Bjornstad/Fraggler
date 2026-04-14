import sys
from pathlib import Path
sys.path.insert(0, "/Users/christian/Desktop/OUS")
from core.runner import run_pipeline_job_collect

files = list(Path("/Volumes/T7 Shield/DATA/flt3/2026").rglob("*.fsa"))[:5]
print("Running pipeline on", len(files), "files")
entries = run_pipeline_job_collect(
    fsa_dir=None,
    base_outdir=Path("/tmp/test_out"),
    out_folder_name="tmp",
    scope="all",
    needle="",
    files=files,
    chunk_files=False,
    update_tracking_workbook=False
)
for e in entries:
    print(e.get("file_name"), "r2:", e.get("ladder_fit_metrics", {}).get("r2"))
    rat = e.get("ratio_resolution", {})
    print("  WT Area:", rat.get("selected_wt_area"), "MUT Area:", rat.get("selected_mutant_area"), "Ratio:", rat.get("ratio"))
