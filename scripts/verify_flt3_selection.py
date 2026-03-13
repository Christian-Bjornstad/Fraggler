import os
from pathlib import Path
import pandas as pd
from core.analyses.flt3.pipeline import run_pipeline

# Test on a subset of the rerun data to verify selection logic
fsa_dir_3s = Path("/Users/christian/Desktop/OUS/data/flt3/data/rerun_all")
fsa_dir_1s = Path("/Users/christian/Desktop/OUS/data/flt3/data/rerun_all_2")

# We want to see if the pipeline correctly picks from both if we give it a parent dir
# Or we can run twice and check the QC report

outdir = Path("/tmp/flt3_test_results")

print("Running pipeline on 3s injections (rerun_all)...")
entries_3s = run_pipeline(fsa_dir_3s, base_outdir=outdir, assay_folder_name="REPORTS_3S", return_entries=True)

print("\nRunning pipeline on 1s injections (rerun_all_2)...")
entries_1s = run_pipeline(fsa_dir_1s, base_outdir=outdir, assay_folder_name="REPORTS_1S", return_entries=True)

# Check QC report
qc_3s = pd.read_csv(outdir / "REPORTS_3S" / "QC_FLT3_Injections.csv")
qc_1s = pd.read_csv(outdir / "REPORTS_1S" / "QC_FLT3_Injections.csv")

print("\nQC Report (3s) - Top 5:")
print(qc_3s.head())

print("\nQC Report (1s) - Top 5:")
print(qc_1s.head())

# Verify D835 peaks in 3s
d835_3s = [e for e in entries_3s if e["assay"] == "FLT3-D835"]
if d835_3s:
    print(f"\nD835 (3s) Sample: {d835_3s[0]['fsa'].file_name}")
    peaks = d835_3s[0]["peaks_by_channel"][d835_3s[0]["primary_peak_channel"]]
    print(peaks[peaks.label.isin(["WT", "MUT"])])
