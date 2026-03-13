
import sys
from pathlib import Path

# Add project root to sys.path at the front
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import APP_SETTINGS
from core.pipeline import run_pipeline
from core.analyses.registry import get_active_analysis_name

def test_flt3_pipeline():
    import fraggler
    import fraggler.ladders
    print(f"Fraggler path: {fraggler.__file__}")
    print(f"LADDERS: {list(fraggler.ladders.LADDERS.keys())}")
    print(f"Current active analysis: {get_active_analysis_name()}")
    
    # 1. Set FLT3 as active
    APP_SETTINGS["active_analysis"] = "flt3"
    print(f"Switched to: {get_active_analysis_name()}")
    
    # 2. Run pipeline on rerun data
    datasets = ["rerun_all", "rerun_all_2"]
    out_dir = Path("TEST_REPORTS_RERUN").absolute()
    
    for ds in datasets:
        fsa_dir = Path(f"data/flt3/data/{ds}").absolute()
        print(f"\n--- Running pipeline on {ds} ---")
        
        entries = run_pipeline(
            fsa_dir=fsa_dir,
            base_outdir=out_dir,
            make_dit_reports=True,
            return_entries=True,
        )
        
        if entries:
            print(f"Success! Processed {len(entries)} entries from {ds}.")
            for e in entries:
                ratio = e.get("ratio", "N/A")
                par = e.get("parallel", "none")
                print(f" - {e['fsa'].file_name} ({par}): Ratio: {ratio}")
        else:
            print(f"Pipeline failed for {ds}.")
        
    if entries:
        print(f"Success! Processed {len(entries)} entries.")
        # Group by DIT manually for printing if needed
        for e in entries:
            ratio = e.get("ratio", "N/A")
            dit = e.get("dit", "NO_DIT")
            print(f" - {e['fsa'].file_name}: {e['assay']} ({e['group']}), DIT: {dit}, Ratio: {ratio}")
            
        # Check if reports exist
        reports = list(out_dir.glob("REPORTS/*.html"))
        print(f"Generated {len(reports)} reports:")
        for r in reports:
            print(f" - {r.name}")
    else:
        print("Pipeline failed to process any files.")

if __name__ == "__main__":
    test_flt3_pipeline()
