"""
Fraggler Diagnostics — Pipeline Dispatcher.
"""
from __future__ import annotations
import inspect
from pathlib import Path
from core.analyses.registry import get_analysis_module

def _scan_files(fsa_dir: Path, mode: str = "all") -> list[Path]:
    """Compatibility wrapper for tests and shared callers."""
    mod = get_analysis_module("pipeline")
    scanner = getattr(mod, "_scan_files", None)
    if scanner is None:
        return []
    params = inspect.signature(scanner).parameters
    if "mode" in params:
        return scanner(fsa_dir, mode=mode)
    return scanner(fsa_dir)

def run_pipeline(
    fsa_dir: Path,
    base_outdir: Path | None = None,
    assay_folder_name: str | None = None,
    return_entries: bool = False,
    make_dit_reports: bool = True,
    mode: str = "all",
) -> list[dict] | None:
    """Delegates pipeline execution to the active analysis module."""
    mod = get_analysis_module("pipeline")
    return mod.run_pipeline(
        fsa_dir=fsa_dir,
        base_outdir=base_outdir,
        assay_folder_name=assay_folder_name,
        return_entries=return_entries,
        make_dit_reports=make_dit_reports,
        mode=mode,
    )
