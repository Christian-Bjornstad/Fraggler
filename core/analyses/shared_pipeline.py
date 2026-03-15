"""
Shared helpers for analysis pipeline implementations.
"""
from __future__ import annotations

from pathlib import Path

from fraggler.fraggler import print_green, print_warning

from core.html_reports import build_dit_html_reports
from core.utils import is_control_file, is_water_file


def normalize_pipeline_paths(
    fsa_dir: Path,
    base_outdir: Path | None,
    assay_folder_name: str | None,
) -> tuple[Path, Path]:
    """Resolve common input/output locations for all analyses."""
    fsa_dir = Path(fsa_dir).expanduser()
    base_outdir = Path(base_outdir or fsa_dir).expanduser()
    assay_dir = base_outdir / (assay_folder_name or "REPORTS").strip()
    return fsa_dir, assay_dir


def scan_fsa_files(
    fsa_dir: Path,
    mode: str = "all",
    *,
    include_controls_only: bool = False,
    recursive: bool = False,
) -> list[Path]:
    """Scan a folder for .fsa files with shared filtering semantics."""
    if not fsa_dir.exists():
        print_warning(f"FSA-katalog finnes ikke: {fsa_dir}")
        return []

    iterator = fsa_dir.rglob("*.fsa") if recursive else fsa_dir.glob("*.fsa")
    fsa_files = [
        p for p in sorted(iterator)
        if not is_water_file(p.name)
    ]

    controls_mode = mode == "controls" or include_controls_only
    if controls_mode:
        fsa_files = [p for p in fsa_files if is_control_file(p.name)]
        print_green(f"[INFO] Controls mode: {len(fsa_files)} control files selected.")

    if not fsa_files:
        print_warning(f"Fant ingen .fsa-filer i {fsa_dir}.")
    return fsa_files


def finalize_pipeline_run(
    entries: list[dict],
    assay_dir: Path,
    *,
    return_entries: bool,
    make_dit_reports: bool,
    mode: str,
) -> list[dict] | None:
    """Apply the common report-generation and return-value contract."""
    if not entries:
        return [] if return_entries else None

    if make_dit_reports and mode != "controls":
        build_dit_html_reports(entries, assay_dir)

    return entries if return_entries else None
