from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

from scripts.run_clonality_yearly import (
    build_arg_parser as _build_generic_arg_parser,
    discover_month_folders as _discover_month_folders,
    normalize_month_keys as _normalize_month_keys,
    _invoke_month_validation as _invoke_month_validation,
    run_yearly_validation,
    write_month_folder_lists as _write_month_folder_lists,
)

DEFAULT_INPUT_ROOT = Path("/Users/christian/Desktop/data/Klonalitet/2025_data")
DEFAULT_OUTPUT_ROOT = Path("/Users/christian/Desktop/Excel_Fraggler/full_2025_runs")
MONTH_KEYS = [f"2025_{month:02d}" for month in range(1, 13)]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = _build_generic_arg_parser(default_year_label="2025")
    parser.description = (
        "Run the clonality validation harness across all 2025 months, "
        "writing a fresh workbook/state/artifact set per month into a new output tree."
    )
    parser.set_defaults(
        input_root=DEFAULT_INPUT_ROOT,
        output_root=DEFAULT_OUTPUT_ROOT,
        year_label="2025",
    )
    return parser


def normalize_month_keys(values: Iterable[str]) -> list[str]:
    return _normalize_month_keys("2025", values)


def discover_month_folders(input_root: Path) -> dict[str, list[Path]]:
    return _discover_month_folders(input_root, "2025")


def write_month_folder_lists(month_map: dict[str, list[Path]], output_dir: Path) -> dict[str, Path]:
    return _write_month_folder_lists(month_map, output_dir, year_label="2025")


def run_full_2025_validation(argv: list[str] | None = None) -> dict[str, object]:
    args = build_arg_parser().parse_args(argv)
    return run_yearly_validation(
        year_label=args.year_label,
        input_root=args.input_root,
        output_root=args.output_root,
        run_name=args.run_name,
        months=args.months,
        max_workers=args.max_workers,
        folder_workers=args.folder_workers,
        refresh_each_folder=args.refresh_each_folder,
        include_sl=args.include_sl,
        cleanup_staging_root=args.cleanup_staging_root,
        resume_existing=args.resume_existing,
        invoke_month_validation=_invoke_month_validation,
    )


def main(argv: list[str] | None = None) -> int:
    run_full_2025_validation(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
