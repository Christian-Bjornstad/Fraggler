from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import validate_clonality_hard_cases as harness


DEFAULT_INPUT_ROOT = Path("/Users/christian/Desktop/data/Klonalitet/2025_data")
DEFAULT_OUTPUT_ROOT = Path("/Users/christian/Desktop/Excel_Fraggler/full_2025_runs")
MONTH_KEYS = [f"2025_{month:02d}" for month in range(1, 13)]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the clonality validation harness across all 2025 months, "
            "writing a fresh workbook/state/artifact set per month into a new output tree."
        )
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=f"Root directory containing the 2025 clonality run folders. Default: {DEFAULT_INPUT_ROOT}",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Fresh output root for this full-year run. Default: {DEFAULT_OUTPUT_ROOT}",
    )
    parser.add_argument(
        "--run-name",
        type=str,
        default=None,
        help="Optional run directory name. Defaults to full_2025_validation_<timestamp>.",
    )
    parser.add_argument(
        "--month",
        dest="months",
        action="append",
        default=[],
        help="Optional month key like 2025_03. Repeat to restrict the run to specific months.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Optional backfill max worker count passed through to the validation harness.",
    )
    parser.add_argument(
        "--folder-workers",
        type=int,
        default=None,
        help="Optional same-month folder concurrency passed through to the validation harness.",
    )
    parser.add_argument(
        "--refresh-each-folder",
        action="store_true",
        help="Refresh the tracking workbook after each folder instead of deferring to month boundaries.",
    )
    parser.add_argument(
        "--include-sl",
        action="store_true",
        help="Include SL rows in exported feature/candidate artifacts.",
    )
    parser.add_argument(
        "--cleanup-staging-root",
        action="store_true",
        help="Delete each month staging root after that month completes.",
    )
    parser.add_argument(
        "--resume-existing",
        action="store_true",
        help="Resume an existing run directory instead of failing when it already exists.",
    )
    return parser


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def normalize_month_keys(values: Iterable[str]) -> list[str]:
    requested = [str(value).strip() for value in values if str(value).strip()]
    if not requested:
        return list(MONTH_KEYS)
    normalized: list[str] = []
    for value in requested:
        month_key = value.replace("-", "_")
        if len(month_key) == 7 and month_key.startswith("2025_") and month_key[5:7].isdigit():
            if month_key not in MONTH_KEYS:
                raise ValueError(f"Unsupported month key: {value}")
            normalized.append(month_key)
            continue
        raise ValueError(f"Month must look like 2025_03, got: {value}")
    deduped: list[str] = []
    for month_key in normalized:
        if month_key not in deduped:
            deduped.append(month_key)
    return deduped


def discover_month_folders(input_root: Path) -> dict[str, list[Path]]:
    root = Path(input_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Input root not found: {root}")
    month_map: dict[str, list[Path]] = {month_key: [] for month_key in MONTH_KEYS}
    for child in sorted(root.iterdir(), key=lambda path: path.name):
        if not child.is_dir():
            continue
        for month_key in MONTH_KEYS:
            if child.name.startswith(f"{month_key}_"):
                month_map[month_key].append(child.resolve())
                break
    return month_map


def write_month_folder_lists(month_map: dict[str, list[Path]], output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for month_key in MONTH_KEYS:
        txt_path = output_dir / f"{month_key}_folders.txt"
        folders = month_map.get(month_key, [])
        txt_path.write_text(
            "".join(f"{folder}\n" for folder in folders),
            encoding="utf-8",
        )
        written[month_key] = txt_path
    return written


def _build_month_argv(args: argparse.Namespace, month_key: str, month_txt: Path, month_runs_dir: Path) -> list[str]:
    argv_for_month = [
        "--folders-file",
        str(month_txt),
        "--output-dir",
        str(month_runs_dir),
        "--run-name",
        month_key,
    ]
    if args.max_workers is not None:
        argv_for_month.extend(["--max-workers", str(args.max_workers)])
    if args.folder_workers is not None:
        argv_for_month.extend(["--folder-workers", str(args.folder_workers)])
    if args.refresh_each_folder:
        argv_for_month.append("--refresh-each-folder")
    if args.include_sl:
        argv_for_month.append("--include-sl")
    if args.cleanup_staging_root:
        argv_for_month.append("--cleanup-staging-root")
    return argv_for_month


def _invoke_month_validation(argv_for_month: list[str]) -> None:
    command = [sys.executable, str(REPO_ROOT / "scripts" / "validate_clonality_hard_cases.py"), *argv_for_month]
    subprocess.run(command, check=True)


def _load_month_summary(month_runs_dir: Path, month_key: str) -> dict[str, object]:
    summary_path = month_runs_dir / month_key / "run_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Missing month summary after validation run: {summary_path}")
    return json.loads(summary_path.read_text(encoding="utf-8"))


def _build_month_record(
    month_key: str,
    folders: list[Path],
    month_txt: Path,
    summary: dict[str, object],
    *,
    resumed: bool,
) -> dict[str, object]:
    return {
        "month": month_key,
        "folder_count": len(folders),
        "folders_file": str(month_txt),
        "status": "done",
        "resumed": resumed,
        "run_dir": str(summary["run_dir"]),
        "summary_json": str(summary["summary_json"]),
        "workbook_path": str(summary["workbook_path"]),
        "state_file": str(summary["state_file"]),
        "artifact_dir": str(summary["artifact_dir"]),
        "candidate_artifact_dir": str(summary["candidate_artifact_dir"]),
        "timing_seconds": summary["timing_seconds"],
    }


def _write_manifest(
    run_dir: Path,
    input_root: Path,
    output_root: Path,
    selected_months: list[str],
    month_lists_dir: Path,
    month_runs_dir: Path,
    month_summaries: list[dict[str, object]],
) -> dict[str, object]:
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(input_root),
        "output_root": str(output_root),
        "run_dir": str(run_dir),
        "selected_months": selected_months,
        "month_folder_lists_dir": str(month_lists_dir),
        "month_runs_dir": str(month_runs_dir),
        "months": month_summaries,
    }
    manifest_path = run_dir / "full_2025_run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def run_full_2025_validation(argv: list[str] | None = None) -> dict[str, object]:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    input_root = Path(args.input_root).expanduser().resolve()
    output_root = Path(args.output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    run_name = args.run_name or f"full_2025_validation_{_timestamp()}"
    run_dir = output_root / run_name
    if run_dir.exists():
        if not args.resume_existing:
            raise FileExistsError(f"Run directory already exists: {run_dir}")
    else:
        run_dir.mkdir(parents=True, exist_ok=False)

    month_map = discover_month_folders(input_root)
    selected_months = normalize_month_keys(args.months)
    month_lists_dir = run_dir / "month_folder_lists"
    month_runs_dir = run_dir / "month_runs"
    month_runs_dir.mkdir(parents=True, exist_ok=True)
    month_list_paths = write_month_folder_lists(month_map, month_lists_dir)

    month_summaries: list[dict[str, object]] = []
    for month_key in selected_months:
        folders = month_map.get(month_key, [])
        month_txt = month_list_paths[month_key]
        month_record: dict[str, object] = {
            "month": month_key,
            "folder_count": len(folders),
            "folders_file": str(month_txt),
            "status": "skipped_empty",
        }
        if not folders:
            month_summaries.append(month_record)
            _write_manifest(
                run_dir,
                input_root,
                output_root,
                selected_months,
                month_lists_dir,
                month_runs_dir,
                month_summaries,
            )
            continue

        existing_summary_path = month_runs_dir / month_key / "run_summary.json"
        if args.resume_existing and existing_summary_path.is_file():
            summary = _load_month_summary(month_runs_dir, month_key)
            month_summaries.append(
                _build_month_record(month_key, folders, month_txt, summary, resumed=True)
            )
            _write_manifest(
                run_dir,
                input_root,
                output_root,
                selected_months,
                month_lists_dir,
                month_runs_dir,
                month_summaries,
            )
            continue

        argv_for_month = _build_month_argv(args, month_key, month_txt, month_runs_dir)
        _invoke_month_validation(argv_for_month)
        summary = _load_month_summary(month_runs_dir, month_key)
        month_summaries.append(
            _build_month_record(month_key, folders, month_txt, summary, resumed=False)
        )
        _write_manifest(
            run_dir,
            input_root,
            output_root,
            selected_months,
            month_lists_dir,
            month_runs_dir,
            month_summaries,
        )

    manifest = _write_manifest(
        run_dir,
        input_root,
        output_root,
        selected_months,
        month_lists_dir,
        month_runs_dir,
        month_summaries,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def main(argv: list[str] | None = None) -> int:
    run_full_2025_validation(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
