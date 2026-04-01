from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

ProgressCallback = Callable[[dict[str, Any]], None]
StatusCallback = Callable[[str], None]

DEFAULT_INPUT_ROOT = Path("/Users/christian/Desktop/data/Klonalitet")
DEFAULT_OUTPUT_ROOT = Path("/Users/christian/Desktop/Excel_Fraggler/full_year_runs")
DEFAULT_YEAR_LABEL = "2025"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _month_keys(year_label: str) -> list[str]:
    return [f"{year_label}_{month:02d}" for month in range(1, 13)]


def build_arg_parser(default_year_label: str = DEFAULT_YEAR_LABEL) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the clonality validation harness across all months for a given year, "
            "writing a fresh workbook/state/artifact set per month into a new output tree."
        )
    )
    parser.add_argument(
        "--year-label",
        default=default_year_label,
        help=f"Year label used to discover folders and name outputs. Default: {default_year_label}",
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help="Root directory containing the clonality run folders for the selected year.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Fresh output root for this yearly run.",
    )
    parser.add_argument(
        "--run-name",
        type=str,
        default=None,
        help="Optional run directory name. Defaults to full_<year>_validation_<timestamp>.",
    )
    parser.add_argument(
        "--month",
        dest="months",
        action="append",
        default=[],
        help=f"Optional month key like {default_year_label}_03. Repeat to restrict the run to specific months.",
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


def normalize_month_keys(year_label: str, values: Iterable[str]) -> list[str]:
    year_label = str(year_label).strip()
    if not year_label or not year_label.isdigit():
        raise ValueError(f"Year label must be numeric, got: {year_label!r}")
    requested = [str(value).strip() for value in values if str(value).strip()]
    month_keys = _month_keys(year_label)
    if not requested:
        return month_keys

    normalized: list[str] = []
    for value in requested:
        month_key = value.replace("-", "_")
        if len(month_key) == 7 and month_key.startswith(f"{year_label}_") and month_key[5:7].isdigit():
            if month_key not in month_keys:
                raise ValueError(f"Unsupported month key: {value}")
            normalized.append(month_key)
            continue
        raise ValueError(f"Month must look like {year_label}_03, got: {value}")

    deduped: list[str] = []
    for month_key in normalized:
        if month_key not in deduped:
            deduped.append(month_key)
    return deduped


def discover_month_folders(input_root: Path, year_label: str) -> dict[str, list[Path]]:
    root = Path(input_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Input root not found: {root}")
    month_keys = _month_keys(str(year_label))
    month_map: dict[str, list[Path]] = {month_key: [] for month_key in month_keys}
    for child in sorted(root.iterdir(), key=lambda path: path.name):
        if not child.is_dir():
            continue
        for month_key in month_keys:
            if child.name.startswith(f"{month_key}_"):
                month_map[month_key].append(child.resolve())
                break
    return month_map


def write_month_folder_lists(
    month_map: dict[str, list[Path]],
    output_dir: Path,
    *,
    year_label: str,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for month_key in _month_keys(year_label):
        txt_path = output_dir / f"{month_key}_folders.txt"
        folders = month_map.get(month_key, [])
        txt_path.write_text("".join(f"{folder}\n" for folder in folders), encoding="utf-8")
        written[month_key] = txt_path
    return written


def _build_month_argv(
    args: argparse.Namespace,
    month_key: str,
    month_txt: Path,
    month_runs_dir: Path,
) -> list[str]:
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
    if getattr(sys, "frozen", False):
        _invoke_month_validation_in_process(argv_for_month)
        return
    command = [sys.executable, str(REPO_ROOT / "scripts" / "validate_clonality_hard_cases.py"), *argv_for_month]
    subprocess.run(command, check=True)


def _invoke_month_validation_in_process(argv_for_month: list[str]) -> None:
    from scripts.validate_clonality_hard_cases import run_validation

    run_validation(argv_for_month)


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
    status: str = "done",
) -> dict[str, object]:
    return {
        "month": month_key,
        "folder_count": len(folders),
        "folders_file": str(month_txt),
        "status": status,
        "resumed": resumed,
        "run_dir": str(summary["run_dir"]),
        "summary_json": str(summary["summary_json"]),
        "workbook_path": str(summary["workbook_path"]),
        "state_file": str(summary["state_file"]),
        "artifact_dir": str(summary["artifact_dir"]),
        "candidate_artifact_dir": str(summary["candidate_artifact_dir"]),
        "timing_seconds": summary["timing_seconds"],
    }


def _emit(callback: ProgressCallback | StatusCallback | None, payload: dict[str, Any] | str) -> None:
    if callback is None:
        return
    callback(payload)  # type: ignore[arg-type]


def _write_manifest(
    run_dir: Path,
    input_root: Path,
    output_root: Path,
    year_label: str,
    selected_months: list[str],
    month_lists_dir: Path,
    month_runs_dir: Path,
    month_summaries: list[dict[str, object]],
) -> dict[str, object]:
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "year_label": year_label,
        "input_root": str(input_root),
        "output_root": str(output_root),
        "run_dir": str(run_dir),
        "selected_months": selected_months,
        "month_folder_lists_dir": str(month_lists_dir),
        "month_runs_dir": str(month_runs_dir),
        "months": month_summaries,
    }
    manifest_path = run_dir / f"full_{year_label}_run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def run_yearly_validation(
    *,
    year_label: str = DEFAULT_YEAR_LABEL,
    input_root: Path,
    output_root: Path,
    run_name: str | None = None,
    months: Iterable[str] | None = None,
    max_workers: int | None = None,
    folder_workers: int | None = None,
    refresh_each_folder: bool = False,
    include_sl: bool = False,
    cleanup_staging_root: bool = False,
    resume_existing: bool = False,
    progress_callback: ProgressCallback | None = None,
    status_callback: StatusCallback | None = None,
    invoke_month_validation: Callable[[list[str]], None] = _invoke_month_validation,
) -> dict[str, object]:
    year_label = str(year_label).strip()
    if not year_label.isdigit():
        raise ValueError(f"Year label must be numeric, got: {year_label!r}")

    input_root = Path(input_root).expanduser().resolve()
    output_root = Path(output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    run_name = run_name or f"full_{year_label}_validation_{_timestamp()}"
    run_dir = output_root / run_name
    if run_dir.exists():
        if not resume_existing:
            raise FileExistsError(f"Run directory already exists: {run_dir}")
    else:
        run_dir.mkdir(parents=True, exist_ok=False)

    selected_months = normalize_month_keys(year_label, months or [])
    month_map = discover_month_folders(input_root, year_label)
    month_lists_dir = run_dir / "month_folder_lists"
    month_runs_dir = run_dir / "month_runs"
    month_runs_dir.mkdir(parents=True, exist_ok=True)
    month_list_paths = write_month_folder_lists(month_map, month_lists_dir, year_label=year_label)

    _emit(progress_callback, {"event": "run_started", "year_label": year_label, "run_dir": str(run_dir)})
    _emit(status_callback, f"Starting yearly run for {year_label}")

    month_summaries: list[dict[str, object]] = []
    for month_key in selected_months:
        folders = month_map.get(month_key, [])
        month_txt = month_list_paths[month_key]
        base_record: dict[str, object] = {
            "month": month_key,
            "folder_count": len(folders),
            "folders_file": str(month_txt),
        }
        if not folders:
            month_record = {**base_record, "status": "skipped_empty", "resumed": False}
            month_summaries.append(month_record)
            manifest = _write_manifest(
                run_dir,
                input_root,
                output_root,
                year_label,
                selected_months,
                month_lists_dir,
                month_runs_dir,
                month_summaries,
            )
            _emit(progress_callback, {"event": "month_skipped_empty", "year_label": year_label, "month": month_key, "run_dir": str(run_dir)})
            _emit(progress_callback, {"event": "manifest_written", "year_label": year_label, "run_dir": str(run_dir), "manifest_path": str(run_dir / f"full_{year_label}_run_manifest.json")})
            continue

        existing_summary_path = month_runs_dir / month_key / "run_summary.json"
        if resume_existing and existing_summary_path.is_file():
            _emit(progress_callback, {"event": "month_resumed", "year_label": year_label, "month": month_key, "run_dir": str(run_dir)})
            summary = _load_month_summary(month_runs_dir, month_key)
            month_summaries.append(_build_month_record(month_key, folders, month_txt, summary, resumed=True))
            manifest = _write_manifest(
                run_dir,
                input_root,
                output_root,
                year_label,
                selected_months,
                month_lists_dir,
                month_runs_dir,
                month_summaries,
            )
            _emit(progress_callback, {"event": "manifest_written", "year_label": year_label, "run_dir": str(run_dir), "manifest_path": str(run_dir / f"full_{year_label}_run_manifest.json")})
            continue

        _emit(progress_callback, {"event": "month_started", "year_label": year_label, "month": month_key, "folder_count": len(folders), "run_dir": str(run_dir)})
        _emit(status_callback, f"Running {month_key} ({len(folders)} folders)")
        argv_for_month = _build_month_argv(
            argparse.Namespace(
                max_workers=max_workers,
                folder_workers=folder_workers,
                refresh_each_folder=refresh_each_folder,
                include_sl=include_sl,
                cleanup_staging_root=cleanup_staging_root,
            ),
            month_key,
            month_txt,
            month_runs_dir,
        )
        try:
            invoke_month_validation(argv_for_month)
        except Exception as exc:
            month_summaries.append({**base_record, "status": "error", "resumed": False, "error": str(exc)})
            manifest = _write_manifest(
                run_dir,
                input_root,
                output_root,
                year_label,
                selected_months,
                month_lists_dir,
                month_runs_dir,
                month_summaries,
            )
            _emit(progress_callback, {"event": "run_failed", "year_label": year_label, "month": month_key, "run_dir": str(run_dir), "error": str(exc), "manifest_path": str(run_dir / f"full_{year_label}_run_manifest.json")})
            _emit(status_callback, f"Yearly run failed in {month_key}")
            raise
        summary = _load_month_summary(month_runs_dir, month_key)
        month_summaries.append(_build_month_record(month_key, folders, month_txt, summary, resumed=False))
        _emit(progress_callback, {"event": "month_finished", "year_label": year_label, "month": month_key, "run_dir": str(run_dir), "summary_json": str(summary["summary_json"])})
        manifest = _write_manifest(
            run_dir,
            input_root,
            output_root,
            year_label,
            selected_months,
            month_lists_dir,
            month_runs_dir,
            month_summaries,
        )
        _emit(progress_callback, {"event": "manifest_written", "year_label": year_label, "run_dir": str(run_dir), "manifest_path": str(run_dir / f"full_{year_label}_run_manifest.json")})

    manifest = _write_manifest(
        run_dir,
        input_root,
        output_root,
        year_label,
        selected_months,
        month_lists_dir,
        month_runs_dir,
        month_summaries,
    )
    manifest_path = run_dir / f"full_{year_label}_run_manifest.json"
    _emit(progress_callback, {"event": "run_finished", "year_label": year_label, "run_dir": str(run_dir), "manifest_path": str(manifest_path)})
    _emit(status_callback, f"Finished yearly run for {year_label}")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def run_yearly_validation_from_argv(argv: list[str] | None = None) -> dict[str, object]:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
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
    )


def main(argv: list[str] | None = None) -> int:
    run_yearly_validation_from_argv(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
