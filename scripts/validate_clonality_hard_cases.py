from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import threading


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from core.analyses.clonality.feature_artifacts import write_clonality_feature_artifacts
from core.analyses.clonality.candidate_artifacts import write_clonality_candidate_artifacts
from core.clonality_backfill import run_clonality_backfill


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run a focused clonality validation against a specified set of 2025 folders, "
            "capture timings/output paths, and export feature artifacts."
        )
    )
    parser.add_argument(
        "--folder",
        dest="folders",
        action="append",
        type=Path,
        default=[],
        help="Input 2025 clonality folder. Repeat for multiple folders.",
    )
    parser.add_argument(
        "--folders-file",
        type=Path,
        default=None,
        help="Optional text file with one folder path per line.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where the validation run outputs will be written.",
    )
    parser.add_argument(
        "--run-name",
        type=str,
        default=None,
        help="Optional name for the run directory. Defaults to a timestamp-based name.",
    )
    parser.add_argument(
        "--workbook-path",
        type=Path,
        default=None,
        help="Optional workbook path. Defaults to <output-dir>/<run-name>/track-clonality.xlsx.",
    )
    parser.add_argument(
        "--state-file",
        type=Path,
        default=None,
        help="Optional backfill state file. Defaults to <output-dir>/<run-name>/backfill_state.json.",
    )
    parser.add_argument(
        "--analysis-output-base",
        type=Path,
        default=None,
        help="Optional analysis output base. Defaults to <output-dir>/<run-name>/analysis.",
    )
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=None,
        help="Optional feature-artifact output directory. Defaults to <output-dir>/<run-name>/feature_artifacts.",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="Optional summary JSON path. Defaults to <output-dir>/<run-name>/run_summary.json.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Optional maximum worker threads for the backfill run.",
    )
    parser.add_argument(
        "--folder-workers",
        type=int,
        default=None,
        help="Optional same-month folder concurrency for the backfill run.",
    )
    parser.add_argument(
        "--refresh-each-folder",
        action="store_true",
        help="Refresh the tracking workbook after each folder instead of deferring to month boundaries.",
    )
    parser.add_argument(
        "--include-sl",
        action="store_true",
        help="Include SL rows when exporting feature artifacts.",
    )
    parser.add_argument(
        "--cleanup-staging-root",
        action="store_true",
        help="Delete the temporary symlink staging root after the run completes.",
    )
    return parser


def _read_folder_list(path: Path) -> list[Path]:
    if not path.is_file():
        raise FileNotFoundError(f"Folders file not found: {path}")
    folders: list[Path] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        folders.append(Path(value))
    return folders


def _normalize_folders(args: argparse.Namespace) -> list[Path]:
    folders = list(args.folders or [])
    if args.folders_file is not None:
        folders.extend(_read_folder_list(args.folders_file))
    normalized = [Path(folder).expanduser().resolve() for folder in folders]
    if not normalized:
        raise ValueError("At least one folder must be supplied via --folder or --folders-file.")
    missing = [folder for folder in normalized if not folder.is_dir()]
    if missing:
        raise FileNotFoundError(f"Missing input folders: {', '.join(str(folder) for folder in missing)}")
    basenames = [folder.name for folder in normalized]
    duplicates = sorted({name for name in basenames if basenames.count(name) > 1})
    if duplicates:
        raise ValueError(
            "Duplicate top-level folder names are not supported in the validation harness: "
            + ", ".join(duplicates)
        )
    return normalized


def _build_run_layout(output_dir: Path, run_name: str | None) -> dict[str, Path]:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    safe_run_name = run_name or f"hard_case_validation_{stamp}_{os.getpid()}_{uuid.uuid4().hex[:8]}"
    run_dir = output_dir / safe_run_name
    return {
        "run_dir": run_dir,
        "staging_root": run_dir / "staging_root",
        "analysis_output_base": run_dir / "analysis",
        "workbook_path": run_dir / "track-clonality.xlsx",
        "state_file": run_dir / "backfill_state.json",
        "artifact_dir": run_dir / "feature_artifacts",
        "candidate_artifact_dir": run_dir / "candidate_artifacts",
        "candidate_entries_pickle": run_dir / "candidate_entries.pkl",
        "summary_json": run_dir / "run_summary.json",
    }


def _create_staging_root(folders: list[Path], staging_root: Path) -> Path:
    staging_root.mkdir(parents=True, exist_ok=True)
    for folder in folders:
        staged = staging_root / folder.name
        if staged.exists() or staged.is_symlink():
            staged.unlink()
        staged.symlink_to(folder, target_is_directory=True)
    return staging_root


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _summarize_state(state: dict[str, Any], elapsed_seconds: float, folders: list[Path]) -> dict[str, Any]:
    folder_rows: list[dict[str, Any]] = []
    folder_state = state.get("folders", {}) if isinstance(state, dict) else {}
    for folder in folders:
        item = folder_state.get(folder.name, {}) if isinstance(folder_state, dict) else {}
        started = str(item.get("last_started_at") or "")
        finished = str(item.get("last_finished_at") or "")
        folder_rows.append(
            {
                "folder_name": folder.name,
                "folder_path": str(folder),
                "status": str(item.get("status") or ""),
                "attempts": int(item.get("attempts") or 0),
                "job_count": int(item.get("job_count") or 0),
                "file_count": int(item.get("file_count") or 0),
                "patient_count": int(item.get("patient_count") or 0),
                "qc_file_count": int(item.get("qc_file_count") or 0),
                "last_started_at": started,
                "last_finished_at": finished,
                "last_note": str(item.get("last_note") or ""),
                "error": str(item.get("error") or ""),
            }
        )
    return {
        "elapsed_seconds": round(elapsed_seconds, 3),
        "folder_results": folder_rows,
    }


def run_validation(argv: list[str] | None = None) -> dict[str, Any]:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    folders = _normalize_folders(args)

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    layout = _build_run_layout(output_dir, args.run_name)
    run_dir = layout["run_dir"]
    run_dir.mkdir(parents=True, exist_ok=True)
    staging_root = _create_staging_root(folders, layout["staging_root"])

    workbook_path = Path(args.workbook_path).expanduser().resolve() if args.workbook_path else layout["workbook_path"]
    state_file = Path(args.state_file).expanduser().resolve() if args.state_file else layout["state_file"]
    analysis_output_base = Path(args.analysis_output_base).expanduser().resolve() if args.analysis_output_base else layout["analysis_output_base"]
    artifact_dir = Path(args.artifact_dir).expanduser().resolve() if args.artifact_dir else layout["artifact_dir"]
    summary_json = Path(args.summary_json).expanduser().resolve() if args.summary_json else layout["summary_json"]

    timing: dict[str, float] = {}
    overall_started = time.monotonic()
    timing["setup_started"] = overall_started
    collected_entries: list[dict[str, Any]] = []
    collected_entries_lock = threading.Lock()

    def _record_entries(_folder_name: str, entries: list[dict[str, Any]]) -> None:
        if not entries:
            return
        with collected_entries_lock:
            collected_entries.extend(entries)

    try:
        timing["backfill_started"] = time.monotonic()
        state = run_clonality_backfill(
            input_root=staging_root,
            output_base=analysis_output_base,
            tracking_excel_path=workbook_path,
            state_file=state_file,
            max_workers=args.max_workers,
            folder_workers=args.folder_workers,
            retry_failed=False,
            defer_tracking_refresh=not args.refresh_each_folder,
            collected_entries_callback=_record_entries,
        )
        timing["backfill_finished"] = time.monotonic()

        if not workbook_path.exists():
            raise FileNotFoundError(f"Expected workbook was not created: {workbook_path}")

        timing["artifact_export_started"] = time.monotonic()
        artifact_paths = write_clonality_feature_artifacts(
            workbook_path,
            artifact_dir,
            include_sl=args.include_sl,
        )

        candidate_artifact_paths: dict[str, Path] = {}
        candidate_entries_pickle = layout["candidate_entries_pickle"]
        if collected_entries:
            candidate_entries_pickle.write_bytes(
                pickle.dumps(collected_entries, protocol=pickle.HIGHEST_PROTOCOL)
            )
            candidate_artifact_paths = write_clonality_candidate_artifacts(
                layout["candidate_artifact_dir"],
                collected_entries,
                include_sl=args.include_sl,
                write_gold_label_template=True,
            )
        timing["artifact_export_finished"] = time.monotonic()

        elapsed = time.monotonic() - overall_started
        state_payload = state
        if state_file.exists():
            try:
                state_payload = json.loads(state_file.read_text(encoding="utf-8"))
            except Exception:
                pass

        summary = {
            "generated_at_utc": _timestamp(),
            "folders": [str(folder) for folder in folders],
            "output_dir": str(output_dir),
            "run_dir": str(run_dir),
            "staging_root": str(staging_root),
            "workbook_path": str(workbook_path),
            "state_file": str(state_file),
            "analysis_output_base": str(analysis_output_base),
            "artifact_dir": str(artifact_dir),
            "candidate_artifact_dir": str(layout["candidate_artifact_dir"]),
            "candidate_entries_pickle": str(candidate_entries_pickle) if collected_entries else "",
            "summary_json": str(summary_json),
            "artifact_outputs": {name: str(path) for name, path in artifact_paths.items()},
            "candidate_artifact_outputs": {name: str(path) for name, path in candidate_artifact_paths.items()},
            "backfill_options": {
                "max_workers": args.max_workers,
                "folder_workers": args.folder_workers,
                "defer_tracking_refresh": not args.refresh_each_folder,
                "include_sl": bool(args.include_sl),
            },
            "timing_seconds": {
                "total": round(elapsed, 3),
                "backfill": round(timing["backfill_finished"] - timing["backfill_started"], 3),
                "artifact_export": round(timing["artifact_export_finished"] - timing["artifact_export_started"], 3),
            },
            "candidate_entry_count": int(len(collected_entries)),
            "state_summary": _summarize_state(state_payload if isinstance(state_payload, dict) else {}, elapsed, folders),
        }
        summary_json.parent.mkdir(parents=True, exist_ok=True)
        summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        print(json.dumps(summary, indent=2, sort_keys=True))
        return summary
    finally:
        if args.cleanup_staging_root:
            import shutil

            shutil.rmtree(staging_root, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    run_validation(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
