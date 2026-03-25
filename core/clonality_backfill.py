from __future__ import annotations

import argparse
import copy
import json
import os
import re
import threading
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from config import APP_SETTINGS
from core.batch import generate_jobs, run_batch_jobs
from core.log import log


DEFAULT_INPUT_ROOT = Path("/Users/christian/Desktop/DATA/Klonalitet/2025_data")
DEFAULT_OUTPUT_BASE = Path("/Users/christian/Desktop/FINAL")
DEFAULT_TRACKING_EXCEL = Path("/Users/christian/Desktop/Excel_Fraggler/track-clonality.xlsx")
DEFAULT_STATE_FILE = Path("/Users/christian/Desktop/Excel_Fraggler/backfill_2025_state.json")
BACKFILL_OUTPUT_DIRNAME = "backfill_2025"
BACKFILL_REPORT_DIRNAME = "reports_backfill"
STATE_VERSION = 1
RUNNING_HEARTBEAT_TTL_SECONDS = 300
TERMINAL_PHASES = {"done", "failed"}
PHASE_RANK = {
    "": 0,
    "pending": 0,
    "folder_start": 1,
    "job_start": 2,
    "stage_files": 3,
    "collect_entries": 4,
    "analyze": 5,
    "build_report": 6,
    "write_tracking_excel": 7,
    "done": 8,
    "failed": 8,
}


def discover_top_level_run_folders(input_root: Path, month: str | None = None) -> list[Path]:
    if not input_root.is_dir():
        raise ValueError(f"Input root not found: {input_root}")
    month_prefix = str(month or "").strip()
    folders = [p for p in input_root.iterdir() if p.is_dir()]
    if month_prefix:
        folders = [p for p in folders if p.name.startswith(month_prefix)]
    return sorted(folders, key=lambda p: p.name)


def run_clonality_backfill(
    input_root: Path = DEFAULT_INPUT_ROOT,
    month: str | None = None,
    output_base: Path = DEFAULT_OUTPUT_BASE,
    tracking_excel_path: Path = DEFAULT_TRACKING_EXCEL,
    state_file: Path = DEFAULT_STATE_FILE,
    max_workers: int = 1,
    retry_failed: bool = False,
) -> dict[str, Any]:
    input_root = Path(input_root).expanduser()
    output_base = Path(output_base).expanduser()
    tracking_excel_path = Path(tracking_excel_path).expanduser()
    state_file = Path(state_file).expanduser()
    tracking_excel_path.parent.mkdir(parents=True, exist_ok=True)
    state_file.parent.mkdir(parents=True, exist_ok=True)

    settings_backup = copy.deepcopy(APP_SETTINGS)
    APP_SETTINGS["active_analysis"] = "clonality"
    analysis_batch = APP_SETTINGS.setdefault("analyses", {}).setdefault("clonality", {}).setdefault("batch", {})
    patient_regex = analysis_batch.get("patient_id_regex", r"\d{2}OUM\d{5}")
    pipeline_settings = APP_SETTINGS.setdefault("analyses", {}).setdefault("clonality", {}).setdefault("pipeline", {})
    previous_disable_mp = os.environ.get("FRAGGLER_DISABLE_MULTIPROCESSING")
    os.environ["FRAGGLER_DISABLE_MULTIPROCESSING"] = "1"

    try:
        folders = discover_top_level_run_folders(input_root, month=month)
        state = _load_state(state_file, input_root, output_base, tracking_excel_path, patient_regex)
        state_lock = threading.Lock()
        with state_lock:
            _sync_state_folders(state, folders, patient_regex)
            _reset_stale_running_items(state)
            _save_state(state_file, state)

        if not folders:
            log(f"[BACKFILL] No run folders found for month={month or 'all'} under {input_root}")
            return state

        total_folders = len(folders)
        done_this_run = 0
        failed_this_run = 0
        touched_months: set[str] = set()

        log(
            f"[BACKFILL] Starting clonality backfill for {month or 'all'} with "
            f"{total_folders} top-level folders, tracking workbook {tracking_excel_path}, "
            f"max_workers={max(1, int(max_workers))}."
        )

        for index, folder in enumerate(folders, start=1):
            folder_name = folder.name
            item = state["folders"][folder_name]
            month_key = _month_key(folder_name)
            touched_months.add(month_key)

            if item["status"] == "done":
                log(f"[BACKFILL] Skipping completed folder {folder_name} ({index}/{total_folders}).")
                continue
            if item["status"] == "failed" and not retry_failed:
                log(f"[BACKFILL] Skipping failed folder {folder_name} ({index}/{total_folders}); use retry_failed to rerun.")
                continue

            log(
                f"[BACKFILL] Folder {index}/{total_folders} start: {folder_name} | "
                f"files={item['file_count']} patients={item['patient_count']} qc_files={item['qc_file_count']}"
            )
            with state_lock:
                _mark_folder_running(item, output_base, month_key, folder_name, tracking_excel_path)
                _save_state(state_file, state)

            try:
                jobs = generate_jobs([folder], aggregate_patients=True, patient_regex=patient_regex)
                if not jobs:
                    raise RuntimeError("No jobs generated for folder.")

                folder_output_base = Path(item["output_base"])
                def _progress_callback(event: dict[str, Any]) -> None:
                    payload = dict(event)
                    payload["folder_name"] = folder_name
                    with state_lock:
                        _record_folder_progress(state_file, state, folder_name, payload)

                result = run_batch_jobs(
                    jobs=jobs,
                    output_base=folder_output_base,
                    out_folder_tmpl="ASSAY_REPORTS",
                    outfile_html_tmpl="QC_REPORT_{name}.html",
                    excel_name_tmpl="Fraggler_QC_Trends.xlsx",
                    pipeline_scope=pipeline_settings.get("mode", "all"),
                    assay_filter=pipeline_settings.get("assay_filter_substring", ""),
                    aggregate_dit_reports=True,
                    continue_on_error=True,
                    update_callback=None,
                    progress_callback=_progress_callback,
                    max_workers=max_workers,
                    tracking_excel_path=tracking_excel_path,
                    aggregate_outdir_name=BACKFILL_REPORT_DIRNAME,
                )
                failed_jobs = list((result or {}).get("failed_jobs", []))
                completed_jobs = list((result or {}).get("completed_jobs", []))

                with state_lock:
                    item["completed_jobs"] = completed_jobs
                    item["failed_jobs"] = failed_jobs
                    item["job_count"] = len(jobs)
                    item["jobs_total"] = len(jobs)
                    item["updated_at"] = _timestamp()
                    item["last_finished_at"] = item["updated_at"]
                    item["owner_pid"] = 0
                    item["job_progress"] = {}

                    if failed_jobs:
                        item["status"] = "failed"
                        item["error"] = f"Failed jobs: {', '.join(failed_jobs)}"
                        item["current_phase"] = "failed"
                        failed_this_run += 1
                        log(
                            f"[BACKFILL] Folder failed: {folder_name} | failed_jobs={len(failed_jobs)} "
                            f"| cumulative_done={done_this_run} cumulative_failed={failed_this_run}"
                        )
                    else:
                        item["status"] = "done"
                        item["error"] = ""
                        item["current_phase"] = "done"
                        done_this_run += 1
                        log(
                            f"[BACKFILL] Folder complete: {folder_name} | reports={len(completed_jobs)} "
                            f"| cumulative_done={done_this_run} cumulative_failed={failed_this_run}"
                        )
            except Exception as exc:
                with state_lock:
                    item["status"] = "failed"
                    item["error"] = str(exc)
                    item["failed_jobs"] = [str(exc)]
                    item["updated_at"] = _timestamp()
                    item["last_finished_at"] = item["updated_at"]
                    item["current_phase"] = "failed"
                    item["last_note"] = str(exc)
                    item["owner_pid"] = 0
                    item["job_progress"] = {}
                    failed_this_run += 1
                log(
                    f"[BACKFILL] Folder failed with exception: {folder_name} | {exc} "
                    f"| cumulative_done={done_this_run} cumulative_failed={failed_this_run}"
                )
            finally:
                with state_lock:
                    _save_state(state_file, state)

        for month_key in sorted(touched_months):
            _write_month_summary(state, tracking_excel_path, state_file.parent, month_key)

        log(
            f"[BACKFILL] Finished {month or 'all'} | done={done_this_run} failed={failed_this_run} "
            f"| state={state_file} workbook={tracking_excel_path}"
        )
        return state
    finally:
        if previous_disable_mp is None:
            os.environ.pop("FRAGGLER_DISABLE_MULTIPROCESSING", None)
        else:
            os.environ["FRAGGLER_DISABLE_MULTIPROCESSING"] = previous_disable_mp
        APP_SETTINGS.clear()
        APP_SETTINGS.update(settings_backup)


def _load_state(
    state_file: Path,
    input_root: Path,
    output_base: Path,
    tracking_excel_path: Path,
    patient_regex: str,
) -> dict[str, Any]:
    if state_file.exists():
        try:
            state = json.loads(state_file.read_text(encoding="utf-8"))
        except Exception:
            state = {}
    else:
        state = {}

    state.setdefault("version", STATE_VERSION)
    state.setdefault("created_at", _timestamp())
    state["updated_at"] = _timestamp()
    state["input_root"] = str(input_root)
    state["output_base"] = str(output_base)
    state["tracking_excel_path"] = str(tracking_excel_path)
    state["patient_regex"] = patient_regex
    state.setdefault("folders", {})
    return state


def _save_state(state_file: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = _timestamp()
    tmp_file = state_file.with_name(f"{state_file.name}.tmp")
    tmp_file.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
    tmp_file.replace(state_file)


def _sync_state_folders(state: dict[str, Any], folders: list[Path], patient_regex: str) -> None:
    for folder in folders:
        state["folders"].setdefault(folder.name, _build_folder_item(folder, patient_regex))


def _reset_stale_running_items(state: dict[str, Any]) -> None:
    for name, item in state.get("folders", {}).items():
        if item.get("status") == "running":
            owner_pid = int(item.get("owner_pid") or 0)
            last_heartbeat = _parse_timestamp(item.get("last_heartbeat_at"))
            heartbeat_stale = (
                last_heartbeat is None
                or (datetime.now() - last_heartbeat).total_seconds() > RUNNING_HEARTBEAT_TTL_SECONDS
            )
            pid_alive = _pid_is_alive(owner_pid)
            if heartbeat_stale or not pid_alive:
                item["status"] = "pending"
                item["error"] = "Reset stale running state after restart."
                item["updated_at"] = _timestamp()
                _clear_live_progress(item, phase="pending")
                log(f"[BACKFILL] Reset stale running folder to pending: {name}")


def _build_folder_item(folder: Path, patient_regex: str) -> dict[str, Any]:
    file_count = 0
    qc_file_count = 0
    patient_ids: set[str] = set()
    patient_re = re.compile(patient_regex) if patient_regex else None
    qc_re = re.compile(r"^(PK|NK|RK)(\d+)?[_-]", re.IGNORECASE)

    for file_path in sorted(folder.rglob("*.fsa")):
        file_count += 1
        name = file_path.name
        if qc_re.match(name):
            qc_file_count += 1
            continue
        if patient_re:
            match = patient_re.search(name)
            if match:
                patient_ids.add(match.group())

    return {
        "folder_name": folder.name,
        "folder_path": str(folder),
        "month": _month_key(folder.name),
        "status": "pending",
        "attempts": 0,
        "file_count": file_count,
        "patient_count": len(patient_ids),
        "qc_file_count": qc_file_count,
        "job_count": 0,
        "completed_jobs": [],
        "failed_jobs": [],
        "report_dir": "",
        "tracking_excel_path": "",
        "output_base": "",
        "current_job": "",
        "current_phase": "",
        "current_file": "",
        "files_done": 0,
        "files_total": 0,
        "jobs_done": 0,
        "jobs_total": 0,
        "last_started_at": "",
        "last_finished_at": "",
        "last_heartbeat_at": "",
        "last_note": "",
        "owner_pid": 0,
        "job_progress": {},
        "updated_at": _timestamp(),
        "error": "",
    }


def _mark_folder_running(item: dict[str, Any], output_base: Path, month_key: str, folder_name: str, tracking_excel_path: Path) -> None:
    report_base = output_base / BACKFILL_OUTPUT_DIRNAME / month_key / folder_name
    report_dir = report_base / BACKFILL_REPORT_DIRNAME
    item["status"] = "running"
    item["attempts"] = int(item.get("attempts", 0)) + 1
    item["last_started_at"] = _timestamp()
    item["updated_at"] = item["last_started_at"]
    item["output_base"] = str(report_base)
    item["report_dir"] = str(report_dir)
    item["tracking_excel_path"] = str(tracking_excel_path)
    item["current_job"] = ""
    item["current_phase"] = "folder_start"
    item["current_file"] = ""
    item["files_done"] = 0
    item["files_total"] = 0
    item["jobs_done"] = 0
    item["jobs_total"] = 0
    item["last_heartbeat_at"] = item["last_started_at"]
    item["last_note"] = "folder_started"
    item["owner_pid"] = os.getpid()
    item["job_progress"] = {}
    item["error"] = ""


def _record_folder_progress(
    state_file: Path,
    state: dict[str, Any],
    folder_name: str,
    event: dict[str, Any],
) -> None:
    item = state["folders"][folder_name]
    phase = str(event.get("phase", "") or "")
    if item.get("status") in TERMINAL_PHASES and phase not in TERMINAL_PHASES:
        return

    heartbeat_at = str(event.get("heartbeat_at") or _timestamp())
    job_name = str(event.get("job_name", "") or "")
    file_name = str(event.get("file_name", "") or "")
    note = str(event.get("note", "") or "")
    files_done = _coerce_int(event.get("files_done"))
    files_total = _coerce_int(event.get("files_total"))

    job_progress = item.setdefault("job_progress", {})
    slot_key = job_name or "__folder__"
    slot = job_progress.get(slot_key, {})
    if not _should_accept_progress_event(slot, heartbeat_at, phase, files_done):
        return

    slot.update(
        {
            "job_name": job_name,
            "phase": phase,
            "file_name": file_name,
            "files_done": files_done,
            "files_total": files_total,
            "heartbeat_at": heartbeat_at,
            "note": note,
        }
    )
    job_progress[slot_key] = slot

    item["jobs_done"] = max(
        int(item.get("jobs_done", 0) or 0),
        _coerce_int(event.get("jobs_done")) or 0,
    )
    item["jobs_total"] = max(
        int(item.get("jobs_total", 0) or 0),
        _coerce_int(event.get("jobs_total")) or 0,
    )
    aggregated_files_done = 0
    aggregated_files_total = 0
    for progress in job_progress.values():
        slot_done = _coerce_int(progress.get("files_done")) or 0
        slot_total = _coerce_int(progress.get("files_total")) or 0
        aggregated_files_done += min(slot_done, slot_total) if slot_total > 0 else slot_done
        aggregated_files_total += slot_total
    item["files_done"] = aggregated_files_done
    item["files_total"] = aggregated_files_total

    latest_progress = max(
        job_progress.values(),
        key=lambda progress: (
            _parse_timestamp(progress.get("heartbeat_at")) or datetime.min,
            _phase_rank(progress.get("phase")),
            progress.get("job_name", ""),
        ),
    )
    item["current_job"] = str(latest_progress.get("job_name", "") or "")
    item["current_phase"] = str(latest_progress.get("phase", "") or "")
    item["current_file"] = str(latest_progress.get("file_name", "") or "")
    item["last_heartbeat_at"] = str(latest_progress.get("heartbeat_at") or heartbeat_at)
    item["last_note"] = str(latest_progress.get("note", "") or note)
    item["updated_at"] = item["last_heartbeat_at"]
    _save_state(state_file, state)


def _clear_live_progress(item: dict[str, Any], *, phase: str = "") -> None:
    item["current_job"] = ""
    item["current_phase"] = phase
    item["current_file"] = ""
    item["files_done"] = 0
    item["files_total"] = 0
    item["jobs_done"] = 0
    item["jobs_total"] = 0
    item["last_note"] = ""
    item["last_heartbeat_at"] = ""
    item["owner_pid"] = 0
    item["job_progress"] = {}


def _coerce_int(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _phase_rank(phase: Any) -> int:
    return PHASE_RANK.get(str(phase or ""), 0)


def _should_accept_progress_event(
    slot: dict[str, Any],
    heartbeat_at: str,
    phase: str,
    files_done: int | None,
) -> bool:
    previous_phase = str(slot.get("phase", "") or "")
    if previous_phase in TERMINAL_PHASES and phase not in TERMINAL_PHASES:
        return False

    previous_heartbeat = _parse_timestamp(slot.get("heartbeat_at"))
    current_heartbeat = _parse_timestamp(heartbeat_at)
    if previous_heartbeat and current_heartbeat and current_heartbeat < previous_heartbeat:
        return False

    previous_rank = _phase_rank(previous_phase)
    current_rank = _phase_rank(phase)
    if previous_phase and current_rank < previous_rank:
        return False

    previous_done = _coerce_int(slot.get("files_done"))
    if current_rank == previous_rank and previous_done is not None and files_done is not None and files_done < previous_done:
        return False

    return True


def _write_month_summary(
    state: dict[str, Any],
    tracking_excel_path: Path,
    out_dir: Path,
    month_key: str,
) -> Path:
    month_prefix = month_key.replace("_", "-")
    month_items = [item for item in state.get("folders", {}).values() if item.get("month") == month_key]
    done_items = [item for item in month_items if item.get("status") == "done"]
    failed_items = [item for item in month_items if item.get("status") == "failed"]
    ladder_review_count = 0
    pk_exception_count = 0

    if tracking_excel_path.exists():
        try:
            patient = pd.read_excel(tracking_excel_path, sheet_name="Patient_Runs", engine="openpyxl")
            control = pd.read_excel(tracking_excel_path, sheet_name="Control_Runs", engine="openpyxl")
            peaks = pd.read_excel(tracking_excel_path, sheet_name="PK_Peaks", engine="openpyxl")
            patient["RunDate"] = patient["RunDate"].astype(str)
            control["RunDate"] = control["RunDate"].astype(str)
            peaks["RunDate"] = peaks["RunDate"].astype(str)
            ladder_review_count = int(
                (
                    patient.loc[
                        patient["RunDate"].str.startswith(month_prefix)
                        & patient["LadderQC"].astype(str).str.strip().ne("")
                        & patient["LadderQC"].astype(str).str.strip().str.lower().ne("ok")
                    ].shape[0]
                    + control.loc[
                        control["RunDate"].str.startswith(month_prefix)
                        & control["LadderQC"].astype(str).str.strip().ne("")
                        & control["LadderQC"].astype(str).str.strip().str.lower().ne("ok")
                    ].shape[0]
                )
            )
            peaks["AbsDeltaBP"] = pd.to_numeric(peaks.get("AbsDeltaBP"), errors="coerce")
            peaks["OK"] = peaks.get("OK").astype(str).str.lower()
            peaks["Kind"] = peaks.get("Kind").astype(str).str.lower()
            pk_exception_count = int(
                peaks.loc[
                    peaks["RunDate"].str.startswith(month_prefix)
                    & (peaks["Kind"] == "sample")
                    & ((peaks["OK"] != "true") | (peaks["AbsDeltaBP"] > 2.0))
                ].shape[0]
            )
        except Exception as exc:
            log(f"[BACKFILL] Could not compute workbook-backed month summary for {month_key}: {exc}")

    summary_lines = [
        f"Backfill month summary: {month_key}",
        f"Generated: {_timestamp()}",
        f"Tracking workbook: {tracking_excel_path}",
        f"Processed folders: {len(done_items)}",
        f"Failed folders: {len(failed_items)}",
        f"Ladder review count: {ladder_review_count}",
        f"PK marker exceptions: {pk_exception_count}",
        "",
        "Failed folders detail:",
    ]
    if failed_items:
        for item in failed_items:
            summary_lines.append(f"- {item['folder_name']}: {item.get('error') or 'failed'}")
    else:
        summary_lines.append("- none")

    out_path = out_dir / f"backfill_2025_summary_{month_key}.txt"
    out_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    log(f"[BACKFILL] Wrote month summary for {month_key} to {out_path}")
    return out_path


def _month_key(folder_name: str) -> str:
    return folder_name[:7] if len(folder_name) >= 7 else "unknown"


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _parse_timestamp(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.strptime(str(value), "%Y-%m-%d %H:%M:%S")
    except ValueError:
        return None


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run resumable clonality backfill over historical data.")
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--month", default="all", help="Month prefix like 2025_01, or 'all'.")
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--tracking-excel-path", type=Path, default=DEFAULT_TRACKING_EXCEL)
    parser.add_argument("--state-file", type=Path, default=DEFAULT_STATE_FILE)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--retry-failed", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    month = None if str(args.month).strip().lower() == "all" else str(args.month).strip()
    run_clonality_backfill(
        input_root=args.input_root,
        month=month,
        output_base=args.output_base,
        tracking_excel_path=args.tracking_excel_path,
        state_file=args.state_file,
        max_workers=args.max_workers,
        retry_failed=args.retry_failed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
