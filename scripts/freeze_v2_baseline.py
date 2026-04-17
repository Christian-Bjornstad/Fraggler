"""Freeze Python reference outputs and timings for the Fraggler v2 rewrite."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config import APP_SETTINGS
from core.analyses.clonality.classification import classify_fsa
from core.analysis import (
    LADDER_FIT_PROFILE_CLONALITY_LIZ500,
    LADDER_FIT_PROFILE_CLONALITY_ROX400HD,
    analyse_fsa_liz,
    analyse_fsa_rox,
)
from core.html_reports import build_dit_html_reports
from core.pipeline import run_pipeline
from core.qc.qc_rules import QCRules
from core.runner import run_qc_job


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Freeze timings and artifact manifests from the Python Fraggler implementation."
    )
    parser.add_argument(
        "--scenario-file",
        type=Path,
        required=True,
        help="JSON file describing baseline scenarios.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where the baseline manifest and copied artifacts are written.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data = json.loads(args.scenario_file.read_text(encoding="utf-8"))
    scenarios = data.get("scenarios", [])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []

    for scenario in scenarios:
        name = str(scenario.get("name") or "unnamed")
        scenario_dir = args.output_dir / name
        scenario_dir.mkdir(parents=True, exist_ok=True)
        result = run_scenario(scenario, scenario_dir)
        results.append(result)
        (scenario_dir / "summary.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    manifest = {
        "generated_at_epoch_s": time.time(),
        "scenario_file": str(args.scenario_file),
        "scenario_count": len(results),
        "results": results,
    }
    (args.output_dir / "baseline_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


def run_scenario(scenario: dict[str, Any], scenario_dir: Path) -> dict[str, Any]:
    kind = str(scenario.get("kind") or "").strip()
    started = time.perf_counter()
    try:
        if kind == "clonality_file_analysis":
            details = run_clonality_file_analysis(scenario)
        elif kind == "combined_qc_dit":
            details = run_combined_qc_dit(scenario, scenario_dir)
        elif kind == "command":
            details = run_command_scenario(scenario, scenario_dir)
        else:
            raise ValueError(f"Unsupported scenario kind: {kind}")
        status = "ok"
    except Exception as exc:  # pragma: no cover - defensive for ad-hoc benchmark use
        details = {"error": f"{type(exc).__name__}: {exc}"}
        status = "error"
    elapsed = time.perf_counter() - started
    return {
        "name": str(scenario.get("name") or "unnamed"),
        "kind": kind,
        "status": status,
        "elapsed_seconds": round(elapsed, 4),
        "details": details,
    }


def run_clonality_file_analysis(scenario: dict[str, Any]) -> dict[str, Any]:
    APP_SETTINGS["active_analysis"] = "clonality"
    input_file = Path(str(scenario["input_file"]))
    classified = classify_fsa(input_file)
    if classified is None:
        raise RuntimeError(f"Could not classify {input_file}")

    assay, group, ladder, trace_channels, peak_channels, primary_peak_channel, bp_min, bp_max = classified
    sample_channel = trace_channels[0]
    started = time.perf_counter()
    if ladder == "LIZ":
        fsa = analyse_fsa_liz(
            input_file,
            sample_channel,
            ladder_name="LIZ500_250",
            ladder_fit_profile=LADDER_FIT_PROFILE_CLONALITY_LIZ500,
        )
    else:
        fsa = analyse_fsa_rox(
            input_file,
            sample_channel,
            ladder_name="ROX400HD",
            ladder_fit_profile=LADDER_FIT_PROFILE_CLONALITY_ROX400HD,
        )
    analysis_elapsed = time.perf_counter() - started
    return {
        "assay": assay,
        "group": group,
        "ladder": ladder,
        "trace_channels": trace_channels,
        "peak_channels": peak_channels,
        "primary_peak_channel": primary_peak_channel,
        "bp_min": bp_min,
        "bp_max": bp_max,
        "analysis_seconds": round(analysis_elapsed, 4),
        "fsa_file_name": getattr(fsa, "file_name", None),
        "ladder_fit_strategy": getattr(fsa, "ladder_fit_strategy", None),
        "ladder_review_required": getattr(fsa, "ladder_review_required", None),
        "ladder_missing_expected_steps": list(
            map(float, getattr(fsa, "ladder_missing_expected_steps", []))
        ),
    }


def run_combined_qc_dit(scenario: dict[str, Any], scenario_dir: Path) -> dict[str, Any]:
    APP_SETTINGS["active_analysis"] = "clonality"
    source_dir = Path(str(scenario["source_dir"]))
    patient_prefixes = tuple(str(x) for x in scenario.get("patient_prefixes", []))
    control_prefixes = tuple(str(x) for x in scenario.get("control_prefixes", []))
    selected = sorted(
        p for p in source_dir.glob("*.fsa")
        if p.name.startswith(patient_prefixes) or p.name.startswith(control_prefixes)
    )
    patient_files = [p for p in selected if p.name.startswith(patient_prefixes)]
    control_files = [p for p in selected if p not in patient_files]
    selected_names = {p.name for p in selected}

    pipeline_start = time.perf_counter()
    all_entries = run_pipeline(
        fsa_dir=source_dir,
        base_outdir=scenario_dir,
        assay_folder_name=str(scenario.get("output_folder_name") or scenario_dir.name),
        return_entries=True,
        make_dit_reports=False,
        mode="all",
        update_tracking_workbook=False,
    ) or []
    pipeline_seconds = time.perf_counter() - pipeline_start

    all_entries = [
        e for e in all_entries
        if getattr(e.get("fsa"), "file_name", "") in selected_names
    ]
    patient_only_entries = [
        e for e in all_entries
        if not getattr(e.get("fsa"), "file_name", "").startswith(control_prefixes)
    ]

    qc_start = time.perf_counter()
    qc_report_path, qc_entries = run_qc_job(
        fsa_dir=None,
        base_outdir=scenario_dir,
        out_html_name="QC_REPORT_v2_baseline.html",
        excel_name="QC_REPORT_v2_baseline.xlsx",
        rules=QCRules(),
        files=control_files,
        update_tracking_workbook=False,
        return_entries=True,
    )
    qc_seconds = time.perf_counter() - qc_start

    report_dir = scenario_dir / "REPORTS"
    report_dir.mkdir(parents=True, exist_ok=True)
    dit_start = time.perf_counter()
    build_dit_html_reports(patient_only_entries + (qc_entries or []), report_dir)
    dit_seconds = time.perf_counter() - dit_start

    report_paths = sorted(report_dir.glob("*.html"))
    return {
        "selected_files": len(selected),
        "patient_files": len(patient_files),
        "control_files": len(control_files),
        "all_entries": len(all_entries),
        "patient_only_entries": len(patient_only_entries),
        "qc_entries": len(qc_entries or []),
        "pipeline_seconds": round(pipeline_seconds, 4),
        "qc_seconds": round(qc_seconds, 4),
        "dit_build_seconds": round(dit_seconds, 4),
        "qc_report_path": str(qc_report_path) if qc_report_path else None,
        "report_paths": [str(p) for p in report_paths],
        "report_sizes": {p.name: p.stat().st_size for p in report_paths},
    }


def run_command_scenario(scenario: dict[str, Any], scenario_dir: Path) -> dict[str, Any]:
    argv = [str(arg) for arg in scenario.get("argv", [])]
    if not argv:
        raise ValueError("command scenario requires argv")

    completed = subprocess.run(
        argv,
        cwd=str(Path.cwd()),
        capture_output=True,
        text=True,
        check=False,
    )
    stdout_path = scenario_dir / "stdout.log"
    stderr_path = scenario_dir / "stderr.log"
    stdout_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    stderr_path.write_text(completed.stderr, encoding="utf-8", errors="replace")

    copied_artifacts: list[str] = []
    for artifact in scenario.get("copy_artifacts", []):
        src = Path(str(artifact))
        if src.exists():
            dst = scenario_dir / src.name
            if src.is_dir():
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dst)
            copied_artifacts.append(str(dst))

    return {
        "argv": argv,
        "returncode": completed.returncode,
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "copied_artifacts": copied_artifacts,
    }


if __name__ == "__main__":
    raise SystemExit(main())
