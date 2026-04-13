from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


from scripts.validate_flt3_ladder_fits import DEFAULT_DATA_ROOT, run_validation  # noqa: E402


DEFAULT_YEARS = ("2024", "2025", "2026")
PRIMARY_RUNS = (
    ("2024_dit", {"years": ["2024"], "dit_only": True}),
    ("2025_dit", {"years": ["2025"], "dit_only": True}),
    ("2026_dit", {"years": ["2026"], "dit_only": True}),
    ("all_years_dit", {"years": list(DEFAULT_YEARS), "dit_only": True}),
)
SECONDARY_RUNS = (
    ("all_years_all_files", {"years": list(DEFAULT_YEARS), "dit_only": False}),
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run standardized FLT3 ladder-fit verification across multiple years and "
            "write a consolidated summary/report bundle."
        )
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help=f"Root containing FLT3 year folders. Default: {DEFAULT_DATA_ROOT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where yearly verification outputs and consolidated summaries are written.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Worker count passed to validate_flt3_ladder_fits.py.",
    )
    parser.add_argument(
        "--years",
        nargs="+",
        default=list(DEFAULT_YEARS),
        help="Years to include. Defaults to 2024 2025 2026.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional limit passed through to each run. Use 0 for full runs.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=45,
        help="Per-file timeout for FLT3 validation.",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=250,
        help="Checkpoint cadence for each validation run.",
    )
    parser.add_argument(
        "--baseline-summary",
        type=Path,
        default=None,
        help="Optional historical FLT3 summary.json to compare against.",
    )
    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Also run the broader all-files secondary verification scope.",
    )
    parser.add_argument(
        "--require-run-name-contains",
        default="",
        help="Optional substring required in ABI run_name.",
    )
    parser.add_argument(
        "--include-npm1",
        action="store_true",
        help="Also include NPM1 in the secondary scope if desired.",
    )
    parser.add_argument(
        "--run-name",
        default="flt3_multiyear_verification",
        help="Folder name under output-dir for this verification pass.",
    )
    return parser


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _ordered_fieldnames(rows: list[dict[str, Any]]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key in seen:
                continue
            ordered.append(str(key))
            seen.add(str(key))
    return ordered


def _copy_summary_bundle(src_dir: Path, dst_dir: Path) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "summary.json",
        "summary.md",
        "flt3_ladder_metrics.csv",
        "flt3_ladder_review_manifest.csv",
        "flt3_ladder_fit_failed_manifest.csv",
        "flt3_ladder_actionable_review_manifest.csv",
        "flt3_ladder_manual_adjustments_applied.csv",
    ):
        src = src_dir / name
        if src.exists():
            shutil.copy2(src, dst_dir / name)


def _top_counter(rows: list[dict[str, Any]], field: str, *, limit: int = 10) -> list[list[Any]]:
    counts = Counter(str(row.get(field) or "") for row in rows if str(row.get(field) or "").strip())
    return [[key, value] for key, value in counts.most_common(limit)]


def _load_metrics_rows(run_dir: Path) -> list[dict[str, Any]]:
    metrics_csv = run_dir / "flt3_ladder_metrics.csv"
    if not metrics_csv.exists():
        return []
    with metrics_csv.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _load_review_rows(run_dir: Path) -> list[dict[str, Any]]:
    review_csv = run_dir / "flt3_ladder_review_manifest.csv"
    if not review_csv.exists():
        return []
    with review_csv.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _subset_breakdown(metrics_rows: list[dict[str, Any]]) -> dict[str, Any]:
    processed = [row for row in metrics_rows if not str(row.get("status") or "").startswith("skipped")]
    status_by_analysis_type = Counter(
        f"{str(row.get('analysis_type') or '')}|{str(row.get('status') or '')}" for row in processed
    )
    status_by_assay = Counter(
        f"{str(row.get('assay') or '')}|{str(row.get('status') or '')}" for row in processed
    )

    return {
        "status_by_analysis_type": dict(sorted(status_by_analysis_type.items())),
        "status_by_assay": dict(sorted(status_by_assay.items())),
        "top_review_reasons": _top_counter(processed, "review_reason"),
        "top_strategies": _top_counter(processed, "ladder_fit_strategy"),
    }


def _summary_row(label: str, summary: dict[str, Any]) -> dict[str, Any]:
    status_counts = summary.get("status_counts", {}) if isinstance(summary.get("status_counts"), dict) else {}
    return {
        "run": label,
        "validated_count": int(summary.get("validated_count", 0) or 0),
        "ok": int(status_counts.get("ok", 0) or 0),
        "review_required": int(status_counts.get("review_required", 0) or 0),
        "ladder_fit_failed": int(status_counts.get("ladder_fit_failed", 0) or 0),
        "manual_adjustment": int(status_counts.get("manual_adjustment", 0) or 0),
        "manual_review_count": int(summary.get("manual_review_count", 0) or 0),
        "median_r2": summary.get("median_r2"),
        "median_max_residual_bp": summary.get("median_max_residual_bp"),
        "max_residual_bp": summary.get("max_residual_bp"),
    }


def _diff_vs_baseline(row: dict[str, Any], baseline: dict[str, Any] | None) -> dict[str, Any]:
    if not baseline:
        return {}
    baseline_row = _summary_row("baseline", baseline)
    fields = (
        "validated_count",
        "ok",
        "review_required",
        "ladder_fit_failed",
        "manual_adjustment",
        "manual_review_count",
    )
    diff = {f"{field}_delta": int(row[field]) - int(baseline_row[field]) for field in fields}
    for field in ("median_r2", "median_max_residual_bp", "max_residual_bp"):
        current = row.get(field)
        base = baseline_row.get(field)
        if current is None or base is None:
            diff[f"{field}_delta"] = None
        else:
            diff[f"{field}_delta"] = float(current) - float(base)
    return diff


def _markdown_report(
    rows: list[dict[str, Any]],
    baseline_name: str,
    per_run_details: dict[str, dict[str, Any]],
) -> str:
    header = (
        "| Run | Validated | OK | Review Required | Fit Failed | Manual Adjustment | "
        "Manual Review Count | Median R2 | Median Max Residual | Worst Max Residual |"
    )
    separator = "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
    lines = [
        "# FLT3 Multiyear Verification",
        "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"Baseline comparison: `{baseline_name}`",
        "",
        header,
        separator,
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["run"]),
                    str(row["validated_count"]),
                    str(row["ok"]),
                    str(row["review_required"]),
                    str(row["ladder_fit_failed"]),
                    str(row["manual_adjustment"]),
                    str(row["manual_review_count"]),
                    str(row["median_r2"]),
                    str(row["median_max_residual_bp"]),
                    str(row["max_residual_bp"]),
                ]
            )
            + " |"
        )

    for label, detail in per_run_details.items():
        breakdown = detail.get("breakdown", {})
        lines.extend(
            [
                "",
                f"## {label}",
                "",
                "Top review reasons:",
            ]
        )
        top_reasons = breakdown.get("top_review_reasons", [])
        if top_reasons:
            lines.extend([f"- {reason}: {count}" for reason, count in top_reasons])
        else:
            lines.append("- none")
        lines.append("")
        lines.append("Top strategies:")
        top_strategies = breakdown.get("top_strategies", [])
        if top_strategies:
            lines.extend([f"- {strategy}: {count}" for strategy, count in top_strategies])
        else:
            lines.append("- none")
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)

    output_root = args.output_dir.expanduser().resolve() / str(args.run_name)
    output_root.mkdir(parents=True, exist_ok=True)
    baseline_summary = (
        _load_json(args.baseline_summary.expanduser().resolve())
        if args.baseline_summary is not None and args.baseline_summary.expanduser().resolve().exists()
        else None
    )
    baseline_name = (
        str(args.baseline_summary.expanduser().resolve())
        if args.baseline_summary is not None
        else ""
    )

    selected_years = [str(year) for year in args.years]
    run_specs = []
    for label, config in PRIMARY_RUNS:
        years = [year for year in config["years"] if year in selected_years]
        if not years:
            continue
        run_specs.append((label, {"years": years, "dit_only": bool(config["dit_only"])}))
    if args.include_secondary:
        for label, config in SECONDARY_RUNS:
            years = [year for year in config["years"] if year in selected_years]
            if not years:
                continue
            run_specs.append((label, {"years": years, "dit_only": bool(config["dit_only"])}))

    consolidated_rows: list[dict[str, Any]] = []
    per_run_details: dict[str, dict[str, Any]] = {}

    for label, config in run_specs:
        run_dir = output_root / label
        summary = run_validation(
            args.data_root,
            run_dir,
            workers=max(1, int(args.workers)),
            selected_years=config["years"],
            limit=max(0, int(args.limit)),
            input_manifest=None,
            include_npm1=bool(args.include_npm1 and not config["dit_only"]),
            suppress_worker_output=True,
            timeout_seconds=max(0, int(args.timeout_seconds)),
            checkpoint_every=max(0, int(args.checkpoint_every)),
            dit_only=bool(config["dit_only"]),
            required_run_name=str(args.require_run_name_contains or ""),
        )

        bundle_dir = output_root / "bundles" / label
        _copy_summary_bundle(run_dir, bundle_dir)
        metrics_rows = _load_metrics_rows(run_dir)
        row = _summary_row(label, summary)
        row.update(_diff_vs_baseline(row, baseline_summary if label == "all_years_dit" else None))
        consolidated_rows.append(row)
        per_run_details[label] = {
            "summary": summary,
            "breakdown": _subset_breakdown(metrics_rows),
        }

    consolidated_payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_root": str(args.data_root.expanduser().resolve()),
        "output_root": str(output_root),
        "baseline_summary": baseline_name,
        "runs": consolidated_rows,
        "details": per_run_details,
    }
    _write_json(output_root / "verification_summary.json", consolidated_payload)
    (output_root / "verification_summary.md").write_text(
        _markdown_report(consolidated_rows, baseline_name or "(none)", per_run_details),
        encoding="utf-8",
    )

    consolidated_csv = output_root / "verification_summary.csv"
    with consolidated_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_ordered_fieldnames(consolidated_rows) if consolidated_rows else ["run"])
        writer.writeheader()
        for row in consolidated_rows:
            writer.writerow(row)

    print(json.dumps(consolidated_payload, indent=2, ensure_ascii=False))
    print(f"Verification summary JSON: {output_root / 'verification_summary.json'}")
    print(f"Verification summary Markdown: {output_root / 'verification_summary.md'}")
    print(f"Verification summary CSV: {consolidated_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
