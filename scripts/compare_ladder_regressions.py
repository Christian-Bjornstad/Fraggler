from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare FLT3 ladder validation output with clonality hard-case ladder summaries."
    )
    parser.add_argument("--flt3-summary", type=Path, required=True, help="Path to FLT3 summary.json.")
    parser.add_argument(
        "--clonality-summary",
        type=Path,
        required=True,
        help="Path to clonality run_summary.json from validate_clonality_hard_cases.py.",
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
        help="Output format.",
    )
    return parser


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.expanduser().resolve().read_text(encoding="utf-8"))


def _flt3_row(summary: dict[str, Any]) -> dict[str, Any]:
    status_counts = summary.get("status_counts", {}) if isinstance(summary.get("status_counts"), dict) else {}
    return {
        "dataset": "FLT3",
        "ok": int(status_counts.get("ok", 0)),
        "review_required": int(status_counts.get("review_required", 0)),
        "ladder_fit_failed": int(status_counts.get("ladder_fit_failed", 0)),
        "manual_review_count": int(summary.get("manual_review_count", 0) or 0),
        "top_review_reason": (
            next(iter(summary.get("worst_cases", []) or []), {}).get("review_reason", "")
            if isinstance(summary.get("worst_cases"), list)
            else ""
        ),
        "top_strategy": max(
            (summary.get("strategy_counts", {}) or {}).items(),
            key=lambda item: item[1],
            default=("", 0),
        )[0],
    }


def _clonality_row(summary: dict[str, Any]) -> dict[str, Any]:
    ladder_validation = summary.get("ladder_validation", {}) if isinstance(summary.get("ladder_validation"), dict) else {}
    status_counts = ladder_validation.get("status_counts", {}) if isinstance(ladder_validation.get("status_counts"), dict) else {}
    return {
        "dataset": "Clonality hard-cases",
        "ok": int(status_counts.get("ok", 0)),
        "review_required": int(status_counts.get("review_required", 0)),
        "ladder_fit_failed": int(status_counts.get("ladder_fit_failed", 0)),
        "manual_review_count": int(ladder_validation.get("review_required_count", 0) or 0),
        "top_review_reason": "",
        "top_strategy": "",
    }


def _markdown_table(rows: list[dict[str, Any]]) -> str:
    headers = [
        "Dataset",
        "OK",
        "Review Required",
        "Ladder Fit Failed",
        "Manual Review Count",
        "Top Review Reason",
        "Top Strategy",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["dataset"]),
                    str(row["ok"]),
                    str(row["review_required"]),
                    str(row["ladder_fit_failed"]),
                    str(row["manual_review_count"]),
                    str(row["top_review_reason"]),
                    str(row["top_strategy"]),
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    flt3_summary = _load_json(args.flt3_summary)
    clonality_summary = _load_json(args.clonality_summary)

    payload = {
        "rows": [
            _flt3_row(flt3_summary),
            _clonality_row(clonality_summary),
        ]
    }
    if args.format == "json":
        print(json.dumps(payload, indent=2, ensure_ascii=False))
    else:
        print(_markdown_table(payload["rows"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
