from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize human-reviewed ladder bundle labels and selected-peak agreement."
    )
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        required=True,
        help="Directory containing ladder_review_cases.csv and ladder_review_candidates.csv.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        help="Optional JSON path. Defaults to <bundle-dir>/ladder_review_eval.json.",
    )
    return parser


def _selected_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"1", "true", "yes"})


def evaluate_ladder_review_bundle(bundle_dir: Path) -> dict[str, object]:
    case_path = bundle_dir / "ladder_review_cases.csv"
    candidate_path = bundle_dir / "ladder_review_candidates.csv"
    cases = pd.read_csv(case_path).fillna("")
    candidates = pd.read_csv(candidate_path).fillna("")

    reviewed_cases = cases[cases["label"].astype(str).str.strip() != ""].copy()
    reviewed_candidates = candidates[candidates["human_label"].astype(str).str.strip() != ""].copy()

    selected = reviewed_candidates[_selected_mask(reviewed_candidates["selected_for_fit"])].copy()
    selected_keep = selected[selected["human_label"] == "keep_peak"]
    selected_reject = selected[selected["human_label"] == "reject_peak"]
    keep = reviewed_candidates[reviewed_candidates["human_label"] == "keep_peak"]

    precision = float(len(selected_keep) / len(selected)) if len(selected) else 0.0
    recall = float(len(selected_keep) / len(keep)) if len(keep) else 0.0
    f1 = float((2 * precision * recall) / (precision + recall)) if (precision + recall) else 0.0

    summary: dict[str, object] = {
        "bundle_dir": str(bundle_dir),
        "reviewed_case_count": int(len(reviewed_cases)),
        "reviewed_candidate_count": int(len(reviewed_candidates)),
        "case_label_counts": reviewed_cases["label"].value_counts().to_dict(),
        "candidate_label_counts": reviewed_candidates["human_label"].value_counts().to_dict(),
        "selected_peak_count": int(len(selected)),
        "selected_keep_count": int(len(selected_keep)),
        "selected_reject_count": int(len(selected_reject)),
        "selected_peak_precision": precision,
        "selected_peak_recall": recall,
        "selected_peak_f1": f1,
    }
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    summary = evaluate_ladder_review_bundle(args.bundle_dir)
    output_json = args.output_json or (args.bundle_dir / "ladder_review_eval.json")
    output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
