from __future__ import annotations

import argparse
import pickle
from pathlib import Path

from core.analyses.clonality.candidate_artifacts import write_clonality_candidate_artifacts


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export candidate-level clonality artifacts from collected live entries.")
    parser.add_argument("--entries-pickle", type=Path, required=True, help="Path to a pickle file containing collected clonality entries.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for CSV/JSON candidate artifacts.")
    parser.add_argument("--include-sl", action="store_true", help="Include SL rows in exported candidate tables.")
    parser.add_argument(
        "--write-gold-label-template",
        action="store_true",
        help="Also write an empty clonality_gold_labels.csv template beside the candidate artifacts.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    entries = pickle.loads(Path(args.entries_pickle).expanduser().read_bytes())
    outputs = write_clonality_candidate_artifacts(
        args.output_dir,
        entries,
        include_sl=args.include_sl,
        write_gold_label_template=args.write_gold_label_template,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
