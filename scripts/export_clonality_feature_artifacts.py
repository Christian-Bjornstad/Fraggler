from __future__ import annotations

import argparse
import json
from pathlib import Path

from core.analyses.clonality.feature_artifacts import write_clonality_feature_artifacts


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export clonality feature artifacts for ladder and PK training.")
    parser.add_argument("--workbook", type=Path, required=True, help="Path to Clonality_Tracking.xlsx")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for CSV/JSON feature artifacts")
    parser.add_argument(
        "--metadata-json",
        type=Path,
        default=None,
        help="Optional JSON file containing entry metadata or an identity-key map.",
    )
    parser.add_argument(
        "--include-sl",
        action="store_true",
        help="Include SL rows in exported feature tables for monitoring.",
    )
    return parser


def _load_metadata(path: Path | None):
    if path is None:
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    metadata = _load_metadata(args.metadata_json)
    outputs = write_clonality_feature_artifacts(
        args.workbook,
        args.output_dir,
        entry_metadata=metadata,
        include_sl=args.include_sl,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
