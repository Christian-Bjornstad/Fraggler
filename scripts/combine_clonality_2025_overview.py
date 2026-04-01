from __future__ import annotations

from scripts import combine_clonality_yearly_overview as yearly


build_arg_parser = yearly.build_arg_parser
combine_run_root = yearly.combine_run_root


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    output_path = args.output or (args.run_root / f"track-clonality-{args.year_label}-overview.xlsx")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = combine_run_root(args.run_root, output_path, year_label=str(args.year_label))
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
