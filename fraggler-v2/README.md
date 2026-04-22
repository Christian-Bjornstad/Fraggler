# HemaFrag Rust Engine

`fraggler-v2` is the Rust-first engine workspace for HemaFrag Diagnostics.

Current scope in this scaffold:
- establish the new workspace layout
- define the first version of the engine/desktop JSON contract
- stand up the CLI and desktop shell entrypoints
- freeze Python reference timings and outputs before engine porting begins
- start Phase 2 engine porting in `fraggler-core`
- wire the Slint desktop shell to real Rust analyze runs

## Phase 2 status

Implemented in Rust so far:
- native ABIF/FSA directory parsing
- raw signal extraction for `DATA*` channels and `DyeN*` metadata
- baseline correction primitives
- peak detection primitives
- first-pass ladder candidate generation and curvature scoring
- first-pass sizing-model fit and ladder QC metrics for the best candidate
- conservative local refinement around the best ladder candidate when QC suggests it
- sample-trace to basepair mapping preview from the fitted sizing model
- sample-peak preview in basepair space from the mapped sample trace
- assay-agnostic sample peak grouping preview in basepair space
- first-pass clonality assay matching using filename, channel compatibility, bp overlap, and group-level dominance signals
- first-pass FLT3 assay preview using filename-aware assay detection, preferred channel selection, and WT/mutant peak preview
- Willros-style monotone spline ladder sizing as the primary Rust sizing strategy, with polynomial fallback
- multi-file `analyze` request handling in the core engine
- desktop shell fields for analysis, input path, output path, and live Rust engine log/status updates
- desktop shell browse buttons for input/output selection, run summary, clear-log action, and direct open-output/open-artifact actions
- persisted analyze artifacts (`analyze_summary.json`, `primitive_result_preview.json`) with desktop open-output/open-artifact actions

Current limitation:
- the engine does not yet reproduce the full Python ladder-fit/refinement/rescue workflow
- QC, FLT3 validation, and report generation remain scaffold-only in Rust

## Workspace layout

- `crates/fraggler-core`
  - shared contract types
  - engine boundary
  - report payload types
- `crates/fraggler-cli`
  - headless commands for `analyze`, `qc`, `validate-flt3`, `build-report`
- `crates/fraggler-desktop`
  - Slint desktop shell scaffold
- `schemas/fraggler-contract-v1.schema.json`
  - v1 JSON contract for desktop ↔ core/CLI communication
- `baselines/scenarios.example.json`
  - example baseline-freeze scenarios for the Python reference implementation

## Product rules

- Offline-first desktop delivery
- Native support for macOS, Windows, and Linux
- Standalone HTML/Plotly reports remain the report model
- Python remains the reference implementation until Rust passes parity + performance gates

## Expected development flow

1. Freeze baseline outputs and timings with the Python implementation.
2. Port the engine into `fraggler-core`.
3. Drive the engine through `fraggler-cli`.
4. Attach the Slint shell after the engine reaches parity.

## Baseline freeze

Use the existing Python application to freeze reference timings and artifacts:

```bash
python3 scripts/freeze_v2_baseline.py \
  --scenario-file fraggler-v2/baselines/scenarios.example.json \
  --output-dir validation_outputs/fraggler_v2_baseline
```

The script writes:
- `baseline_manifest.json`
- per-scenario summaries
- copied artifacts when configured by the scenario

## Example requests

Example v1 request payloads live in `fraggler-v2/examples/`:
- `analyze_request.json`
- `qc_request.json`
- `build_report_request.json`

Use them with:

```bash
cd fraggler-v2
cargo run -p fraggler-cli -- analyze --json-request examples/analyze_request.json
```

## Build intent

Current local verification commands:

```bash
cd fraggler-v2
cargo fmt
cargo check
cargo test -p fraggler-core
cargo run -p fraggler-cli -- analyze --help
cargo run -p fraggler-desktop
```
