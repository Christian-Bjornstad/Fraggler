# Third-Party Notices

Fraggler Diagnostics contains and builds on MIT-licensed code and concepts from the upstream
[`willros/fraggler`](https://github.com/willros/fraggler) project.

## Upstream project

- Project: `fraggler`
- Upstream repository: <https://github.com/willros/fraggler>
- Authors listed upstream: William Rosenbaum and Pär Larsson
- Upstream license: MIT
- Upstream copyright notice: `Copyright (c) 2023 Clinical Genomic Umea`

The MIT license text used by the upstream project is included verbatim in
[`LICENSES/fraggler_MIT.txt`](LICENSES/fraggler_MIT.txt).

## Local areas derived from or closely based on upstream

The following areas in this repository should be treated as likely upstream-derived or vendored
foundations that preserve the upstream MIT notice:

- [`fraggler/`](fraggler)
  - Local embedded runtime module used by the application
  - Includes ladder definitions and FSA/baseline-related logic derived from the upstream package
- [`core/analysis.py`](core/analysis.py)
  - Contains explicit references to upstream Fraggler logic and builds on the local `fraggler`
    runtime module

These areas have been substantially modified and extended for the Fraggler Diagnostics application,
including clinical workflow, GUI, QC tracking, reporting, archive validation, and assay-specific
logic that are not part of the upstream package.

## Repository licensing status

This repository does **not** currently declare a single root open-source license for the whole
project. The application contains third-party MIT-licensed components and derived code, while the
repository as a whole remains otherwise unpublished under a separate open-source license.

This notice is intended to preserve attribution and license visibility for upstream-derived code.
