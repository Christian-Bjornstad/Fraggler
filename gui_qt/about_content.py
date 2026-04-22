from __future__ import annotations

from pathlib import Path

from app_meta import APP_NAME, APP_VERSION


REPO_ROOT = Path(__file__).resolve().parents[1]
THIRD_PARTY_NOTICE_PATH = REPO_ROOT / "THIRD_PARTY_NOTICES.md"
UPSTREAM_LICENSE_PATH = REPO_ROOT / "LICENSES" / "fraggler_MIT.txt"

APP_OVERVIEW = {
    "title": f"About {APP_NAME}",
    "subtitle": (
        "Clinical fragment-analysis application for assay processing, QC review, "
        "tracking workbooks, and interactive HTML reporting."
    ),
    "version_label": f"Version {APP_VERSION}",
    "owner_context": (
        "This application is maintained as a local clinical diagnostics workflow tool and "
        "includes both custom code and selected upstream-derived components."
    ),
    "repo_license_status": (
        "The repository does not currently publish a single root open-source license for the "
        "whole project. Third-party MIT notices are preserved for the upstream-derived parts."
    ),
}

THIRD_PARTY_SOFTWARE = [
    {
        "name": "fraggler",
        "homepage": "https://github.com/willros/fraggler",
        "authors": "William Rosenbaum and Pär Larsson",
        "copyright": "Clinical Genomic Umea",
        "license_name": "MIT",
        "summary": (
            "HemaFrag Diagnostics includes and builds on local embedded code derived from the "
            "upstream fraggler package, especially in the bundled fraggler runtime module."
        ),
        "derived_paths": [
            "fraggler/",
            "core/analysis.py",
        ],
    },
]


def load_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return f"Missing file: {path}"
