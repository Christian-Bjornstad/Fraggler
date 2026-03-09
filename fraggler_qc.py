"""
Backward-compatible facade — imports from core.qc.* submodules.

This file is kept for backward compatibility so that existing code
that does ``import fraggler_qc`` or
``from fraggler_qc import QCRules, build_qc_html`` still works.
"""
from __future__ import annotations

# Re-export everything from the QC subpackage
from core.qc.qc_rules import *  # noqa: F401,F403
from core.qc.qc_markers import *  # noqa: F401,F403
from core.qc.qc_plots import *  # noqa: F401,F403
from core.qc.qc_excel import *  # noqa: F401,F403
from core.qc.qc_html import *  # noqa: F401,F403
from core.qc.qc_main import main  # noqa: F401


if __name__ == "__main__":
    main()