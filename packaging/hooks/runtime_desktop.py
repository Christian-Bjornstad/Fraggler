"""
Runtime defaults for packaged Fraggler desktop builds.
"""
from __future__ import annotations

import os
import sys

# Packaged desktop builds should not start the legacy embedded Panel server by default.
os.environ.setdefault("FRAGGLER_ENABLE_LEGACY_PANEL", "0")
# Packaged desktop builds should stay single-process to avoid GUI child launches.
os.environ.setdefault("FRAGGLER_DISABLE_MULTIPROCESSING", "1")

# Fedora 35 / offline Linux target: prefer X11/xcb for stable Qt startup.
if sys.platform.startswith("linux"):
    os.environ.setdefault("QT_QPA_PLATFORM", "xcb")
