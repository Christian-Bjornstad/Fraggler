#!/usr/bin/env python3
"""
Legacy launcher for the old browser-first packaging model.

Current desktop release builds use `qt_app.py` directly and do not use this file
as the packaged entrypoint.
"""
from __future__ import annotations

import os
import sys
import signal
import threading
import webbrowser

# ── Resolve paths for both frozen (PyInstaller) and dev mode ──
if getattr(sys, "frozen", False):
    # Running as a PyInstaller bundle
    BUNDLE_DIR = sys._MEIPASS          # temp dir with extracted data
    APP_DIR = os.path.dirname(sys.executable)  # where the .exe lives
else:
    # Running from source
    BUNDLE_DIR = os.path.dirname(os.path.abspath(__file__))
    APP_DIR = os.path.dirname(BUNDLE_DIR)  # packaging/ → OUS/

# The actual application code lives in APP_DIR
os.chdir(APP_DIR)
sys.path.insert(0, APP_DIR)

PORT = 5078
HOST = "localhost"


def open_browser_delayed():
    """Open the browser after a short delay so the server can start."""
    import time
    time.sleep(2.5)
    url = f"http://{HOST}:{PORT}/app"
    print(f"\n✦ Opening browser at {url}")
    webbrowser.open(url)


def main():
    print("=" * 60)
    print("  Fraggler Diagnostics")
    print("=" * 60)
    print(f"  App directory : {APP_DIR}")
    print(f"  Bundle dir    : {BUNDLE_DIR}")
    print(f"  Server        : http://{HOST}:{PORT}/app")
    print("=" * 60)
    print()

    # ── Graceful shutdown on Ctrl+C ──
    def handle_sigint(sig, frame):
        print("\n\n✦ Shutting down Fraggler Diagnostics …")
        sys.exit(0)

    signal.signal(signal.SIGINT, handle_sigint)

    # ── Start browser in background ──
    browser_thread = threading.Thread(target=open_browser_delayed, daemon=True)
    browser_thread.start()

    # ── Import Panel and serve ──
    import panel as pn

    pn.extension(
        "plotly",
        "tabulator",
        sizing_mode="stretch_width",
        notifications=True,
    )

    # Serve the app
    pn.serve(
        {"app": os.path.join(APP_DIR, "app.py")},
        port=PORT,
        address=HOST,
        show=False,
        title="Fraggler Diagnostics",
        verbose=False,
    )


if __name__ == "__main__":
    main()
