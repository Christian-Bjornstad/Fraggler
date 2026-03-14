"""
Fraggler Diagnostics — Main Entry Point for PyQt6 UI
"""
import sys
import os
from pathlib import Path

# Force X11 (xcb) on Linux to avoid Wayland symbol mismatches (e.g., wl_proxy_marshal_flags)
if sys.platform == "linux":
    os.environ["QT_QPA_PLATFORM"] = "xcb"

from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QIcon

from gui_qt.main_window import MainWindow
from core.log import log

LEGACY_PANEL_HOST = "localhost"
LEGACY_PANEL_PORT = 5078
LEGACY_PANEL_ENABLED = os.environ.get("FRAGGLER_ENABLE_LEGACY_PANEL", "").lower() in {
    "1",
    "true",
    "yes",
}

def start_panel_server():
    """Start the legacy Panel server when explicitly requested."""
    try:
        import panel as pn
        # Resolve path to app.py relative to this file
        bundle_dir = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
        app_path = os.path.join(bundle_dir, "app.py")
        
        if not os.path.exists(app_path):
            log(f"[WARN] Could not find legacy app.py at {app_path}. Server not started.")
            return

        log(
            f"[INFO] Starting legacy Panel server at "
            f"http://{LEGACY_PANEL_HOST}:{LEGACY_PANEL_PORT}/app"
        )

        pn.serve(
            {"app": app_path},
            port=LEGACY_PANEL_PORT,
            address=LEGACY_PANEL_HOST,
            show=False,
            title="Fraggler Diagnostics",
            verbose=False,
        )
    except Exception as e:
        log(f"[ERROR] Failed to start legacy web server: {e}")

def exception_hook(exctype, value, tb):
    """Global exception handler to prevent silent crashes in slots."""
    import traceback
    from PyQt6.QtWidgets import QMessageBox
    
    err_msg = "".join(traceback.format_exception(exctype, value, tb))
    print(err_msg, file=sys.stderr)
    
    # Try to show a message box
    try:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Critical)
        msg.setText("An unexpected error occurred.")
        msg.setInformativeText(str(value))
        msg.setDetailedText(err_msg)
        msg.setWindowTitle("Error")
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()
    except:
        pass
    
    sys.__excepthook__(exctype, value, tb)

sys.excepthook = exception_hook

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("Fraggler")
    app.setOrganizationName("OUS")
    
    icon_path = Path(__file__).parent / "assets" / "app_icon.png"
    if icon_path.exists():
        app_icon = QIcon(str(icon_path))
        app.setWindowIcon(app_icon)
    
    if LEGACY_PANEL_ENABLED:
        import threading
        server_thread = threading.Thread(target=start_panel_server, daemon=True)
        server_thread.start()

    window = MainWindow()
    window.resize(1200, 800)
    if icon_path.exists():
        window.setWindowIcon(app_icon)
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
