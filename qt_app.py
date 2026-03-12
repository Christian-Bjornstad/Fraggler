"""
Fraggler Diagnostics — Main Entry Point for PyQt6 UI
"""
import sys
import os

# Force X11 (xcb) on Linux to avoid Wayland symbol mismatches (e.g., wl_proxy_marshal_flags)
if sys.platform == "linux":
    os.environ["QT_QPA_PLATFORM"] = "xcb"

from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QIcon
from pathlib import Path

from gui_qt.main_window import MainWindow

# ── Web Server Integration ──
import threading
import time

def start_panel_server():
    """Start the Panel server in a separate thread."""
    try:
        import panel as pn
        # Resolve path to app.py relative to this file
        bundle_dir = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
        app_path = os.path.join(bundle_dir, "app.py")
        
        if not os.path.exists(app_path):
            print(f"WARN: Could not find app.py at {app_path}. Server not started.")
            return

        print(f"✦ Starting background web server at http://{HOST}:{PORT}/app")

        pn.serve(
            {"app": app_path},
            port=PORT,
            address=HOST,
            show=False,
            title="Fraggler Diagnostics",
            verbose=False,
        )
    except Exception as e:
        print(f"ERROR: Failed to start web server: {e}")

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
    
    # Check if Inter font is available, else we fallback automatically defined in styles.py
    
    # Start web server in background thread
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
