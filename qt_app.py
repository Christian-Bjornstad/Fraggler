"""
Fraggler Diagnostics — Main Entry Point for PyQt6 UI
"""
import multiprocessing
import sys
import os
import locale
from pathlib import Path

from app_meta import APP_VERSION

# Force X11 (xcb) on Linux to avoid Wayland symbol mismatches (e.g., wl_proxy_marshal_flags)
if sys.platform == "linux":
    os.environ["QT_QPA_PLATFORM"] = "xcb"
    os.environ.setdefault("LANG", "C.UTF-8")
    os.environ.setdefault("LC_ALL", "C.UTF-8")
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")

from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QIcon

from core.log import log

LEGACY_PANEL_HOST = "localhost"
LEGACY_PANEL_PORT = 5078
LEGACY_PANEL_ENABLED = os.environ.get("FRAGGLER_ENABLE_LEGACY_PANEL", "").lower() in {
    "1",
    "true",
    "yes",
}


def _remove_macos_metadata_files(bundle_dir: Path) -> None:
    """Delete AppleDouble/Finder metadata files that can break Linux runtime imports."""
    patterns = ("._*", ".DS_Store")
    removed = 0
    for pattern in patterns:
        for path in bundle_dir.rglob(pattern):
            if not path.is_file():
                continue
            try:
                path.unlink()
                removed += 1
            except OSError:
                continue
    if removed:
        log(f"[INFO] Removed {removed} macOS metadata file(s) from bundle: {bundle_dir}")


def _prepare_runtime_bundle() -> Path:
    bundle_dir = Path(getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__))))
    _remove_macos_metadata_files(bundle_dir)
    return bundle_dir


_BUNDLE_DIR = _prepare_runtime_bundle()

from gui_qt.main_window import MainWindow

def start_panel_server():
    """Start the legacy Panel server when explicitly requested."""
    try:
        import panel as pn
        # Resolve path to app.py relative to this file
        app_path = os.path.join(_BUNDLE_DIR, "app.py")
        
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
    if sys.platform == "linux":
        try:
            locale.setlocale(locale.LC_ALL, "")
        except locale.Error:
            try:
                locale.setlocale(locale.LC_ALL, "C.UTF-8")
            except locale.Error:
                pass

    app = QApplication(sys.argv)
    app.setApplicationName("Fraggler")
    app.setOrganizationName("OUS")
    app.setApplicationVersion(APP_VERSION)
    
    icon_path = _BUNDLE_DIR / "assets" / "app_icon.png"
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
    # Ensure stdout/stderr use UTF-8 regardless of environment locales
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')

    multiprocessing.freeze_support()
    
    try:
        main()
    except Exception as e:
        import traceback
        err_msg = traceback.format_exc()
        print(f"CRITICAL STARTUP ERROR:\n{err_msg}", file=sys.stderr)
        
        # If QApplication was already created (unlikely here but for safety)
        # we try to show a message box.
        app = QApplication.instance()
        if not app:
            app = QApplication(sys.argv)
            
        from PyQt6.QtWidgets import QMessageBox
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Critical)
        msg.setWindowTitle("Startup Error")
        msg.setText("Fraggler Diagnostics failed to start.")
        msg.setInformativeText(str(e))
        msg.setDetailedText(err_msg)
        msg.exec()
        sys.exit(1)
