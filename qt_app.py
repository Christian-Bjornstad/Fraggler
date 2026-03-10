"""
Fraggler Diagnostics — Main Entry Point for PyQt6 UI
"""
import sys
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QIcon
from pathlib import Path

from gui_qt.main_window import MainWindow

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
    
    window = MainWindow()
    window.resize(1200, 800)
    if icon_path.exists():
        window.setWindowIcon(app_icon)
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
