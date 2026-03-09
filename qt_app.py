"""
Fraggler Diagnostics — Main Entry Point for PyQt6 UI
"""
import sys
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QIcon
from pathlib import Path

from gui_qt.main_window import MainWindow

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
