from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QStackedWidget, QPushButton, QLabel, QFrame
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QIcon

from gui_qt.styles import VIBRANT_PRO_QSS
from gui_qt.tabs.tab_batch import TabBatch
from gui_qt.tabs.tab_log import TabLog
from gui_qt.tabs.tab_settings import TabSettings

class SidebarButton(QPushButton):
    def __init__(self, text, icon_name=None, parent=None):
        super().__init__(text, parent)
        self.setObjectName("SidebarButton")
        self.setCheckable(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Fraggler Diagnostics")
        self.setStyleSheet(VIBRANT_PRO_QSS)
        
        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        
        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # --- Sidebar ---
        self.sidebar_container = QWidget()
        self.sidebar_container.setObjectName("Sidebar")
        self.sidebar_container.setFixedWidth(220)
        
        sidebar_layout = QVBoxLayout(self.sidebar_container)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        sidebar_layout.setSpacing(4)
        
        brand = QLabel("FRAGGLER")
        brand.setObjectName("SidebarBrand")
        sidebar_layout.addWidget(brand)
        
        # Nav Buttons
        self.btn_run = SidebarButton("Run")
        self.btn_log = SidebarButton("Log")
        self.btn_settings = SidebarButton("Settings")
        
        self.nav_buttons = [self.btn_run, self.btn_log, self.btn_settings]
        
        for btn in self.nav_buttons:
            sidebar_layout.addWidget(btn)
            btn.clicked.connect(self.on_nav_clicked)
            
        sidebar_layout.addStretch()
        
        # --- Stacked Widget (Content) ---
        self.stacked_widget = QStackedWidget()
        
        # Tabs
        self.tab_run = TabBatch()
        self.tab_log = TabLog()
        self.tab_settings = TabSettings()
        
        # Connect settings saved to reload in batch tab
        self.tab_settings.settings_saved.connect(self.tab_run.load_from_settings)
        
        # Connect global core logging to this tab
        from gui_qt.log_handler import qt_log_handler
        qt_log_handler.emitter.log_signal.connect(self.tab_log.append_log)
        
        # Add to stack
        self.stacked_widget.addWidget(self.tab_run)
        self.stacked_widget.addWidget(self.tab_log)
        self.stacked_widget.addWidget(self.tab_settings)
        
        # Content Container (for padding)
        content_container = QWidget()
        content_layout = QVBoxLayout(content_container)
        content_layout.setContentsMargins(32, 32, 32, 32)
        content_layout.addWidget(self.stacked_widget)
        
        # Add to main
        main_layout.addWidget(self.sidebar_container)
        main_layout.addWidget(content_container, stretch=1)
        
        # Initialize
        self.btn_run.setChecked(True)
        self.stacked_widget.setCurrentIndex(0)
        
    def on_nav_clicked(self):
        sender = self.sender()
        # Enforce single selection
        for i, btn in enumerate(self.nav_buttons):
            if btn == sender:
                btn.setChecked(True)
                self.stacked_widget.setCurrentIndex(i)
            else:
                btn.setChecked(False)
