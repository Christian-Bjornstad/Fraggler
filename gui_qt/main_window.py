from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QStackedWidget, QPushButton, QLabel, QFrame, QComboBox
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QIcon

from gui_qt.styles import VIBRANT_PRO_QSS
from gui_qt.tabs.tab_batch import TabBatch
from gui_qt.tabs.tab_log import TabLog
from gui_qt.tabs.tab_settings import TabSettings
from config import APP_SETTINGS, save_settings

class SidebarButton(QPushButton):
    def __init__(self, text, icon_name=None, parent=None):
        super().__init__(text, parent)
        self.setObjectName("SidebarButton")
        self.setCheckable(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

class AnalysisGroupHeader(QPushButton):
    """The main button for an analysis type (e.g. Klonalitet)."""
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setObjectName("AnalysisGroupHeader")
        self.setCheckable(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

class AnalysisSubButton(QPushButton):
    """The sub-buttons that appear when a group is expanded."""
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setObjectName("AnalysisSubButton")
        self.setCheckable(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)

class AnalysisGroup(QWidget):
    """Container for an analysis header and its sub-buttons."""
    def __init__(self, name, internal_id, on_sub_clicked, parent=None):
        super().__init__(parent)
        self.internal_id = internal_id
        self.on_sub_clicked = on_sub_clicked
        
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)
        
        self.header = AnalysisGroupHeader(name)
        self.layout.addWidget(self.header)
        
        self.sub_container = QWidget()
        self.sub_layout = QVBoxLayout(self.sub_container)
        self.sub_layout.setContentsMargins(0, 0, 0, 0)
        self.sub_layout.setSpacing(0)
        
        self.btn_run = AnalysisSubButton("•  Run")
        self.btn_log = AnalysisSubButton("•  Log")
        self.btn_settings = AnalysisSubButton("•  Settings")
        
        self.sub_buttons = [self.btn_run, self.btn_log, self.btn_settings]
        for i, btn in enumerate(self.sub_buttons):
            self.sub_layout.addWidget(btn)
            btn.clicked.connect(lambda _, b=btn, idx=i: self._handle_sub_click(b, idx))
            
        self.layout.addWidget(self.sub_container)
        self.sub_container.setVisible(False)
        
        self.header.clicked.connect(self.toggle_expansion)
        
    def _handle_sub_click(self, clicked_btn, tab_idx):
        for btn in self.sub_buttons:
            btn.setChecked(btn == clicked_btn)
        self.on_sub_clicked(self.internal_id, tab_idx)

    def toggle_expansion(self):
        # We handle expansion control from MainWindow to ensure only one is open
        pass
        
    def set_expanded(self, expanded):
        self.header.setChecked(expanded)
        self.sub_container.setVisible(expanded)
        if not expanded:
            for btn in self.sub_buttons:
                btn.setChecked(False)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Fraggler Diagnostics v1.1.1")
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
        
        sidebar_layout.addSpacing(10)
        
        # --- Analysis Groups ---
        self.groups = []
        
        self.group_clonality = AnalysisGroup("Klonalitet", "clonality", self.on_sub_tab_clicked)
        self.group_flt3 = AnalysisGroup("FLT3 Analysis", "flt3", self.on_sub_tab_clicked)
        
        self.groups = [self.group_clonality, self.group_flt3]
        for g in self.groups:
            sidebar_layout.addWidget(g)
            g.header.clicked.connect(lambda _, grp=g: self.on_group_clicked(grp))
            
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
        active_ana = APP_SETTINGS.get("active_analysis", "clonality")
        start_group = self.group_clonality if active_ana == "clonality" else self.group_flt3
        self.on_group_clicked(start_group)
        start_group.btn_run.setChecked(True)
        self.stacked_widget.setCurrentIndex(0)
        
    def on_group_clicked(self, group):
        # Update active analysis in core
        new_ana = group.internal_id
        if APP_SETTINGS.get("active_analysis") != new_ana:
            APP_SETTINGS["active_analysis"] = new_ana
            save_settings(APP_SETTINGS)
            print(f"[UI] Analysis switched to: {new_ana}")
            # Refresh tabs if needed
            self.tab_run._detected_jobs = []
            self.tab_run._rebuild_table()

        # Update Sidebar expansion
        for g in self.groups:
            expanded = (g == group)
            g.set_expanded(expanded)
            if expanded:
                # Automatically select the first sub-tab (Run) when expanding a new group
                g.btn_run.setChecked(True)
                self.on_sub_tab_clicked(g.internal_id, 0)
            
    def on_sub_tab_clicked(self, analysis_id, tab_idx):
        # Ensure we are on the right analysis
        if APP_SETTINGS.get("active_analysis") != analysis_id:
            # This shouldn't happen with our set_expanded logic but good to have
            APP_SETTINGS["active_analysis"] = analysis_id
            save_settings(APP_SETTINGS)
            
        self.stacked_widget.setCurrentIndex(tab_idx)
