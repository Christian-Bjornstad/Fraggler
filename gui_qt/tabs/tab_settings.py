from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QDoubleSpinBox, QCheckBox, QFileDialog
)
from pathlib import Path
from PyQt6.QtCore import pyqtSignal
from config import APP_SETTINGS, save_settings

class TabSettings(QWidget):
    settings_saved = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(16)
        
        # Header
        header = QVBoxLayout()
        title = QLabel("Settings")
        title.setObjectName("PageTitle")
        sub = QLabel("Configure standard properties.")
        sub.setObjectName("PageSubtitle")
        header.addWidget(title)
        header.addWidget(sub)
        
        # Card
        card = QWidget()
        card.setObjectName("Card")
        c_layout = QFormLayout(card)
        
        s_gen = APP_SETTINGS.get("general", {})
        s_qc = APP_SETTINGS.get("qc", {})
        s_batch = APP_SETTINGS.get("batch", {})
        
        # Paths
        c_layout.addRow(QLabel("<b>Path Defaults</b>"))
        
        row_in = QHBoxLayout()
        self.default_input = QLineEdit(s_batch.get("base_input_dir", str(Path.home())))
        btn_browse_in = QPushButton("Browse...")
        btn_browse_in.clicked.connect(lambda: self._browse_dir(self.default_input))
        row_in.addWidget(self.default_input, stretch=1)
        row_in.addWidget(btn_browse_in)
        c_layout.addRow("Default Input Folder:", row_in)
        
        row_out = QHBoxLayout()
        self.default_output = QLineEdit(s_batch.get("output_base", str(Path.home())))
        btn_browse_out = QPushButton("Browse...")
        btn_browse_out.clicked.connect(lambda: self._browse_dir(self.default_output))
        row_out.addWidget(self.default_output, stretch=1)
        row_out.addWidget(btn_browse_out)
        c_layout.addRow("Default Output Folder:", row_out)

        self.author = QLineEdit(s_gen.get("author", "OUS"))
        
        c_layout.addRow("Author (for PDF templates):", self.author)
        
        # QC Parameters
        c_layout.addRow(QLabel("<b>QC Parameters</b>"))
        self.d_min_r2_ok = QDoubleSpinBox(); self.d_min_r2_ok.setRange(0, 1); self.d_min_r2_ok.setSingleStep(0.001); self.d_min_r2_ok.setDecimals(3); self.d_min_r2_ok.setValue(float(s_qc.get("min_r2_ok", 0.995)))
        self.d_min_r2_warn = QDoubleSpinBox(); self.d_min_r2_warn.setRange(0, 1); self.d_min_r2_warn.setSingleStep(0.001); self.d_min_r2_warn.setDecimals(3); self.d_min_r2_warn.setValue(float(s_qc.get("min_r2_warn", 0.990)))
        
        c_layout.addRow("Min R² (OK):", self.d_min_r2_ok)
        c_layout.addRow("Min R² (WARN):", self.d_min_r2_warn)
        
        # Pipeline/Batch Parameters
        c_layout.addRow(QLabel("<b>Batch Options</b>"))
        self.chk_agg_dit = QCheckBox("Aggregate Multiple Folders to 1 DIT Report")
        self.chk_agg_dit.setChecked(bool(s_batch.get("aggregate_dit_reports", True)))
        self.chk_agg_pat = QCheckBox("Group scans by Patient ID (Regex)")
        self.chk_agg_pat.setChecked(bool(s_batch.get("aggregate_by_patient", True)))
        
        c_layout.addRow("", self.chk_agg_dit)
        c_layout.addRow("", self.chk_agg_pat)
        
        btn_save = QPushButton("Save Settings")
        btn_save.setObjectName("PrimaryButton")
        btn_save.clicked.connect(self.save)
        c_layout.addRow("", btn_save)
        
        self.status_lbl = QLabel("")
        self.status_lbl.setStyleSheet("color: #22c55e; font-weight: 500;")
        c_layout.addRow("", self.status_lbl)
        
        main_layout.addLayout(header)
        main_layout.addWidget(card)
        main_layout.addStretch()
        
    def save(self):
        s_gen = APP_SETTINGS.setdefault("general", {})
        s_gen["author"] = self.author.text()
        
        s_qc = APP_SETTINGS.setdefault("qc", {})
        s_qc["min_r2_ok"] = self.d_min_r2_ok.value()
        s_qc["min_r2_warn"] = self.d_min_r2_warn.value()
        
        s_batch = APP_SETTINGS.setdefault("batch", {})
        s_batch["base_input_dir"] = self.default_input.text()
        s_batch["output_base"] = self.default_output.text()
        s_batch["aggregate_dit_reports"] = self.chk_agg_dit.isChecked()
        s_batch["aggregate_by_patient"] = self.chk_agg_pat.isChecked()
        
        save_settings(APP_SETTINGS)
        self.settings_saved.emit()
        self.status_lbl.setText("Settings saved.")

    def _browse_dir(self, line_edit):
        folder = QFileDialog.getExistingDirectory(self, "Select Directory", line_edit.text() or str(Path.home()))
        if folder:
            line_edit.setText(folder)
