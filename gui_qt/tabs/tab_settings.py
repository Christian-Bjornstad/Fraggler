from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFormLayout, QDoubleSpinBox, QCheckBox
)
from config import APP_SETTINGS, save_settings

class TabSettings(QWidget):
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
        s_batch["aggregate_dit_reports"] = self.chk_agg_dit.isChecked()
        s_batch["aggregate_by_patient"] = self.chk_agg_pat.isChecked()
        
        save_settings(APP_SETTINGS)
        self.status_lbl.setText("Settings saved.")
