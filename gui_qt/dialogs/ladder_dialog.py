import numpy as np
import pandas as pd
import copy
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QMessageBox, QSplitter, QWidget
)
from PyQt6.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class LadderAdjustmentDialog(QDialog):
    def __init__(self, fsa, parent=None):
        super().__init__(parent)
        self.fsa = fsa
        self.setWindowTitle(f"Ladder Adjustment - {fsa.file_name}")
        self.resize(1000, 700)
        
        self.ladder_steps = fsa.ladder_steps
        self.candidates = self._get_candidates()
        self.mapping = {} # {step_index: candidate_index}
        
        self._init_ui()
        self._suggest_auto()
        self._plot_ladder()
        self._update_table()

    def _get_candidates(self):
        from core.analysis import get_ladder_candidates
        return get_ladder_candidates(self.fsa)

    def _init_ui(self):
        layout = QVBoxLayout(self)
        
        # Splitter for plot and table
        splitter = QSplitter(Qt.Orientation.Vertical)
        
        # Plot area
        plot_container = QWidget()
        plot_layout = QVBoxLayout(plot_container)
        self.figure, self.ax = plt.subplots(figsize=(10, 4))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        splitter.addWidget(plot_container)
        
        # Table area
        table_container = QWidget()
        table_layout = QVBoxLayout(table_container)
        
        info_label = QLabel("Map Ladder BP to detected peaks. Select a row, then click a peak (X) on the plot.")
        table_layout.addWidget(info_label)
        
        self.table = QTableWidget(len(self.ladder_steps), 4)
        self.table.setHorizontalHeaderLabels(["Ladder BP", "Peak Time", "Intensity", "Status"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.cellDoubleClicked.connect(self._on_cell_double_clicked)
        
        table_layout.addWidget(self.table)
        splitter.addWidget(table_container)
        
        layout.addWidget(splitter)
        
        # Stats & Buttons
        bottom_layout = QHBoxLayout()
        self.stats_label = QLabel("R²: Unknown")
        bottom_layout.addWidget(self.stats_label)
        bottom_layout.addStretch()
        
        btn_auto = QPushButton("Suggest Auto")
        btn_auto.clicked.connect(self._suggest_auto)
        bottom_layout.addWidget(btn_auto)
        
        btn_preview = QPushButton("Preview Fit")
        btn_preview.clicked.connect(self._preview_fit)
        bottom_layout.addWidget(btn_preview)
        
        btn_cancel = QPushButton("Cancel")
        btn_cancel.clicked.connect(self.reject)
        bottom_layout.addWidget(btn_cancel)
        
        btn_apply = QPushButton("Apply")
        btn_apply.setObjectName("PrimaryButton")
        btn_apply.clicked.connect(self.accept)
        bottom_layout.addWidget(btn_apply)
        
        layout.addLayout(bottom_layout)
        
        # Connect plot clicks
        self.canvas.mpl_connect('button_press_event', self._on_plot_click)

    def _plot_ladder(self):
        self.ax.clear()
        trace = self.fsa.size_standard
        self.ax.plot(trace, color='#94a3b8', alpha=0.7, label="Size Standard")
        
        # Plot candidates
        if not self.candidates.empty:
            times = self.candidates['time'].values
            intensities = self.candidates['intensity'].values
            self.ax.plot(times, intensities, 'x', color='#ef4444', label="Candidates", picker=5)
            
        # Plot currently mapped peaks
        for step_idx, cand_idx in self.mapping.items():
            peak_time = self.candidates.iloc[cand_idx]['time']
            peak_intensity = self.candidates.iloc[cand_idx]['intensity']
            bp = self.ladder_steps[step_idx]
            self.ax.plot(peak_time, peak_intensity, 'o', color='#22c55e')
            self.ax.annotate(f"{bp}", (peak_time, peak_intensity), 
                             textcoords="offset points", xytext=(0,10), ha='center',
                             fontsize=8, color='#15803d', fontweight='bold')
            
        self.ax.set_title(f"Ladder Fitting: {self.fsa.ladder}")
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel("Intensity")
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(loc='upper right')
        self.canvas.draw()

    def _update_table(self):
        for i, bp in enumerate(self.ladder_steps):
            self.table.setItem(i, 0, QTableWidgetItem(f"{bp} bp"))
            
            if i in self.mapping:
                cand_idx = self.mapping[i]
                cand = self.candidates.iloc[cand_idx]
                self.table.setItem(i, 1, QTableWidgetItem(f"{cand['time']:.0f}"))
                self.table.setItem(i, 2, QTableWidgetItem(f"{cand['intensity']:.0f}"))
                self.table.setItem(i, 3, QTableWidgetItem("Mapped"))
                self.table.item(i, 3).setForeground(Qt.GlobalColor.darkGreen)
            else:
                self.table.setItem(i, 1, QTableWidgetItem("-"))
                self.table.setItem(i, 2, QTableWidgetItem("-"))
                self.table.setItem(i, 3, QTableWidgetItem("Missing"))
                self.table.item(i, 3).setForeground(Qt.GlobalColor.red)

    def _on_plot_click(self, event):
        if event.inaxes != self.ax:
            return
            
        if self.candidates.empty:
            return
            
        x = event.xdata
        # Find nearest candidate index
        nearest_idx = (self.candidates['time'] - x).abs().idxmin()
        
        # Which row is selected in the table?
        selected_rows = self.table.selectionModel().selectedRows()
        if not selected_rows:
            QMessageBox.information(self, "No Row Selected", "Select a ladder step in the table first.")
            return
            
        row = selected_rows[0].row()
        self.mapping[row] = nearest_idx
        
        # Auto-select next row for convenience
        if row + 1 < len(self.ladder_steps):
            self.table.selectRow(row + 1)
            
        self._update_table()
        self._plot_ladder()

    def _suggest_auto(self):
        best = self.fsa.best_size_standard
        if best is None or len(best) == 0:
            return
            
        ss_peaks = list(self.fsa.size_standard_peaks)
        for i, peak_time in enumerate(best):
            if peak_time > 0 and peak_time in ss_peaks:
                cand_idx = ss_peaks.index(peak_time)
                self.mapping[i] = cand_idx
                
        self._update_table()
        self._plot_ladder()

    def _on_cell_double_clicked(self, row, column):
        if row in self.mapping:
            del self.mapping[row]
            self._update_table()
            self._plot_ladder()

    def _preview_fit(self):
        if len(self.mapping) < 3:
            QMessageBox.warning(self, "Invalid Fit", "Select at least 3 peaks to preview fit.")
            return
            
        from core.analysis import apply_manual_ladder_mapping, compute_ladder_qc_metrics
        try:
            preview_fsa = copy.deepcopy(self.fsa)
            preview_fsa = apply_manual_ladder_mapping(preview_fsa, self.mapping)
            metrics = compute_ladder_qc_metrics(preview_fsa)
            r2 = metrics['r2']
            self.stats_label.setText(f"Preview R²: {r2:.6f}")
            if r2 > 0.999:
                self.stats_label.setStyleSheet("color: #22c55e;")
            else:
                self.stats_label.setStyleSheet("color: #f59e0b;")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not fit: {e}")

    def get_mapping(self):
        return self.mapping
