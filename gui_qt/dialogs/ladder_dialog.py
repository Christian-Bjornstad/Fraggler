from __future__ import annotations

import copy

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QGridLayout,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


class LadderAdjustmentDialog(QDialog):
    def __init__(self, fsa, parent=None):
        super().__init__(parent)
        self.fsa = fsa
        self.setWindowTitle(f"Ladder Adjustment - {fsa.file_name}")
        self.resize(1380, 860)

        self.ladder_steps = np.asarray(fsa.ladder_steps, dtype=float)
        self.candidates = self._get_candidates().reset_index(drop=True)
        self.mapping: dict[int, int] = {}
        self._initial_mapping: dict[int, int] = {}
        self._preview_metrics: dict | None = None

        self._init_ui()
        self._suggest_auto(store_initial=True)
        self._refresh_all()

    def _get_candidates(self):
        from core.analysis import get_ladder_candidates

        return get_ladder_candidates(self.fsa)

    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(12)

        summary_card = QWidget()
        summary_card.setObjectName("Card")
        summary_layout = QVBoxLayout(summary_card)

        summary_title = QLabel("Ladder Adjustment Studio")
        summary_title.setStyleSheet("font-size: 18px; font-weight: 800; color: #0f172a;")
        summary_layout.addWidget(summary_title)

        info_grid = QGridLayout()
        info_grid.setHorizontalSpacing(18)
        info_grid.setVerticalSpacing(8)
        self.meta_labels: dict[str, QLabel] = {}
        meta_rows = [
            ("file", "File"),
            ("ladder", "Ladder"),
            ("candidate_count", "Detected Ladder Peaks"),
            ("mapped_count", "Mapped Steps"),
            ("preview", "Preview"),
        ]
        for row, (key, label) in enumerate(meta_rows):
            left = QLabel(f"{label}:")
            left.setStyleSheet("color: #64748b; font-weight: 700;")
            right = QLabel("—")
            right.setWordWrap(True)
            self.meta_labels[key] = right
            info_grid.addWidget(left, row, 0, alignment=Qt.AlignmentFlag.AlignTop)
            info_grid.addWidget(right, row, 1)
        summary_layout.addLayout(info_grid)

        help_label = QLabel(
            "Workflow: select a ladder step, then either double-click a candidate peak, click a red X in the plot, "
            "or use the assign button. Double-click a mapped ladder row to clear it."
        )
        help_label.setWordWrap(True)
        help_label.setStyleSheet("color: #64748b;")
        summary_layout.addWidget(help_label)
        layout.addWidget(summary_card)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        plot_container = QWidget()
        plot_layout = QVBoxLayout(plot_container)
        self.figure, self.ax = plt.subplots(figsize=(11, 5))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        splitter.addWidget(plot_container)

        side_container = QWidget()
        side_layout = QVBoxLayout(side_container)
        side_layout.setSpacing(10)

        steps_card = QWidget()
        steps_card.setObjectName("Card")
        steps_layout = QVBoxLayout(steps_card)
        steps_title = QLabel("LADDER STEPS")
        steps_title.setObjectName("CardTitle")
        steps_layout.addWidget(steps_title)

        self.table = QTableWidget(len(self.ladder_steps), 5)
        self.table.setHorizontalHeaderLabels(["Ladder BP", "Peak Time", "Intensity", "Candidate", "Status"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.cellDoubleClicked.connect(self._on_step_double_clicked)
        self.table.itemSelectionChanged.connect(self._sync_selection_from_step_table)
        steps_layout.addWidget(self.table)
        side_layout.addWidget(steps_card, stretch=1)

        candidates_card = QWidget()
        candidates_card.setObjectName("Card")
        candidates_layout = QVBoxLayout(candidates_card)
        candidates_title = QLabel("CANDIDATE PEAKS")
        candidates_title.setObjectName("CardTitle")
        candidates_layout.addWidget(candidates_title)

        self.candidate_table = QTableWidget(0 if self.candidates.empty else len(self.candidates), 4)
        self.candidate_table.setHorizontalHeaderLabels(["#", "Time", "Intensity", "Assigned"])
        self.candidate_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.candidate_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.candidate_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.candidate_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self.candidate_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.candidate_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.candidate_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.candidate_table.cellDoubleClicked.connect(self._assign_selected_candidate)
        candidates_layout.addWidget(self.candidate_table)

        candidate_btns = QHBoxLayout()
        self.btn_assign_candidate = QPushButton("Assign Selected Candidate")
        self.btn_assign_candidate.clicked.connect(self._assign_selected_candidate)
        self.btn_clear_step = QPushButton("Clear Selected Step")
        self.btn_clear_step.clicked.connect(self._clear_selected_step)
        candidate_btns.addWidget(self.btn_assign_candidate)
        candidate_btns.addWidget(self.btn_clear_step)
        candidates_layout.addLayout(candidate_btns)
        side_layout.addWidget(candidates_card, stretch=1)

        splitter.addWidget(side_container)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter, stretch=1)

        bottom_layout = QHBoxLayout()
        self.stats_label = QLabel("Preview: not run")
        self.stats_label.setStyleSheet("color: #64748b; font-weight: 600;")
        bottom_layout.addWidget(self.stats_label)
        bottom_layout.addStretch()

        btn_auto = QPushButton("Suggest Auto")
        btn_auto.clicked.connect(lambda: self._suggest_auto(store_initial=False))
        bottom_layout.addWidget(btn_auto)

        btn_reset = QPushButton("Reset To Initial")
        btn_reset.clicked.connect(self._reset_to_initial)
        bottom_layout.addWidget(btn_reset)

        btn_clear_all = QPushButton("Clear All")
        btn_clear_all.clicked.connect(self._clear_all)
        bottom_layout.addWidget(btn_clear_all)

        btn_preview = QPushButton("Preview Fit")
        btn_preview.clicked.connect(self._preview_fit)
        bottom_layout.addWidget(btn_preview)

        btn_cancel = QPushButton("Cancel")
        btn_cancel.clicked.connect(self.reject)
        bottom_layout.addWidget(btn_cancel)

        btn_apply = QPushButton("Apply")
        btn_apply.setObjectName("PrimaryButton")
        btn_apply.clicked.connect(self._on_apply)
        bottom_layout.addWidget(btn_apply)

        layout.addLayout(bottom_layout)
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)

    def _refresh_all(self):
        self._update_meta()
        self._update_step_table()
        self._update_candidate_table()
        self._plot_ladder()

    def _update_meta(self):
        self.meta_labels["file"].setText(self.fsa.file_name)
        self.meta_labels["ladder"].setText(str(self.fsa.ladder))
        self.meta_labels["candidate_count"].setText(str(len(self.candidates)))
        self.meta_labels["mapped_count"].setText(f"{len(self.mapping)} / {len(self.ladder_steps)}")

        if self._preview_metrics:
            r2 = self._preview_metrics.get("r2", float("nan"))
            n = self._preview_metrics.get("n_ladder_steps", 0)
            txt = f"R² {r2:.6f} with {n} ladder steps"
        else:
            txt = "Not previewed yet"
        self.meta_labels["preview"].setText(txt)

    def _candidate_used_by(self, cand_idx: int) -> int | None:
        for step_idx, mapped_idx in self.mapping.items():
            if mapped_idx == cand_idx:
                return step_idx
        return None

    def _update_step_table(self):
        selected_step = self._selected_step_row()
        for row, bp in enumerate(self.ladder_steps):
            self.table.setItem(row, 0, QTableWidgetItem(f"{bp:.0f} bp"))
            if row in self.mapping and not self.candidates.empty:
                cand_idx = self.mapping[row]
                cand = self.candidates.iloc[cand_idx]
                self.table.setItem(row, 1, QTableWidgetItem(f"{float(cand['time']):.0f}"))
                self.table.setItem(row, 2, QTableWidgetItem(f"{float(cand['intensity']):.0f}"))
                self.table.setItem(row, 3, QTableWidgetItem(f"#{cand_idx}"))
                status_text = "Mapped"
                status_item = QTableWidgetItem(status_text)
                status_item.setForeground(Qt.GlobalColor.darkGreen)
                self.table.setItem(row, 4, status_item)
            else:
                self.table.setItem(row, 1, QTableWidgetItem("-"))
                self.table.setItem(row, 2, QTableWidgetItem("-"))
                self.table.setItem(row, 3, QTableWidgetItem("-"))
                status_item = QTableWidgetItem("Missing")
                status_item.setForeground(Qt.GlobalColor.red)
                self.table.setItem(row, 4, status_item)
        if selected_step is not None and 0 <= selected_step < self.table.rowCount():
            self.table.selectRow(selected_step)

    def _update_candidate_table(self):
        selected_candidate = self._selected_candidate_row()
        self.candidate_table.setRowCount(len(self.candidates))
        for row in range(len(self.candidates)):
            cand = self.candidates.iloc[row]
            assigned_step = self._candidate_used_by(row)
            assigned_text = f"{self.ladder_steps[assigned_step]:.0f} bp" if assigned_step is not None else "Free"

            items = [
                QTableWidgetItem(str(row)),
                QTableWidgetItem(f"{float(cand['time']):.0f}"),
                QTableWidgetItem(f"{float(cand['intensity']):.0f}"),
                QTableWidgetItem(assigned_text),
            ]
            if assigned_step is not None:
                items[3].setForeground(Qt.GlobalColor.darkGreen)
            for col, item in enumerate(items):
                self.candidate_table.setItem(row, col, item)
        if selected_candidate is not None and 0 <= selected_candidate < self.candidate_table.rowCount():
            self.candidate_table.selectRow(selected_candidate)

    def _plot_ladder(self):
        self.ax.clear()
        trace = self.fsa.size_standard
        self.ax.plot(trace, color="#94a3b8", alpha=0.8, linewidth=1.2, label="Size Standard")

        selected_candidate = self._selected_candidate_row()
        selected_step = self._selected_step_row()

        if not self.candidates.empty:
            times = self.candidates["time"].to_numpy(dtype=float)
            intensities = self.candidates["intensity"].to_numpy(dtype=float)
            self.ax.scatter(times, intensities, marker="x", color="#ef4444", s=48, linewidths=1.4, label="Candidates")

            if selected_candidate is not None and 0 <= selected_candidate < len(self.candidates):
                c = self.candidates.iloc[selected_candidate]
                self.ax.scatter(
                    [float(c["time"])],
                    [float(c["intensity"])],
                    s=120,
                    facecolors="none",
                    edgecolors="#2563eb",
                    linewidths=2,
                    label="Selected Candidate",
                )

        for step_idx, cand_idx in self.mapping.items():
            cand = self.candidates.iloc[cand_idx]
            peak_time = float(cand["time"])
            peak_intensity = float(cand["intensity"])
            bp = self.ladder_steps[step_idx]
            marker_color = "#2563eb" if step_idx == selected_step else "#22c55e"
            self.ax.scatter([peak_time], [peak_intensity], s=70, color=marker_color, zorder=3)
            self.ax.annotate(
                f"{bp:.0f}",
                (peak_time, peak_intensity),
                textcoords="offset points",
                xytext=(0, 10),
                ha="center",
                fontsize=8,
                color=marker_color,
                fontweight="bold",
            )

        self.ax.set_title(f"Ladder Fitting: {self.fsa.ladder}")
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel("Intensity")
        self.ax.grid(True, alpha=0.25)
        self.ax.legend(loc="upper right")
        self.canvas.draw_idle()

    def _selected_step_row(self) -> int | None:
        selected_rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        return selected_rows[0].row() if selected_rows else None

    def _selected_candidate_row(self) -> int | None:
        selected_rows = self.candidate_table.selectionModel().selectedRows() if self.candidate_table.selectionModel() else []
        return selected_rows[0].row() if selected_rows else None

    def _sync_selection_from_step_table(self):
        step = self._selected_step_row()
        if step is not None and step in self.mapping:
            self.candidate_table.selectRow(self.mapping[step])
        self._plot_ladder()

    def _assign_candidate_to_step(self, step_idx: int, cand_idx: int):
        if step_idx < 0 or step_idx >= len(self.ladder_steps):
            return
        if cand_idx < 0 or cand_idx >= len(self.candidates):
            return

        # Enforce one candidate per ladder step and one ladder step per candidate.
        for other_step, other_cand in list(self.mapping.items()):
            if other_step == step_idx:
                continue
            if other_cand == cand_idx:
                del self.mapping[other_step]
        self.mapping[step_idx] = cand_idx
        self._preview_metrics = None
        self._refresh_all()

        if step_idx + 1 < len(self.ladder_steps):
            self.table.selectRow(step_idx + 1)

    def _clear_selected_step(self):
        step_idx = self._selected_step_row()
        if step_idx is None:
            QMessageBox.information(self, "No Step Selected", "Select a ladder step to clear first.")
            return
        if step_idx in self.mapping:
            del self.mapping[step_idx]
            self._preview_metrics = None
            self._refresh_all()
            self.table.selectRow(step_idx)

    def _clear_all(self):
        self.mapping = {}
        self._preview_metrics = None
        self._refresh_all()
        if self.table.rowCount():
            self.table.selectRow(0)

    def _reset_to_initial(self):
        self.mapping = dict(self._initial_mapping)
        self._preview_metrics = None
        self._refresh_all()
        if self.table.rowCount():
            self.table.selectRow(0)

    def _on_plot_click(self, event):
        if event.inaxes != self.ax or self.candidates.empty or event.xdata is None:
            return

        step_idx = self._selected_step_row()
        if step_idx is None:
            QMessageBox.information(self, "No Step Selected", "Select a ladder step first, then click a candidate peak.")
            return

        x = float(event.xdata)
        nearest_idx = int((self.candidates["time"] - x).abs().idxmin())
        self._assign_candidate_to_step(step_idx, nearest_idx)

    def _on_step_double_clicked(self, row, _column):
        if row in self.mapping:
            del self.mapping[row]
            self._preview_metrics = None
            self._refresh_all()
            self.table.selectRow(row)

    def _assign_selected_candidate(self, *_args):
        step_idx = self._selected_step_row()
        cand_idx = self._selected_candidate_row()
        if step_idx is None:
            QMessageBox.information(self, "No Step Selected", "Select a ladder step first.")
            return
        if cand_idx is None:
            QMessageBox.information(self, "No Candidate Selected", "Select a candidate peak first.")
            return
        self._assign_candidate_to_step(step_idx, cand_idx)

    def _suggest_auto(self, store_initial: bool):
        best = getattr(self.fsa, "best_size_standard", None)
        ss_peaks_raw = getattr(self.fsa, "size_standard_peaks", None)
        ss_peaks = list(ss_peaks_raw) if ss_peaks_raw is not None else []
        auto_mapping: dict[int, int] = {}
        if best is not None and len(best) > 0 and ss_peaks:
            for step_idx, peak_time in enumerate(best):
                if peak_time > 0 and peak_time in ss_peaks:
                    auto_mapping[step_idx] = ss_peaks.index(peak_time)

        self.mapping = auto_mapping
        if store_initial:
            self._initial_mapping = dict(auto_mapping)
        self._preview_metrics = None

    def _preview_fit(self):
        if len(self.mapping) < 3:
            QMessageBox.warning(self, "Invalid Fit", "Select at least 3 peaks to preview fit.")
            return

        from core.analysis import apply_manual_ladder_mapping, compute_ladder_qc_metrics

        try:
            preview_fsa = copy.deepcopy(self.fsa)
            preview_fsa = apply_manual_ladder_mapping(preview_fsa, self.mapping)
            metrics = compute_ladder_qc_metrics(preview_fsa)
            self._preview_metrics = metrics
            r2 = metrics["r2"]
            n_steps = metrics.get("n_ladder_steps", 0)
            n_found = metrics.get("n_size_standard_peaks", 0)
            self.stats_label.setText(f"Preview fit OK: R² {r2:.6f} | mapped {len(self.mapping)} | fitted peaks {n_found}/{n_steps}")
            if r2 > 0.999:
                self.stats_label.setStyleSheet("color: #22c55e; font-weight: 700;")
            else:
                self.stats_label.setStyleSheet("color: #f59e0b; font-weight: 700;")
            self._update_meta()
        except Exception as exc:
            self._preview_metrics = None
            self.stats_label.setText(f"Preview failed: {exc}")
            self.stats_label.setStyleSheet("color: #ef4444; font-weight: 700;")
            self._update_meta()
            QMessageBox.critical(self, "Preview Failed", f"Could not fit this mapping:\n{exc}")

    def _on_apply(self):
        if not self.mapping:
            QMessageBox.warning(self, "No Mapping", "Map at least one ladder step before applying.")
            return
        self.accept()

    def get_mapping(self):
        return dict(self.mapping)
