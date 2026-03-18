from __future__ import annotations

import copy
import math

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QDialog,
    QGridLayout,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QScrollArea,
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


PASS_R2 = 0.9995
CHECK_R2 = 0.9990
PASS_MAX_ABS_RESIDUAL = 0.5
CHECK_MAX_ABS_RESIDUAL = 1.5


class LadderAdjustmentDialog(QDialog):
    def __init__(self, fsa, parent=None):
        super().__init__(parent)
        self.fsa = fsa
        self.setWindowTitle(f"Ladder Adjustment - {fsa.file_name}")
        screen = QApplication.primaryScreen()
        available = screen.availableGeometry() if screen is not None else None
        if available is not None:
            target_w = min(available.width() - 64, 1700)
            target_h = min(available.height() - 64, 1040)
            self.resize(max(target_w, 1180), max(target_h, 780))
            self.setMinimumSize(min(1180, available.width()), min(780, available.height()))
        else:
            self.resize(1660, 980)
            self.setMinimumSize(1180, 780)

        self.fitted_ladder_steps = np.asarray(fsa.ladder_steps, dtype=float)
        self.ladder_steps = np.asarray(
            getattr(fsa, "expected_ladder_steps", self.fitted_ladder_steps),
            dtype=float,
        )
        self.candidates = self._get_candidates().reset_index(drop=True)
        self.mapping: dict[int, int] = {}
        self._initial_mapping: dict[int, int] = {}
        self._manual_candidate_times: list[float] = []
        self._add_peak_mode = False
        self._preview_fsa = None
        self._preview_metrics: dict | None = None
        self._fit_rows: list[dict] = []
        self._fit_grade = "unknown"
        self._fit_reason = "Preview not run"
        self._missing_order = "ascending"

        self._init_ui()
        self._suggest_auto(store_initial=True)
        self._refresh_preview_state(show_errors=False)
        self._refresh_all()
        self._focus_initial_step()

    def _get_candidates(self):
        from core.analysis import get_ladder_candidates

        df = get_ladder_candidates(self.fsa).copy()
        if "source" not in df.columns:
            df["source"] = "auto"
        return df

    def _init_ui(self):
        self.setObjectName("LadderDialog")
        self.setStyleSheet(
            """
            QDialog#LadderDialog {
                background: #eef4fb;
            }
            QWidget#WorkspaceCard {
                background: #ffffff;
                border: 1px solid #d7e4f3;
                border-radius: 18px;
            }
            QLabel#WorkspaceTitle {
                font-size: 19px;
                font-weight: 800;
                color: #0f172a;
            }
            QLabel#WorkspaceSubtitle {
                color: #53657f;
                font-size: 12px;
                font-weight: 500;
            }
            QLabel#WorkspaceEyebrow {
                color: #6b7b95;
                font-size: 11px;
                font-weight: 700;
                letter-spacing: 1px;
            }
            QWidget#MetaChip {
                background: #f7fbff;
                border: 1px solid #dce8f5;
                border-radius: 14px;
            }
            QLabel#MetaValue {
                color: #10233d;
                font-weight: 700;
            }
            QLabel#MetaChipValue {
                color: #10233d;
                font-weight: 800;
                font-size: 14px;
            }
            QLabel#MetaChipLabel {
                color: #69809d;
                font-size: 10px;
                font-weight: 700;
                letter-spacing: 1px;
            }
            QLabel#CardTitle {
                color: #10233d;
                font-size: 13px;
                font-weight: 800;
                letter-spacing: 0.6px;
            }
            QTableWidget {
                background: #fbfdff;
                border: 1px solid #e1ebf5;
                border-radius: 14px;
                gridline-color: #e8eef6;
                alternate-background-color: #f4f8fc;
                selection-background-color: #d9ebff;
                selection-color: #10233d;
            }
            QHeaderView::section {
                background: #f3f7fc;
                color: #516784;
                border: none;
                border-right: 1px solid #e3ebf5;
                border-bottom: 1px solid #e3ebf5;
                padding: 8px 10px;
                font-weight: 700;
            }
            QSplitter::handle {
                background: #dbe6f2;
                border-radius: 3px;
            }
            QSplitter::handle:horizontal {
                width: 8px;
            }
            QSplitter::handle:vertical {
                height: 8px;
            }
            QPushButton {
                background: #f7faff;
                color: #17314f;
                border: 1px solid #d4e2f1;
                border-radius: 12px;
                padding: 10px 16px;
                font-weight: 700;
                min-height: 40px;
            }
            QPushButton:hover {
                background: #edf5ff;
                border-color: #bdd3ea;
            }
            QPushButton:pressed {
                background: #e3eefb;
            }
            QPushButton#PrimaryButton {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #0891b2, stop:1 #4f46e5);
                color: white;
                border: none;
                padding: 11px 18px;
            }
            QPushButton#PrimaryButton:hover {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #0ea5c7, stop:1 #5b54f0);
            }
            """
        )
        outer_layout = QVBoxLayout(self)
        outer_layout.setSpacing(0)
        outer_layout.setContentsMargins(0, 0, 0, 0)

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setFrameShape(QScrollArea.Shape.NoFrame)
        outer_layout.addWidget(scroll_area)

        content = QWidget()
        scroll_area.setWidget(content)

        layout = QVBoxLayout(content)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 16, 16, 16)

        summary_card = QWidget()
        summary_card.setObjectName("WorkspaceCard")
        summary_layout = QVBoxLayout(summary_card)
        summary_layout.setContentsMargins(18, 14, 18, 12)
        summary_layout.setSpacing(8)

        summary_title = QLabel("Ladder Adjustment Studio")
        summary_title.setObjectName("WorkspaceTitle")
        summary_layout.addWidget(summary_title)

        summary_subtitle = QLabel("Manual correction workspace for ladder matching, fit review, and safe save/apply.")
        summary_subtitle.setObjectName("WorkspaceSubtitle")
        summary_layout.addWidget(summary_subtitle)

        info_row = QHBoxLayout()
        info_row.setSpacing(10)
        self.meta_labels: dict[str, QLabel] = {}
        meta_rows = [
            ("file", "File"),
            ("ladder", "Ladder"),
            ("expected_count", "Expected Ladder Sizes"),
            ("candidate_count", "Detected Ladder Peaks"),
            ("mapped_count", "Mapped Steps"),
            ("preview", "Preview"),
        ]
        for key, label in meta_rows:
            chip = QWidget()
            chip.setObjectName("MetaChip")
            chip_layout = QVBoxLayout(chip)
            chip_layout.setContentsMargins(12, 9, 12, 9)
            chip_layout.setSpacing(2)
            left = QLabel(label.upper())
            left.setObjectName("MetaChipLabel")
            right = QLabel("—")
            right.setObjectName("MetaChipValue")
            right.setWordWrap(True)
            self.meta_labels[key] = right
            chip_layout.addWidget(left)
            chip_layout.addWidget(right)
            info_row.addWidget(chip, 1)
        summary_layout.addLayout(info_row)

        help_label = QLabel(
            "Select a ladder row, then click a candidate peak or use the tables. Double-click a mapped row to clear it."
        )
        help_label.setWordWrap(True)
        help_label.setStyleSheet("color: #64748b;")
        summary_layout.addWidget(help_label)
        layout.addWidget(summary_card)

        main_splitter = QSplitter(Qt.Orientation.Vertical)
        main_splitter.setChildrenCollapsible(False)
        main_splitter.setHandleWidth(8)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.setHandleWidth(8)

        plot_container = QWidget()
        plot_container.setObjectName("WorkspaceCard")
        plot_layout = QVBoxLayout(plot_container)
        plot_layout.setContentsMargins(14, 14, 14, 14)
        plot_layout.setSpacing(10)
        plot_title = QLabel("LADDER TRACE")
        plot_title.setObjectName("WorkspaceEyebrow")
        plot_layout.addWidget(plot_title)
        self.figure, self.ax = plt.subplots(figsize=(11, 5))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setFixedHeight(38)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        splitter.addWidget(plot_container)

        side_container = QWidget()
        side_container.setObjectName("WorkspaceCard")
        side_container.setMinimumWidth(420)
        side_layout = QVBoxLayout(side_container)
        side_layout.setContentsMargins(14, 14, 14, 14)
        side_layout.setSpacing(10)
        side_title = QLabel("EDITOR")
        side_title.setObjectName("WorkspaceEyebrow")
        side_layout.addWidget(side_title)

        side_splitter = QSplitter(Qt.Orientation.Vertical)
        side_splitter.setChildrenCollapsible(False)
        side_splitter.setHandleWidth(8)

        steps_card = QWidget()
        steps_card.setObjectName("WorkspaceCard")
        steps_layout = QVBoxLayout(steps_card)
        steps_layout.setContentsMargins(12, 12, 12, 12)
        steps_title = QLabel("LADDER MATCHES")
        steps_title.setObjectName("CardTitle")
        steps_layout.addWidget(steps_title)

        self.missing_steps_label = QLabel("Missing ladder sizes: none")
        self.missing_steps_label.setWordWrap(True)
        self.missing_steps_label.setStyleSheet("color: #64748b; font-weight: 600;")
        steps_layout.addWidget(self.missing_steps_label)

        self.missing_list = QListWidget()
        self.missing_list.setMaximumHeight(110)
        self.missing_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.missing_list.itemSelectionChanged.connect(self._sync_selection_from_missing_list)
        steps_layout.addWidget(self.missing_list)

        self.table = QTableWidget(len(self.ladder_steps), 6)
        self.table.setHorizontalHeaderLabels(["Expected bp", "Observed pos", "Assignment", "Residual", "Confidence", "Status"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(30)
        self.table.setWordWrap(False)
        self.table.setMinimumHeight(220)
        self.table.cellDoubleClicked.connect(self._on_step_double_clicked)
        self.table.itemSelectionChanged.connect(self._sync_selection_from_match_table)
        steps_layout.addWidget(self.table)
        side_splitter.addWidget(steps_card)

        candidates_card = QWidget()
        candidates_card.setObjectName("WorkspaceCard")
        candidates_layout = QVBoxLayout(candidates_card)
        candidates_layout.setContentsMargins(12, 12, 12, 12)
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
        self.candidate_table.setAlternatingRowColors(True)
        self.candidate_table.verticalHeader().setVisible(False)
        self.candidate_table.verticalHeader().setDefaultSectionSize(30)
        self.candidate_table.setMinimumHeight(170)
        self.candidate_table.cellDoubleClicked.connect(self._assign_selected_candidate)
        self.candidate_table.itemSelectionChanged.connect(self._sync_selection_from_candidate_table)
        candidates_layout.addWidget(self.candidate_table)

        candidate_btns_top = QHBoxLayout()
        self.btn_add_peak = QPushButton("Add Missing From Plot")
        self.btn_add_peak.setCheckable(True)
        self.btn_add_peak.toggled.connect(self._toggle_add_peak_mode)
        self.btn_next_missing = QPushButton("Next Missing")
        self.btn_next_missing.clicked.connect(self._select_next_missing_step)
        self.btn_missing_order = QPushButton()
        self.btn_missing_order.setCheckable(True)
        self.btn_missing_order.toggled.connect(self._toggle_missing_order)
        candidate_btns_top.addWidget(self.btn_add_peak)
        candidate_btns_top.addWidget(self.btn_next_missing)
        candidate_btns_top.addWidget(self.btn_missing_order)
        candidates_layout.addLayout(candidate_btns_top)

        candidate_btns_bottom = QHBoxLayout()
        self.btn_assign_candidate = QPushButton("Assign Selected Candidate")
        self.btn_assign_candidate.clicked.connect(self._assign_selected_candidate)
        self.btn_clear_step = QPushButton("Clear Selected Step")
        self.btn_clear_step.clicked.connect(self._clear_selected_step)
        candidate_btns_bottom.addWidget(self.btn_assign_candidate)
        candidate_btns_bottom.addWidget(self.btn_clear_step)
        candidates_layout.addLayout(candidate_btns_bottom)
        side_splitter.addWidget(candidates_card)

        side_layout.addWidget(side_splitter, stretch=1)

        splitter.addWidget(side_container)
        splitter.setStretchFactor(0, 5)
        splitter.setStretchFactor(1, 3)
        splitter.setSizes([1080, 520])
        side_splitter.setSizes([360, 240])

        qc_card = QWidget()
        qc_card.setObjectName("WorkspaceCard")
        qc_layout = QVBoxLayout(qc_card)
        qc_layout.setContentsMargins(16, 14, 16, 14)
        qc_layout.setSpacing(10)
        qc_title = QLabel("SIZING QC")
        qc_title.setObjectName("CardTitle")
        qc_layout.addWidget(qc_title)

        qc_header = QHBoxLayout()
        self.qc_grade_label = QLabel("UNKNOWN")
        self.qc_grade_label.setStyleSheet("font-size: 16px; font-weight: 800; color: #64748b;")
        self.qc_summary_label = QLabel("Preview not run")
        self.qc_summary_label.setWordWrap(True)
        self.qc_summary_label.setStyleSheet("color: #475569; font-weight: 600;")
        qc_header.addWidget(self.qc_grade_label)
        qc_header.addSpacing(12)
        qc_header.addWidget(self.qc_summary_label, stretch=1)
        qc_layout.addLayout(qc_header)

        self.qc_reason_label = QLabel("Map all ladder steps to inspect residuals and sizing quality.")
        self.qc_reason_label.setWordWrap(True)
        self.qc_reason_label.setStyleSheet("color: #64748b;")
        qc_layout.addWidget(self.qc_reason_label)

        self.residual_figure, self.residual_ax = plt.subplots(figsize=(11, 1.85))
        self.residual_canvas = FigureCanvas(self.residual_figure)
        qc_layout.addWidget(self.residual_canvas)
        main_splitter.addWidget(splitter)
        main_splitter.addWidget(qc_card)
        main_splitter.setStretchFactor(0, 5)
        main_splitter.setStretchFactor(1, 2)
        main_splitter.setSizes([640, 260])
        layout.addWidget(main_splitter, stretch=1)

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
        self._sync_missing_order_button()

    def _refresh_all(self):
        self._update_meta()
        self._update_match_table()
        self._update_candidate_table()
        self._update_missing_steps_label()
        self._plot_ladder()
        self._update_qc_panel()
        self._plot_residuals()

    def _update_meta(self):
        self.meta_labels["file"].setText(self.fsa.file_name)
        self.meta_labels["ladder"].setText(str(self.fsa.ladder))
        self.meta_labels["expected_count"].setText(str(len(self.ladder_steps)))
        manual_count = 0
        if "source" in self.candidates.columns:
            manual_count = int(self.candidates["source"].astype(str).eq("manual").sum())
        candidate_text = str(len(self.candidates))
        if manual_count:
            candidate_text += f" ({manual_count} manual)"
        self.meta_labels["candidate_count"].setText(candidate_text)
        self.meta_labels["mapped_count"].setText(f"{len(self.mapping)} / {len(self.ladder_steps)}")

        if self._preview_metrics:
            r2 = self._preview_metrics.get("r2", float("nan"))
            n = self._preview_metrics.get("n_ladder_steps", 0)
            fit_method = self._fit_method_name()
            txt = f"{self._fit_grade.upper()} · {fit_method} · R² {r2:.6f} with {n} ladder steps"
        else:
            txt = "Not previewed yet"
        self.meta_labels["preview"].setText(txt)

    def _candidate_used_by(self, cand_idx: int) -> int | None:
        for step_idx, mapped_idx in self.mapping.items():
            if mapped_idx == cand_idx:
                return step_idx
        return None

    def _row_fit_state(self, row: int) -> dict:
        base = {
            "expected_bp": float(self.ladder_steps[row]),
            "observed_pos": None,
            "assignment": "Missing",
            "residual": None,
            "confidence": "None",
            "status": "Missing",
        }
        if row >= len(self._fit_rows):
            return base
        return {**base, **self._fit_rows[row]}

    def _missing_step_indices(self) -> list[int]:
        missing = [idx for idx in range(len(self.ladder_steps)) if idx not in self.mapping]
        if self._missing_order == "descending":
            missing.reverse()
        return missing

    def _focus_initial_step(self):
        missing = self._missing_step_indices()
        if missing:
            self.table.selectRow(missing[0])
            return
        if self.table.rowCount():
            self.table.selectRow(0)

    def _update_missing_steps_label(self):
        missing = self._missing_step_indices()
        self.missing_list.clear()
        if not missing:
            self.missing_steps_label.setText("Missing ladder sizes: none")
            self.missing_steps_label.setStyleSheet("color: #16a34a; font-weight: 700;")
            return
        for idx in missing:
            item = QListWidgetItem(f"{self.ladder_steps[idx]:.0f} bp")
            item.setData(Qt.ItemDataRole.UserRole, idx)
            self.missing_list.addItem(item)
        bp_text = ", ".join(f"{self.ladder_steps[idx]:.0f} bp" for idx in missing)
        self.missing_steps_label.setText(f"Missing ladder sizes ({len(missing)} remaining): {bp_text}")
        self.missing_steps_label.setStyleSheet("color: #dc2626; font-weight: 700;")

    def _sync_selection_from_missing_list(self):
        items = self.missing_list.selectedItems()
        if not items:
            return
        step_idx = items[0].data(Qt.ItemDataRole.UserRole)
        if step_idx is not None:
            self.table.selectRow(int(step_idx))

    def _next_missing_step(self, current_step: int | None = None) -> int | None:
        missing = self._missing_step_indices()
        if not missing:
            return None
        if current_step is None:
            return missing[0]
        if self._missing_order == "descending":
            for step in missing:
                if step < current_step:
                    return step
        else:
            for step in missing:
                if step > current_step:
                    return step
        return missing[0]

    def _sync_missing_order_button(self):
        if self._missing_order == "descending":
            self.btn_missing_order.setText("Order: High → Low")
            self.btn_missing_order.setChecked(True)
        else:
            self.btn_missing_order.setText("Order: Low → High")
            self.btn_missing_order.setChecked(False)

    def _recommended_missing_order(self) -> str:
        missing = [idx for idx in range(len(self.ladder_steps)) if idx not in self.mapping]
        if not missing:
            return "ascending"

        ladder_len = len(self.ladder_steps)
        low_end_missing = sum(1 for idx in missing if idx < max(3, ladder_len // 3))
        high_end_mapped = sum(1 for idx in self.mapping if idx >= ladder_len // 2)
        if low_end_missing and high_end_mapped:
            return "descending"
        return "ascending"

    def _toggle_missing_order(self, checked: bool):
        self._missing_order = "descending" if checked else "ascending"
        self._sync_missing_order_button()
        self._update_missing_steps_label()
        step_idx = self._selected_step_row()
        next_missing = self._next_missing_step(step_idx)
        if next_missing is not None:
            self.table.selectRow(next_missing)

    def _select_next_missing_step(self):
        step_idx = self._next_missing_step(self._selected_step_row())
        if step_idx is None:
            QMessageBox.information(self, "No Missing Steps", "All ladder steps are currently assigned.")
            return
        self.table.selectRow(step_idx)
        self.stats_label.setText(
            f"Selected missing ladder step {self.ladder_steps[step_idx]:.0f} bp. Click the trace to add its peak."
        )
        self.stats_label.setStyleSheet("color: #0f766e; font-weight: 700;")

    def _update_match_table(self):
        selected_step = self._selected_step_row()
        for row in range(len(self.ladder_steps)):
            row_state = self._row_fit_state(row)
            items = [
                QTableWidgetItem(f"{row_state['expected_bp']:.0f} bp"),
                QTableWidgetItem("—" if row_state["observed_pos"] is None else f"{row_state['observed_pos']:.0f}"),
                QTableWidgetItem(str(row_state["assignment"])),
                QTableWidgetItem("—" if row_state["residual"] is None else f"{row_state['residual']:+.2f} bp"),
                QTableWidgetItem(str(row_state["confidence"])),
                QTableWidgetItem(str(row_state["status"])),
            ]
            status = str(row_state["status"])
            if status == "Missing":
                items[5].setForeground(Qt.GlobalColor.red)
            elif status == "Outlier":
                items[5].setForeground(Qt.GlobalColor.darkYellow)
            elif status == "Weak":
                items[5].setForeground(Qt.GlobalColor.darkYellow)
            else:
                items[5].setForeground(Qt.GlobalColor.darkGreen)
            if str(row_state["assignment"]).startswith("Manual"):
                items[2].setForeground(Qt.GlobalColor.darkBlue)
            for col, item in enumerate(items):
                self.table.setItem(row, col, item)
        if selected_step is not None and 0 <= selected_step < self.table.rowCount():
            self.table.selectRow(selected_step)

    def _update_candidate_table(self):
        selected_candidate = self._selected_candidate_row()
        self.candidate_table.setRowCount(len(self.candidates))
        for row in range(len(self.candidates)):
            cand = self.candidates.iloc[row]
            assigned_step = self._candidate_used_by(row)
            assigned_text = f"{self.ladder_steps[assigned_step]:.0f} bp" if assigned_step is not None else "Free"
            source = str(cand.get("source", "auto"))
            row_label = str(row)
            if source == "manual":
                row_label += " *"

            items = [
                QTableWidgetItem(row_label),
                QTableWidgetItem(f"{float(cand['time']):.0f}"),
                QTableWidgetItem(f"{float(cand['intensity']):.0f}"),
                QTableWidgetItem(assigned_text),
            ]
            if assigned_step is not None:
                items[3].setForeground(Qt.GlobalColor.darkGreen)
            if source == "manual":
                items[0].setForeground(Qt.GlobalColor.darkBlue)
                items[1].setForeground(Qt.GlobalColor.darkBlue)
            for col, item in enumerate(items):
                self.candidate_table.setItem(row, col, item)
        if selected_candidate is not None and 0 <= selected_candidate < self.candidate_table.rowCount():
            self.candidate_table.selectRow(selected_candidate)

    def _build_adjustment_payload(self) -> dict:
        mapping_times: dict[int, float] = {}
        for step_idx, cand_idx in self.mapping.items():
            if 0 <= cand_idx < len(self.candidates):
                mapping_times[int(step_idx)] = float(self.candidates.iloc[cand_idx]["time"])
        return {
            "mapping": dict(self.mapping),
            "mapping_times": mapping_times,
            "manual_candidates": list(self._manual_candidate_times),
        }

    def _candidate_time_exists(self, peak_time: float, tolerance: float = 2.0) -> int | None:
        if self.candidates.empty:
            return None
        diff = (self.candidates["time"].astype(float) - float(peak_time)).abs()
        matches = diff[diff <= tolerance]
        if matches.empty:
            return None
        return int(matches.index[0])

    def _find_local_peak_time(self, x_value: float, search_radius: int = 18) -> tuple[float, float]:
        trace = np.asarray(self.fsa.size_standard, dtype=float)
        if trace.size == 0:
            raise ValueError("No size-standard trace available.")
        center = int(round(float(x_value)))
        lo = max(center - search_radius, 0)
        hi = min(center + search_radius + 1, trace.size)
        if lo >= hi:
            raise ValueError("Could not inspect the selected ladder region.")
        window = trace[lo:hi]
        local_index = int(np.argmax(window))
        peak_index = lo + local_index
        return float(peak_index), float(trace[peak_index])

    def _insert_manual_candidate(self, peak_time: float, intensity: float) -> int:
        existing_idx = self._candidate_time_exists(peak_time)
        if existing_idx is not None:
            return existing_idx

        if not any(math.isclose(float(existing), float(peak_time), abs_tol=1e-6) for existing in self._manual_candidate_times):
            self._manual_candidate_times.append(float(peak_time))
            self._manual_candidate_times.sort()

        manual_row = pd.DataFrame(
            [
                {
                    "index": len(self.candidates),
                    "time": float(peak_time),
                    "intensity": float(intensity),
                    "source": "manual",
                }
            ]
        )
        self.candidates = pd.concat([self.candidates, manual_row], ignore_index=True)
        return int(self.candidates.index[-1])

    def _add_manual_peak_from_plot(self, x_value: float, assign_to_step: int | None = None) -> None:
        peak_time, intensity = self._find_local_peak_time(x_value)
        cand_idx = self._insert_manual_candidate(peak_time, intensity)
        if assign_to_step is not None:
            self._assign_candidate_to_step(assign_to_step, cand_idx)
            return
        self._refresh_preview_state(show_errors=False)
        self._refresh_all()
        self.candidate_table.selectRow(cand_idx)

    def _fit_method_name(self) -> str:
        model = getattr(self._preview_fsa or self.fsa, "ladder_model", None)
        if model is None:
            return "unknown"
        name = model.__class__.__name__.lower()
        if "spline" in name:
            return "spline"
        if "poly" in name:
            return "polynomial"
        return name.replace("model", "")

    def _lookup_fitted_bp(self, peak_time: float) -> float | None:
        preview_fsa = self._preview_fsa
        if preview_fsa is None:
            return None
        df = getattr(preview_fsa, "sample_data_with_basepairs", None)
        if df is not None and {"time", "basepairs"}.issubset(df.columns):
            row = df.loc[df["time"] == int(peak_time)]
            if not row.empty:
                return float(row["basepairs"].iloc[0])
        ladder_model = getattr(preview_fsa, "ladder_model", None)
        if ladder_model is not None:
            try:
                return float(ladder_model.predict(np.array([[peak_time]], dtype=float))[0])
            except Exception:
                return None
        return None

    def _candidate_intensity_median(self) -> float:
        if self.candidates.empty:
            return 0.0
        return float(self.candidates["intensity"].median())

    def _build_fit_rows(self) -> list[dict]:
        rows: list[dict] = []
        intensity_median = self._candidate_intensity_median()
        for step_idx, bp in enumerate(self.ladder_steps):
            if step_idx not in self.mapping or self.candidates.empty:
                rows.append(
                    {
                        "expected_bp": float(bp),
                        "observed_pos": None,
                        "assignment": "Missing",
                        "residual": None,
                        "confidence": "None",
                        "status": "Missing",
                    }
                )
                continue

            cand_idx = self.mapping[step_idx]
            cand = self.candidates.iloc[cand_idx]
            peak_time = float(cand["time"])
            intensity = float(cand["intensity"])
            fitted_bp = self._lookup_fitted_bp(peak_time)
            residual = None if fitted_bp is None else float(fitted_bp - bp)

            assignment_prefix = "Auto" if self._initial_mapping.get(step_idx) == cand_idx else "Manual"
            status = "Mapped"
            if residual is not None and abs(residual) > CHECK_MAX_ABS_RESIDUAL:
                status = "Outlier"
            elif intensity_median > 0 and intensity < intensity_median * 0.35:
                status = "Weak"
            elif assignment_prefix == "Manual":
                status = "Manual"

            if residual is None:
                confidence = "Low"
            elif abs(residual) <= 0.35 and (intensity_median <= 0 or intensity >= intensity_median * 0.6):
                confidence = "High"
            elif abs(residual) <= 1.0:
                confidence = "Medium"
            else:
                confidence = "Low"

            rows.append(
                {
                    "expected_bp": float(bp),
                    "observed_pos": peak_time,
                    "assignment": f"{assignment_prefix} #{cand_idx}",
                    "residual": residual,
                    "confidence": confidence,
                    "status": status,
                }
            )
        return rows

    def _grade_preview_state(self) -> tuple[str, str]:
        missing_count = sum(1 for row in self._fit_rows if row["status"] == "Missing")
        outlier_count = sum(1 for row in self._fit_rows if row["status"] == "Outlier")
        if self._preview_metrics is None:
            if len(self.mapping) < 3:
                return "check", "Map at least 3 ladder steps to preview the fit."
            if missing_count:
                return "check", f"{missing_count} ladder step(s) are still missing from the current edit."
            return "unknown", "Preview not run"

        r2 = float(self._preview_metrics.get("r2", float("nan")))
        max_abs = float(self._preview_metrics.get("max_abs_error_bp", float("inf")))
        if missing_count or outlier_count or r2 < CHECK_R2 or max_abs > CHECK_MAX_ABS_RESIDUAL:
            return "fail", "Fit needs attention: missing steps, low R², or high residual outlier detected."
        if r2 < PASS_R2 or max_abs > PASS_MAX_ABS_RESIDUAL:
            return "check", "Fit is usable, but one or more residuals still need review."
        return "pass", "Stable ladder fit with low residuals across mapped steps."

    def _refresh_preview_state(self, show_errors: bool) -> None:
        self._preview_fsa = None
        self._preview_metrics = None
        self._fit_rows = []
        self._fit_grade = "unknown"
        self._fit_reason = "Preview not run"

        if len(self.mapping) < 3:
            self._fit_rows = self._build_fit_rows()
            self._fit_grade, self._fit_reason = self._grade_preview_state()
            return

        missing_steps = [idx for idx in range(len(self.ladder_steps)) if idx not in self.mapping]
        if missing_steps:
            self._fit_rows = self._build_fit_rows()
            self._fit_grade, self._fit_reason = self._grade_preview_state()
            return

        from core.analysis import apply_manual_ladder_mapping, compute_ladder_qc_metrics

        try:
            preview_fsa = copy.deepcopy(self.fsa)
            preview_fsa.expected_ladder_steps = np.array(self.ladder_steps, dtype=float).copy()
            preview_fsa.ladder_steps = np.array(self.ladder_steps, dtype=float).copy()
            preview_fsa = apply_manual_ladder_mapping(preview_fsa, self._build_adjustment_payload())
            self._preview_fsa = preview_fsa
            self._preview_metrics = compute_ladder_qc_metrics(preview_fsa)
        except Exception as exc:
            self._preview_fsa = None
            self._preview_metrics = None
            self._fit_reason = str(exc)
            if show_errors:
                QMessageBox.critical(self, "Preview Failed", f"Could not fit this mapping:\n{exc}")
        self._fit_rows = self._build_fit_rows()
        self._fit_grade, self._fit_reason = self._grade_preview_state()

    def _update_qc_panel(self):
        color_map = {
            "pass": "#16a34a",
            "check": "#d97706",
            "fail": "#dc2626",
            "unknown": "#64748b",
        }
        label = self._fit_grade.upper()
        self.qc_grade_label.setText(label)
        self.qc_grade_label.setStyleSheet(f"font-size: 16px; font-weight: 800; color: {color_map.get(self._fit_grade, '#64748b')};")

        missing_count = sum(1 for row in self._fit_rows if row["status"] == "Missing")
        extra_count = max(len(self.candidates) - len(self.mapping), 0)
        if self._preview_metrics is None:
            self.qc_summary_label.setText(
                f"{self._fit_method_name()} · mapped {len(self.mapping)}/{len(self.ladder_steps)} · missing {missing_count} · extra {extra_count}"
            )
            self.qc_reason_label.setText(self._fit_reason)
            self.stats_label.setText(f"Preview pending: {self._fit_reason}")
            self.stats_label.setStyleSheet("color: #d97706; font-weight: 700;")
            return

        r2 = float(self._preview_metrics.get("r2", float("nan")))
        mean_abs = float(self._preview_metrics.get("mean_abs_error_bp", float("inf")))
        max_abs = float(self._preview_metrics.get("max_abs_error_bp", float("inf")))
        self.qc_summary_label.setText(
            f"{self._fit_method_name()} · R² {r2:.6f} · mean {mean_abs:.2f} bp · max {max_abs:.2f} bp · missing {missing_count} · extra {extra_count}"
        )
        self.qc_reason_label.setText(self._fit_reason)
        self.stats_label.setText(
            f"Preview fit {label}: R² {r2:.6f} | mean {mean_abs:.2f} bp | max {max_abs:.2f} bp"
        )
        self.stats_label.setStyleSheet(f"color: {color_map.get(self._fit_grade, '#64748b')}; font-weight: 700;")

    def _plot_residuals(self):
        self.residual_ax.clear()
        xs = []
        ys = []
        colors = []
        for row in self._fit_rows:
            if row["residual"] is None:
                continue
            xs.append(row["expected_bp"])
            ys.append(row["residual"])
            if abs(row["residual"]) <= PASS_MAX_ABS_RESIDUAL:
                colors.append("#16a34a")
            elif abs(row["residual"]) <= CHECK_MAX_ABS_RESIDUAL:
                colors.append("#d97706")
            else:
                colors.append("#dc2626")

        self.residual_ax.axhline(0.0, color="#94a3b8", linestyle="--", linewidth=1.0)
        if xs:
            self.residual_ax.scatter(xs, ys, c=colors, s=42, zorder=3)
            self.residual_ax.plot(xs, ys, color="#cbd5e1", linewidth=1.0, zorder=2)
            self.residual_ax.set_ylabel("Residual (bp)")
        else:
            self.residual_ax.text(
                0.5,
                0.5,
                "Residuals will appear after a valid fit preview.",
                transform=self.residual_ax.transAxes,
                ha="center",
                va="center",
                color="#64748b",
            )
        self.residual_ax.set_xlabel("Expected ladder step (bp)")
        self.residual_ax.grid(True, alpha=0.2)
        self.residual_figure.tight_layout()
        self.residual_canvas.draw_idle()

    def _plot_ladder(self):
        self.ax.clear()
        trace = self.fsa.size_standard
        self.ax.plot(trace, color="#8fa6c1", alpha=0.95, linewidth=1.35, label="Size Standard")

        selected_candidate = self._selected_candidate_row()
        selected_step = self._selected_step_row()

        if not self.candidates.empty:
            times = self.candidates["time"].to_numpy(dtype=float)
            intensities = self.candidates["intensity"].to_numpy(dtype=float)
            manual_mask = self.candidates["source"].astype(str).eq("manual").to_numpy(dtype=bool) if "source" in self.candidates.columns else np.zeros(len(self.candidates), dtype=bool)
            auto_mask = ~manual_mask
            if np.any(auto_mask):
                self.ax.scatter(
                    times[auto_mask],
                    intensities[auto_mask],
                    marker="x",
                    color="#ef4444",
                    s=48,
                    linewidths=1.4,
                    label="Candidates",
                )
            if np.any(manual_mask):
                self.ax.scatter(
                    times[manual_mask],
                    intensities[manual_mask],
                    marker="D",
                    color="#0f766e",
                    s=46,
                    linewidths=1.0,
                    edgecolors="#0f766e",
                    label="Manual peaks",
                )

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
            offset_y = 10 if step_idx % 2 == 0 else 22
            self.ax.annotate(
                f"{bp:.0f}",
                (peak_time, peak_intensity),
                textcoords="offset points",
                xytext=(0, offset_y),
                ha="center",
                fontsize=8,
                color=marker_color,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec=marker_color, lw=0.8, alpha=0.92),
            )

        if not self.candidates.empty:
            x_min = max(float(self.candidates["time"].min()) - 180.0, 0.0)
            x_max = min(float(self.candidates["time"].max()) + 180.0, float(len(trace)))
            y_max = max(float(self.candidates["intensity"].max()) * 1.28, float(np.max(trace)) * 0.95, 1.0)
            self.ax.set_xlim(x_min, x_max)
            self.ax.set_ylim(min(-150.0, float(np.min(trace)) * 1.05), y_max)

        self.ax.set_title(f"Ladder Trace · {self.fsa.ladder}", fontsize=14, fontweight="bold")
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel("Intensity")
        self.ax.grid(True, alpha=0.25)
        self.ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#dbe4ef")
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _selected_step_row(self) -> int | None:
        selected_rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        return selected_rows[0].row() if selected_rows else None

    def _selected_candidate_row(self) -> int | None:
        selected_rows = self.candidate_table.selectionModel().selectedRows() if self.candidate_table.selectionModel() else []
        return selected_rows[0].row() if selected_rows else None

    def _sync_selection_from_match_table(self):
        step = self._selected_step_row()
        if step is not None and step in self.mapping:
            self.candidate_table.selectRow(self.mapping[step])
        self._plot_ladder()

    def _sync_selection_from_candidate_table(self):
        cand_idx = self._selected_candidate_row()
        if cand_idx is not None:
            step_idx = self._candidate_used_by(cand_idx)
            if step_idx is not None:
                self.table.selectRow(step_idx)
        self._plot_ladder()

    def _toggle_add_peak_mode(self, checked: bool):
        self._add_peak_mode = checked
        if checked:
            step_idx = self._selected_step_row()
            if step_idx is None or step_idx in self.mapping:
                next_missing = self._next_missing_step(step_idx)
                if next_missing is not None:
                    self.table.selectRow(next_missing)
                    step_idx = next_missing
            if step_idx is None:
                self.stats_label.setText("Add-missing mode: all ladder steps are already assigned.")
                self.stats_label.setStyleSheet("color: #64748b; font-weight: 700;")
                return
            direction = "high → low" if self._missing_order == "descending" else "low → high"
            self.stats_label.setText(
                f"Add-missing mode ({direction}): click the trace to place {self.ladder_steps[step_idx]:.0f} bp at the local maximum."
            )
            self.stats_label.setStyleSheet("color: #0f766e; font-weight: 700;")
        else:
            self._update_qc_panel()

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
        self._refresh_preview_state(show_errors=False)
        self._refresh_all()

        if self._add_peak_mode:
            next_missing = self._next_missing_step(step_idx)
            if next_missing is not None:
                self.table.selectRow(next_missing)
                self.stats_label.setText(
                    f"Added {self.ladder_steps[step_idx]:.0f} bp. Click to place the next missing step: {self.ladder_steps[next_missing]:.0f} bp."
                )
                self.stats_label.setStyleSheet("color: #0f766e; font-weight: 700;")
            else:
                self.btn_add_peak.setChecked(False)
                self.stats_label.setText("All ladder steps are now assigned. Review the fit and save if it looks good.")
                self.stats_label.setStyleSheet("color: #16a34a; font-weight: 700;")
        elif step_idx + 1 < len(self.ladder_steps):
            self.table.selectRow(step_idx + 1)

    def _clear_selected_step(self):
        step_idx = self._selected_step_row()
        if step_idx is None:
            QMessageBox.information(self, "No Step Selected", "Select a ladder step to clear first.")
            return
        if step_idx in self.mapping:
            del self.mapping[step_idx]
            self._refresh_preview_state(show_errors=False)
            self._refresh_all()
            self.table.selectRow(step_idx)

    def _clear_all(self):
        self.mapping = {}
        self._refresh_preview_state(show_errors=False)
        self._refresh_all()
        if self.table.rowCount():
            self.table.selectRow(0)

    def _reset_to_initial(self):
        self.mapping = dict(self._initial_mapping)
        self._refresh_preview_state(show_errors=False)
        self._refresh_all()
        if self.table.rowCount():
            self.table.selectRow(0)

    def _on_plot_click(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return

        step_idx = self._selected_step_row()
        if self._add_peak_mode:
            if step_idx is None:
                QMessageBox.information(self, "No Step Selected", "Select a ladder step first, then add the missing peak from the plot.")
                return
            self._add_manual_peak_from_plot(float(event.xdata), assign_to_step=step_idx)
            self.btn_add_peak.setChecked(False)
            return
        if self.candidates.empty:
            return
        if step_idx is None:
            QMessageBox.information(self, "No Step Selected", "Select a ladder step first, then click a candidate peak.")
            return

        x = float(event.xdata)
        peak_time, _intensity = self._find_local_peak_time(x)
        existing_idx = self._candidate_time_exists(peak_time, tolerance=2.0)
        if existing_idx is not None:
            self._assign_candidate_to_step(step_idx, existing_idx)
            return

        self._add_manual_peak_from_plot(x, assign_to_step=step_idx)

    def _on_step_double_clicked(self, row, _column):
        if row in self.mapping:
            del self.mapping[row]
            self._refresh_preview_state(show_errors=False)
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
            fitted_steps = np.asarray(getattr(self.fsa, "ladder_steps", self.ladder_steps), dtype=float)
            for fitted_idx, peak_time in enumerate(best):
                if peak_time > 0 and peak_time in ss_peaks:
                    matches = np.where(np.isclose(self.ladder_steps, fitted_steps[fitted_idx], atol=1e-6))[0]
                    if matches.size == 0:
                        continue
                    step_idx = int(matches[0])
                    auto_mapping[step_idx] = ss_peaks.index(peak_time)

        self.mapping = auto_mapping
        if not store_initial:
            self._manual_candidate_times = []
            self.candidates = self._get_candidates().reset_index(drop=True)
        if store_initial:
            self._initial_mapping = dict(auto_mapping)
        self._missing_order = self._recommended_missing_order()
        self._sync_missing_order_button()
        self._refresh_preview_state(show_errors=False)

    def _preview_fit(self):
        if len(self.mapping) < 3:
            QMessageBox.warning(self, "Invalid Fit", "Select at least 3 peaks to preview fit.")
            return
        self._refresh_preview_state(show_errors=True)
        self._refresh_all()

    def _on_apply(self):
        if not self.mapping:
            QMessageBox.warning(self, "No Mapping", "Map at least one ladder step before applying.")
            return
        missing_steps = self._missing_step_indices()
        if missing_steps:
            missing_text = ", ".join(f"{self.ladder_steps[idx]:.0f} bp" for idx in missing_steps[:8])
            if len(missing_steps) > 8:
                missing_text += ", ..."
            QMessageBox.warning(
                self,
                "Incomplete Ladder Mapping",
                "All expected ladder steps must be assigned before saving this adjustment.\n\n"
                f"Missing: {missing_text}",
            )
            return

        self._refresh_preview_state(show_errors=True)
        self._refresh_all()
        if self._preview_metrics is None:
            QMessageBox.warning(
                self,
                "Preview Required",
                "This ladder correction could not be previewed successfully yet. Fix the fit before saving.",
            )
            return
        self.accept()

    def get_mapping(self):
        return dict(self.mapping)

    def get_adjustment_payload(self):
        return self._build_adjustment_payload()
