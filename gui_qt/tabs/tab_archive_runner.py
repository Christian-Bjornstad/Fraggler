from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from PyQt6.QtCore import Qt, QThreadPool
from PyQt6.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from config import APP_SETTINGS, get_analysis_settings, save_settings
from gui_qt.worker import Worker
from scripts.combine_clonality_yearly_overview import combine_run_root
from scripts.run_clonality_yearly import discover_month_folders, normalize_month_keys, run_yearly_validation


def _open_path(path: Path) -> None:
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    elif sys.platform == "win32":
        subprocess.Popen(["explorer", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


class TabArchiveRunner(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.threadpool = QThreadPool.globalInstance()
        self._workflow_state = "ready"
        self._current_analysis_id = APP_SETTINGS.get("active_analysis", "clonality")
        self._active_worker: Worker | None = None
        self._current_run_root: Path | None = None
        self._current_manifest_path: Path | None = None
        self._current_workbook_path: Path | None = None
        self._month_checkboxes: dict[str, QCheckBox] = {}
        self.month_checkboxes = self._month_checkboxes
        self._month_row_map: dict[str, int] = {}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(18)

        header = QVBoxLayout()
        title = QLabel("Archive Runner")
        title.setObjectName("PageTitle")
        subtitle = QLabel(
            "Run year-scale clonality backfills with safe fresh output folders, explicit resume support, and an optional combined yearly workbook."
        )
        subtitle.setObjectName("PageSubtitle")
        subtitle.setWordWrap(True)
        header.addWidget(title)
        header.addWidget(subtitle)
        layout.addLayout(header)

        layout.addWidget(self._build_settings_card())
        layout.addWidget(self._build_month_card())
        layout.addWidget(self._build_dashboard_card())
        layout.addWidget(self._build_output_card())

        self.refresh_from_settings()
        self.set_analysis(self._current_analysis_id)

    def _build_settings_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        form = QFormLayout(card)
        form.addRow(QLabel("<b>Run Settings</b>"))

        self.year_input = QLineEdit()
        self.year_input.setPlaceholderText("2025")
        form.addRow("Year:", self.year_input)

        input_row = QHBoxLayout()
        self.input_root = QLineEdit()
        self.input_root.setPlaceholderText("/path/to/Klonalitet/2025_data")
        btn_input = QPushButton("Browse...")
        btn_input.clicked.connect(lambda: self._browse_directory(self.input_root))
        input_row.addWidget(self.input_root, stretch=1)
        input_row.addWidget(btn_input)
        form.addRow("Input Root:", input_row)

        output_row = QHBoxLayout()
        self.output_root = QLineEdit()
        self.output_root.setPlaceholderText("/path/to/output/full_year_runs")
        btn_output = QPushButton("Browse...")
        btn_output.clicked.connect(lambda: self._browse_directory(self.output_root))
        output_row.addWidget(self.output_root, stretch=1)
        output_row.addWidget(btn_output)
        form.addRow("Output Root:", output_row)

        self.run_name = QLineEdit()
        self.run_name.setPlaceholderText("Optional run folder name")
        form.addRow("Run Name:", self.run_name)

        self.max_workers = QSpinBox()
        self.max_workers.setRange(1, 64)
        form.addRow("Max Workers:", self.max_workers)

        self.folder_workers = QSpinBox()
        self.folder_workers.setRange(1, 64)
        form.addRow("Folder Workers:", self.folder_workers)

        self.chk_resume = QCheckBox("Resume existing run folder")
        self.chk_include_sl = QCheckBox("Include SL in exported artifacts")
        self.chk_refresh_each_folder = QCheckBox("Refresh workbook after each folder")
        self.chk_cleanup_staging = QCheckBox("Delete month staging roots after completion")
        form.addRow("", self.chk_resume)
        form.addRow("", self.chk_include_sl)
        form.addRow("", self.chk_refresh_each_folder)
        form.addRow("", self.chk_cleanup_staging)
        return card

    def _build_month_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        layout = QVBoxLayout(card)

        title = QLabel("MONTHS")
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        quick_row = QHBoxLayout()
        btn_all = QPushButton("Select All")
        btn_all.clicked.connect(lambda: self._set_all_months(True))
        btn_none = QPushButton("Clear")
        btn_none.clicked.connect(lambda: self._set_all_months(False))
        quick_row.addWidget(btn_all)
        quick_row.addWidget(btn_none)
        quick_row.addStretch()
        layout.addLayout(quick_row)

        grid = QGridLayout()
        for idx in range(12):
            month = f"{idx + 1:02d}"
            checkbox = QCheckBox(month)
            checkbox.setChecked(True)
            checkbox.toggled.connect(self._rebuild_month_table)
            self._month_checkboxes[month] = checkbox
            grid.addWidget(checkbox, idx // 6, idx % 6)
        layout.addLayout(grid)

        note = QLabel("Month selection is per-session only. New runs default to all months selected.")
        note.setObjectName("MutedText")
        note.setWordWrap(True)
        layout.addWidget(note)
        return card

    def _build_dashboard_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("DashboardCard")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(14)

        header_row = QHBoxLayout()
        title = QLabel("Yearly Workflow")
        title.setObjectName("DashboardTitle")
        header_row.addWidget(title)
        header_row.addStretch()

        self.status_badge = QLabel("READY")
        self.status_badge.setObjectName("WorkflowStatusBadge")
        self.status_badge.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.status_badge.setMinimumWidth(140)
        header_row.addWidget(self.status_badge)
        layout.addLayout(header_row)

        action_row = QHBoxLayout()
        self.btn_run = QPushButton("Run Yearly Backfill")
        self.btn_run.setObjectName("PrimaryButton")
        self.btn_run.clicked.connect(self.on_run_yearly)
        self.btn_combine = QPushButton("Build Combined Workbook")
        self.btn_combine.clicked.connect(self.on_build_combined_workbook)
        self.btn_open_run = QPushButton("Open Run Folder")
        self.btn_open_run.clicked.connect(self.on_open_run_folder)
        self.btn_open_workbook = QPushButton("Open Combined Workbook")
        self.btn_open_workbook.clicked.connect(self.on_open_combined_workbook)
        for button in (self.btn_run, self.btn_combine, self.btn_open_run, self.btn_open_workbook):
            action_row.addWidget(button)
        action_row.addStretch()
        layout.addLayout(action_row)

        self.status_lbl = QLabel("Ready")
        self.status_lbl.setWordWrap(True)
        layout.addWidget(self.status_lbl)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self.month_table = QTableWidget(0, 4)
        self.month_table.setHorizontalHeaderLabels(["Month", "Status", "Folders", "Run Dir"])
        self.month_table.verticalHeader().setVisible(False)
        self.month_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.month_table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self.month_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.month_table)
        return card

    def _build_output_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        form = QFormLayout(card)
        form.addRow(QLabel("<b>Current Outputs</b>"))

        run_root_row = QHBoxLayout()
        self.selected_run_root = QLineEdit()
        self.selected_run_root.setReadOnly(True)
        btn_select_run_root = QPushButton("Choose Run Folder...")
        btn_select_run_root.clicked.connect(self.on_choose_run_root)
        run_root_row.addWidget(self.selected_run_root, stretch=1)
        run_root_row.addWidget(btn_select_run_root)
        form.addRow("Selected Run:", run_root_row)

        self.run_root_lbl = QLabel("—")
        self.run_root_lbl.setWordWrap(True)
        form.addRow("Run Root:", self.run_root_lbl)

        self.manifest_lbl = QLabel("—")
        self.manifest_lbl.setWordWrap(True)
        form.addRow("Manifest:", self.manifest_lbl)

        self.workbook_lbl = QLabel("—")
        self.workbook_lbl.setWordWrap(True)
        form.addRow("Combined Workbook:", self.workbook_lbl)
        return card

    def _profile(self) -> dict:
        return get_analysis_settings("clonality")

    def set_analysis(self, analysis_id: str) -> None:
        self._current_analysis_id = analysis_id
        self.setEnabled(analysis_id == "clonality")

    def refresh_from_settings(self) -> None:
        archive = self._profile().get("archive_runner", {})
        self.year_input.setText(str(archive.get("year_label", "2025")))
        self.input_root.setText(str(archive.get("input_root", "")))
        self.output_root.setText(str(archive.get("output_root", "")))
        self.run_name.setText(str(archive.get("run_name", "")))
        self.max_workers.setValue(int(archive.get("max_workers", 1) or 1))
        self.folder_workers.setValue(int(archive.get("folder_workers", 1) or 1))
        # Fresh output folders stay the default even if the user resumed a prior run earlier.
        self.chk_resume.setChecked(False)
        self.chk_include_sl.setChecked(bool(archive.get("include_sl", False)))
        self.chk_refresh_each_folder.setChecked(bool(archive.get("refresh_each_folder", False)))
        self.chk_cleanup_staging.setChecked(bool(archive.get("cleanup_staging_root", False)))

        last_selected_run_root = str(
            archive.get("last_selected_run_root", "") or archive.get("last_run_root", "") or ""
        ).strip()
        self._current_run_root = Path(last_selected_run_root).expanduser() if last_selected_run_root else None
        self._current_manifest_path = self._guess_manifest_path()
        self._current_workbook_path = self._guess_workbook_path()
        self._set_all_months(True)
        self._refresh_output_labels()
        self._refresh_action_buttons()

    def _persist_settings(self) -> None:
        archive = APP_SETTINGS.setdefault("analyses", {}).setdefault("clonality", {}).setdefault("archive_runner", {})
        archive.update(self._collect_settings())
        archive["last_selected_run_root"] = str(self._current_run_root or "")
        archive["last_run_root"] = str(self._current_run_root or "")
        archive["combined_workbook_path"] = str(self._current_workbook_path or "")
        save_settings(APP_SETTINGS)

    def save_defaults(self) -> None:
        self._persist_settings()

    def _collect_settings(self) -> dict[str, object]:
        return {
            "input_root": self.input_root.text().strip(),
            "output_root": self.output_root.text().strip(),
            "year_label": self.year_input.text().strip(),
            "run_name": self.run_name.text().strip(),
            "max_workers": self.max_workers.value(),
            "folder_workers": self.folder_workers.value(),
            "last_run_root": str(self._current_run_root or ""),
            "last_selected_run_root": str(self._current_run_root or ""),
            "combined_workbook_path": str(self._current_workbook_path or ""),
            "resume_existing": self.chk_resume.isChecked(),
            "include_sl": self.chk_include_sl.isChecked(),
            "refresh_each_folder": self.chk_refresh_each_folder.isChecked(),
            "cleanup_staging_root": self.chk_cleanup_staging.isChecked(),
        }

    def _browse_directory(self, target: QLineEdit) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Select Directory", target.text() or str(Path.home()))
        if folder:
            target.setText(folder)

    def _selected_months(self) -> list[str]:
        year_label = self.year_input.text().strip()
        raw_months = [f"{year_label}_{month}" for month, checkbox in self._month_checkboxes.items() if checkbox.isChecked()]
        return normalize_month_keys(year_label, raw_months)

    def _selected_month_keys(self) -> list[str]:
        return self._selected_months()

    def _set_all_months(self, checked: bool) -> None:
        for checkbox in self._month_checkboxes.values():
            checkbox.setChecked(checked)
        self._rebuild_month_table()

    def _set_workflow_status(self, message: str, state: str) -> None:
        self._workflow_state = state
        self.status_lbl.setText(message)
        self.status_badge.setText(state.replace("_", " ").upper())
        self.status_badge.setProperty("state", state)
        self.status_lbl.setProperty("state", state)
        self._restyle_widget(self.status_badge)
        self._restyle_widget(self.status_lbl)

    def _restyle_widget(self, widget: QWidget) -> None:
        widget.style().unpolish(widget)
        widget.style().polish(widget)
        widget.update()

    def _refresh_output_labels(self) -> None:
        self.selected_run_root.setText(str(self._current_run_root) if self._current_run_root else "")
        self.run_root_lbl.setText(str(self._current_run_root) if self._current_run_root else "—")
        self.manifest_lbl.setText(str(self._current_manifest_path) if self._current_manifest_path else "—")
        self.workbook_lbl.setText(str(self._current_workbook_path) if self._current_workbook_path else "—")

    def _refresh_action_buttons(self) -> None:
        has_run_root = self._current_run_root is not None and self._current_run_root.exists()
        has_workbook = self._current_workbook_path is not None and self._current_workbook_path.exists()
        self.btn_open_run.setEnabled(has_run_root)
        self.btn_combine.setEnabled(has_run_root)
        self.btn_open_workbook.setEnabled(has_workbook)

    def _set_busy(self, busy: bool) -> None:
        self.btn_run.setEnabled(not busy)
        self.btn_combine.setEnabled(not busy and self._current_run_root is not None and self._current_run_root.exists())
        self.btn_open_run.setEnabled(not busy and self._current_run_root is not None and self._current_run_root.exists())
        self.btn_open_workbook.setEnabled(not busy and self._current_workbook_path is not None and self._current_workbook_path.exists())

    def _guess_manifest_path(self) -> Path | None:
        if self._current_run_root is None:
            return None
        year_label = self.year_input.text().strip()
        path = self._current_run_root / f"full_{year_label}_run_manifest.json"
        return path if path.exists() else None

    def _guess_workbook_path(self) -> Path | None:
        if self._current_run_root is None:
            return None
        return self._current_run_root / f"track-clonality-{self.year_input.text().strip()}-overview.xlsx"

    def _month_counts(self) -> dict[str, int]:
        year_label = self.year_input.text().strip()
        input_root = self.input_root.text().strip()
        if len(year_label) != 4 or not year_label.isdigit() or not input_root:
            return {}
        try:
            month_map = discover_month_folders(Path(input_root), year_label)
        except Exception:
            return {}
        return {month: len(paths) for month, paths in month_map.items()}

    def _rebuild_month_table(self) -> None:
        counts = self._month_counts()
        selected_months = self._selected_months() if self.year_input.text().strip().isdigit() else []
        old_state: dict[str, tuple[str, str, str]] = {}
        for month, row in self._month_row_map.items():
            old_state[month] = (
                self.month_table.item(row, 1).text() if self.month_table.item(row, 1) else "pending",
                self.month_table.item(row, 2).text() if self.month_table.item(row, 2) else "",
                self.month_table.item(row, 3).text() if self.month_table.item(row, 3) else "",
            )

        self.month_table.setRowCount(0)
        self._month_row_map = {}
        for row, month_key in enumerate(selected_months):
            self.month_table.insertRow(row)
            self._month_row_map[month_key] = row
            status, folders, run_dir = old_state.get(month_key, ("pending", str(counts.get(month_key, "")), ""))
            self.month_table.setItem(row, 0, QTableWidgetItem(month_key))
            self.month_table.setItem(row, 1, QTableWidgetItem(status))
            self.month_table.setItem(row, 2, QTableWidgetItem(folders))
            self.month_table.setItem(row, 3, QTableWidgetItem(run_dir))
        self.month_table.resizeColumnsToContents()

    def _update_month_row(self, month_key: str, *, status: str, folder_count: int | None = None, run_dir: str = "") -> None:
        row = self._month_row_map.get(month_key)
        if row is None:
            return
        self.month_table.setItem(row, 1, QTableWidgetItem(status))
        if folder_count is not None:
            self.month_table.setItem(row, 2, QTableWidgetItem(str(folder_count)))
        if run_dir:
            self.month_table.setItem(row, 3, QTableWidgetItem(run_dir))

    def _validated_inputs(self) -> tuple[str, Path, Path, list[str]]:
        year_label = self.year_input.text().strip()
        if len(year_label) != 4 or not year_label.isdigit():
            raise ValueError("Year must be four digits, for example 2025.")
        input_root = Path(self.input_root.text().strip()).expanduser()
        if not input_root.is_dir():
            raise FileNotFoundError(f"Input root not found: {input_root}")
        output_root = Path(self.output_root.text().strip()).expanduser()
        output_root.mkdir(parents=True, exist_ok=True)
        months = self._selected_months()
        if not months:
            raise ValueError("Select at least one month.")
        return year_label, input_root, output_root, months

    def on_run_yearly(self) -> None:
        try:
            year_label, input_root, output_root, months = self._validated_inputs()
        except Exception as exc:
            QMessageBox.warning(self, "Archive Runner", str(exc))
            return

        self._current_run_root = None
        self._current_manifest_path = None
        self._current_workbook_path = None
        self._refresh_output_labels()
        self._rebuild_month_table()
        self._persist_settings()
        self._set_busy(True)
        self.progress.setRange(0, max(len(months), 1))
        self.progress.setValue(0)
        self._set_workflow_status(f"Starting yearly run for {year_label}", "running")

        worker = Worker(
            run_yearly_validation,
            year_label=year_label,
            input_root=input_root,
            output_root=output_root,
            run_name=self.run_name.text().strip() or None,
            months=months,
            max_workers=self.max_workers.value(),
            folder_workers=self.folder_workers.value(),
            refresh_each_folder=self.chk_refresh_each_folder.isChecked(),
            include_sl=self.chk_include_sl.isChecked(),
            cleanup_staging_root=self.chk_cleanup_staging.isChecked(),
            resume_existing=self.chk_resume.isChecked(),
        )
        worker.kwargs["progress_callback"] = worker.signals.event.emit
        worker.kwargs["status_callback"] = worker.signals.status.emit
        worker.signals.event.connect(self._on_runner_event)
        worker.signals.status.connect(self._on_runner_status)
        worker.signals.result.connect(self._on_runner_finished)
        worker.signals.error.connect(self._on_worker_error)
        worker.signals.finished.connect(self._on_worker_finished)
        self._active_worker = worker
        self.threadpool.start(worker)

    def _run_yearly_job(
        self,
        *,
        year_label: str,
        input_root: Path,
        output_root: Path,
        run_name: str | None,
        months: list[str],
        max_workers: int,
        folder_workers: int,
        resume_existing: bool,
        include_sl: bool,
        refresh_each_folder: bool,
        cleanup_staging_root: bool,
        bridge,
    ) -> dict[str, object]:
        return run_yearly_validation(
            year_label=year_label,
            input_root=input_root,
            output_root=output_root,
            run_name=run_name,
            months=months,
            max_workers=max_workers,
            folder_workers=folder_workers,
            resume_existing=resume_existing,
            include_sl=include_sl,
            refresh_each_folder=refresh_each_folder,
            cleanup_staging_root=cleanup_staging_root,
            progress_callback=lambda payload: bridge.progress.emit(payload),
            status_callback=lambda message: bridge.status.emit(message),
        )

    def on_build_combined_workbook(self) -> None:
        run_root = self._current_run_root
        if run_root is None or not run_root.exists():
            QMessageBox.warning(self, "Archive Runner", "No run root is available yet.")
            return
        year_label = self.year_input.text().strip()
        self._set_busy(True)
        self._set_workflow_status(f"Building combined workbook for {year_label}", "running")

        worker = Worker(combine_run_root, run_root, run_root / f"track-clonality-{year_label}-overview.xlsx", year_label=year_label)
        worker.signals.result.connect(self._on_combine_finished)
        worker.signals.error.connect(self._on_worker_error)
        worker.signals.finished.connect(self._on_worker_finished)
        self._active_worker = worker
        self.threadpool.start(worker)

    def on_open_run_folder(self) -> None:
        if self._current_run_root and self._current_run_root.exists():
            _open_path(self._current_run_root)

    def on_open_combined_workbook(self) -> None:
        if self._current_workbook_path and self._current_workbook_path.exists():
            _open_path(self._current_workbook_path)

    def on_choose_run_root(self) -> None:
        start_dir = str(self._current_run_root) if self._current_run_root else (self.output_root.text().strip() or str(Path.home()))
        folder = QFileDialog.getExistingDirectory(self, "Select Existing Run Folder", start_dir)
        if not folder:
            return
        self._current_run_root = Path(folder).expanduser()
        self._current_manifest_path = self._guess_manifest_path()
        self._current_workbook_path = self._guess_workbook_path()
        self._persist_settings()
        self._refresh_output_labels()
        self._refresh_action_buttons()

    def _on_runner_event(self, payload: object) -> None:
        if not isinstance(payload, dict):
            return

        event = str(payload.get("event") or "")
        month = str(payload.get("month") or "")
        run_dir = str(payload.get("run_dir") or "")
        folder_count = payload.get("folder_count")

        if event == "run_started" and run_dir:
            self._current_run_root = Path(run_dir)
        elif event == "month_started" and month:
            self._update_month_row(month, status="running", folder_count=int(folder_count or 0), run_dir=run_dir)
            self.progress.setValue(min(self.progress.value() + 1, self.progress.maximum()))
        elif event == "month_resumed" and month:
            self._update_month_row(month, status="resumed", run_dir=run_dir)
        elif event == "month_skipped_empty" and month:
            self._update_month_row(month, status="skipped_empty", folder_count=0, run_dir=run_dir)
        elif event == "month_finished" and month:
            self._update_month_row(month, status="done", run_dir=run_dir)
        elif event == "manifest_written":
            manifest_path = payload.get("manifest_path")
            if manifest_path:
                self._current_manifest_path = Path(str(manifest_path))
        elif event == "run_finished":
            manifest_path = payload.get("manifest_path")
            if run_dir:
                self._current_run_root = Path(run_dir)
            if manifest_path:
                self._current_manifest_path = Path(str(manifest_path))

        self._current_workbook_path = self._guess_workbook_path()
        self._refresh_output_labels()
        self._refresh_action_buttons()

    def _on_runner_status(self, message: str) -> None:
        self._set_workflow_status(message, "running")

    def _on_runner_finished(self, manifest: object) -> None:
        if isinstance(manifest, dict):
            run_dir = manifest.get("run_dir")
            if run_dir:
                self._current_run_root = Path(str(run_dir))
            self._current_manifest_path = self._guess_manifest_path()
            self._current_workbook_path = self._guess_workbook_path()
            self._persist_settings()
        self.progress.setValue(self.progress.maximum())
        self._set_workflow_status("Yearly backfill finished.", "success")
        self._refresh_output_labels()
        self._refresh_action_buttons()

    def _on_combine_finished(self, workbook_path: object) -> None:
        if workbook_path:
            self._current_workbook_path = Path(str(workbook_path))
            self._persist_settings()
        self._set_workflow_status("Combined workbook created.", "success")
        self._refresh_output_labels()
        self._refresh_action_buttons()

    def _on_worker_error(self, err_tuple) -> None:
        message = str(err_tuple[1]) if isinstance(err_tuple, tuple) and len(err_tuple) > 1 else "Archive Runner failed."
        self._set_workflow_status(message, "error")
        QMessageBox.critical(self, "Archive Runner", message)

    def _on_worker_finished(self) -> None:
        self._active_worker = None
        self._set_busy(False)
