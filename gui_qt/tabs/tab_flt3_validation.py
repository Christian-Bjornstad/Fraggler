from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from PyQt6.QtCore import Qt, QThreadPool
from PyQt6.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QProgressBar,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from config import APP_SETTINGS, get_analysis_settings, save_settings
from gui_qt.worker import Worker
from scripts.run_flt3_backfill_validation import (
    DEFAULT_EXCLUDED_BASENAMES,
    run_backfill_validation,
)


def _open_path(path: Path) -> None:
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    elif sys.platform == "win32":
        subprocess.Popen(["explorer", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


class TabFlt3Validation(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.threadpool = QThreadPool.globalInstance()
        self._current_analysis_id = APP_SETTINGS.get("active_analysis", "clonality")
        self._active_worker: Worker | None = None
        self._current_run_dir: Path | None = None
        self._current_workbook_path: Path | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(18)

        header = QVBoxLayout()
        title = QLabel("Archive Runner")
        title.setObjectName("PageTitle")
        subtitle = QLabel(
            "Run a 3730DNA-filtered FLT3 backfill validation with residual exports so you can validate ladder-fit quality across the FLT3 archive."
        )
        subtitle.setObjectName("PageSubtitle")
        subtitle.setWordWrap(True)
        header.addWidget(title)
        header.addWidget(subtitle)
        layout.addLayout(header)

        layout.addWidget(self._build_settings_card())
        layout.addWidget(self._build_dashboard_card())
        layout.addWidget(self._build_output_card())

        self.refresh_from_settings()
        self.set_analysis(self._current_analysis_id)

    def _build_settings_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        form = QFormLayout(card)
        form.addRow(QLabel("<b>Archive Runner Settings</b>"))

        data_row = QHBoxLayout()
        self.data_root = QLineEdit()
        self.data_root.setPlaceholderText("/Volumes/T7 Shield/DATA/flt3")
        btn_data = QPushButton("Browse...")
        btn_data.clicked.connect(lambda: self._browse_directory(self.data_root))
        data_row.addWidget(self.data_root, stretch=1)
        data_row.addWidget(btn_data)
        form.addRow("Data Root:", data_row)

        output_row = QHBoxLayout()
        self.output_root = QLineEdit()
        self.output_root.setPlaceholderText("/path/to/validation_outputs")
        btn_output = QPushButton("Browse...")
        btn_output.clicked.connect(lambda: self._browse_directory(self.output_root))
        output_row.addWidget(self.output_root, stretch=1)
        output_row.addWidget(btn_output)
        form.addRow("Output Root:", output_row)

        self.run_name = QLineEdit()
        self.run_name.setPlaceholderText("Optional run folder name")
        form.addRow("Run Name:", self.run_name)

        self.require_run_name_contains = QLineEdit()
        self.require_run_name_contains.setPlaceholderText("3730DNA")
        form.addRow("Run Name Filter:", self.require_run_name_contains)

        self.years = QLineEdit()
        self.years.setPlaceholderText("2025,2026")
        form.addRow("Years:", self.years)

        self.workers = QSpinBox()
        self.workers.setRange(1, 64)
        form.addRow("Workers:", self.workers)

        self.limit = QSpinBox()
        self.limit.setRange(0, 500000)
        self.limit.setSingleStep(100)
        form.addRow("Limit (0 = all):", self.limit)

        self.timeout_seconds = QSpinBox()
        self.timeout_seconds.setRange(1, 3600)
        self.timeout_seconds.setSingleStep(5)
        form.addRow("Timeout / file (s):", self.timeout_seconds)

        self.checkpoint_every = QSpinBox()
        self.checkpoint_every.setRange(0, 50000)
        self.checkpoint_every.setSingleStep(25)
        form.addRow("Checkpoint Every:", self.checkpoint_every)

        self.exclude_basenames = QPlainTextEdit()
        self.exclude_basenames.setPlaceholderText("One basename per line for known human/machine-error files")
        self.exclude_basenames.setMinimumHeight(90)
        form.addRow("Excluded Files:", self.exclude_basenames)

        self.chk_include_npm1 = QCheckBox("Include NPM1")
        self.chk_dit_only = QCheckBox("DIT only")
        form.addRow("", self.chk_include_npm1)
        form.addRow("", self.chk_dit_only)
        return card

    def _build_dashboard_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("DashboardCard")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(14)

        header_row = QHBoxLayout()
        title = QLabel("FLT3 Archive Workflow")
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
        self.btn_run = QPushButton("Run FLT3 Archive Runner")
        self.btn_run.setObjectName("PrimaryButton")
        self.btn_run.clicked.connect(self.on_run_validation)
        self.btn_open_run = QPushButton("Open Run Folder")
        self.btn_open_run.clicked.connect(self.on_open_run_folder)
        self.btn_open_workbook = QPushButton("Open Workbook")
        self.btn_open_workbook.clicked.connect(self.on_open_workbook)
        for button in (self.btn_run, self.btn_open_run, self.btn_open_workbook):
            action_row.addWidget(button)
        action_row.addStretch()
        layout.addLayout(action_row)

        self.status_lbl = QLabel("Ready")
        self.status_lbl.setWordWrap(True)
        layout.addWidget(self.status_lbl)

        self.remaining_lbl = QLabel("Remaining: —")
        self.remaining_lbl.setWordWrap(True)
        layout.addWidget(self.remaining_lbl)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        self.command_preview = QPlainTextEdit()
        self.command_preview.setReadOnly(True)
        self.command_preview.setMaximumBlockCount(200)
        self.command_preview.setMinimumHeight(120)
        layout.addWidget(self.command_preview)

        self.summary_preview = QPlainTextEdit()
        self.summary_preview.setReadOnly(True)
        self.summary_preview.setMaximumBlockCount(2000)
        self.summary_preview.setMinimumHeight(220)
        layout.addWidget(self.summary_preview)
        return card

    def _build_output_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        form = QFormLayout(card)
        form.addRow(QLabel("<b>Current Outputs</b>"))

        self.run_root_lbl = QLabel("—")
        self.run_root_lbl.setWordWrap(True)
        form.addRow("Run Root:", self.run_root_lbl)

        self.summary_json_lbl = QLabel("—")
        self.summary_json_lbl.setWordWrap(True)
        form.addRow("Summary JSON:", self.summary_json_lbl)

        self.residual_json_lbl = QLabel("—")
        self.residual_json_lbl.setWordWrap(True)
        form.addRow("Residual JSON:", self.residual_json_lbl)

        self.workbook_lbl = QLabel("—")
        self.workbook_lbl.setWordWrap(True)
        form.addRow("Workbook:", self.workbook_lbl)
        return card

    def _profile(self) -> dict:
        return get_analysis_settings("flt3")

    def set_analysis(self, analysis_id: str) -> None:
        self._current_analysis_id = analysis_id
        self.setEnabled(analysis_id == "flt3")

    def _browse_directory(self, target: QLineEdit) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Select Directory", target.text() or str(Path.home()))
        if folder:
            target.setText(folder)
            self._persist_settings()
            self._refresh_command_preview()

    def refresh_from_settings(self) -> None:
        validation = self._profile().get("validation", {})
        self.data_root.setText(str(validation.get("data_root", "")))
        self.output_root.setText(str(validation.get("output_root", "")))
        self.run_name.setText(str(validation.get("run_name", "")))
        self.require_run_name_contains.setText(str(validation.get("require_run_name_contains", "3730DNA")))
        self.years.setText(",".join(str(year) for year in validation.get("years", ["2025", "2026"])))
        self.workers.setValue(int(validation.get("workers", 6) or 6))
        self.limit.setValue(int(validation.get("limit", 0) or 0))
        self.timeout_seconds.setValue(int(validation.get("timeout_seconds", 45) or 45))
        self.checkpoint_every.setValue(int(validation.get("checkpoint_every", 100) or 100))
        excluded = validation.get("excluded_basenames", DEFAULT_EXCLUDED_BASENAMES)
        self.exclude_basenames.setPlainText("\n".join(excluded))
        self.chk_include_npm1.setChecked(bool(validation.get("include_npm1", False)))
        self.chk_dit_only.setChecked(bool(validation.get("dit_only", False)))
        last_run_dir = str(validation.get("last_run_dir", "") or "").strip()
        self._current_run_dir = Path(last_run_dir).expanduser() if last_run_dir else None
        last_workbook_path = str(validation.get("last_workbook_path", "") or "").strip()
        self._current_workbook_path = Path(last_workbook_path).expanduser() if last_workbook_path else None
        self._refresh_output_labels()
        self._refresh_action_buttons()
        self._refresh_command_preview()

    def _collect_settings(self) -> dict[str, object]:
        years = [part.strip() for part in self.years.text().split(",") if part.strip()]
        return {
            "data_root": self.data_root.text().strip(),
            "output_root": self.output_root.text().strip(),
            "run_name": self.run_name.text().strip(),
            "require_run_name_contains": self.require_run_name_contains.text().strip() or "3730DNA",
            "years": years or ["2025", "2026"],
            "workers": self.workers.value(),
            "limit": self.limit.value(),
            "timeout_seconds": self.timeout_seconds.value(),
            "checkpoint_every": self.checkpoint_every.value(),
            "excluded_basenames": [
                line.strip()
                for line in self.exclude_basenames.toPlainText().splitlines()
                if line.strip()
            ],
            "include_npm1": self.chk_include_npm1.isChecked(),
            "dit_only": self.chk_dit_only.isChecked(),
            "last_run_dir": str(self._current_run_dir or ""),
            "last_workbook_path": str(self._current_workbook_path or ""),
        }

    def _persist_settings(self) -> None:
        validation = APP_SETTINGS.setdefault("analyses", {}).setdefault("flt3", {}).setdefault("validation", {})
        validation.update(self._collect_settings())
        save_settings(APP_SETTINGS)

    def save_defaults(self) -> None:
        self._persist_settings()

    def _selected_years(self) -> list[str]:
        years = [part.strip() for part in self.years.text().split(",") if part.strip()]
        return years or ["2025", "2026"]

    def _refresh_command_preview(self) -> None:
        lines = [
            "python3 scripts/run_flt3_backfill_validation.py",
            f"  --data-root '{self.data_root.text().strip() or '/Volumes/T7 Shield/DATA/flt3'}'",
            f"  --output-root '{self.output_root.text().strip() or str(Path.home())}'",
            f"  --require-run-name-contains '{self.require_run_name_contains.text().strip() or '3730DNA'}'",
            f"  --workers {self.workers.value()}",
            f"  --timeout-seconds {self.timeout_seconds.value()}",
            f"  --checkpoint-every {self.checkpoint_every.value()}",
        ]
        if self.run_name.text().strip():
            lines.append(f"  --run-name '{self.run_name.text().strip()}'")
        for year in self._selected_years():
            lines.append(f"  --year {year}")
        if self.limit.value() > 0:
            lines.append(f"  --limit {self.limit.value()}")
        if self.chk_include_npm1.isChecked():
            lines.append("  --include-npm1")
        if self.chk_dit_only.isChecked():
            lines.append("  --dit-only")
        for basename in [
            line.strip() for line in self.exclude_basenames.toPlainText().splitlines() if line.strip()
        ]:
            lines.append(f"  --exclude-basename '{basename}'")
        self.command_preview.setPlainText(" \\\n".join(lines))

    def _set_workflow_status(self, message: str, state: str) -> None:
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
        self.run_root_lbl.setText(str(self._current_run_dir) if self._current_run_dir else "—")
        summary_json = (self._current_run_dir / "summary.json") if self._current_run_dir else None
        residual_json = (self._current_run_dir / "residual_summary.json") if self._current_run_dir else None
        self.summary_json_lbl.setText(str(summary_json) if summary_json and summary_json.exists() else "—")
        self.residual_json_lbl.setText(str(residual_json) if residual_json and residual_json.exists() else "—")
        self.workbook_lbl.setText(
            str(self._current_workbook_path)
            if self._current_workbook_path and self._current_workbook_path.exists()
            else "—"
        )

        summary_preview = ""
        if self._current_run_dir:
            summary_md = self._current_run_dir / "summary.md"
            residual_md = self._current_run_dir / "residual_summary.md"
            if summary_md.exists():
                summary_preview += summary_md.read_text(encoding="utf-8", errors="replace")
            if residual_md.exists():
                if summary_preview:
                    summary_preview += "\n\n---\n\n"
                summary_preview += residual_md.read_text(encoding="utf-8", errors="replace")
        self.summary_preview.setPlainText(summary_preview or "No FLT3 validation run loaded yet.")

    def _refresh_action_buttons(self) -> None:
        self.btn_open_run.setEnabled(self._current_run_dir is not None and self._current_run_dir.exists())
        self.btn_open_workbook.setEnabled(
            self._current_workbook_path is not None and self._current_workbook_path.exists()
        )

    def _set_busy(self, busy: bool) -> None:
        self.btn_run.setEnabled(not busy)
        self.btn_open_run.setEnabled(not busy and self._current_run_dir is not None and self._current_run_dir.exists())
        self.progress.setVisible(busy)

    def _validated_inputs(self) -> tuple[Path, Path, list[str]]:
        data_root = Path(self.data_root.text().strip()).expanduser()
        if not data_root.is_dir():
            raise FileNotFoundError(f"Data root not found: {data_root}")
        output_root = Path(self.output_root.text().strip()).expanduser()
        output_root.mkdir(parents=True, exist_ok=True)
        years = self._selected_years()
        return data_root, output_root, years

    def on_run_validation(self) -> None:
        try:
            data_root, output_root, years = self._validated_inputs()
        except Exception as exc:
            QMessageBox.warning(self, "FLT3 Archive Runner", str(exc))
            return

        self._persist_settings()
        self._set_busy(True)
        self.progress.setRange(0, 0)
        self.progress.setValue(0)
        self._set_workflow_status("Running FLT3 3730 archive backfill", "running")
        self.remaining_lbl.setText("Remaining: calculating...")
        self.summary_preview.setPlainText("Validation running...")

        worker = Worker(
            run_backfill_validation,
            data_root=data_root,
            output_root=output_root,
            run_name=self.run_name.text().strip() or None,
            years=years,
            workers=self.workers.value(),
            limit=self.limit.value(),
            input_manifest=None,
            include_npm1=self.chk_include_npm1.isChecked(),
            dit_only=self.chk_dit_only.isChecked(),
            timeout_seconds=self.timeout_seconds.value(),
            checkpoint_every=self.checkpoint_every.value(),
            required_run_name_contains=self.require_run_name_contains.text().strip() or "3730DNA",
            excluded_basenames=[
                line.strip() for line in self.exclude_basenames.toPlainText().splitlines() if line.strip()
            ],
        )
        worker.kwargs["progress_callback"] = worker.signals.progress.emit
        worker.kwargs["progress_max_callback"] = worker.signals.progress_max.emit
        worker.kwargs["status_callback"] = worker.signals.status.emit
        worker.signals.progress.connect(self._on_progress)
        worker.signals.progress_max.connect(self._on_progress_max)
        worker.signals.status.connect(self._on_status)
        worker.signals.result.connect(self._on_run_finished)
        worker.signals.error.connect(self._on_worker_error)
        worker.signals.finished.connect(self._on_worker_finished)
        self._active_worker = worker
        self.threadpool.start(worker)

    def on_open_run_folder(self) -> None:
        if self._current_run_dir and self._current_run_dir.exists():
            _open_path(self._current_run_dir)

    def on_open_workbook(self) -> None:
        if self._current_workbook_path and self._current_workbook_path.exists():
            _open_path(self._current_workbook_path)

    def _on_run_finished(self, payload: object) -> None:
        if isinstance(payload, dict):
            run_dir = payload.get("run_dir")
            if run_dir:
                self._current_run_dir = Path(str(run_dir))
            workbook_path = payload.get("workbook_path")
            if workbook_path:
                self._current_workbook_path = Path(str(workbook_path))
                self._persist_settings()
                self._refresh_output_labels()
                self._refresh_action_buttons()
            status_counts = ((payload.get("validator_summary") or {}).get("status_counts") or {})
            self._set_workflow_status(f"FLT3 archive run finished. Status counts: {status_counts}", "success")
            self.remaining_lbl.setText("Remaining: 0")
            self.summary_preview.setPlainText(json.dumps(payload, indent=2, sort_keys=True))

    def _on_progress_max(self, maximum: int) -> None:
        maximum = max(int(maximum or 0), 0)
        if maximum > 0:
            self.progress.setRange(0, maximum)
            self.progress.setValue(0)
            self.remaining_lbl.setText(f"Remaining: {maximum}")
        else:
            self.progress.setRange(0, 0)
            self.remaining_lbl.setText("Remaining: calculating...")

    def _on_progress(self, value: int) -> None:
        value = max(int(value or 0), 0)
        if self.progress.maximum() > 0:
            self.progress.setValue(min(value, self.progress.maximum()))
            remaining = max(self.progress.maximum() - value, 0)
            self.remaining_lbl.setText(f"Remaining: {remaining} of {self.progress.maximum()}")

    def _on_status(self, message: str) -> None:
        self._set_workflow_status(message, "running")

    def _on_worker_error(self, err_tuple) -> None:
        message = str(err_tuple[1]) if isinstance(err_tuple, tuple) and len(err_tuple) > 1 else "FLT3 archive runner failed."
        self._set_workflow_status(message, "error")
        QMessageBox.critical(self, "FLT3 Archive Runner", message)

    def _on_worker_finished(self) -> None:
        self._active_worker = None
        self._set_busy(False)
