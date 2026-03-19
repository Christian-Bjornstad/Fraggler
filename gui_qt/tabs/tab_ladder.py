from __future__ import annotations

from pathlib import Path
import copy
import subprocess
import sys

from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QListWidget,
    QListWidgetItem,
    QGridLayout,
    QMessageBox,
    QAbstractItemView,
)
from PyQt6.QtCore import Qt, QThreadPool

from config import APP_SETTINGS, get_analysis_settings
from core.analysis import load_ladder_adjustment, save_ladder_adjustment
from core.html_reports import extract_dit_from_name
from gui_qt.dialogs.ladder_dialog import LadderAdjustmentDialog
from gui_qt.ladder_utils import detect_fsa_for_ladder, load_adjustable_fsa
from gui_qt.worker import Worker


def _open_path(path: Path) -> None:
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    elif sys.platform == "win32":
        subprocess.Popen(["explorer", str(path)])
    else:
        subprocess.Popen(["xdg-open", str(path)])


class TabLadder(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.threadpool = QThreadPool.globalInstance()
        self._all_files: list[Path] = []
        self._current_file: Path | None = None
        self._current_meta: dict | None = None
        self._current_fsa = None
        self._report_matches: list[Path] = []
        self._current_analysis_id = APP_SETTINGS.get("active_analysis", "clonality")
        self._scan_request_id = 0
        self._metadata_request_id = 0
        self._report_request_id = 0
        self._metadata_loading = False

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(16)

        header = QVBoxLayout()
        title = QLabel("Ladder Studio")
        title.setObjectName("PageTitle")
        sub = QLabel("Pick one .fsa file, inspect its ladder metadata, and open a focused ladder-adjustment workflow.")
        sub.setObjectName("PageSubtitle")
        header.addWidget(title)
        header.addWidget(sub)
        main_layout.addLayout(header)

        main_layout.addWidget(self._build_source_card(), stretch=1)
        main_layout.addWidget(self._build_details_card())
        main_layout.addWidget(self._build_report_card(), stretch=1)

        self.status_lbl = QLabel("Ready — scan a folder or browse directly to a single .fsa file.")
        self.status_lbl.setStyleSheet("color: #64748b; font-weight: 500;")
        main_layout.addWidget(self.status_lbl)

        self._load_defaults()

    def _build_source_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        layout = QVBoxLayout(card)

        title = QLabel("SOURCE FILES")
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        row1 = QHBoxLayout()
        self.source_dir = QLineEdit()
        self.source_dir.setPlaceholderText("/path/to/folder with .fsa files")
        btn_browse_dir = QPushButton("Browse Folder...")
        btn_browse_dir.clicked.connect(self._choose_source_dir)
        self.btn_scan = QPushButton("Scan .fsa Files")
        self.btn_scan.clicked.connect(self._scan_files)
        btn_browse_file = QPushButton("Open Single File...")
        btn_browse_file.clicked.connect(self._choose_single_file)
        row1.addWidget(QLabel("Input Folder:"))
        row1.addWidget(self.source_dir, stretch=1)
        row1.addWidget(btn_browse_dir)
        row1.addWidget(self.btn_scan)
        row1.addWidget(btn_browse_file)
        layout.addLayout(row1)

        row2 = QHBoxLayout()
        self.file_filter = QLineEdit()
        self.file_filter.setPlaceholderText("Filter by filename, DIT, assay, plate position...")
        self.file_filter.textChanged.connect(self._rebuild_file_list)
        row2.addWidget(QLabel("Filter:"))
        row2.addWidget(self.file_filter, stretch=1)
        layout.addLayout(row2)

        self.file_list = QListWidget()
        self.file_list.setMinimumHeight(220)
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.file_list.itemSelectionChanged.connect(self._on_file_selected)
        self.file_list.itemDoubleClicked.connect(lambda _: self._open_ladder_editor())
        layout.addWidget(self.file_list)

        return card

    def _build_details_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        layout = QVBoxLayout(card)

        title = QLabel("SELECTED FILE")
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        details = QGridLayout()
        details.setHorizontalSpacing(18)
        details.setVerticalSpacing(10)
        self.detail_labels: dict[str, QLabel] = {}

        fields = [
            ("file", "File"),
            ("assay", "Assay"),
            ("ladder", "Ladder"),
            ("fit_strategy", "Fit Strategy"),
            ("fit_counts", "Expected / Fitted"),
            ("review_state", "Review State"),
            ("missing_steps", "Missing Steps"),
            ("adjustment", "Saved Adjustment"),
        ]

        for row, (key, label) in enumerate(fields):
            lbl_key = QLabel(f"{label}:")
            lbl_key.setStyleSheet("color: #64748b; font-weight: 700;")
            lbl_val = QLabel("—")
            lbl_val.setWordWrap(True)
            self.detail_labels[key] = lbl_val
            details.addWidget(lbl_key, row, 0, alignment=Qt.AlignmentFlag.AlignTop)
            details.addWidget(lbl_val, row, 1)

        layout.addLayout(details)

        actions = QHBoxLayout()
        self.btn_refresh_meta = QPushButton("Refresh Metadata")
        self.btn_refresh_meta.clicked.connect(self._refresh_current_metadata)
        self.btn_open_editor = QPushButton("Open Ladder Editor")
        self.btn_open_editor.setObjectName("PrimaryButton")
        self.btn_open_editor.clicked.connect(self._open_ladder_editor)
        self.btn_remove_adjustment = QPushButton("Remove Saved Adjustment")
        self.btn_remove_adjustment.clicked.connect(self._remove_saved_adjustment)
        self.btn_open_file_folder = QPushButton("Open File Folder")
        self.btn_open_file_folder.clicked.connect(self._open_file_folder)

        for btn in [
            self.btn_refresh_meta,
            self.btn_open_editor,
            self.btn_remove_adjustment,
            self.btn_open_file_folder,
        ]:
            btn.setEnabled(False)
            actions.addWidget(btn)
        actions.addStretch()
        layout.addLayout(actions)

        hint = QLabel(
            "Tip: double-click a file to jump straight into the ladder editor. "
            "Inside the editor you can re-map peaks, preview the fit, and save the adjustment for re-runs."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #64748b;")
        layout.addWidget(hint)
        return card

    def _build_report_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        layout = QVBoxLayout(card)

        title = QLabel("MATCHING REPORTS")
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        row1 = QHBoxLayout()
        self.report_root = QLineEdit()
        self.report_root.setPlaceholderText("/optional/path/to/report root")
        btn_browse = QPushButton("Browse Reports...")
        btn_browse.clicked.connect(self._choose_report_root)
        self.btn_find_reports = QPushButton("Find Matching Reports")
        self.btn_find_reports.clicked.connect(self._refresh_report_matches)
        row1.addWidget(QLabel("Report Root:"))
        row1.addWidget(self.report_root, stretch=1)
        row1.addWidget(btn_browse)
        row1.addWidget(self.btn_find_reports)
        layout.addLayout(row1)

        self.report_list = QListWidget()
        self.report_list.itemDoubleClicked.connect(self._open_selected_report)
        layout.addWidget(self.report_list)

        row2 = QHBoxLayout()
        self.btn_open_report = QPushButton("Open Selected Report")
        self.btn_open_report.clicked.connect(self._open_selected_report)
        self.btn_open_report_folder = QPushButton("Open Report Folder")
        self.btn_open_report_folder.clicked.connect(self._open_selected_report_folder)
        self.btn_open_report.setEnabled(False)
        self.btn_open_report_folder.setEnabled(False)
        row2.addWidget(self.btn_open_report)
        row2.addWidget(self.btn_open_report_folder)
        row2.addStretch()
        layout.addLayout(row2)

        self.report_list.itemSelectionChanged.connect(self._update_report_buttons)
        return card

    def _load_defaults(self) -> None:
        profile = get_analysis_settings(self._current_analysis_id)
        input_dir = profile.get("batch", {}).get("base_input_dir", "")
        output_dir = profile.get("batch", {}).get("output_base", "")

        if input_dir:
            self.source_dir.setText(input_dir)
        elif self._current_analysis_id == "clonality":
            default_source = Path("data/euroclonality")
            if default_source.exists():
                self.source_dir.setText(str(default_source))

        if output_dir:
            self.report_root.setText(output_dir)
        elif Path("final").exists():
            self.report_root.setText("final")

    def set_analysis(self, analysis_id: str) -> None:
        previous_profile = get_analysis_settings(self._current_analysis_id)
        next_profile = get_analysis_settings(analysis_id)

        previous_input = previous_profile.get("batch", {}).get("base_input_dir", "")
        previous_output = previous_profile.get("batch", {}).get("output_base", "")

        if not self.source_dir.text().strip() or self.source_dir.text().strip() == previous_input:
            self.source_dir.setText(next_profile.get("batch", {}).get("base_input_dir", ""))
        if not self.report_root.text().strip() or self.report_root.text().strip() == previous_output:
            self.report_root.setText(next_profile.get("batch", {}).get("output_base", ""))

        self._current_analysis_id = analysis_id
        self._current_meta = None
        self._current_fsa = None
        self._clear_details()
        if self._current_file:
            self._set_status(
                f"Analysis switched to {analysis_id}. Refresh metadata to re-evaluate the current file."
            )

    def _choose_source_dir(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Folder With .fsa Files",
            self.source_dir.text() or str(Path.home()),
        )
        if folder:
            self.source_dir.setText(folder)

    def _choose_report_root(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Report Root",
            self.report_root.text() or str(Path.home()),
        )
        if folder:
            self.report_root.setText(folder)
            self._refresh_report_matches()

    def _choose_single_file(self) -> None:
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            "Open .fsa File",
            self.source_dir.text() or str(Path.home()),
            "FSA files (*.fsa)",
        )
        if file_name:
            file_path = Path(file_name)
            if file_path.parent.exists():
                self.source_dir.setText(str(file_path.parent))
            if file_path not in self._all_files:
                self._all_files.append(file_path)
                self._all_files.sort(key=lambda p: p.name.lower())
            self._rebuild_file_list()
            self._select_file(file_path)

    def _scan_files(self) -> None:
        source = Path(self.source_dir.text().strip()).expanduser()
        if not source.exists() or not source.is_dir():
            self._set_status("Input folder does not exist.", error=True)
            return

        self._scan_request_id += 1
        request_id = self._scan_request_id
        self.btn_scan.setEnabled(False)
        self._set_status(f"Scanning {source} for .fsa files...")

        worker = Worker(self._scan_fsa_files_worker, source)
        worker.signals.result.connect(lambda files, rid=request_id, src=source: self._on_scan_result(rid, src, files))
        worker.signals.error.connect(lambda err, rid=request_id: self._on_scan_error(rid, err))
        self.threadpool.start(worker)

    def _rebuild_file_list(self) -> None:
        active_path = self._current_file
        text = self.file_filter.text().strip().lower()
        self.file_list.clear()

        matches = []
        for path in self._all_files:
            haystack = str(path).lower()
            if text and text not in haystack:
                continue
            item = QListWidgetItem(path.name)
            item.setData(Qt.ItemDataRole.UserRole, str(path))
            self.file_list.addItem(item)
            matches.append(path)

        if matches and active_path in matches:
            self._select_file(active_path)
        elif matches and self.file_list.currentRow() < 0:
            self.file_list.setCurrentRow(0)
        else:
            self._update_current_file(None)

    def _select_file(self, file_path: Path) -> None:
        file_str = str(file_path)
        for idx in range(self.file_list.count()):
            item = self.file_list.item(idx)
            if item.data(Qt.ItemDataRole.UserRole) == file_str:
                self.file_list.setCurrentItem(item)
                return

    def _on_file_selected(self) -> None:
        items = self.file_list.selectedItems()
        if not items:
            self._update_current_file(None)
            return
        self._update_current_file(Path(items[0].data(Qt.ItemDataRole.UserRole)))

    def _update_current_file(self, file_path: Path | None) -> None:
        self._current_file = file_path
        self._current_meta = None
        self._current_fsa = None
        self._clear_details()

        enabled = file_path is not None
        for btn in [
            self.btn_refresh_meta,
            self.btn_open_editor,
            self.btn_remove_adjustment,
            self.btn_open_file_folder,
        ]:
            btn.setEnabled(enabled)

        if not file_path:
            self.report_list.clear()
            self._report_matches = []
            self._update_report_buttons()
            return

        self.detail_labels["file"].setText(str(file_path))
        self._refresh_current_metadata()
        self._refresh_report_matches()

    def _refresh_current_metadata(self) -> None:
        if not self._current_file:
            return

        self.detail_labels["file"].setText(str(self._current_file))
        self.detail_labels["assay"].setText("Loading...")
        self.detail_labels["ladder"].setText("Loading...")
        self.detail_labels["adjustment"].setText("Loading...")
        self.detail_labels["fit_strategy"].setText("Loading...")
        self.detail_labels["fit_counts"].setText("—")
        self.detail_labels["review_state"].setText("—")
        self.detail_labels["missing_steps"].setText("—")
        self._start_metadata_load(self._current_file)

    def _clear_details(self) -> None:
        for label in self.detail_labels.values():
            label.setText("—")

    def _open_ladder_editor(self) -> None:
        if not self._current_file:
            return

        if self._metadata_loading:
            QMessageBox.information(self, "Metadata Loading", "Metadata is still loading for the selected file.")
            return
        if self._current_meta is None or self._current_fsa is None:
            self._refresh_current_metadata()
            QMessageBox.information(self, "Metadata Required", "Refresh metadata before opening the ladder editor.")
            return

        fsa = copy.deepcopy(self._current_fsa)
        dialog = LadderAdjustmentDialog(fsa, self)
        if dialog.exec():
            adjustment = dialog.get_adjustment_payload()
            save_ladder_adjustment(fsa, adjustment)
            self._refresh_current_metadata()
            self._set_status(
                f"Saved ladder adjustment for {self._current_file.name}. Re-run the analysis to use the new fit."
            )
            QMessageBox.information(
                self,
                "Adjustment Saved",
                f"Ladder adjustment saved for {self._current_file.name}.\n\n"
                "Re-run the relevant report job to apply it in the outputs.",
            )
        else:
            self._set_status(f"Closed ladder editor for {self._current_file.name}.")

    def _remove_saved_adjustment(self) -> None:
        if not self._current_file:
            return

        adj_path = self._current_file.with_suffix(".ladder_adj.json")
        if not adj_path.exists():
            QMessageBox.information(self, "No Adjustment", "There is no saved ladder adjustment for this file.")
            return

        reply = QMessageBox.question(
            self,
            "Remove Adjustment",
            f"Delete the saved ladder adjustment for {self._current_file.name}?",
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        adj_path.unlink(missing_ok=True)
        self._refresh_current_metadata()
        self._set_status(f"Removed saved ladder adjustment for {self._current_file.name}.")

    def _open_file_folder(self) -> None:
        if self._current_file:
            _open_path(self._current_file.parent)

    def _refresh_report_matches(self) -> None:
        root_text = self.report_root.text().strip()
        if not self._current_file or not root_text:
            self.report_list.clear()
            self._report_matches = []
            self._update_report_buttons()
            return

        self._report_request_id += 1
        request_id = self._report_request_id
        self.btn_find_reports.setEnabled(False)
        self.report_list.clear()
        self._report_matches = []
        self._update_report_buttons()

        worker = Worker(self._find_report_matches_worker, self._current_file, root_text)
        worker.signals.result.connect(lambda result, rid=request_id: self._on_report_matches_result(rid, result))
        worker.signals.error.connect(lambda err, rid=request_id: self._on_report_matches_error(rid, err))
        self.threadpool.start(worker)

    def _update_report_buttons(self) -> None:
        has_selection = bool(self.report_list.selectedItems())
        self.btn_open_report.setEnabled(has_selection)
        self.btn_open_report_folder.setEnabled(has_selection)

    def _open_selected_report(self) -> None:
        items = self.report_list.selectedItems()
        if not items:
            return
        _open_path(Path(items[0].data(Qt.ItemDataRole.UserRole)))

    def _open_selected_report_folder(self) -> None:
        items = self.report_list.selectedItems()
        if not items:
            return
        _open_path(Path(items[0].data(Qt.ItemDataRole.UserRole)).parent)

    def _set_status(self, text: str, error: bool = False) -> None:
        color = "#ef4444" if error else "#64748b"
        self.status_lbl.setText(text)
        self.status_lbl.setStyleSheet(f"color: {color}; font-weight: 500;")

    @staticmethod
    def _adjustment_status_for(file_path: Path) -> str:
        payload = load_ladder_adjustment(type("Dummy", (), {"file": file_path})())
        return "Saved" if payload else "None"

    def _start_metadata_load(self, file_path: Path) -> None:
        self._metadata_request_id += 1
        request_id = self._metadata_request_id
        analysis_id = APP_SETTINGS.get("active_analysis")
        self._metadata_loading = True
        self.btn_open_editor.setEnabled(False)
        self._set_status(f"Loading ladder metadata for {file_path.name}...")

        worker = Worker(self._load_metadata_worker, file_path, analysis_id)
        worker.signals.result.connect(lambda result, rid=request_id: self._on_metadata_result(rid, result))
        worker.signals.error.connect(lambda err, rid=request_id: self._on_metadata_error(rid, err))
        self.threadpool.start(worker)

    @staticmethod
    def _scan_fsa_files_worker(source: Path) -> list[Path]:
        return sorted(source.rglob("*.fsa"), key=lambda p: p.name.lower())

    @staticmethod
    def _load_metadata_worker(file_path: Path, analysis_id: str | None) -> dict:
        meta = detect_fsa_for_ladder(file_path, preferred_analysis=analysis_id)
        if not meta:
            return {"file_path": file_path, "meta": None, "fsa": None}
        fsa, refreshed_meta = load_adjustable_fsa(file_path, preferred_analysis=analysis_id, metadata=meta)
        return {"file_path": file_path, "meta": refreshed_meta, "fsa": fsa}

    @staticmethod
    def _find_report_matches_worker(file_path: Path, root_text: str) -> dict:
        root = Path(root_text).expanduser()
        if not root.exists():
            raise FileNotFoundError("Report root does not exist.")

        dit = extract_dit_from_name(file_path.name)
        stem = file_path.stem.lower()
        tokens = [token for token in [dit, stem] if token]

        html_matches = []
        for path in root.rglob("*.html"):
            lower = path.name.lower()
            if any(token.lower() in lower for token in tokens):
                html_matches.append(path)

        return {"root": root, "matches": sorted(set(html_matches))}

    def _on_scan_result(self, request_id: int, source: Path, files: list[Path]) -> None:
        if request_id != self._scan_request_id:
            return
        self.btn_scan.setEnabled(True)
        self._all_files = files
        self._rebuild_file_list()
        self._set_status(f"Found {len(self._all_files)} .fsa files in {source}.")

    def _on_scan_error(self, request_id: int, err_tuple) -> None:
        if request_id != self._scan_request_id:
            return
        self.btn_scan.setEnabled(True)
        self._set_status(f"Could not scan .fsa files: {err_tuple[1]}", error=True)

    def _on_metadata_result(self, request_id: int, result: dict) -> None:
        if request_id != self._metadata_request_id:
            return
        self._metadata_loading = False

        file_path = result["file_path"]
        if file_path != self._current_file:
            return

        meta = result["meta"]
        if not meta:
            for key in [
                "assay",
                "ladder",
                "fit_strategy",
                "fit_counts",
                "review_state",
                "missing_steps",
                "adjustment",
            ]:
                self.detail_labels[key].setText("Could not classify")
            self._set_status(f"Could not classify {file_path.name}.", error=True)
            return

        self._current_meta = meta
        self._current_fsa = result["fsa"]
        adj_status = self._adjustment_status_for(file_path)
        assay_label = meta.get("assay") or meta.get("analysis", "").capitalize() or "—"
        self.detail_labels["assay"].setText(assay_label)
        self.detail_labels["ladder"].setText(meta["ladder"])
        self.detail_labels["adjustment"].setText(adj_status)

        fsa = self._current_fsa
        fit_strategy = str(getattr(fsa, "ladder_fit_strategy", "auto_full")).replace("_", " ")
        expected_steps = list(map(float, getattr(fsa, "expected_ladder_steps", getattr(fsa, "ladder_steps", []))))
        fitted_steps = list(map(float, getattr(fsa, "ladder_steps", [])))
        missing_steps = list(map(float, getattr(fsa, "ladder_missing_expected_steps", [])))
        fit_note = str(getattr(fsa, "ladder_fit_note", ""))
        review_required = bool(getattr(fsa, "ladder_review_required", bool(missing_steps)))

        if getattr(fsa, "ladder_fit_strategy", "") == "manual_adjustment":
            review_state = "Manual correction active"
        elif review_required:
            review_state = "Usable but incomplete"
        else:
            review_state = "Full fit"

        self.detail_labels["fit_strategy"].setText(fit_strategy)
        self.detail_labels["fit_counts"].setText(f"{len(expected_steps)} / {len(fitted_steps)}")
        self.detail_labels["review_state"].setText(review_state)
        self.detail_labels["missing_steps"].setText(
            ", ".join(f"{bp:.0f}" for bp in missing_steps) if missing_steps else "None"
        )
        self.btn_open_editor.setEnabled(True)
        self._set_status(fit_note or f"Loaded metadata for {file_path.name}.")

    def _on_metadata_error(self, request_id: int, err_tuple) -> None:
        if request_id != self._metadata_request_id:
            return
        self._metadata_loading = False
        self.btn_open_editor.setEnabled(self._current_file is not None)
        self.detail_labels["fit_strategy"].setText("Could not load")
        self.detail_labels["fit_counts"].setText("—")
        self.detail_labels["review_state"].setText("Unknown")
        self.detail_labels["missing_steps"].setText("—")
        self._set_status(f"Loaded metadata, but not ladder state: {err_tuple[1]}", error=True)

    def _on_report_matches_result(self, request_id: int, result: dict) -> None:
        if request_id != self._report_request_id:
            return
        self.btn_find_reports.setEnabled(True)
        root = result["root"]
        self._report_matches = result["matches"]
        self.report_list.clear()
        for path in self._report_matches:
            item = QListWidgetItem(str(path.relative_to(root)))
            item.setData(Qt.ItemDataRole.UserRole, str(path))
            self.report_list.addItem(item)

        if self._report_matches:
            self.report_list.setCurrentRow(0)
            self._set_status(f"Found {len(self._report_matches)} matching reports.")
        else:
            self._set_status("No matching reports found under the selected report root.")
        self._update_report_buttons()

    def _on_report_matches_error(self, request_id: int, err_tuple) -> None:
        if request_id != self._report_request_id:
            return
        self.btn_find_reports.setEnabled(True)
        self.report_list.clear()
        self._report_matches = []
        self._update_report_buttons()
        self._set_status(str(err_tuple[1]), error=True)
