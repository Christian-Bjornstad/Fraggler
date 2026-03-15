from __future__ import annotations

from pathlib import Path
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
from PyQt6.QtCore import Qt

from config import APP_SETTINGS
from core.analysis import save_ladder_adjustment
from core.html_reports import extract_dit_from_name
from gui_qt.dialogs.ladder_dialog import LadderAdjustmentDialog
from gui_qt.ladder_utils import detect_fsa_for_ladder, load_adjustable_fsa


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
        self._all_files: list[Path] = []
        self._current_file: Path | None = None
        self._current_meta: dict | None = None
        self._report_matches: list[Path] = []

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

        main_layout.addWidget(self._build_source_card())
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
        btn_scan = QPushButton("Scan .fsa Files")
        btn_scan.clicked.connect(self._scan_files)
        btn_browse_file = QPushButton("Open Single File...")
        btn_browse_file.clicked.connect(self._choose_single_file)
        row1.addWidget(QLabel("Input Folder:"))
        row1.addWidget(self.source_dir, stretch=1)
        row1.addWidget(btn_browse_dir)
        row1.addWidget(btn_scan)
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
            ("analysis", "Analysis"),
            ("assay", "Assay"),
            ("group", "Group"),
            ("ladder", "Ladder"),
            ("sample_channel", "Sample Channel"),
            ("trace_channels", "Trace Channels"),
            ("bp_range", "bp Range"),
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
        btn_find = QPushButton("Find Matching Reports")
        btn_find.clicked.connect(self._refresh_report_matches)
        row1.addWidget(QLabel("Report Root:"))
        row1.addWidget(self.report_root, stretch=1)
        row1.addWidget(btn_browse)
        row1.addWidget(btn_find)
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
        default_source = Path("data/euroclonality")
        if default_source.exists():
            self.source_dir.setText(str(default_source))
            self._scan_files()

        example = default_source / "00011_87685fba_24OUM10035_SL__060824_A01_C990GXRS.fsa"
        if example.exists():
            self._select_file(example)

        report_root = Path("final")
        if report_root.exists():
            self.report_root.setText(str(report_root))

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

        self._all_files = sorted(source.rglob("*.fsa"), key=lambda p: p.name.lower())
        self._rebuild_file_list()
        self._set_status(f"Found {len(self._all_files)} .fsa files in {source}.")

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

        self._refresh_current_metadata()
        self._refresh_report_matches()

    def _refresh_current_metadata(self) -> None:
        if not self._current_file:
            return

        meta = detect_fsa_for_ladder(
            self._current_file,
            preferred_analysis=APP_SETTINGS.get("active_analysis"),
        )
        self._current_meta = meta

        self.detail_labels["file"].setText(str(self._current_file))
        if not meta:
            for key in [
                "analysis",
                "assay",
                "group",
                "ladder",
                "sample_channel",
                "trace_channels",
                "bp_range",
                "adjustment",
            ]:
                self.detail_labels[key].setText("Could not classify")
            self._set_status(f"Could not classify {self._current_file.name}.", error=True)
            return

        adj_status = "Saved" if self._current_file.with_suffix(".ladder_adj.json").exists() else "None"
        self.detail_labels["analysis"].setText(meta["analysis"])
        self.detail_labels["assay"].setText(meta["assay"])
        self.detail_labels["group"].setText(meta["group"])
        self.detail_labels["ladder"].setText(meta["ladder"])
        self.detail_labels["sample_channel"].setText(meta.get("sample_channel") or "—")
        self.detail_labels["trace_channels"].setText(", ".join(meta.get("trace_channels", [])) or "—")
        self.detail_labels["bp_range"].setText(f"{meta['bp_min']:.1f} - {meta['bp_max']:.1f}")
        self.detail_labels["adjustment"].setText(adj_status)
        self._set_status(f"Loaded metadata for {self._current_file.name}.")

    def _clear_details(self) -> None:
        for label in self.detail_labels.values():
            label.setText("—")

    def _open_ladder_editor(self) -> None:
        if not self._current_file:
            return

        try:
            fsa, meta = load_adjustable_fsa(
                self._current_file,
                preferred_analysis=APP_SETTINGS.get("active_analysis"),
            )
        except Exception as exc:
            QMessageBox.critical(self, "Load Failed", f"Could not load file for ladder adjustment:\n{exc}")
            self._set_status(f"Could not load {self._current_file.name}.", error=True)
            return

        dialog = LadderAdjustmentDialog(fsa, self)
        if dialog.exec():
            mapping = dialog.get_mapping()
            save_ladder_adjustment(fsa, mapping)
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
        self.report_list.clear()
        self._report_matches = []

        root_text = self.report_root.text().strip()
        if not self._current_file or not root_text:
            self._update_report_buttons()
            return

        root = Path(root_text).expanduser()
        if not root.exists():
            self._set_status("Report root does not exist.", error=True)
            self._update_report_buttons()
            return

        dit = extract_dit_from_name(self._current_file.name)
        stem = self._current_file.stem.lower()
        tokens = [token for token in [dit, stem] if token]

        html_matches = []
        for path in root.rglob("*.html"):
            lower = path.name.lower()
            if any(token.lower() in lower for token in tokens):
                html_matches.append(path)

        self._report_matches = sorted(set(html_matches))
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
