from pathlib import Path
import subprocess
import sys
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QTableWidget, QTableWidgetItem, QProgressBar,
    QHeaderView, QAbstractItemView, QFileDialog, QListWidget, QSizePolicy, QScroller,
    QCheckBox, QComboBox, QFrame, QGridLayout, QLayout,
)
from PyQt6.QtCore import Qt, QThreadPool, QEvent, QSignalBlocker, QPoint, QRect, QSize
from gui_qt.worker import Worker
from config import (
    APP_SETTINGS,
    get_analysis_settings,
    save_settings,
    GENERAL_DEFAULT_LADDER,
    GENERAL_DEFAULT_TRACE_CHANNELS,
    GENERAL_DEFAULT_PRIMARY_CHANNEL,
)


ANALYSIS_LABELS = {
    "clonality": "Klonalitet",
    "flt3": "FLT3 Analysis",
    "general": "General",
}

GENERAL_LADDER_OPTIONS = [
    ("LIZ500", "LIZ500_250"),
    ("ROX400HD", "ROX400HD"),
    ("GS500ROX", "GS500ROX"),
]
GENERAL_TRACE_OPTIONS = [
    ("DATA1", "Blue trace"),
    ("DATA2", "Green trace"),
    ("DATA3", "Yellow / Black trace"),
]


class FlowLayout(QLayout):
    """Simple wrapping layout for compact option cards."""

    def __init__(self, parent=None, margin: int = 0, h_spacing: int = 10, v_spacing: int = 10):
        super().__init__(parent)
        self._items = []
        self._h_spacing = h_spacing
        self._v_spacing = v_spacing
        self.setContentsMargins(margin, margin, margin, margin)

    def addItem(self, item):
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index: int):
        if 0 <= index < len(self._items):
            return self._items[index]
        return None

    def takeAt(self, index: int):
        if 0 <= index < len(self._items):
            return self._items.pop(index)
        return None

    def expandingDirections(self):
        return Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._do_layout(QRect(0, 0, width, 0), test_only=True)

    def setGeometry(self, rect: QRect) -> None:
        super().setGeometry(rect)
        self._do_layout(rect, test_only=False)

    def sizeHint(self) -> QSize:
        return self.minimumSize()

    def minimumSize(self) -> QSize:
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        margins = self.contentsMargins()
        size += QSize(margins.left() + margins.right(), margins.top() + margins.bottom())
        return size

    def _do_layout(self, rect: QRect, *, test_only: bool) -> int:
        margins = self.contentsMargins()
        effective = rect.adjusted(margins.left(), margins.top(), -margins.right(), -margins.bottom())
        x = effective.x()
        y = effective.y()
        line_height = 0

        for item in self._items:
            hint = item.sizeHint()
            next_x = x + hint.width()
            if line_height > 0 and next_x > effective.right() + 1:
                x = effective.x()
                y += line_height + self._v_spacing
                next_x = x + hint.width()
                line_height = 0

            if not test_only:
                item.setGeometry(QRect(QPoint(x, y), hint))

            x = next_x + self._h_spacing
            line_height = max(line_height, hint.height())

        return y + line_height - rect.y() + margins.bottom()


class GeneralTraceCard(QFrame):
    def __init__(self, channel_id: str, subtitle: str, parent=None):
        super().__init__(parent)
        self.setObjectName("GeneralTraceCard")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumWidth(170)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(6)

        self.checkbox = QCheckBox(channel_id)
        self.checkbox.setObjectName("GeneralTraceCheckbox")
        layout.addWidget(self.checkbox)

        subtitle_lbl = QLabel(subtitle)
        subtitle_lbl.setObjectName("MutedText")
        subtitle_lbl.setWordWrap(True)
        layout.addWidget(subtitle_lbl)

        self.checkbox.toggled.connect(self._sync_checked_state)
        self._sync_checked_state(self.checkbox.isChecked())

    def _sync_checked_state(self, checked: bool) -> None:
        self.setProperty("checked", checked)
        self.style().unpolish(self)
        self.style().polish(self)
        self.update()

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton and self.isEnabled():
            self.checkbox.toggle()
            event.accept()
            return
        super().mousePressEvent(event)


class JobsTableWidget(QTableWidget):
    """QTableWidget with reliable wheel scrolling on the viewport."""

    def __init__(self, rows: int, columns: int, parent=None):
        super().__init__(rows, columns, parent)
        self.viewport().installEventFilter(self)
        QScroller.grabGesture(self.viewport(), QScroller.ScrollerGestureType.TouchGesture)

    def eventFilter(self, source, event):
        if source is self.viewport() and event.type() == QEvent.Type.Wheel:
            self.wheelEvent(event)
            return event.isAccepted()
        return super().eventFilter(source, event)

    def wheelEvent(self, event):
        delta = event.pixelDelta().y()
        if not delta:
            delta = event.angleDelta().y()
            if delta:
                delta = int(delta / 120) * self.verticalScrollBar().singleStep()

        if delta:
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                bar = self.horizontalScrollBar()
            else:
                bar = self.verticalScrollBar()
            bar.setValue(bar.value() - delta)
            event.accept()
            return

        super().wheelEvent(event)

class TabBatch(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.threadpool = QThreadPool.globalInstance()
        self._detected_jobs = []
        self._job_states = {}
        self._current_analysis_id = APP_SETTINGS.get("active_analysis", "clonality")
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(18)
        
        # Header
        header = QVBoxLayout()
        self.title_lbl = QLabel("Run Fraggler")
        self.title_lbl.setObjectName("PageTitle")
        self.subtitle_lbl = QLabel("")
        self.subtitle_lbl.setObjectName("PageSubtitle")
        self.subtitle_lbl.setWordWrap(True)
        header.addWidget(self.title_lbl)
        header.addWidget(self.subtitle_lbl)

        self._general_controls_ready = False
        self._general_trace_checkboxes: dict[str, QCheckBox] = {}

        # 1. General runtime card
        self.general_card = self._build_general_card()

        # 2. Folders / Files Card
        f_card = QWidget()
        f_card.setObjectName("Card")
        f_layout = QVBoxLayout(f_card)
        f_layout.setSpacing(12)
        
        l_ftitle = QLabel("SAMPLES")
        l_ftitle.setObjectName("CardTitle")

        row1 = QHBoxLayout()
        row1.setSpacing(10)
        self.folder_list = QListWidget()
        self.folder_list.setMaximumHeight(100)
        self.folder_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.folder_list.setAcceptDrops(True)
        self.folder_list.setAlternatingRowColors(True)
        
        # Inject Drag & Drop support
        def _dragEnterEvent(e):
            if e.mimeData().hasUrls():
                e.acceptProposedAction()
        
        def _dropEvent(e):
            for url in e.mimeData().urls():
                if url.isLocalFile():
                    path = url.toLocalFile()
                    self._add_source_item(path)
            e.acceptProposedAction()
                
        self.folder_list.dragEnterEvent = _dragEnterEvent
        self.folder_list.dragMoveEvent = _dragEnterEvent
        self.folder_list.dropEvent = _dropEvent
        
        btn_layout = QVBoxLayout()
        self.btn_add_folders = QPushButton("Add Folders...")
        self.btn_add_files = QPushButton("Add Files...")
        self.btn_remove_sources = QPushButton("Remove Selected")
        self.btn_add_folders.clicked.connect(self._add_folders)
        self.btn_add_files.clicked.connect(self._add_files)
        self.btn_remove_sources.clicked.connect(self._remove_sources)
        btn_layout.addWidget(self.btn_add_folders)
        btn_layout.addWidget(self.btn_add_files)
        btn_layout.addWidget(self.btn_remove_sources)
        btn_layout.addStretch()
        
        self.input_label = QLabel("Samples:")
        row1.addWidget(self.input_label)
        row1.addWidget(self.folder_list, stretch=1)
        row1.addLayout(btn_layout)
        
        row2 = QHBoxLayout()
        row2.setSpacing(10)
        self.output_base = QLineEdit("")
        self.output_base.setClearButtonEnabled(True)
        btn_browse_out = QPushButton("Browse...")
        btn_browse_out.clicked.connect(lambda: self._ask_dir(self.output_base))
        self.output_label = QLabel("Save To:")
        row2.addWidget(self.output_label)
        row2.addWidget(self.output_base, stretch=1)
        row2.addWidget(btn_browse_out)
        
        f_layout.addWidget(l_ftitle)
        f_layout.addLayout(row1)
        f_layout.addLayout(row2)
        
        # 2. Actions & Progress
        a_layout = QHBoxLayout()
        a_layout.setSpacing(10)
        self.btn_scan = QPushButton("Find Jobs")
        self.btn_run = QPushButton("Run Batch")
        self.btn_run.setObjectName("PrimaryButton")
        self.btn_run.setEnabled(False)
        self.btn_open = QPushButton("Open Output")
        
        self.progress = QProgressBar()
        self.progress.setValue(0)
        
        self.status_lbl = QLabel("Ready — review the loaded folders and click Find Jobs.")
        self.status_lbl.setStyleSheet("color: #64748b; font-weight: 500;")
        
        a_layout.addWidget(self.btn_scan)
        a_layout.addWidget(self.btn_run)
        a_layout.addWidget(self.btn_open)
        a_layout.addStretch()
        
        # 3. Jobs Table
        t_card = QWidget()
        t_card.setObjectName("Card")
        t_layout = QVBoxLayout(t_card)
        t_title = QLabel("DETECTED JOBS")
        t_title.setObjectName("CardTitle")
        
        t_btns = QHBoxLayout()
        self.btn_sel_all = QPushButton("Select All")
        self.btn_sel_none = QPushButton("Select None")
        t_btns.addWidget(self.btn_sel_all)
        t_btns.addWidget(self.btn_sel_none)
        t_btns.addStretch()
        
        self.table = JobsTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Name", "Type", "Source", "Files", "Status"])
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setMinimumSectionSize(80)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setWordWrap(False)
        self.table.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        self.table.setShowGrid(False)
        self.table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setAutoScroll(False)
        self.table.setMinimumHeight(280)
        self.table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.table.setFocusPolicy(Qt.FocusPolicy.WheelFocus)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(32)
        self.table.verticalHeader().setMinimumSectionSize(28)
        self.table.verticalScrollBar().setSingleStep(24)
        
        self.btn_sel_all.clicked.connect(self.table.selectAll)
        self.btn_sel_none.clicked.connect(self.table.clearSelection)
        self.btn_scan.clicked.connect(self.on_scan)
        self.btn_run.clicked.connect(self.on_run)
        self.btn_open.clicked.connect(self.on_open_output)
        
        t_layout.addWidget(t_title)
        t_layout.addLayout(t_btns)
        t_layout.addWidget(self.table, stretch=1)
        
        # Add to main
        main_layout.addLayout(header)
        main_layout.addWidget(self.general_card)
        main_layout.addWidget(f_card)
        main_layout.addLayout(a_layout)
        main_layout.addWidget(self.progress)
        main_layout.addWidget(self.status_lbl)
        main_layout.addWidget(t_card, stretch=1)

        self.set_analysis(self._current_analysis_id, force_replace_inputs=True)
        
    def _build_general_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        layout = QVBoxLayout(card)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(14)

        title = QLabel("General Analysis Runtime")
        title.setStyleSheet("font-size: 16px; font-weight: 700; color: #0f172a;")
        subtitle = QLabel(
            "Choose the ladder, the trace channels to show, and which channel should be used as the primary peak view."
        )
        subtitle.setObjectName("MutedText")
        subtitle.setWordWrap(True)
        layout.addWidget(title)
        layout.addWidget(subtitle)

        selector_grid = QGridLayout()
        selector_grid.setHorizontalSpacing(12)
        selector_grid.setVerticalSpacing(12)
        selector_grid.setColumnStretch(0, 1)
        selector_grid.setColumnStretch(1, 1)

        self.general_ladder_combo = QComboBox()
        for label, value in GENERAL_LADDER_OPTIONS:
            self.general_ladder_combo.addItem(label, value)
        self.general_ladder_combo.setMinimumContentsLength(12)
        self.general_ladder_combo.currentIndexChanged.connect(self._on_general_runtime_changed)
        selector_grid.addWidget(
            self._build_general_selector_card(
                "Ladder",
                "Supported ladders for the general workflow.",
                self.general_ladder_combo,
            ),
            0,
            0,
        )

        self.general_primary_combo = QComboBox()
        self.general_primary_combo.setMinimumContentsLength(12)
        self.general_primary_combo.currentIndexChanged.connect(self._on_general_runtime_changed)
        selector_grid.addWidget(
            self._build_general_selector_card(
                "Primary Peak Channel",
                "Used as the default sample trace for peak editing.",
                self.general_primary_combo,
            ),
            0,
            1,
        )
        layout.addLayout(selector_grid)

        trace_title = QLabel("Trace Channels")
        trace_title.setStyleSheet("font-size: 14px; font-weight: 700; color: #0f172a;")
        trace_note = QLabel(
            "Pick one or more trace channels. The cards wrap automatically when the window gets narrower."
        )
        trace_note.setObjectName("MutedText")
        trace_note.setWordWrap(True)
        layout.addWidget(trace_title)
        layout.addWidget(trace_note)

        trace_box = QWidget()
        trace_layout = FlowLayout(trace_box, h_spacing=12, v_spacing=12)
        for key, subtitle_text in GENERAL_TRACE_OPTIONS:
            option_card = GeneralTraceCard(key, subtitle_text)
            option_card.checkbox.toggled.connect(self._on_general_trace_toggled)
            self._general_trace_checkboxes[key] = option_card.checkbox
            trace_layout.addWidget(option_card)
        layout.addWidget(trace_box)

        note = QLabel(
            "General keeps the selected ladder and channel choices in the per-analysis settings so the backend can reuse them."
        )
        note.setObjectName("MutedText")
        note.setWordWrap(True)
        layout.addWidget(note)

        return card

    def _build_general_selector_card(self, title: str, subtitle: str, field: QWidget) -> QWidget:
        card = QFrame()
        card.setObjectName("GeneralSelectorCard")
        card.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        layout = QVBoxLayout(card)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(8)

        heading = QLabel(title)
        heading.setStyleSheet("font-size: 13px; font-weight: 700; color: #0f172a;")
        helper = QLabel(subtitle)
        helper.setObjectName("MutedText")
        helper.setWordWrap(True)

        layout.addWidget(heading)
        layout.addWidget(helper)
        layout.addWidget(field)
        return card

    def _profile_for(self, analysis_id: str | None = None) -> dict:
        return get_analysis_settings(analysis_id or self._current_analysis_id)

    def set_analysis(self, analysis_id: str, force_replace_inputs: bool = False) -> None:
        previous_profile = self._profile_for(self._current_analysis_id)
        previous_default = previous_profile.get("batch", {}).get("base_input_dir", "")
        current_items = [self.folder_list.item(i).text() for i in range(self.folder_list.count())]
        should_replace_inputs = force_replace_inputs or not current_items or current_items == [previous_default]

        self._current_analysis_id = analysis_id
        self.load_from_settings(replace_inputs=should_replace_inputs)
        pretty_name = ANALYSIS_LABELS.get(analysis_id, analysis_id.capitalize())
        self.title_lbl.setText(f"Run {pretty_name}")
        self.subtitle_lbl.setText(
            f"Saved defaults for {pretty_name.lower()} load automatically. Change them only if needed, then find and run jobs."
        )
        is_general = self._is_general_analysis()
        self.general_card.setVisible(is_general)
        self.btn_add_files.setVisible(is_general)
        self.input_label.setText("Files / Folders:" if is_general else "Samples:")
        self.status_lbl.setText(
            "Ready — review the selected files/folders and click Find Jobs."
            if is_general
            else "Ready — review the loaded folders and click Find Jobs."
        )

    def load_from_settings(self, replace_inputs: bool = False):
        """Reload analysis-specific defaults from APP_SETTINGS."""
        profile = self._profile_for()
        batch_settings = profile.get("batch", {})
        pipeline_settings = profile.get("pipeline", {})

        saved_output = batch_settings.get("output_base", "")
        self.output_base.setText(saved_output)
        self.output_base.setPlaceholderText(
            saved_output or "/path/to/output (leave empty to use the saved output or the first sample folder)"
        )

        default_dir = batch_settings.get("base_input_dir", "")
        if replace_inputs:
            self.folder_list.clear()
        if default_dir and self.folder_list.count() == 0:
            self.folder_list.addItem(default_dir)

        if self._is_general_analysis():
            self._load_general_runtime_controls(pipeline_settings)
        else:
            self._general_controls_ready = False

    def _ask_dir(self, widget: QLineEdit):
        folder = QFileDialog.getExistingDirectory(self, "Select Directory", widget.text() or str(Path.home()))
        if folder:
            widget.setText(folder)

    def _is_general_analysis(self) -> bool:
        return self._current_analysis_id == "general"

    def _load_general_runtime_controls(self, pipeline_settings: dict) -> None:
        self._general_controls_ready = False
        try:
            ladder = str(pipeline_settings.get("ladder", GENERAL_DEFAULT_LADDER))
            idx = self.general_ladder_combo.findData(ladder)
            if idx < 0 and ladder.upper() == "LIZ500":
                idx = self.general_ladder_combo.findData("LIZ500_250")
            if idx < 0:
                idx = 0
            ladder_blocker = QSignalBlocker(self.general_ladder_combo)
            self.general_ladder_combo.setCurrentIndex(idx)
            del ladder_blocker

            selected_traces = pipeline_settings.get("trace_channels", list(GENERAL_DEFAULT_TRACE_CHANNELS))
            if not isinstance(selected_traces, list):
                selected_traces = list(GENERAL_DEFAULT_TRACE_CHANNELS)
            selected_traces = [ch for ch in selected_traces if ch in dict(GENERAL_TRACE_OPTIONS)]
            if not selected_traces:
                selected_traces = list(GENERAL_DEFAULT_TRACE_CHANNELS)

            for key, checkbox in self._general_trace_checkboxes.items():
                trace_blocker = QSignalBlocker(checkbox)
                checkbox.setChecked(key in selected_traces)
                del trace_blocker

            self._refresh_general_primary_combo(
                preferred=str(pipeline_settings.get("primary_peak_channel", GENERAL_DEFAULT_PRIMARY_CHANNEL))
            )
        finally:
            self._general_controls_ready = True

    def _selected_general_trace_channels(self, *, fallback: bool = True) -> list[str]:
        selected = [key for key, checkbox in self._general_trace_checkboxes.items() if checkbox.isChecked()]
        if selected or not fallback:
            return selected
        return list(GENERAL_DEFAULT_TRACE_CHANNELS)

    def _refresh_general_primary_combo(self, preferred: str | None = None) -> None:
        selected = self._selected_general_trace_channels()
        target = preferred if preferred in selected else selected[0]
        combo_blocker = QSignalBlocker(self.general_primary_combo)
        self.general_primary_combo.clear()
        for channel in selected:
            self.general_primary_combo.addItem(channel, channel)
        index = self.general_primary_combo.findData(target)
        if index < 0:
            index = 0
        self.general_primary_combo.setCurrentIndex(index)
        del combo_blocker

    def _persist_general_runtime_settings(self) -> None:
        if not self._is_general_analysis():
            return
        trace_channels = self._selected_general_trace_channels()
        primary_channel = self.general_primary_combo.currentData() or self.general_primary_combo.currentText() or trace_channels[0]
        if primary_channel not in trace_channels:
            primary_channel = trace_channels[0]

        profile = APP_SETTINGS.setdefault("analyses", {}).setdefault("general", {})
        pipeline_settings = profile.setdefault("pipeline", {})
        pipeline_settings["ladder"] = self.general_ladder_combo.currentData() or GENERAL_DEFAULT_LADDER
        pipeline_settings["trace_channels"] = trace_channels
        pipeline_settings["peak_channels"] = list(trace_channels)
        pipeline_settings["primary_peak_channel"] = primary_channel

        if APP_SETTINGS.get("active_analysis") == "general":
            APP_SETTINGS.setdefault("pipeline", {}).update(pipeline_settings)

        save_settings(APP_SETTINGS)

    def _on_general_runtime_changed(self, *_args) -> None:
        if not self._is_general_analysis() or not self._general_controls_ready:
            return
        preferred = self.general_primary_combo.currentData() or self.general_primary_combo.currentText()
        self._refresh_general_primary_combo(preferred=str(preferred) if preferred else None)
        self._persist_general_runtime_settings()

    def _on_general_trace_toggled(self, *_args) -> None:
        if not self._is_general_analysis() or not self._general_controls_ready:
            return
        selected = self._selected_general_trace_channels(fallback=False)
        if not selected:
            blocker = QSignalBlocker(self._general_trace_checkboxes["DATA1"])
            self._general_trace_checkboxes["DATA1"].setChecked(True)
            del blocker
        preferred = self.general_primary_combo.currentData() or self.general_primary_combo.currentText()
        self._refresh_general_primary_combo(preferred=str(preferred) if preferred else None)
        self._persist_general_runtime_settings()

    def _add_folders(self):
        dialog = QFileDialog(self, "Add Folders", str(Path.home()))
        dialog.setFileMode(QFileDialog.FileMode.Directory)
        dialog.setOption(QFileDialog.Option.DontUseNativeDialog, True)
        
        # Enable multiple selection in the dialog's views
        for view in dialog.findChildren(QAbstractItemView):
            view.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
            
        if dialog.exec():
            folders = dialog.selectedFiles()
            existing = {self.folder_list.item(i).text() for i in range(self.folder_list.count())}
            for folder in folders:
                if folder not in existing:
                    self._add_source_item(folder)

    def _add_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Add .fsa Files",
            str(Path.home()),
            "FSA files (*.fsa)",
        )
        for file_name in files:
            self._add_source_item(file_name)

    def _remove_sources(self):
        for item in self.folder_list.selectedItems():
            self.folder_list.takeItem(self.folder_list.row(item))

    def _add_source_item(self, path_text: str) -> None:
        path = Path(path_text).expanduser()
        if not path.exists():
            return
        if path.is_file() and path.suffix.lower() != ".fsa":
            return
        if not self._is_general_analysis() and path.is_file():
            return
        existing = {self.folder_list.item(i).text() for i in range(self.folder_list.count())}
        normalized = str(path)
        if normalized not in existing:
            self.folder_list.addItem(normalized)

    def _expand_general_path(self, path: Path) -> list[Path]:
        if path.is_file():
            return [path] if path.suffix.lower() == ".fsa" else []
        if path.is_dir():
            return sorted(
                {
                    candidate
                    for candidate in path.rglob("*")
                    if candidate.is_file() and candidate.suffix.lower() == ".fsa"
                },
                key=lambda p: p.as_posix().lower(),
            )
        return []

    def _build_general_jobs_worker(self, input_paths: list[Path]) -> list[dict]:
        jobs: list[dict] = []
        seen_names: set[str] = set()
        seen_files: set[str] = set()

        for source in input_paths:
            files = self._expand_general_path(source)
            if not files:
                continue
            deduped_files: list[Path] = []
            for file_path in files:
                key = file_path.resolve().as_posix() if file_path.exists() else str(file_path)
                if key in seen_files:
                    continue
                seen_files.add(key)
                deduped_files.append(file_path)
            if not deduped_files:
                continue

            base_name = source.stem if source.is_file() else source.name
            job_name = base_name or source.name
            suffix = 2
            while job_name in seen_names:
                job_name = f"{base_name}_{suffix}"
                suffix += 1
            seen_names.add(job_name)

            jobs.append(
                {
                    "name": job_name,
                    "type": "pipeline",
                    "path": source,
                    "files": deduped_files,
                }
            )

        return jobs

    def _general_selected_paths(self) -> list[Path]:
        paths: list[Path] = []
        for i in range(self.folder_list.count()):
            p_str = self.folder_list.item(i).text().strip()
            if p_str:
                paths.append(Path(p_str).expanduser())
        return paths

    def _build_general_scan_worker(self, paths: list[Path]):
        return Worker(self._build_general_jobs_worker, paths)

    def _rebuild_table(self):
        selected_names = {
            self.table.item(index.row(), 0).text()
            for index in self.table.selectionModel().selectedRows()
            if self.table.item(index.row(), 0) is not None
        } if self.table.selectionModel() else set()
        v_scroll = self.table.verticalScrollBar().value()
        h_scroll = self.table.horizontalScrollBar().value()

        self.table.setUpdatesEnabled(False)
        self.table.setRowCount(0)
        for row_idx, j in enumerate(self._detected_jobs):
            self.table.insertRow(row_idx)
            
            state = self._job_states.get(j["name"], "pending")
            
            item_name = QTableWidgetItem(j["name"])
            
            jtype = j.get("type", "unknown")
            item_type = QTableWidgetItem(jtype.upper())
            if jtype == "qc":
                item_type.setForeground(Qt.GlobalColor.darkMagenta)
            else:
                item_type.setForeground(Qt.GlobalColor.darkCyan)
            
            src = str(j["path"]) if j.get("path") else "[Aggregated]"
            item_src = QTableWidgetItem(src)
            
            files = str(len(j.get("files", []))) if j.get("files") else "auto"
            item_files = QTableWidgetItem(files)
            
            display_state = state.upper()
            if ":" in display_state:
                display_state = display_state.split(":", 1)[0]
            item_state = QTableWidgetItem(display_state)
            if state == "success" or state == "done":
                item_state.setForeground(Qt.GlobalColor.darkGreen)
            elif state == "error" or state.startswith("error"):
                item_state.setForeground(Qt.GlobalColor.red)
            elif state == "running":
                item_state.setForeground(Qt.GlobalColor.blue)
            else:
                item_state.setForeground(Qt.GlobalColor.darkGray)
                
            self.table.setItem(row_idx, 0, item_name)
            self.table.setItem(row_idx, 1, item_type)
            self.table.setItem(row_idx, 2, item_src)
            self.table.setItem(row_idx, 3, item_files)
            self.table.setItem(row_idx, 4, item_state)
            if j["name"] in selected_names:
                self.table.selectRow(row_idx)

        self.table.setUpdatesEnabled(True)
        self.table.verticalScrollBar().setValue(v_scroll)
        self.table.horizontalScrollBar().setValue(h_scroll)
    def on_scan(self):
        paths = self._general_selected_paths()
        if not paths:
            self.status_lbl.setText(
                "No files or folders selected." if self._is_general_analysis() else "No input folders selected."
            )
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
            return
            
        self.btn_scan.setEnabled(False)
        self.btn_run.setEnabled(False)
        self.progress.setRange(0, 0) # Indeterminate spinner
        self.status_lbl.setText("Finding jobs...")
        self.status_lbl.setStyleSheet("color: #f59e0b; font-weight: 500;")

        if self._is_general_analysis():
            worker = self._build_general_scan_worker(paths)
        else:
            from core.batch import generate_jobs

            batch_settings = self._profile_for().get("batch", {})
            agg_pat = bool(batch_settings.get("aggregate_by_patient", True))
            regex = batch_settings.get("patient_id_regex", r"\d{2}OUM\d{5}")

            worker = Worker(
                generate_jobs,
                input_paths=paths,
                aggregate_patients=agg_pat,
                patient_regex=regex
            )
        worker.signals.result.connect(self._on_scan_result)
        worker.signals.error.connect(self._on_scan_error)
        
        self.threadpool.start(worker)
        
    def _on_scan_result(self, jobs):
        self._detected_jobs = jobs
        self._job_states = {j["name"]: "pending" for j in jobs}
        
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        
        if not jobs:
            self.status_lbl.setText(
                "No jobs found — check files/folders."
                if self._is_general_analysis()
                else "No jobs found — check input folders."
            )
            self.status_lbl.setStyleSheet("color: #f59e0b; font-weight: 500;")
        else:
            self.status_lbl.setText(f"Found {len(jobs)} jobs — ready to run.")
            self.status_lbl.setStyleSheet("color: #22c55e; font-weight: 500;")
            self.btn_run.setEnabled(True)
            
        self._rebuild_table()
        self.btn_scan.setEnabled(True)
        
    def _on_scan_error(self, err_tuple):
        self.status_lbl.setText(f"Scan error: {err_tuple[1]}")
        self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.btn_scan.setEnabled(True)
        self.btn_run.setEnabled(bool(self._detected_jobs))

    def _on_run_error(self, err_tuple):
        self.status_lbl.setText(f"Run error: {err_tuple[1]}")
        self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
        self.progress.setRange(0, max(len(self._detected_jobs), 1))
        self.btn_scan.setEnabled(True)
        self.btn_run.setEnabled(True)
        
    def on_run(self):
        from core.batch import run_batch_jobs
        
        selected_rows = [index.row() for index in self.table.selectionModel().selectedRows()]
        if not selected_rows:
            self.status_lbl.setText(
                "No files or folders selected — check rows in the table."
                if self._is_general_analysis()
                else "No jobs selected — check rows in the table."
            )
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
            return
            
        jobs_to_run = [self._detected_jobs[i] for i in selected_rows]
        
        out_path_str = self._resolve_output_path_str()
            
        out_path_obj = Path(out_path_str).expanduser() if out_path_str else None
        
        if not out_path_obj or not out_path_obj.exists():
            self.status_lbl.setText("Output folder does not exist — set it before running.")
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
            return
            
        for j in jobs_to_run:
            self._job_states[j["name"]] = "running"
        self._rebuild_table()
        
        self.btn_scan.setEnabled(False)
        self.btn_run.setEnabled(False)
        self.progress.setRange(0, len(jobs_to_run))
        self.progress.setValue(0)
        
        self.status_lbl.setText(f"Running {len(jobs_to_run)} jobs...")
        self.status_lbl.setStyleSheet("color: #f59e0b; font-weight: 500;")
        
        profile = self._profile_for()
        s_pipe = profile.get("pipeline", {})
        s_batch = profile.get("batch", {})
        p_scope = s_pipe.get("mode", "all")
        a_filter = s_pipe.get("assay_filter_substring", "")
        aggregate_dit_reports = bool(s_batch.get("aggregate_dit_reports", True))
        if self._is_general_analysis():
            self._persist_general_runtime_settings()
        
        worker = Worker(
            run_batch_jobs,
            jobs=jobs_to_run,
            output_base=out_path_obj,
            out_folder_tmpl="ASSAY_REPORTS",
            outfile_html_tmpl="QC_REPORT_{name}.html",
            excel_name_tmpl="Fraggler_QC_Trends.xlsx",
            pipeline_scope=p_scope,
            assay_filter=a_filter,
            aggregate_dit_reports=aggregate_dit_reports,
            continue_on_error=True,
            update_callback=None, # Passed explicitly as kwarg below
        )
        # Assign the emit method of our new progress_ext signal as the callback
        worker.kwargs['update_callback'] = worker.signals.progress_ext.emit
        
        worker.signals.result.connect(self._on_run_finished)
        worker.signals.progress_ext.connect(self._update_progress_from_thread)
        worker.signals.error.connect(self._on_run_error)
        
        self.threadpool.start(worker)
        
    def _update_progress_from_thread(self, idx, total, name, state):
        self._job_states[name] = state
        self._rebuild_table()
        self.progress.setValue(idx)
        if state.startswith("error"):
            self.status_lbl.setText(f"Run error in {name} ({idx}/{total})")
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
        elif state == "success":
            self.status_lbl.setText(f"Completed: {name} ({idx}/{total})")
            self.status_lbl.setStyleSheet("color: #22c55e; font-weight: 500;")
        elif state == "done":
            pass
        else:
            self.status_lbl.setText(f"Running: {name} ({idx}/{total})")
            self.status_lbl.setStyleSheet("color: #f59e0b; font-weight: 500;")
        
    def _on_run_finished(self, result):
        self.progress.setRange(0, 100)
        self.progress.setValue(100)
        failed_jobs = (result or {}).get("failed_jobs", [])
        if failed_jobs:
            self.status_lbl.setText(f"Batch finished with {len(failed_jobs)} failed job(s).")
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
        else:
            self.status_lbl.setText("Batch complete.")
            self.status_lbl.setStyleSheet("color: #22c55e; font-weight: 500;")
        self.btn_scan.setEnabled(True)
        self.btn_run.setEnabled(True)
        
    def on_open_output(self):
        p_str = self._resolve_output_path_str()
            
        p = Path(p_str).expanduser() if p_str else None
        if p and p.exists():
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(p)])
            elif sys.platform == "win32":
                subprocess.Popen(["explorer", str(p)])
            else:
                subprocess.Popen(["xdg-open", str(p)])

    def _resolve_output_path_str(self) -> str:
        explicit_output = self.output_base.text().strip()
        if explicit_output:
            return explicit_output

        saved_output = self._profile_for().get("batch", {}).get("output_base", "").strip()
        if saved_output:
            return saved_output

        if self.folder_list.count() > 0:
            first_item = Path(self.folder_list.item(0).text().strip()).expanduser()
            if first_item.is_file():
                return str(first_item.parent)
            return str(first_item)
        return ""
