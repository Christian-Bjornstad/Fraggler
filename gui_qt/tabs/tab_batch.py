from pathlib import Path
import subprocess
import sys
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, 
    QTableWidget, QTableWidgetItem, QProgressBar, 
    QHeaderView, QAbstractItemView, QFileDialog, QListWidget
)
from PyQt6.QtCore import Qt, QThreadPool
from gui_qt.worker import Worker
from config import APP_SETTINGS, get_analysis_settings


ANALYSIS_LABELS = {
    "clonality": "Klonalitet",
    "flt3": "FLT3 Analysis",
}

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
        header.addWidget(self.title_lbl)
        header.addWidget(self.subtitle_lbl)
        
        # 1. Folders Card
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
                    if Path(path).is_dir():
                        existing = [self.folder_list.item(i).text() for i in range(self.folder_list.count())]
                        if path not in existing:
                            self.folder_list.addItem(path)
            e.acceptProposedAction()
                
        self.folder_list.dragEnterEvent = _dragEnterEvent
        self.folder_list.dragMoveEvent = _dragEnterEvent
        self.folder_list.dropEvent = _dropEvent
        
        btn_layout = QVBoxLayout()
        btn_add = QPushButton("Add Samples...")
        btn_add.clicked.connect(self._add_folders)
        btn_remove = QPushButton("Remove Selected")
        btn_remove.clicked.connect(self._remove_folders)
        btn_layout.addWidget(btn_add)
        btn_layout.addWidget(btn_remove)
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
        
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Name", "Type", "Source", "Files", "Status"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        
        self.btn_sel_all.clicked.connect(self.table.selectAll)
        self.btn_sel_none.clicked.connect(self.table.clearSelection)
        self.btn_scan.clicked.connect(self.on_scan)
        self.btn_run.clicked.connect(self.on_run)
        self.btn_open.clicked.connect(self.on_open_output)
        
        t_layout.addWidget(t_title)
        t_layout.addLayout(t_btns)
        t_layout.addWidget(self.table)
        
        # Add to main
        main_layout.addLayout(header)
        main_layout.addWidget(f_card)
        main_layout.addLayout(a_layout)
        main_layout.addWidget(self.progress)
        main_layout.addWidget(self.status_lbl)
        main_layout.addWidget(t_card, stretch=1)

        self.set_analysis(self._current_analysis_id, force_replace_inputs=True)
        
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

    def load_from_settings(self, replace_inputs: bool = False):
        """Reload analysis-specific defaults from APP_SETTINGS."""
        profile = self._profile_for()
        batch_settings = profile.get("batch", {})

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
            
    def _ask_dir(self, widget: QLineEdit):
        folder = QFileDialog.getExistingDirectory(self, "Select Directory", widget.text() or str(Path.home()))
        if folder:
            widget.setText(folder)
            
    def _add_folders(self):
        dialog = QFileDialog(self, "Add Folders", str(Path.home()))
        dialog.setFileMode(QFileDialog.FileMode.Directory)
        dialog.setOption(QFileDialog.Option.DontUseNativeDialog, True)
        
        # Enable multiple selection in the dialog's views
        for view in dialog.findChildren(QAbstractItemView):
            view.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
            
        if dialog.exec():
            folders = dialog.selectedFiles()
            existing = [self.folder_list.item(i).text() for i in range(self.folder_list.count())]
            for folder in folders:
                if folder not in existing:
                    self.folder_list.addItem(folder)

    def _remove_folders(self):
        for item in self.folder_list.selectedItems():
            self.folder_list.takeItem(self.folder_list.row(item))
            
    def _rebuild_table(self):
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
    def on_scan(self):
        from core.batch import generate_jobs
        
        paths = []
        for i in range(self.folder_list.count()):
            p_str = self.folder_list.item(i).text().strip()
            if p_str:
                paths.append(Path(p_str).expanduser())
                
        if not paths:
            self.status_lbl.setText("No input folders selected.")
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
            return
            
        self.btn_scan.setEnabled(False)
        self.btn_run.setEnabled(False)
        self.progress.setRange(0, 0) # Indeterminate spinner
        self.status_lbl.setText("Finding jobs...")
        self.status_lbl.setStyleSheet("color: #f59e0b; font-weight: 500;")
        
        # Build Worker for the scan
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
            self.status_lbl.setText("No jobs found — check input folders.")
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
            self.status_lbl.setText("No jobs selected — check rows in the table.")
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
            return self.folder_list.item(0).text().strip()
        return ""
