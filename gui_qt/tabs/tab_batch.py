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
from config import APP_SETTINGS

class TabBatch(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.threadpool = QThreadPool.globalInstance()
        self._detected_jobs = []
        self._job_states = {}
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(16)
        
        # Header
        header = QVBoxLayout()
        title = QLabel("Run Fraggler")
        title.setObjectName("PageTitle")
        sub = QLabel("Select folders containing patient or QC data, then scan and run.")
        sub.setObjectName("PageSubtitle")
        header.addWidget(title)
        header.addWidget(sub)
        
        s_batch = APP_SETTINGS.get("batch", {})
        
        # 1. Folders Card
        f_card = QWidget()
        f_card.setObjectName("Card")
        f_layout = QVBoxLayout(f_card)
        
        l_ftitle = QLabel("FOLDERS")
        l_ftitle.setObjectName("CardTitle")
        
        row1 = QHBoxLayout()
        self.folder_list = QListWidget()
        self.folder_list.setMaximumHeight(100)
        self.folder_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.folder_list.setAcceptDrops(True)
        
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
        
        # Seed with default if set
        default_dir = s_batch.get("base_input_dir", "")
        if default_dir:
            self.folder_list.addItem(default_dir)
            
        btn_layout = QVBoxLayout()
        btn_add = QPushButton("Add Folders...")
        btn_add.clicked.connect(self._add_folders)
        btn_remove = QPushButton("Remove Selected")
        btn_remove.clicked.connect(self._remove_folders)
        btn_layout.addWidget(btn_add)
        btn_layout.addWidget(btn_remove)
        btn_layout.addStretch()
        
        row1.addWidget(QLabel("Input Folders:"))
        row1.addWidget(self.folder_list, stretch=1)
        row1.addLayout(btn_layout)
        
        row2 = QHBoxLayout()
        self.output_base = QLineEdit(s_batch.get("output_base", ""))
        self.output_base.setPlaceholderText("/path/to/output (leave empty = same as input)")
        btn_browse_out = QPushButton("Browse...")
        btn_browse_out.clicked.connect(lambda: self._ask_dir(self.output_base))
        row2.addWidget(QLabel("Output Folder:"))
        row2.addWidget(self.output_base, stretch=1)
        row2.addWidget(btn_browse_out)
        
        f_layout.addWidget(l_ftitle)
        f_layout.addLayout(row1)
        f_layout.addLayout(row2)
        
        # 2. Actions & Progress
        a_layout = QHBoxLayout()
        self.btn_scan = QPushButton("Scan Jobs")
        self.btn_run = QPushButton("Run Batch")
        self.btn_run.setObjectName("PrimaryButton")
        self.btn_run.setEnabled(False)
        self.btn_open = QPushButton("Open Output")
        
        self.progress = QProgressBar()
        self.progress.setValue(0)
        
        self.status_lbl = QLabel("Ready — add folders and click Scan Jobs.")
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
        
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.on_context_menu)
        
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
        
    def load_from_settings(self):
        """Reload default folders and output base from APP_SETTINGS."""
        from config import APP_SETTINGS
        s_batch = APP_SETTINGS.get("batch", {})
        
        # Output base
        self.output_base.setText(s_batch.get("output_base", ""))
        
        # Base input dir - only add if list is empty or we specifically want to sync
        # For now, let's just update the list if it's currently empty
        default_dir = s_batch.get("base_input_dir", "")
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
            
            item_state = QTableWidgetItem(state.upper())
            if state == "success" or state == "done":
                item_state.setForeground(Qt.GlobalColor.darkGreen)
            elif state == "error":
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
        self.status_lbl.setText("Scanning for jobs...")
        self.status_lbl.setStyleSheet("color: #f59e0b; font-weight: 500;")
        
        # Build Worker for the scan
        from config import APP_SETTINGS
        s_batch = APP_SETTINGS.get("batch", {})
        agg_pat = True # User requested mandatory patient ID aggregation
        regex = s_batch.get("patient_id_regex", r"\d{2}OUM\d{5}")
        
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
        
    def on_run(self):
        from core.batch import run_batch_jobs
        
        selected_rows = [index.row() for index in self.table.selectionModel().selectedRows()]
        if not selected_rows:
            self.status_lbl.setText("No jobs selected — check rows in the table.")
            self.status_lbl.setStyleSheet("color: #ef4444; font-weight: 500;")
            return
            
        jobs_to_run = [self._detected_jobs[i] for i in selected_rows]
        
        out_path_str = self.output_base.text().strip()
        # Fallback to the first input folder if output is empty
        if not out_path_str and self.folder_list.count() > 0:
            out_path_str = self.folder_list.item(0).text().strip()
            
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
        
        # Read pipeline_scope and assay_filter from settings
        from config import APP_SETTINGS
        s_pipe = APP_SETTINGS.get("pipeline", {})
        p_scope = s_pipe.get("mode", "all")
        a_filter = s_pipe.get("assay_filter_substring", "")
        
        worker = Worker(
            run_batch_jobs,
            jobs=jobs_to_run,
            output_base=out_path_obj,
            out_folder_tmpl="ASSAY_REPORTS",
            outfile_html_tmpl="QC_REPORT_{name}.html",
            excel_name_tmpl="Fraggler_QC_Trends.xlsx",
            pipeline_scope=p_scope,
            assay_filter=a_filter,
            continue_on_error=True,
            update_callback=None, # Passed explicitly as kwarg below
        )
        # Assign the emit method of our new progress_ext signal as the callback
        worker.kwargs['update_callback'] = worker.signals.progress_ext.emit
        
        worker.signals.result.connect(self._on_run_finished)
        worker.signals.progress_ext.connect(self._update_progress_from_thread)
        worker.signals.error.connect(self._on_scan_error)
        
        self.threadpool.start(worker)
        
    def _update_progress_from_thread(self, idx, total, name, state):
        self._job_states[name] = state
        self._rebuild_table()
        self.progress.setValue(idx)
        self.status_lbl.setText(f"Running: {name} ({idx}/{total})")
        
    def _on_run_finished(self, result):
        self.progress.setRange(0, 100)
        self.progress.setValue(100)
        self.status_lbl.setText("Batch complete.")
        self.status_lbl.setStyleSheet("color: #22c55e; font-weight: 500;")
        self.btn_scan.setEnabled(True)
        self.btn_run.setEnabled(True)
        
    def on_open_output(self):
        p_str = self.output_base.text().strip()
        if not p_str and self.folder_list.count() > 0:
            p_str = self.folder_list.item(0).text().strip()
            
        p = Path(p_str).expanduser() if p_str else None
        if p and p.exists():
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(p)])
            elif sys.platform == "win32":
                subprocess.Popen(["explorer", str(p)])
            else:
                subprocess.Popen(["xdg-open", str(p)])

    def on_context_menu(self, pos):
        from PyQt6.QtWidgets import QMenu
        from PyQt6.QtGui import QAction
        menu = QMenu()
        adj_action = menu.addAction("Adjust Ladder...")
        adj_action.setEnabled(len(self.table.selectionModel().selectedRows()) == 1)
        
        action = menu.exec(self.table.viewport().mapToGlobal(pos))
        if action == adj_action:
            self._on_adjust_ladder()

    def _on_adjust_ladder(self):
        selected_rows = [index.row() for index in self.table.selectionModel().selectedRows()]
        if not selected_rows:
            return
            
        job = self._detected_jobs[selected_rows[0]]
        files = job.get("files", [])
        if not files:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "No Files", "This job has no files to adjust.")
            return
            
        # Adjust first file in the job's list
        fsa_path = files[0]
        
        from gui_qt.dialogs.ladder_dialog import LadderAdjustmentDialog
        from core.classification import classify_fsa
        from fraggler.fraggler import FsaFile
        from core.analysis import analyse_fsa_liz, analyse_fsa_rox
        from PyQt6.QtWidgets import QMessageBox

        self.status_lbl.setText(f"Loading {fsa_path.name} for adjustment...")
        
        # Classification to get params
        classified = classify_fsa(fsa_path)
        if not classified:
            QMessageBox.warning(self, "Error", "Could not classify file.")
            return
            
        (assay, group, ladder, trace_channels, peak_channels, primary, bp_min, bp_max) = classified
        sample_channel = trace_channels[0]
        
        # Initial fit attempt
        if ladder == "LIZ":
            fsa = analyse_fsa_liz(fsa_path, sample_channel)
        else:
            fsa = analyse_fsa_rox(fsa_path, sample_channel)
            
        if not fsa:
            from core.assay_config import (
                LIZ_LADDER, ROX_LADDER, 
                MIN_DISTANCE_BETWEEN_PEAKS_LIZ, MIN_SIZE_STANDARD_HEIGHT_LIZ,
                MIN_DISTANCE_BETWEEN_PEAKS_ROX, MIN_SIZE_STANDARD_HEIGHT_ROX
            )
            if ladder == "LIZ":
                fsa = FsaFile(str(fsa_path), LIZ_LADDER, sample_channel, 
                             MIN_DISTANCE_BETWEEN_PEAKS_LIZ, MIN_SIZE_STANDARD_HEIGHT_LIZ,
                             size_standard_channel="DATA105")
            else:
                fsa = FsaFile(str(fsa_path), ROX_LADDER, sample_channel,
                             MIN_DISTANCE_BETWEEN_PEAKS_ROX, MIN_SIZE_STANDARD_HEIGHT_ROX,
                             size_standard_channel="DATA4")
                             
        dialog = LadderAdjustmentDialog(fsa, self)
        if dialog.exec():
            mapping = dialog.get_mapping()
            from core.analysis import apply_manual_ladder_mapping, save_ladder_adjustment
            apply_manual_ladder_mapping(fsa, mapping)
            
            # SAVE TO DISK for persistence
            save_ladder_adjustment(fsa, mapping)
            
            # Message to user
            QMessageBox.information(self, "Success", f"Ladder for {fsa.file_name} adjusted manually and saved. \n"
                                                      "To apply the correction to the reports, RE-RUN the job.")
            self.status_lbl.setText(f"Manual adjustment for {fsa.file_name} saved.")
