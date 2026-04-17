from __future__ import annotations

import copy
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication

from config import APP_SETTINGS
from gui_qt.main_window import MainWindow
from gui_qt.tabs.tab_archive_runner import TabArchiveRunner
from scripts.combine_clonality_yearly_overview import combine_run_root


_APP = QApplication.instance() or QApplication([])


class TestTabArchiveRunner(unittest.TestCase):
    def setUp(self):
        self._settings_backup = copy.deepcopy(APP_SETTINGS)
        APP_SETTINGS["active_analysis"] = "clonality"
        APP_SETTINGS.setdefault("analyses", {})
        APP_SETTINGS["analyses"]["clonality"] = {
            "batch": {
                "base_input_dir": "/tmp/clonality-input",
                "output_base": "/tmp/clonality-output",
                "tracking_excel_path": "",
                "aggregate_by_patient": True,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": True,
            },
            "pipeline": {
                "mode": "all",
                "assay_filter_substring": "",
            },
            "archive_runner": {
                "input_root": "/tmp/archive-input",
                "output_root": "/tmp/archive-output",
                "year_label": "2024",
                "run_name": "saved_run",
                "max_workers": 2,
                "folder_workers": 3,
                "include_sl": True,
                "refresh_each_folder": True,
                "cleanup_staging_root": True,
                "last_selected_run_root": "",
            },
        }

    def tearDown(self):
        APP_SETTINGS.clear()
        APP_SETTINGS.update(self._settings_backup)

    def test_tab_loads_saved_defaults(self):
        widget = TabArchiveRunner()

        self.assertEqual(widget.year_input.text(), "2024")
        self.assertEqual(widget.input_root.text(), "/tmp/archive-input")
        self.assertEqual(widget.output_root.text(), "/tmp/archive-output")
        self.assertEqual(widget.run_name.text(), "saved_run")
        self.assertEqual(widget.max_workers.value(), 2)
        self.assertEqual(widget.folder_workers.value(), 3)
        self.assertTrue(widget.chk_include_sl.isChecked())
        self.assertTrue(widget.chk_refresh_each_folder.isChecked())
        self.assertTrue(widget.chk_cleanup_staging.isChecked())
        self.assertFalse(widget.chk_resume.isChecked())
        self.assertEqual(len(widget._selected_months()), 12)

    def test_main_window_exposes_archive_runner_only_for_clonality(self):
        window = MainWindow()

        self.assertIsNotNone(window.group_clonality.btn_archive)
        self.assertEqual(len(window.group_clonality.sub_buttons), 5)
        self.assertIsNotNone(window.group_flt3.btn_archive)
        self.assertEqual(len(window.group_flt3.sub_buttons), 5)
        self.assertIsNone(window.group_general.btn_archive)
        self.assertEqual(len(window.group_general.sub_buttons), 4)
        self.assertEqual(window.btn_about.text(), "About")

    def test_run_yearly_launches_worker_with_expected_args(self):
        with TemporaryDirectory() as tmp:
            input_root = Path(tmp) / "input"
            input_root.mkdir()
            (input_root / "2024_01_example").mkdir()
            output_root = Path(tmp) / "out"
            output_root.mkdir()

            widget = TabArchiveRunner()
            widget.year_input.setText("2024")
            widget.input_root.setText(str(input_root))
            widget.output_root.setText(str(output_root))
            widget.run_name.setText("archive_run")
            widget.max_workers.setValue(4)
            widget.folder_workers.setValue(2)
            widget.chk_include_sl.setChecked(True)
            widget.chk_refresh_each_folder.setChecked(True)
            widget.chk_cleanup_staging.setChecked(True)
            widget._set_all_months(False)
            widget._month_checkboxes["01"].setChecked(True)

            captured = {}

            def fake_start(worker):
                captured["worker"] = worker

            widget.threadpool.start = fake_start  # type: ignore[method-assign]
            widget.on_run_yearly()

            worker = captured["worker"]
            self.assertEqual(worker.fn.__name__, "run_yearly_validation")
            self.assertEqual(worker.kwargs["year_label"], "2024")
            self.assertEqual(worker.kwargs["input_root"], input_root)
            self.assertEqual(worker.kwargs["output_root"], output_root)
            self.assertEqual(worker.kwargs["run_name"], "archive_run")
            self.assertEqual(worker.kwargs["months"], ["2024_01"])
            self.assertEqual(worker.kwargs["max_workers"], 4)
            self.assertEqual(worker.kwargs["folder_workers"], 2)
            self.assertTrue(worker.kwargs["include_sl"])
            self.assertTrue(worker.kwargs["refresh_each_folder"])
            self.assertTrue(worker.kwargs["cleanup_staging_root"])
            self.assertTrue(callable(worker.kwargs["progress_callback"]))
            self.assertTrue(callable(worker.kwargs["status_callback"]))

    def test_build_combined_workbook_launches_worker_for_current_run_root(self):
        with TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "full_2024_validation"
            run_root.mkdir()
            widget = TabArchiveRunner()
            widget.year_input.setText("2024")
            widget._current_run_root = run_root
            widget._refresh_action_buttons()

            captured = {}

            def fake_start(worker):
                captured["worker"] = worker

            widget.threadpool.start = fake_start  # type: ignore[method-assign]
            widget.on_build_combined_workbook()

            worker = captured["worker"]
            self.assertEqual(worker.fn, combine_run_root)
            self.assertEqual(worker.args[0], run_root)
            self.assertEqual(worker.args[1], run_root / "track-clonality-2024-overview.xlsx")
            self.assertEqual(worker.kwargs["year_label"], "2024")

    def test_runner_events_update_month_rows_and_paths(self):
        with TemporaryDirectory() as tmp:
            widget = TabArchiveRunner()
            widget.year_input.setText("2024")
            widget._set_all_months(False)
            widget._month_checkboxes["01"].setChecked(True)
            widget._rebuild_month_table()

            run_root = Path(tmp) / "year_2024"
            manifest_path = run_root / "full_2024_run_manifest.json"
            widget._on_runner_event({"event": "run_started", "run_dir": str(run_root)})
            widget._on_runner_event({"event": "month_started", "month": "2024_01", "folder_count": 7, "run_dir": str(run_root)})
            widget._on_runner_event({"event": "manifest_written", "manifest_path": str(manifest_path)})

            row = widget._month_row_map["2024_01"]
            self.assertEqual(widget.month_table.item(row, 1).text(), "running")
            self.assertEqual(widget.month_table.item(row, 2).text(), "7")
            self.assertEqual(widget.run_root_lbl.text(), str(run_root))
            self.assertEqual(widget.manifest_lbl.text(), str(manifest_path))


if __name__ == "__main__":
    unittest.main()
