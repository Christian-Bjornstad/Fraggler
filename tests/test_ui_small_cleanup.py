import copy
import os
import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication, QScrollArea, QTextBrowser

from config import APP_SETTINGS
from gui_qt.main_window import MainWindow
from gui_qt.tabs.tab_archive_runner import TabArchiveRunner
from gui_qt.tabs.tab_about import TabAbout
from gui_qt.tabs.tab_batch import TabBatch
from gui_qt.tabs.tab_settings import TabAnalysisSettings


_APP = QApplication.instance() or QApplication([])


class TestUiSmallCleanup(unittest.TestCase):
    def setUp(self):
        self._settings_backup = copy.deepcopy(APP_SETTINGS)
        APP_SETTINGS["active_analysis"] = "clonality"
        APP_SETTINGS.setdefault("analyses", {})
        APP_SETTINGS["analyses"]["clonality"] = {
            "batch": {
                "base_input_dir": "/tmp/clonality-input",
                "output_base": "/tmp/clonality-output",
                "tracking_excel_path": "/tmp/clonality-output/clonality.xlsx",
                "aggregate_by_patient": True,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": True,
            },
            "archive_runner": {
                "year_label": "2025",
                "input_root": "/tmp/clonality-input",
                "output_root": "/tmp/clonality-output",
                "run_name": "full_2025_validation_test",
                "last_run_root": "/tmp/clonality-output/full_2025_validation_test",
                "combined_workbook_path": "/tmp/clonality-output/full_2025_validation_test/track-clonality-2025-overview.xlsx",
                "max_workers": 1,
                "folder_workers": 1,
                "resume_existing": True,
                "include_sl": False,
                "refresh_each_folder": False,
                "cleanup_staging_root": False,
            },
            "pipeline": {
                "mode": "all",
                "assay_filter_substring": "",
            },
        }
        APP_SETTINGS["analyses"]["general"] = {
            "batch": {
                "base_input_dir": "/tmp/general-input",
                "output_base": "/tmp/general-output",
                "tracking_excel_path": "",
                "aggregate_by_patient": False,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": False,
            },
            "pipeline": {
                "mode": "all",
                "assay_filter_substring": "",
                "ladder": "ROX400HD",
                "trace_channels": ["DATA1"],
                "peak_channels": ["DATA1"],
                "primary_peak_channel": "DATA1",
                "bp_min": 50.0,
                "bp_max": 1000.0,
            },
        }

    def tearDown(self):
        APP_SETTINGS.clear()
        APP_SETTINGS.update(self._settings_backup)

    def test_batch_resolves_saved_output_before_input_folder(self):
        widget = TabBatch()
        widget.output_base.setText("")
        widget.folder_list.clear()
        widget.folder_list.addItem("/tmp/manual-input")

        self.assertEqual(widget._resolve_output_path_str(), "/tmp/clonality-output")

    def test_batch_uses_input_folder_if_no_saved_output_exists(self):
        APP_SETTINGS["analyses"]["clonality"]["batch"]["output_base"] = ""
        widget = TabBatch()
        widget.output_base.setText("")
        widget.folder_list.clear()
        widget.folder_list.addItem("/tmp/manual-input")

        self.assertEqual(widget._resolve_output_path_str(), "/tmp/manual-input")

    def test_batch_input_change_invalidates_scanned_queue(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            original_run = base / "run-a"
            added_run = base / "run-b"
            original_run.mkdir()
            added_run.mkdir()
            first = original_run / "first.fsa"
            second = added_run / "second.fsa"
            first.write_text("ok", encoding="utf-8")
            second.write_text("ok", encoding="utf-8")

            widget = TabBatch()
            widget.folder_list.clear()
            widget._detected_jobs = [
                {"name": "patient_01", "type": "pipeline", "path": original_run, "files": [first]},
            ]
            widget._job_states = {"patient_01": "pending"}
            widget.btn_run.setEnabled(True)
            widget._rebuild_table()

            widget._add_source_item(str(added_run))

        self.assertEqual(widget._detected_jobs, [])
        self.assertEqual(widget._job_states, {})
        self.assertEqual(widget.table.rowCount(), 0)
        self.assertFalse(widget.btn_run.isEnabled())
        self.assertEqual(widget.status_lbl.text(), "Ready")
        self.assertEqual(widget.status_badge.property("state"), "ready")
        self.assertEqual(widget.progress.maximum(), 100)
        self.assertEqual(widget.progress.value(), 0)

    def test_batch_stale_scan_result_is_ignored_while_newer_scan_is_active(self):
        widget = TabBatch()
        widget._detected_jobs = [
            {"name": "older_job", "type": "pipeline", "path": None, "files": [Path("/tmp/a.fsa")]},
        ]
        widget._job_states = {"older_job": "pending"}
        widget._rebuild_table()
        widget.btn_scan.setEnabled(False)
        widget._scan_request_counter = 2
        widget._active_scan_request_id = 2

        widget._on_scan_result(
            [{"name": "stale_job", "type": "pipeline", "path": None, "files": [Path("/tmp/b.fsa")]}],
            request_id=1,
        )

        self.assertEqual([job["name"] for job in widget._detected_jobs], ["older_job"])
        self.assertFalse(widget.btn_scan.isEnabled())

    def test_batch_stale_scan_error_is_ignored_while_newer_scan_is_active(self):
        widget = TabBatch()
        widget.btn_scan.setEnabled(False)
        widget._set_workflow_status("Finding jobs...", "running")
        widget._scan_request_counter = 3
        widget._active_scan_request_id = 3

        widget._on_scan_error((RuntimeError, RuntimeError("stale failure"), None), request_id=2)

        self.assertFalse(widget.btn_scan.isEnabled())
        self.assertEqual(widget.status_lbl.text(), "Finding jobs...")
        self.assertEqual(widget.status_badge.property("state"), "running")

    def test_settings_disables_assay_filter_outside_custom_scope(self):
        widget = TabAnalysisSettings("clonality")
        widget.mode_combo.setCurrentText("all")
        widget._sync_scope_controls()
        self.assertFalse(widget.assay_filter.isEnabled())

        widget.mode_combo.setCurrentText("custom")
        widget._sync_scope_controls()
        self.assertTrue(widget.assay_filter.isEnabled())

    def test_settings_loads_and_saves_tracking_excel_path(self):
        widget = TabAnalysisSettings("clonality")
        self.assertEqual(widget.tracking_excel_path.text(), "/tmp/clonality-output/clonality.xlsx")

        with unittest.mock.patch("gui_qt.tabs.tab_settings.save_settings") as mock_save:
            widget.tracking_excel_path.setText("/tmp/custom/patient-tracking.xlsx")
            widget.save()

        self.assertEqual(
            APP_SETTINGS["analyses"]["clonality"]["batch"]["tracking_excel_path"],
            "/tmp/custom/patient-tracking.xlsx",
        )
        self.assertGreaterEqual(mock_save.call_count, 1)

    def test_batch_general_resolves_file_parent_for_output_path(self):
        APP_SETTINGS["active_analysis"] = "general"
        APP_SETTINGS["analyses"]["general"]["batch"]["output_base"] = ""
        with TemporaryDirectory() as tmp:
            sample = Path(tmp) / "sample.fsa"
            sample.write_text("ok", encoding="utf-8")
            widget = TabBatch()
            widget.output_base.setText("")
            widget.folder_list.clear()
            widget.folder_list.addItem(str(sample))

            self.assertEqual(widget._resolve_output_path_str(), tmp)

    def test_batch_general_persists_runtime_ladder_and_channels(self):
        APP_SETTINGS["active_analysis"] = "general"
        with unittest.mock.patch("gui_qt.tabs.tab_batch.save_settings") as mock_save:
            widget = TabBatch()
            widget.general_ladder_combo.setCurrentIndex(widget.general_ladder_combo.findData("GS500ROX"))
            widget._general_trace_checkboxes["DATA1"].setChecked(True)
            widget._general_trace_checkboxes["DATA2"].setChecked(False)
            widget._general_trace_checkboxes["DATA3"].setChecked(True)
            widget._refresh_general_primary_combo(preferred="DATA3")
            widget._persist_general_runtime_settings()

        profile = APP_SETTINGS["analyses"]["general"]["pipeline"]
        self.assertEqual(profile["ladder"], "GS500ROX")
        self.assertEqual(profile["trace_channels"], ["DATA1", "DATA3"])
        self.assertEqual(profile["peak_channels"], ["DATA1", "DATA3"])
        self.assertEqual(profile["primary_peak_channel"], "DATA3")
        self.assertGreaterEqual(mock_save.call_count, 1)

    def test_batch_general_keeps_one_trace_channel_checked(self):
        APP_SETTINGS["active_analysis"] = "general"
        with unittest.mock.patch("gui_qt.tabs.tab_batch.save_settings"):
            widget = TabBatch()
            widget._general_trace_checkboxes["DATA1"].setChecked(False)

        self.assertTrue(widget._general_trace_checkboxes["DATA1"].isChecked())
        self.assertEqual(widget._selected_general_trace_channels(), ["DATA1"])

    def test_batch_general_uses_compact_trace_checkbox_labels(self):
        APP_SETTINGS["active_analysis"] = "general"
        widget = TabBatch()

        self.assertTrue(widget.subtitle_lbl.wordWrap())
        self.assertEqual(widget._general_trace_checkboxes["DATA1"].text(), "DATA1")
        self.assertEqual(widget._general_trace_checkboxes["DATA2"].text(), "DATA2")
        self.assertEqual(widget._general_trace_checkboxes["DATA3"].text(), "DATA3")

    def test_batch_general_builds_jobs_from_files_and_folders(self):
        APP_SETTINGS["active_analysis"] = "general"
        from core.batch import generate_jobs

        with TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run-a"
            run_dir.mkdir()
            (run_dir / "a.fsa").write_text("ok", encoding="utf-8")
            sample = Path(tmp) / "sample.fsa"
            sample.write_text("ok", encoding="utf-8")

            jobs = generate_jobs([run_dir, sample], aggregate_patients=False, patient_regex=r"\d{2}OUM\d{5}")

        self.assertEqual(len(jobs), 2)
        self.assertEqual(jobs[0]["type"], "pipeline")
        self.assertTrue(jobs[0]["files"])
        self.assertEqual(jobs[1]["files"][0].name, "sample.fsa")

    def test_main_window_wraps_pages_in_scroll_areas(self):
        window = MainWindow()

        for index in range(window.stacked_widget.count()):
            page = window.stacked_widget.widget(index)
            self.assertIsInstance(page, QScrollArea)
            self.assertTrue(page.widgetResizable())

    def test_main_window_exposes_archive_runner_only_for_clonality(self):
        window = MainWindow()

        self.assertEqual(len(window.group_clonality.sub_buttons), 5)
        self.assertEqual(window.group_clonality.sub_buttons[2].text(), "•  Archive Runner")
        self.assertEqual(len(window.group_flt3.sub_buttons), 5)
        self.assertEqual(window.group_flt3.sub_buttons[2].text(), "•  Archive Runner")
        self.assertEqual(len(window.group_general.sub_buttons), 4)
        self.assertEqual(window.stacked_widget.count(), 9)
        self.assertEqual(window.stacked_widget.widget(2).widget().__class__.__name__, "TabArchiveRunner")
        self.assertEqual(window.stacked_widget.widget(5).widget().__class__.__name__, "TabAbout")

    def test_main_window_opens_global_about_page(self):
        window = MainWindow()

        window.on_about_clicked()

        self.assertTrue(window.btn_about.isChecked())
        self.assertEqual(window.stacked_widget.currentIndex(), 5)
        self.assertIsInstance(window.stacked_widget.widget(5).widget(), TabAbout)

    def test_about_tab_surfaces_upstream_license_and_notice(self):
        tab = TabAbout()
        text_dump = "\n".join(widget.toPlainText() for widget in tab.findChildren(QTextBrowser))
        self.assertIn("willros/fraggler", text_dump)
        self.assertIn("MIT License", text_dump)
        self.assertIn("THIRD_PARTY_NOTICES.md", text_dump)

    def test_archive_runner_builds_yearly_config_and_persists_defaults(self):
        widget = TabArchiveRunner()
        widget.year_input.setText("2024")
        widget.input_root.setText("/tmp/archive-input")
        widget.output_root.setText("/tmp/archive-output")
        widget.run_name.setText("full_2024_validation_unit")
        widget.max_workers.setValue(3)
        widget.folder_workers.setValue(2)
        widget.chk_resume.setChecked(False)
        widget.chk_include_sl.setChecked(True)
        widget.chk_refresh_each_folder.setChecked(True)
        widget.chk_cleanup_staging.setChecked(True)
        widget.month_checkboxes["01"].setChecked(True)
        widget.month_checkboxes["02"].setChecked(False)

        config = widget._collect_settings()
        self.assertEqual(config["year_label"], "2024")
        self.assertEqual(widget._selected_month_keys()[0], "2024_01")
        self.assertNotIn("2024_02", widget._selected_month_keys())

        with unittest.mock.patch("gui_qt.tabs.tab_archive_runner.save_settings") as mock_save:
            widget.save_defaults()

        archive = APP_SETTINGS["analyses"]["clonality"]["archive_runner"]
        self.assertEqual(archive["year_label"], "2024")
        self.assertEqual(archive["input_root"], "/tmp/archive-input")
        self.assertEqual(archive["output_root"], "/tmp/archive-output")
        self.assertEqual(archive["run_name"], "full_2024_validation_unit")
        self.assertEqual(archive["max_workers"], 3)
        self.assertEqual(archive["folder_workers"], 2)
        self.assertFalse(archive["resume_existing"])
        self.assertTrue(archive["include_sl"])
        self.assertTrue(archive["refresh_each_folder"])
        self.assertTrue(archive["cleanup_staging_root"])
        self.assertGreaterEqual(mock_save.call_count, 1)

    def test_archive_runner_calls_yearly_backend_with_selected_months(self):
        widget = TabArchiveRunner()
        widget.year_input.setText("2025")
        widget.input_root.setText("/tmp/archive-input")
        widget.output_root.setText("/tmp/archive-output")
        widget.run_name.setText("full_2025_validation_test")
        widget._month_checkboxes["02"].setChecked(False)
        widget._month_checkboxes["12"].setChecked(False)

        captured = {"progress": [], "status": []}

        bridge = type(
            "Bridge",
            (),
            {
                "progress": type("Progress", (), {"emit": lambda self, payload: captured["progress"].append(payload)})(),
                "status": type("Status", (), {"emit": lambda self, payload: captured["status"].append(payload)})(),
            },
        )()

        with unittest.mock.patch("gui_qt.tabs.tab_archive_runner.run_yearly_validation") as mock_run:
            mock_run.return_value = {"run_dir": "/tmp/archive-output/full_2025_validation_test"}
            result = widget._run_yearly_job(
                year_label="2025",
                input_root=Path("/tmp/archive-input"),
                output_root=Path("/tmp/archive-output"),
                run_name="full_2025_validation_test",
                months=widget._selected_month_keys(),
                max_workers=1,
                folder_workers=1,
                resume_existing=True,
                include_sl=False,
                refresh_each_folder=False,
                cleanup_staging_root=False,
                bridge=bridge,
            )

        self.assertEqual(result["run_dir"], "/tmp/archive-output/full_2025_validation_test")
        self.assertEqual(mock_run.call_args.kwargs["year_label"], "2025")
        self.assertEqual(mock_run.call_args.kwargs["months"], widget._selected_month_keys())
        self.assertEqual(mock_run.call_args.kwargs["input_root"], Path("/tmp/archive-input"))
        self.assertEqual(mock_run.call_args.kwargs["output_root"], Path("/tmp/archive-output"))
        self.assertGreaterEqual(len(captured["status"]), 0)

    def test_main_window_activates_analysis_through_shared_helper(self):
        with unittest.mock.patch("gui_qt.main_window.save_settings") as mock_save:
            APP_SETTINGS["batch"] = {"output_base": "/tmp/original-output"}
            APP_SETTINGS["pipeline"] = {
                "mode": "all",
                "assay_filter_substring": "",
                "ladder": "GS500ROX",
            }
            window = MainWindow()

            window.on_group_clicked(window.group_flt3)
            self.assertEqual(APP_SETTINGS["active_analysis"], "flt3")
            self.assertEqual(mock_save.call_count, 1)

            window.on_sub_tab_clicked("general", 1)
            self.assertEqual(APP_SETTINGS["active_analysis"], "general")
            self.assertEqual(mock_save.call_count, 2)
            self.assertEqual(APP_SETTINGS["batch"]["output_base"], "/tmp/general-output")
            self.assertEqual(APP_SETTINGS["pipeline"]["ladder"], "ROX400HD")
            self.assertEqual(window.tab_run._current_analysis_id, "general")
            self.assertEqual(window.tab_ladder._current_analysis_id, "general")

    def test_batch_analysis_switch_resets_queue_state(self):
        widget = TabBatch()
        widget._detected_jobs = [
            {"name": "patient_01", "type": "pipeline", "path": None, "files": [Path("/tmp/a.fsa")]},
        ]
        widget._job_states = {"patient_01": "running"}
        widget.btn_run.setEnabled(True)
        widget._rebuild_table()

        widget.set_analysis("general")

        self.assertEqual(widget._detected_jobs, [])
        self.assertEqual(widget._job_states, {})
        self.assertEqual(widget.table.rowCount(), 0)
        self.assertFalse(widget.btn_run.isEnabled())
        self.assertEqual(widget.status_lbl.text(), "Ready")
        self.assertEqual(widget.status_badge.property("state"), "ready")
        self.assertEqual(widget.progress.maximum(), 100)
        self.assertEqual(widget.progress.value(), 0)


if __name__ == "__main__":
    unittest.main()
