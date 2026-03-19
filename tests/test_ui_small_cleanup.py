import copy
import os
import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication, QScrollArea

from config import APP_SETTINGS
from gui_qt.main_window import MainWindow
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
                "aggregate_by_patient": True,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": True,
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

    def test_settings_disables_assay_filter_outside_custom_scope(self):
        widget = TabAnalysisSettings("clonality")
        widget.mode_combo.setCurrentText("all")
        widget._sync_scope_controls()
        self.assertFalse(widget.assay_filter.isEnabled())

        widget.mode_combo.setCurrentText("custom")
        widget._sync_scope_controls()
        self.assertTrue(widget.assay_filter.isEnabled())

    def test_batch_general_resolves_file_parent_for_output_path(self):
        APP_SETTINGS["active_analysis"] = "general"
        APP_SETTINGS["analyses"]["general"]["batch"]["output_base"] = ""
        with TemporaryDirectory() as tmp:
            sample = Path(tmp) / "sample.fsa"
            sample.write_text("", encoding="utf-8")
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
        widget = TabBatch()
        with TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run-a"
            run_dir.mkdir()
            (run_dir / "a.fsa").write_text("", encoding="utf-8")
            sample = Path(tmp) / "sample.fsa"
            sample.write_text("", encoding="utf-8")

            jobs = widget._build_general_jobs_worker([run_dir, sample])

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


if __name__ == "__main__":
    unittest.main()
