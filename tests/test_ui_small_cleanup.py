import copy
import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication

from config import APP_SETTINGS
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


if __name__ == "__main__":
    unittest.main()
