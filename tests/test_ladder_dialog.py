import copy
import os
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import pandas as pd
from PyQt6.QtWidgets import QApplication

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from gui_qt.dialogs.ladder_dialog import LadderAdjustmentDialog


_APP = QApplication.instance() or QApplication([])


class _DummyModel:
    def predict(self, values):
        arr = np.asarray(values, dtype=float).reshape(-1)
        return arr


class _DummyFsa:
    def __init__(self):
        self.file_name = "dummy.fsa"
        self.ladder = "ROX"
        self.ladder_steps = np.array([100.0, 200.0, 300.0], dtype=float)
        self.expected_ladder_steps = np.array([100.0, 200.0, 300.0], dtype=float)
        self.best_size_standard = np.array([10, 20, 30], dtype=float)
        self.size_standard_peaks = np.array([10, 20, 30], dtype=int)
        self.size_standard = np.zeros(64, dtype=float)
        self.size_standard[10] = 900.0
        self.size_standard[20] = 1250.0
        self.size_standard[30] = 1100.0
        self.size_standard[42] = 1800.0
        self.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": [10, 20, 30, 42],
                "basepairs": [100.1, 199.9, 300.2, 399.8],
            }
        )
        self.ladder_model = _DummyModel()


class TestLadderAdjustmentDialog(unittest.TestCase):
    def setUp(self):
        self.fsa = _DummyFsa()
        self.candidates = pd.DataFrame(
            {
                "index": [0, 1, 2],
                "time": [10.0, 20.0, 30.0],
                "intensity": [900.0, 1250.0, 1100.0],
            }
        )

    def _build_dialog(self):
        preview_metrics = {
            "r2": 0.9997,
            "mean_abs_error_bp": 0.13,
            "max_abs_error_bp": 0.2,
            "n_ladder_steps": 3,
            "n_size_standard_peaks": 3,
        }
        with patch("core.analysis.get_ladder_candidates", return_value=self.candidates.copy()), \
            patch("core.analysis.apply_manual_ladder_mapping", side_effect=lambda fsa, mapping: fsa), \
            patch("core.analysis.compute_ladder_qc_metrics", return_value=preview_metrics):
            return LadderAdjustmentDialog(copy.deepcopy(self.fsa))

    def test_dialog_builds_match_table_and_qc_summary(self):
        dialog = self._build_dialog()
        self.assertEqual(dialog.table.columnCount(), 6)
        self.assertEqual(dialog.table.horizontalHeaderItem(3).text(), "Residual")
        self.assertEqual(dialog.qc_grade_label.text(), "PASS")
        self.assertIn("R² 0.999700", dialog.qc_summary_label.text())
        self.assertEqual(dialog.table.item(0, 2).text(), "Auto #0")

    def test_clearing_step_marks_row_missing_and_downgrades_preview(self):
        dialog = self._build_dialog()
        dialog.table.selectRow(1)
        dialog._clear_selected_step()

        self.assertEqual(dialog.table.item(1, 5).text(), "Missing")
        self.assertEqual(dialog.qc_grade_label.text(), "CHECK")
        self.assertIn("missing", dialog.qc_summary_label.text().lower())

    def test_add_manual_peak_from_plot_builds_payload_and_candidate(self):
        dialog = self._build_dialog()
        dialog.table.selectRow(2)
        dialog._add_manual_peak_from_plot(42.0, assign_to_step=2)

        payload = dialog.get_adjustment_payload()
        self.assertTrue(payload["manual_candidates"])
        self.assertEqual(dialog.candidates.iloc[-1]["source"], "manual")
        self.assertEqual(dialog.table.item(2, 2).text().split()[0], "Manual")

    def test_add_missing_mode_selects_next_missing_step(self):
        dialog = self._build_dialog()
        dialog.mapping = {0: 0}
        dialog._refresh_preview_state(show_errors=False)
        dialog._refresh_all()

        dialog.btn_add_peak.setChecked(True)

        self.assertEqual(dialog._selected_step_row(), 1)
        self.assertIn("200 bp", dialog.stats_label.text())
        self.assertIn("200 bp", dialog.missing_steps_label.text())
        self.assertIn("300 bp", dialog.missing_steps_label.text())
        self.assertEqual(dialog.missing_list.count(), 2)

    def test_dialog_uses_full_expected_ladder_when_fit_was_trimmed(self):
        trimmed_fsa = _DummyFsa()
        trimmed_fsa.ladder_steps = np.array([100.0, 200.0, 300.0], dtype=float)
        trimmed_fsa.expected_ladder_steps = np.array([100.0, 150.0, 200.0, 250.0, 300.0], dtype=float)
        trimmed_fsa.best_size_standard = np.array([10, 20, 30], dtype=float)
        candidates = pd.DataFrame(
            {
                "index": [0, 1, 2],
                "time": [10.0, 20.0, 30.0],
                "intensity": [900.0, 1250.0, 1100.0],
            }
        )
        preview_metrics = {
            "r2": 0.9997,
            "mean_abs_error_bp": 0.13,
            "max_abs_error_bp": 0.2,
            "n_ladder_steps": 5,
            "n_size_standard_peaks": 5,
        }

        with patch("core.analysis.get_ladder_candidates", return_value=candidates.copy()), \
            patch("core.analysis.apply_manual_ladder_mapping", side_effect=lambda fsa, mapping: fsa), \
            patch("core.analysis.compute_ladder_qc_metrics", return_value=preview_metrics):
            dialog = LadderAdjustmentDialog(copy.deepcopy(trimmed_fsa))

        self.assertEqual(dialog.meta_labels["expected_count"].text(), "5")
        self.assertEqual(dialog.table.rowCount(), 5)
        self.assertEqual(dialog.missing_list.count(), 2)
        self.assertIn("150 bp", dialog.missing_steps_label.text())
        self.assertIn("250 bp", dialog.missing_steps_label.text())
        self.assertEqual(dialog.table.item(0, 2).text(), "Auto #0")
        self.assertEqual(dialog.table.item(2, 2).text(), "Auto #1")
        self.assertEqual(dialog.table.item(4, 2).text(), "Auto #2")

    def test_plot_click_creates_manual_candidate_when_peak_is_not_detected(self):
        dialog = self._build_dialog()
        dialog.table.selectRow(2)

        dialog._on_plot_click(SimpleNamespace(inaxes=dialog.ax, xdata=42.0))

        self.assertEqual(dialog.candidates.iloc[-1]["source"], "manual")
        self.assertIn(2, dialog.mapping)
        mapped_idx = dialog.mapping[2]
        self.assertEqual(dialog.candidates.iloc[mapped_idx]["time"], 42.0)

    def test_apply_requires_complete_previewable_mapping(self):
        dialog = self._build_dialog()
        dialog.mapping = {0: 0, 1: 1}
        dialog._refresh_preview_state(show_errors=False)
        dialog._refresh_all()

        with patch("gui_qt.dialogs.ladder_dialog.QMessageBox.warning") as warning:
            dialog._on_apply()

        self.assertFalse(dialog.result())
        warning.assert_called()

    def test_missing_order_toggle_supports_descending_workflow(self):
        dialog = self._build_dialog()
        dialog.mapping = {0: 0}
        dialog._refresh_preview_state(show_errors=False)
        dialog._refresh_all()

        dialog.btn_missing_order.setChecked(True)

        self.assertEqual(dialog._missing_step_indices(), [2, 1])
        self.assertEqual(dialog.missing_list.item(0).text(), "300 bp")


if __name__ == "__main__":
    unittest.main()
