import copy
import unittest
from unittest.mock import patch
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from core.analysis import (
    _select_best_ladder_candidate,
    _candidate_fit_score,
    _rescue_fit_score,
    _clean_rox_size_standard_peaks,
    apply_manual_ladder_mapping,
    compute_ladder_qc_metrics,
    get_ladder_candidates,
    load_ladder_adjustment,
    save_ladder_adjustment,
)
from gui_qt.ladder_utils import detect_fsa_for_ladder, load_adjustable_fsa
from config import APP_SETTINGS


class DummyModel:
    def predict(self, values):
        arr = np.asarray(values, dtype=float).reshape(-1)
        return arr


class DummyFsa:
    def __init__(self):
        self.ladder_steps = np.array([100.0, 200.0, 300.0, 400.0], dtype=float)
        self.expected_ladder_steps = np.array([100.0, 200.0, 300.0, 400.0], dtype=float)
        self.best_size_standard = np.array([10.0, 20.0, 30.0, 40.0], dtype=float)
        self.best_size_standard_combinations = pd.DataFrame(
            {
                "combinations": [
                    np.array([10.0, 20.0, 30.0, 40.0], dtype=float),
                    np.array([11.0, 21.0, 31.0, 41.0], dtype=float),
                ]
            }
        )
        self.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": [10, 20, 30, 40],
                "basepairs": [100.0, 200.0, 300.0, 400.0],
            }
        )
        self.ladder_model = DummyModel()
        self.fitted_to_model = False


class TestLadderSelection(unittest.TestCase):
    def test_select_best_ladder_candidate_prefers_lower_fit_error(self):
        first = DummyFsa()
        second = DummyFsa()

        def fake_fit(trial):
            trial.fitted_to_model = True
            return trial

        fit_scores = [
            {"r2": 0.9997, "mean_abs_error_bp": 0.7, "max_abs_error_bp": 1.2},
            {"r2": 0.9996, "mean_abs_error_bp": 0.15, "max_abs_error_bp": 0.3},
        ]

        with patch("core.analysis.copy.deepcopy", side_effect=[first, second]), \
            patch("core.analysis.fit_size_standard_to_ladder", side_effect=fake_fit), \
            patch("core.analysis.compute_ladder_qc_metrics", side_effect=fit_scores):
            selected = _select_best_ladder_candidate(DummyFsa())

        self.assertIs(selected, second)

    def test_compute_ladder_qc_metrics_returns_error_metrics(self):
        fsa = DummyFsa()
        metrics = compute_ladder_qc_metrics(fsa)
        self.assertIn("r2", metrics)
        self.assertIn("mean_abs_error_bp", metrics)
        self.assertIn("max_abs_error_bp", metrics)
        self.assertEqual(metrics["mean_abs_error_bp"], 0.0)
        self.assertEqual(metrics["max_abs_error_bp"], 0.0)

    def test_load_ladder_adjustment_reads_json_mapping(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp) / "example.fsa"
            base.write_text("")
            adj = base.with_suffix(".ladder_adj.json")
            adj.write_text('{"mapping": {"0": 1, "2": 3}, "mapping_times": {"0": 10.0, "2": 30.0}, "manual_candidates": [55.0]}', encoding="utf-8")

            dummy = type("Dummy", (), {"file": base})()
            mapping = load_ladder_adjustment(dummy)

        self.assertEqual(mapping["mapping"], {0: 1, 2: 3})
        self.assertEqual(mapping["mapping_times"], {0: 10.0, 2: 30.0})
        self.assertEqual(mapping["manual_candidates"], [55.0])

    def test_apply_manual_ladder_mapping_accepts_mapping_times(self):
        fsa = DummyFsa()
        fsa.file = Path("dummy.fsa")
        fsa.size_standard_peaks = np.array([10.0, 20.0], dtype=float)
        payload = {
            "mapping": {0: 0, 1: 1, 2: 2, 3: 3},
            "mapping_times": {0: 10.0, 1: 20.0, 2: 33.0, 3: 44.0},
            "manual_candidates": [33.0, 44.0],
        }

        def fake_fit(trial):
            trial.fitted_to_model = True
            return trial

        with patch("core.analysis.fit_size_standard_to_ladder", side_effect=fake_fit):
            adjusted = apply_manual_ladder_mapping(copy.deepcopy(fsa), payload)

        np.testing.assert_allclose(adjusted.best_size_standard, np.array([10.0, 20.0, 33.0, 44.0]))
        self.assertTrue(adjusted.fitted_to_model)

    def test_apply_manual_ladder_mapping_expands_trimmed_fit_to_expected_ladder(self):
        fsa = DummyFsa()
        fsa.file = Path("dummy.fsa")
        fsa.expected_ladder_steps = np.array([100.0, 150.0, 200.0, 250.0, 300.0], dtype=float)
        fsa.ladder_steps = np.array([100.0, 200.0, 300.0], dtype=float)
        fsa.best_size_standard = np.array([10.0, 20.0, 30.0], dtype=float)
        fsa.size_standard_peaks = np.array([10.0, 20.0, 30.0], dtype=float)
        payload = {
            "mapping": {},
            "mapping_times": {1: 15.0, 3: 25.0},
            "manual_candidates": [15.0, 25.0],
        }

        def fake_fit(trial):
            trial.fitted_to_model = True
            return trial

        with patch("core.analysis.fit_size_standard_to_ladder", side_effect=fake_fit):
            adjusted = apply_manual_ladder_mapping(copy.deepcopy(fsa), payload)

        np.testing.assert_allclose(adjusted.ladder_steps, np.array([100.0, 150.0, 200.0, 250.0, 300.0]))
        np.testing.assert_allclose(adjusted.best_size_standard, np.array([10.0, 15.0, 20.0, 25.0, 30.0]))
        self.assertTrue(adjusted.fitted_to_model)

    def test_get_ladder_candidates_handles_float_peak_times_from_saved_adjustments(self):
        fsa = DummyFsa()
        fsa.size_standard = np.zeros(64, dtype=float)
        fsa.size_standard[10] = 900.0
        fsa.size_standard[20] = 1250.0
        fsa.size_standard[30] = 1100.0
        fsa.size_standard_peaks = np.array([10.0, 20.0, 30.0], dtype=float)
        fsa.manual_ladder_candidates = [20.0]

        df = get_ladder_candidates(fsa)

        self.assertEqual(df["intensity"].tolist(), [900.0, 1250.0, 1100.0])
        self.assertEqual(df["source"].tolist(), ["auto", "manual", "auto"])

    def test_candidate_fit_score_penalizes_weak_selected_ladder_peaks(self):
        strong = DummyFsa()
        strong.size_standard = np.zeros(64, dtype=float)
        strong.size_standard[[10, 20, 30, 40]] = [1000.0, 980.0, 1020.0, 995.0]

        weak = DummyFsa()
        weak.size_standard = np.zeros(64, dtype=float)
        weak.size_standard[[10, 20, 30, 40]] = [1000.0, 980.0, 250.0, 995.0]

        strong_score = _candidate_fit_score(strong)
        weak_score = _candidate_fit_score(weak)

        self.assertLess(strong_score, weak_score)

    def test_candidate_fit_score_penalizes_weak_early_peaks_more_than_late_peaks(self):
        weak_early = DummyFsa()
        weak_early.size_standard = np.zeros(64, dtype=float)
        weak_early.size_standard[[10, 20, 30, 40]] = [260.0, 980.0, 1020.0, 995.0]

        weak_late = DummyFsa()
        weak_late.size_standard = np.zeros(64, dtype=float)
        weak_late.size_standard[[10, 20, 30, 40]] = [1000.0, 980.0, 1020.0, 260.0]

        self.assertLess(_candidate_fit_score(weak_late), _candidate_fit_score(weak_early))

    def test_rescue_fit_score_penalizes_missing_expected_ladder_steps(self):
        full = DummyFsa()
        full.size_standard = np.zeros(64, dtype=float)
        full.size_standard[[10, 20, 30, 40]] = [1000.0, 980.0, 1020.0, 995.0]

        partial = DummyFsa()
        partial.expected_ladder_steps = np.array([100.0, 200.0, 300.0, 400.0, 500.0], dtype=float)
        partial.ladder_steps = np.array([200.0, 300.0, 400.0, 500.0], dtype=float)
        partial.best_size_standard = np.array([10.0, 20.0, 30.0, 40.0], dtype=float)
        partial.size_standard = np.zeros(64, dtype=float)
        partial.size_standard[[10, 20, 30, 40]] = [1000.0, 980.0, 1020.0, 995.0]

        self.assertLess(_rescue_fit_score(full), _rescue_fit_score(partial))

    def test_clean_rox_size_standard_peaks_keeps_isolated_tall_ladder_peaks(self):
        peaks = np.array([1379, 1395, 1422, 1819, 1855, 1876, 2391, 2444, 2568, 2987, 3042, 3169])
        rox = np.zeros(4000, dtype=float)
        heights = {
            1379: 300.0,
            1395: 25351.0,
            1422: 260.0,
            1819: 1028.0,
            1855: 11536.0,
            1876: 942.0,
            2391: 1089.0,
            2444: 6333.0,
            2568: 1098.0,
            2987: 1092.0,
            3042: 2396.0,
            3169: 1000.0,
        }
        for idx, value in heights.items():
            rox[idx] = value

        cleaned = _clean_rox_size_standard_peaks(peaks, rox)

        self.assertNotIn(1395, cleaned.tolist())
        self.assertNotIn(1855, cleaned.tolist())
        self.assertIn(2444, cleaned.tolist())
        self.assertIn(3042, cleaned.tolist())

    def test_load_ladder_adjustment_does_not_fall_back_to_sibling_project_copy(self):
        with TemporaryDirectory() as tmp:
            desktop = Path(tmp) / "Desktop"
            root_a = desktop / "OUS"
            root_b = desktop / "OUS-kopi"
            rel = Path("data/Euroclonality/run/sample.fsa")
            file_a = root_a / rel
            file_b = root_b / rel
            file_a.parent.mkdir(parents=True, exist_ok=True)
            file_b.parent.mkdir(parents=True, exist_ok=True)
            file_a.write_text("")
            file_b.write_text("")
            file_a.with_suffix(".ladder_adj.json").write_text(
                '{"mapping": {"0": 1}, "mapping_times": {"0": 10.0}, "manual_candidates": [12.0]}',
                encoding="utf-8",
            )

            dummy = type("Dummy", (), {"file": file_b})()
            payload = load_ladder_adjustment(dummy)

        self.assertIsNone(payload)

    def test_ladder_utils_prefer_clonality_runtime_even_if_active_analysis_is_flt3(self):
        path = Path("/tmp/sample.fsa")
        clonality_classification = (
            "IGH",
            "sample",
            "ROX",
            ["DATA1"],
            ["DATA1"],
            "DATA1",
            50.0,
            500.0,
        )
        flt3_classification = {
            "assay": "FLT3-D835",
            "group": "sample",
            "ladder": "ROX",
            "trace_channels": ["DATA2"],
            "peak_channels": ["DATA2"],
            "primary_peak_channel": "DATA2",
            "bp_min": 50.0,
            "bp_max": 250.0,
        }
        sentinel_fsa = object()

        original_analysis = APP_SETTINGS.get("active_analysis")
        APP_SETTINGS["active_analysis"] = "flt3"
        try:
            with patch("gui_qt.ladder_utils.classify_clonality_fsa", return_value=clonality_classification), \
                 patch("gui_qt.ladder_utils.classify_flt3_fsa", return_value=flt3_classification) as mock_flt3_classify, \
                 patch("gui_qt.ladder_utils.analyse_fsa_rox", return_value=sentinel_fsa) as mock_analyse:
                meta = detect_fsa_for_ladder(path, preferred_analysis="clonality")
                fsa, refreshed_meta = load_adjustable_fsa(path, preferred_analysis="clonality")

        finally:
            if original_analysis is None:
                APP_SETTINGS.pop("active_analysis", None)
            else:
                APP_SETTINGS["active_analysis"] = original_analysis

        self.assertEqual(meta["analysis"], "clonality")
        self.assertEqual(meta["ladder"], "ROX400HD")
        self.assertIs(fsa, sentinel_fsa)
        self.assertEqual(refreshed_meta["analysis"], "clonality")
        self.assertEqual(refreshed_meta["ladder"], "ROX400HD")
        mock_flt3_classify.assert_not_called()
        mock_analyse.assert_called_once_with(
            path,
            "DATA1",
            ladder_name="ROX400HD",
            min_distance_between_peaks=15.0,
            min_size_standard_height=200.0,
        )

    def test_load_adjustable_fsa_reuses_precomputed_metadata(self):
        path = Path("/tmp/sample.fsa")
        metadata = {
            "analysis": "clonality",
            "assay": "IGH",
            "group": "sample",
            "ladder": "ROX400HD",
            "trace_channels": ["DATA1"],
            "peak_channels": ["DATA1"],
            "primary_peak_channel": "DATA1",
            "bp_min": 50.0,
            "bp_max": 500.0,
            "sample_channel": "DATA1",
        }
        sentinel_fsa = object()

        with patch("gui_qt.ladder_utils.detect_fsa_for_ladder") as mock_detect, \
             patch("gui_qt.ladder_utils.analyse_fsa_rox", return_value=sentinel_fsa):
            fsa, refreshed_meta = load_adjustable_fsa(path, preferred_analysis="clonality", metadata=metadata.copy())

        self.assertIs(fsa, sentinel_fsa)
        self.assertEqual(refreshed_meta["ladder"], "ROX400HD")
        mock_detect.assert_not_called()

    def test_save_ladder_adjustment_stays_with_original_file(self):
        with TemporaryDirectory() as tmp:
            desktop = Path(tmp) / "Desktop"
            root_a = desktop / "OUS"
            root_b = desktop / "OUS-kopi"
            rel = Path("data/Euroclonality/run/sample.fsa")
            file_a = root_a / rel
            file_b = root_b / rel
            file_a.parent.mkdir(parents=True, exist_ok=True)
            file_b.parent.mkdir(parents=True, exist_ok=True)
            file_a.write_text("")
            file_b.write_text("")

            dummy = type("Dummy", (), {"file": file_a})()
            save_ladder_adjustment(dummy, {"mapping": {0: 1}, "mapping_times": {0: 10.0}, "manual_candidates": [12.0]})

            mirrored = file_b.with_suffix(".ladder_adj.json")
            self.assertFalse(mirrored.exists())
            self.assertTrue(file_a.with_suffix(".ladder_adj.json").exists())


if __name__ == "__main__":
    unittest.main()
