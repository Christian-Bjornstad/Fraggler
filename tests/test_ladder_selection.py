import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from core.analysis import _select_best_ladder_candidate


class DummyModel:
    def predict(self, values):
        arr = np.asarray(values, dtype=float).reshape(-1)
        return arr


class DummyFsa:
    def __init__(self):
        self.ladder_steps = np.array([100.0, 200.0, 300.0, 400.0], dtype=float)
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


if __name__ == "__main__":
    unittest.main()
