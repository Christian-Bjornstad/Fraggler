import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

from core.analysis import analyse_fsa_rox, compute_ladder_qc_metrics, _try_apply_saved_ladder_adjustment


class _DummyFsa:
    def __init__(self):
        self.file = Path("dummy.fsa")
        self.file_name = "dummy.fsa"


class TestRoxLadderRegression(unittest.TestCase):
    def test_invalid_saved_ladder_adjustment_falls_back(self):
        with patch("core.analysis.apply_manual_ladder_mapping", side_effect=ValueError("bad mapping")), \
             patch("core.analysis.print_warning") as mock_warning:
            result = _try_apply_saved_ladder_adjustment(_DummyFsa(), {0: 1}, "ROX")
        self.assertIsNone(result)
        mock_warning.assert_called_once()

    def test_bad_rox_sl_examples_improve(self):
        cases = [
            ("./data/Euroclonality/kjøring 3/24OUM10035_SL__060824_A01_C990GXRS.fsa", 0.9980),
            ("./data/Euroclonality/kjøring 3/24OUM10061_SL__060824_C01_C990GXRS.fsa", 0.9980),
        ]

        for path_str, min_r2 in cases:
            path = Path(path_str)
            if not path.exists():
                self.skipTest(f"Real-data fixture not available: {path}")
            fsa = analyse_fsa_rox(path, "DATA1")
            self.assertIsNotNone(fsa, path.name)
            qc = compute_ladder_qc_metrics(fsa)
            self.assertGreaterEqual(qc["r2"], min_r2, path.name)
            self.assertLess(qc["mean_abs_error_bp"], 5.0, path.name)

    def test_problematic_rox_sample_prefers_full_fit_over_shifted_partial_rescue(self):
        path = Path("./data/Euroclonality/kjøring 3/24OUM10061_SL__060824_C01_C990GXRS.fsa")
        if not path.exists():
            self.skipTest(f"Real-data fixture not available: {path}")

        with patch("core.analysis.load_ladder_adjustment", return_value=None):
            fsa = analyse_fsa_rox(path, "DATA1")
        self.assertIsNotNone(fsa, path.name)
        self.assertEqual(getattr(fsa, "ladder_fit_strategy", None), "auto_full")
        self.assertFalse(bool(getattr(fsa, "ladder_review_required", False)))
        self.assertEqual(list(map(float, getattr(fsa, "ladder_missing_expected_steps", []))), [])
        self.assertEqual(len(getattr(fsa, "expected_ladder_steps", [])), 21)
        self.assertEqual(len(getattr(fsa, "ladder_steps", [])), 21)
        np.testing.assert_allclose(
            np.asarray(getattr(fsa, "best_size_standard", []), dtype=float),
            np.array(
                [
                    1602.0,
                    1656.0,
                    1819.0,
                    1876.0,
                    1986.0,
                    2160.0,
                    2219.0,
                    2334.0,
                    2391.0,
                    2444.0,
                    2568.0,
                    2684.0,
                    2804.0,
                    2926.0,
                    2987.0,
                    3042.0,
                    3169.0,
                    3291.0,
                    3414.0,
                    3536.0,
                    3656.0,
                ],
                dtype=float,
            ),
        )
        qc = compute_ladder_qc_metrics(fsa)
        self.assertGreaterEqual(qc["r2"], 0.99998)
        self.assertLessEqual(qc["max_abs_error_bp"], 1.0)

    def test_good_rox_sl_examples_stay_stable(self):
        cases = [
            "./data/Euroclonality/kjøring 2/25OUM10166_SL__30062025_A09_H9C0ZIZO.fsa",
            "./data/Euroclonality/kjøring 3/24OUM10037_SL__060824_B01_C990GXRS.fsa",
        ]

        for path_str in cases:
            path = Path(path_str)
            if not path.exists():
                self.skipTest(f"Real-data fixture not available: {path}")
            fsa = analyse_fsa_rox(path, "DATA1")
            self.assertIsNotNone(fsa, path.name)
            qc = compute_ladder_qc_metrics(fsa)
            self.assertGreaterEqual(qc["r2"], 0.9999, path.name)
            self.assertEqual(
                qc["n_size_standard_peaks"],
                len(getattr(fsa, "ladder_steps", [])),
                path.name,
            )


if __name__ == "__main__":
    unittest.main()
