import unittest
from pathlib import Path
from unittest.mock import patch

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
            self.assertEqual(qc["n_size_standard_peaks"], 20, path.name)


if __name__ == "__main__":
    unittest.main()
