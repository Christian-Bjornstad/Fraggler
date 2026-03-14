import unittest

import numpy as np
import pandas as pd

from core.qc.qc_markers import find_peak_near_bp_with_fallback


class DummyFsa:
    def __init__(self, bp_values, trace_values):
        self.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": np.arange(len(bp_values)),
                "basepairs": np.asarray(bp_values, dtype=float),
            }
        )
        self.fsa = {
            "DATA1": np.asarray(trace_values, dtype=float),
        }


class TestQcMarkers(unittest.TestCase):
    def test_fallback_window_recovers_peak_outside_primary_window(self):
        bp_values = [257.0, 258.0, 259.0, 260.0, 261.0, 262.4, 263.0, 264.0]
        trace_values = [10.0, 12.0, 15.0, 18.0, 16.0, 250.0, 40.0, 30.0]
        fsa = DummyFsa(bp_values, trace_values)

        result = find_peak_near_bp_with_fallback(
            fsa=fsa,
            channel="DATA1",
            target_bp=260.0,
            window_bp=2.0,
            fallback_window_bp=4.0,
            baseline_correct=False,
        )

        self.assertTrue(result["ok"])
        self.assertEqual(result["search_mode"], "fallback")
        self.assertAlmostEqual(result["found_bp"], 262.4)
        self.assertAlmostEqual(result["search_window_bp"], 4.0)


if __name__ == "__main__":
    unittest.main()
