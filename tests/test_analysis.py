import unittest
import numpy as np
from core.analysis import estimate_running_baseline, BASELINE_BIN_SIZE, BASELINE_QUANTILE

class TestAnalysis(unittest.TestCase):
    def test_estimate_running_baseline_zeros(self):
        trace = np.zeros(1000)
        baseline = estimate_running_baseline(trace)
        self.assertTrue(np.allclose(baseline, 0.0))

    def test_estimate_running_baseline_constant(self):
        trace = np.ones(1000) * 100.0
        baseline = estimate_running_baseline(trace)
        self.assertTrue(np.allclose(baseline, 100.0))

    def test_estimate_running_baseline_simple_slope(self):
        trace = np.linspace(0, 100, 1000)
        baseline = estimate_running_baseline(trace, bin_size=100, quantile=0.0)
        self.assertEqual(baseline[0], 0.0)
        self.assertGreater(baseline[-1], 80.0)

    def test_estimate_running_baseline_empty(self):
        trace = np.array([])
        baseline = estimate_running_baseline(trace)
        self.assertEqual(baseline.size, 0)

if __name__ == "__main__":
    unittest.main()
