import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
from core.analysis import (
    estimate_running_baseline,
    BASELINE_BIN_SIZE,
    BASELINE_QUANTILE,
    load_ladder_adjustment,
)

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

    def test_load_ladder_adjustment_finds_sidecar_for_staged_symlink(self):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_dir = tmp_path / "source"
            staged_dir = tmp_path / "stage"
            source_dir.mkdir()
            staged_dir.mkdir()

            source_fsa = source_dir / "26OUM04232_ITD_X25__130326_A05_H9H1DI0C.fsa"
            source_fsa.write_text("stub", encoding="utf-8")
            source_adj = source_fsa.with_suffix(".ladder_adj.json")
            source_adj.write_text(json.dumps({"mapping": {"0": 1}}), encoding="utf-8")

            staged_fsa = staged_dir / "00003_e8aea652_26OUM04232_ITD_X25__130326_A05_H9H1DI0C.fsa"
            staged_fsa.symlink_to(source_fsa)

            payload = load_ladder_adjustment(SimpleNamespace(file=staged_fsa))

        self.assertIsNotNone(payload)
        self.assertEqual(payload["mapping"], {0: 1})

if __name__ == "__main__":
    unittest.main()
