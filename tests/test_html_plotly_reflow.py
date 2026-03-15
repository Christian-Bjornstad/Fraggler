import unittest

import pandas as pd

from core.html_reports import _build_plotly_reflow_script
from core.plotting_plotly import build_interactive_peak_plot_for_entry


class _DummyFsa:
    def __init__(self):
        self.file_name = "dummy_A01_test.fsa"
        self.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": [0, 1, 2, 3, 4],
                "basepairs": [78.0, 79.0, 80.0, 81.0, 82.0],
            }
        )
        self.fsa = {
            "DATA1": [0.0, 120.0, 250.0, 110.0, 0.0],
        }


class TestHtmlPlotlyReflow(unittest.TestCase):
    def test_reflow_script_contains_visibility_relayout_hooks(self):
        script = _build_plotly_reflow_script()
        self.assertIn("window.ReportPlotManager", script)
        self.assertIn("pageshow", script)
        self.assertIn("visibilitychange", script)
        self.assertIn("ResizeObserver", script)
        self.assertIn("Plotly.Plots.resize", script)
        self.assertIn("getInitialStateForPlot", script)
        self.assertIn("getAllStates", script)

    def test_interactive_plot_registers_with_report_plot_manager(self):
        fsa = _DummyFsa()
        entry = {
            "fsa": fsa,
            "primary_peak_channel": "DATA1",
            "trace_channels": ["DATA1"],
            "bp_min": 70.0,
            "bp_max": 90.0,
            "assay": "FLT3-D835",
            "peaks_by_channel": {"DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "label"])},
            "wt_bp": 80.0,
            "mut_bp": 129.0,
        }

        html = build_interactive_peak_plot_for_entry(entry)
        self.assertIsNotNone(html)
        self.assertIn("window.ReportPlotManager.register(g)", html)
        self.assertIn("window.ReportPlotManager.getInitialStateForPlot(divId)", html)


if __name__ == "__main__":
    unittest.main()
