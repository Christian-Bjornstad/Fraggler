import unittest
from pathlib import Path

import pandas as pd

from core.html_reports import _build_plotly_reflow_script, _create_html_header
from core.plotting_plotly import _create_plotly_figure, _prepare_plot_data, build_interactive_peak_plot_for_entry


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

    def test_flt3_itd_plot_exposes_manual_ratio_panel(self):
        fsa = _DummyFsa()
        fsa.fsa["DATA2"] = [0.0, 80.0, 140.0, 90.0, 0.0]
        entry = {
            "fsa": fsa,
            "primary_peak_channel": "DATA1",
            "trace_channels": ["DATA1", "DATA2"],
            "bp_min": 300.0,
            "bp_max": 360.0,
            "assay": "FLT3-ITD",
            "peaks_by_channel": {"DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "label"])},
            "wt_bp": 330.0,
            "mut_bp": 335.0,
        }

        html = build_interactive_peak_plot_for_entry(entry)

        self.assertIsNotNone(html)
        self.assertIn("_flt3_panel", html)
        self.assertIn("_flt3_table", html)
        self.assertIn("flt3_manual_ratio_selection", html)
        self.assertIn("getPeakData", html)
        self.assertIn("Blue area", html)
        self.assertIn("Green area", html)
        self.assertIn("Neste klikk", html)
        self.assertIn("Ratio</div>", html)
        self.assertIn("Nullstill valg", html)
        self.assertIn("radio-knappen", html)
        self.assertIn("data-role='wt'", html)
        self.assertNotIn("Auto ratio", html)
        # Status card should NOT be present
        self.assertNotIn("flt3_ratio_status", html)

    def test_flt3_itd_plot_tracks_manual_peak_ids_and_channel_specific_overlap_clicks(self):
        fsa = _DummyFsa()
        fsa.fsa["DATA2"] = [0.0, 80.0, 140.0, 90.0, 0.0]
        entry = {
            "fsa": fsa,
            "primary_peak_channel": "DATA1",
            "trace_channels": ["DATA1", "DATA2"],
            "bp_min": 300.0,
            "bp_max": 360.0,
            "assay": "FLT3-ITD",
            "peaks_by_channel": {"DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "label"])},
            "wt_bp": 330.0,
            "mut_bp": 335.0,
        }

        html = build_interactive_peak_plot_for_entry(entry)

        self.assertIsNotNone(html)
        self.assertIn("function peakIdFor(peak, idx)", html)
        self.assertIn("ensurePeakIds()", html)
        self.assertIn("manualSelection.mutant_peak_ids.indexOf(peakIdFor(peak, idx)) >= 0", html)
        self.assertIn("if (preferredChannel && peak.source_channel !== preferredChannel) continue;", html)
        self.assertIn("newPeak.peak_id = makePeakId(newPeak, peaks.length);", html)
        # WT selection state is persisted and restored
        self.assertIn("wt_peak_ids", html)
        self.assertIn("isManualWt", html)

    def test_flt3_d835_plot_exposes_manual_ratio_panel(self):
        fsa = _DummyFsa()
        fsa.fsa = {
            "DATA3": [0.0, 120.0, 260.0, 150.0, 0.0],
        }
        entry = {
            "fsa": fsa,
            "primary_peak_channel": "DATA3",
            "trace_channels": ["DATA3"],
            "bp_min": 70.0,
            "bp_max": 160.0,
            "assay": "FLT3-D835",
            "peaks_by_channel": {"DATA3": pd.DataFrame(columns=["basepairs", "peaks", "area", "label"])},
            "wt_bp": 80.0,
            "mut_bp": 129.0,
        }

        html = build_interactive_peak_plot_for_entry(entry)

        self.assertIsNotNone(html)
        self.assertIn("_flt3_panel", html)
        self.assertIn("Nullstill valg", html)
        self.assertIn("WT area", html)
        self.assertIn("Mut area", html)
        self.assertIn("Area</th><th>Status</th><th>Fjern</th>", html)
        self.assertIn("data-role='wt'", html)
        self.assertIn("flt3_manual_ratio_selection", html)

    def test_peak_manager_script_accepts_old_and_new_peak_payloads(self):
        html_lines = []
        _create_html_header("25OUM10166", 2025, 1, Path("/tmp"), html_lines, display_name="Flt3")
        html = "\n".join(html_lines)

        self.assertIn("getInitialPeakDataForPlot", html)
        self.assertIn("_normalizePeakPayload", html)
        self.assertIn("getAllPeakData", html)
        self.assertIn("flt3_manual_ratio_selection", html)
        self.assertIn("Array.isArray(payload)", html)

    def test_flt3_negative_control_plot_uses_minimum_ymax_of_250(self):
        fsa = _DummyFsa()
        entry = {
            "fsa": fsa,
            "primary_peak_channel": "DATA1",
            "trace_channels": ["DATA1"],
            "bp_min": 70.0,
            "bp_max": 90.0,
            "assay": "FLT3-D835",
            "group": "negative_control",
            "peaks_by_channel": {"DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "label"])},
            "wt_bp": 80.0,
            "mut_bp": 129.0,
        }

        data = _prepare_plot_data(entry)
        self.assertIsNotNone(data)

        _, ymax, _ = _create_plotly_figure(data)

        self.assertEqual(ymax, 250.0)

    def test_flt3_itd_plot_uses_all_channels_for_auto_ymax(self):
        fsa = _DummyFsa()
        fsa.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": [0, 1, 2, 3, 4],
                "basepairs": [328.0, 329.0, 330.0, 331.0, 332.0],
            }
        )
        fsa.fsa["DATA1"] = [0.0, 800.0, 1200.0, 700.0, 0.0]
        fsa.fsa["DATA2"] = [0.0, 1500.0, 4000.0, 1800.0, 0.0]
        entry = {
            "fsa": fsa,
            "primary_peak_channel": "DATA1",
            "trace_channels": ["DATA1", "DATA2"],
            "bp_min": 300.0,
            "bp_max": 360.0,
            "assay": "FLT3-ITD",
            "peaks_by_channel": {"DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "label"])},
            "wt_bp": 330.0,
            "mut_bp": 335.0,
        }

        data = _prepare_plot_data(entry)
        self.assertIsNotNone(data)

        _, ymax, _ = _create_plotly_figure(data)

        self.assertGreaterEqual(ymax, 4000.0)



if __name__ == "__main__":
    unittest.main()
