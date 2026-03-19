import copy
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pandas as pd

from config import APP_SETTINGS
from core.analyses.general.classification import classify_fsa
from core.analyses.general.config import GENERAL_ASSAY_NAME, normalize_ladder_name, normalize_trace_channels, resolve_runtime_config
from core.analyses.general.pipeline import run_pipeline
from core.analyses.general.reporting import build_general_html_report


class DummyModel:
    def predict(self, values):
        arr = np.asarray(values, dtype=float).reshape(-1)
        return arr


class DummyFsa:
    def __init__(self, file_name: str = "sample.fsa"):
        self.file = Path(f"/tmp/{file_name}")
        self.file_name = file_name
        self.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": [0, 1, 2, 3],
                "basepairs": [100.0, 200.0, 300.0, 400.0],
            }
        )
        self.fsa = {
            "DATA1": np.array([10.0, 20.0, 30.0, 40.0]),
            "DATA2": np.array([5.0, 15.0, 25.0, 35.0]),
            "DATA3": np.array([1.0, 2.0, 3.0, 4.0]),
            "DATA4": np.array([100.0, 100.0, 100.0, 100.0]),
            "DATA105": np.array([100.0, 100.0, 100.0, 100.0]),
        }
        self.ladder_steps = np.array([100.0, 200.0, 300.0], dtype=float)
        self.expected_ladder_steps = np.array([100.0, 200.0, 300.0], dtype=float)
        self.best_size_standard = np.array([0.0, 1.0, 2.0], dtype=float)
        self.ladder_model = DummyModel()
        self.fitted_to_model = True
        self.ladder_fit_strategy = "auto_full"
        self.ladder_missing_expected_steps = []
        self.ladder_review_required = False
        self.ladder_fit_note = "All expected ladder steps were fitted."


class TestGeneralAnalysis(unittest.TestCase):
    def setUp(self):
        self._settings_backup = copy.deepcopy(APP_SETTINGS)
        APP_SETTINGS.setdefault("analyses", {})
        APP_SETTINGS["analyses"]["general"] = {
            "pipeline": {
                "ladder": "GS500ROX",
                "trace_channels": ["DATA2", "DATA3"],
                "primary_peak_channel": "DATA3",
                "bp_min": 75.0,
                "bp_max": 425.0,
            }
        }

    def tearDown(self):
        APP_SETTINGS.clear()
        APP_SETTINGS.update(self._settings_backup)

    def test_general_runtime_config_normalizes_values(self):
        APP_SETTINGS["analyses"]["general"]["pipeline"]["ladder"] = "liz500"
        APP_SETTINGS["analyses"]["general"]["pipeline"]["trace_channels"] = "data2"
        APP_SETTINGS["analyses"]["general"]["pipeline"]["primary_peak_channel"] = "DATA2"

        runtime = resolve_runtime_config()
        self.assertEqual(runtime["ladder"], "LIZ500_250")
        self.assertEqual(runtime["trace_channels"], ["DATA2"])
        self.assertEqual(runtime["primary_peak_channel"], "DATA2")

    def test_general_classification_uses_runtime_settings(self):
        meta = classify_fsa(Path("/tmp/example.fsa"))
        self.assertIsNotNone(meta)
        self.assertEqual(meta["assay"], GENERAL_ASSAY_NAME)
        self.assertEqual(meta["ladder"], "GS500ROX")
        self.assertEqual(meta["trace_channels"], ["DATA2", "DATA3"])
        self.assertEqual(meta["primary_peak_channel"], "DATA3")

    def test_general_pipeline_builds_standalone_report(self):
        dummy = DummyFsa()

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            with patch("core.analyses.general.pipeline._scan_files", return_value=[Path("/tmp/example.fsa")]), \
                 patch(
                     "core.analyses.general.pipeline.classify_fsa",
                     return_value={
                         "analysis": "general",
                         "assay": GENERAL_ASSAY_NAME,
                         "group": "sample",
                         "ladder": "GS500ROX",
                         "trace_channels": ["DATA2", "DATA3"],
                         "peak_channels": ["DATA2", "DATA3"],
                         "primary_peak_channel": "DATA3",
                         "sample_channel": "DATA3",
                         "bp_min": 75.0,
                         "bp_max": 425.0,
                         "source_run_dir": "run-a",
                     },
                 ), \
                 patch("core.analyses.general.pipeline.analyse_fsa_rox", return_value=dummy) as mock_rox, \
                 patch("core.analyses.general.pipeline.compute_ladder_qc_metrics", return_value={"r2": 0.9992}), \
                 patch("core.analyses.general.pipeline.build_general_html_report") as mock_report:
                result = run_pipeline(tmp_path / "input", base_outdir=tmp_path / "out", assay_folder_name="general", return_entries=True)

        self.assertEqual(len(result or []), 1)
        mock_rox.assert_called_once()
        mock_report.assert_called_once()

    def test_general_report_contains_plotly_editor_and_comments(self):
        entry = {
            "analysis": "general",
            "assay": GENERAL_ASSAY_NAME,
            "group": "sample",
            "fsa": DummyFsa("demo.fsa"),
            "trace_channels": ["DATA1", "DATA2", "DATA3"],
            "peak_channels": ["DATA1", "DATA2", "DATA3"],
            "primary_peak_channel": "DATA1",
            "sample_channel": "DATA1",
            "bp_min": 75.0,
            "bp_max": 425.0,
            "ladder": "GS500ROX",
            "ladder_qc_status": "ok",
            "ladder_r2": 0.9994,
            "peaks_by_channel": {
                "DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "keep"]),
                "DATA2": pd.DataFrame(columns=["basepairs", "peaks", "area", "keep"]),
                "DATA3": pd.DataFrame(columns=["basepairs", "peaks", "area", "keep"]),
            },
        }

        with TemporaryDirectory() as tmp:
            out = build_general_html_report([entry], Path(tmp) / "reports", run_label="job-1")
            self.assertIsNotNone(out)
            html = out.read_text(encoding="utf-8")

        self.assertIn("job-1_General_Report", html)
        self.assertIn("PeakManager.downloadUpdatedHtml()", html)
        self.assertIn("report-comment", html)
        self.assertIn("demo.fsa", html)

    def test_normalizers_handle_supported_options(self):
        self.assertEqual(normalize_ladder_name("GS500ROX"), "GS500ROX")
        self.assertEqual(normalize_trace_channels(["data1", "DATA3", "DATA9"]), ["DATA1", "DATA3"])
