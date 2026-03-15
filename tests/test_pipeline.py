import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch
from math import isnan

from core.pipeline import _scan_files
from core.analyses.flt3.config import PREFERRED_INJECTION_TIME
from core.utils import is_control_file

class TestPipeline(unittest.TestCase):
    def test_scan_files_basic(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "25OUM10166_sample_FR1.fsa").write_text("")
            (base / "NK_sample_IGK.fsa").write_text("")
            (base / "water_blank.fsa").write_text("")
            files = _scan_files(base)
            self.assertEqual([p.name for p in files], ["25OUM10166_sample_FR1.fsa", "NK_sample_IGK.fsa"])

    def test_scan_files_controls_mode(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "25OUM10166_sample_FR1.fsa").write_text("")
            (base / "PK_sample_FR1.fsa").write_text("")
            files = _scan_files(base, mode="controls")
            self.assertEqual([p.name for p in files], ["PK_sample_FR1.fsa"])

    def test_control_detection_includes_flt3_controls(self):
        self.assertTrue(is_control_file("IVS-0000_ITD_ufort__310725_D01.fsa"))
        self.assertTrue(is_control_file("IVS-P001_D8365_kutting__310725_F05.fsa"))
        self.assertTrue(is_control_file("NTC_RATIO__310725_E04.fsa"))

    def test_dit_extraction(self):
        from core.html_reports import extract_dit_from_name
        self.assertEqual(extract_dit_from_name("25OUM10166_some_data.fsa"), "25OUM10166")
        self.assertEqual(extract_dit_from_name("no_dit_here.fsa"), None)
        self.assertEqual(extract_dit_from_name("26OUM00042_ABC.fsa"), "26OUM00042")

    def test_pipeline_dispatches_to_active_analysis_module(self):
        from core.pipeline import run_pipeline

        with patch("core.pipeline.get_analysis_module") as mock_get_module:
            pipeline_mod = type("PipelineMod", (), {"run_pipeline": staticmethod(lambda **kwargs: [{"ok": True, **kwargs}])})
            mock_get_module.return_value = pipeline_mod
            result = run_pipeline(Path("/tmp/input"), return_entries=True)
            self.assertEqual(result[0]["ok"], True)
            self.assertEqual(result[0]["fsa_dir"], Path("/tmp/input"))

    def test_clonality_report_generation_respects_controls_mode(self):
        from core.analyses.clonality import pipeline as clonality_pipeline

        fake_entry = {"fsa": object()}
        with patch.object(clonality_pipeline, "_scan_files", return_value=[Path("sample.fsa")]), \
             patch.object(clonality_pipeline, "_analyze_files", return_value=([fake_entry], 0)), \
             patch("core.analyses.shared_pipeline.build_dit_html_reports") as mock_reports:
            result = clonality_pipeline.run_pipeline(Path("/tmp/in"), return_entries=True, make_dit_reports=True, mode="controls")
            self.assertEqual(result, [fake_entry])
            mock_reports.assert_not_called()

    def test_clonality_report_generation_runs_for_normal_mode(self):
        from core.analyses.clonality import pipeline as clonality_pipeline

        fake_entry = {"fsa": object()}
        with patch.object(clonality_pipeline, "_scan_files", return_value=[Path("sample.fsa")]), \
             patch.object(clonality_pipeline, "_analyze_files", return_value=([fake_entry], 0)), \
             patch("core.analyses.shared_pipeline.build_dit_html_reports") as mock_reports:
            result = clonality_pipeline.run_pipeline(Path("/tmp/in"), return_entries=True, make_dit_reports=True, mode="all")
            self.assertEqual(result, [fake_entry])
            mock_reports.assert_called_once()

    def test_flt3_preferred_injection_time_matches_current_workbook(self):
        self.assertEqual(PREFERRED_INJECTION_TIME["FLT3-ITD"], 1)
        self.assertEqual(PREFERRED_INJECTION_TIME["FLT3-D835"], 3)
        self.assertEqual(PREFERRED_INJECTION_TIME["NPM1"], 3)

    def test_flt3_ratio_thresholds_match_current_workbook_notes(self):
        from core.analyses.flt3.config import ASSAY_CONFIG
        self.assertEqual(ASSAY_CONFIG["FLT3-ITD"]["positive_ratio"], 0.02)
        self.assertEqual(ASSAY_CONFIG["FLT3-D835"]["positive_ratio"], 0.05)

    def test_flt3_d835_real_samples_are_called_positive(self):
        from config import APP_SETTINGS
        from core.pipeline import run_pipeline

        data_dir = Path("/Users/christian/Desktop/OUS/data/flt3/data/rerun_all")
        if not data_dir.exists():
            self.skipTest("FLT3 real-data folder is not available")

        previous = APP_SETTINGS.get("active_analysis", "clonality")
        APP_SETTINGS["active_analysis"] = "flt3"
        try:
            entries = run_pipeline(data_dir, return_entries=True, make_dit_reports=False)
        finally:
            APP_SETTINGS["active_analysis"] = previous

        by_name = {e["fsa"].file_name: e for e in entries or []}
        for filename in [
            "25OUM11316_p1_D8365_kutting__310725_C05_H9C0ZIZJ.fsa",
            "25OUM11316_p2_D8365_kutting__310725_C06_H9C0ZIZJ.fsa",
        ]:
            self.assertIn(filename, by_name)
            ratio = float(by_name[filename].get("ratio", 0.0))
            self.assertGreater(ratio, 0.1, filename)

    def test_flt3_html_validation_strings_are_rendered(self):
        from config import APP_SETTINGS
        from core.pipeline import run_pipeline

        data_dir = Path("/Users/christian/Desktop/OUS/data/flt3/data/rerun_all")
        if not data_dir.exists():
            self.skipTest("FLT3 real-data folder is not available")

        outdir = Path("/Users/christian/Desktop/OUS/tmp/test_flt3_validation_html")
        previous = APP_SETTINGS.get("active_analysis", "clonality")
        APP_SETTINGS["active_analysis"] = "flt3"
        try:
            run_pipeline(data_dir, base_outdir=outdir, assay_folder_name="REPORTS", return_entries=True, make_dit_reports=True)
        finally:
            APP_SETTINGS["active_analysis"] = previous

        html_path = outdir / "REPORTS" / "25OUM11316_Flt3_Resultater.html"
        self.assertTrue(html_path.exists())
        html = html_path.read_text(encoding="utf-8")
        self.assertIn("Validering:", html)
        self.assertIn("150 bp kontroll", html)
        self.assertIn("Mut/(Mut+WT)", html)

if __name__ == "__main__":
    unittest.main()
