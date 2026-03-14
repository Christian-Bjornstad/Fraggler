import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch

from core.pipeline import _scan_files

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

if __name__ == "__main__":
    unittest.main()
