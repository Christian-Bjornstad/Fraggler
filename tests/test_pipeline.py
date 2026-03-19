import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch

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
            (base / "Vann_kontroll.fsa").write_text("")
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

    def test_registry_only_falls_back_when_target_module_is_missing(self):
        from core.analyses.registry import get_analysis_module

        with patch("core.analyses.registry.get_active_analysis_name", return_value="flt3"), \
             patch("core.analyses.registry.importlib.import_module") as mock_import:
            mock_import.side_effect = ModuleNotFoundError("inner dependency missing")
            mock_import.side_effect.name = "numpy"
            with self.assertRaises(ModuleNotFoundError):
                get_analysis_module("pipeline")

    def test_registry_falls_back_when_analysis_module_is_missing(self):
        from core.analyses.registry import get_analysis_module

        with patch("core.analyses.registry.get_active_analysis_name", return_value="missing"), \
             patch("core.analyses.registry.importlib.import_module") as mock_import:
            fallback_mod = object()

            def _side_effect(name):
                if name == "core.analyses.missing.pipeline":
                    err = ModuleNotFoundError(name)
                    err.name = name
                    raise err
                if name == "core.analyses.clonality.pipeline":
                    return fallback_mod
                raise AssertionError(f"Unexpected import: {name}")

            mock_import.side_effect = _side_effect
            self.assertIs(get_analysis_module("pipeline"), fallback_mod)

    def test_registry_falls_back_when_analysis_package_is_missing(self):
        from core.analyses.registry import get_analysis_module

        with patch("core.analyses.registry.get_active_analysis_name", return_value="missing"), \
             patch("core.analyses.registry.importlib.import_module") as mock_import:
            fallback_mod = object()

            def _side_effect(name):
                if name == "core.analyses.missing.pipeline":
                    err = ModuleNotFoundError("No module named 'core.analyses.missing'")
                    err.name = "core.analyses.missing"
                    raise err
                if name == "core.analyses.clonality.pipeline":
                    return fallback_mod
                raise AssertionError(f"Unexpected import: {name}")

            mock_import.side_effect = _side_effect
            self.assertIs(get_analysis_module("pipeline"), fallback_mod)

    def test_run_pipeline_job_collect_chunks_explicit_files(self):
        from core.runner import CHUNK_SIZE, run_pipeline_job_collect

        files = [Path(f"/tmp/sample_{i}.fsa") for i in range(CHUNK_SIZE + 2)]
        staged_dirs = []

        def _fake_stage(chunk):
            tmp = Path(f"/tmp/staged_{len(staged_dirs)}")
            staged_dirs.append((tmp, list(chunk)))
            return tmp

        with patch("core.runner.stage_files", side_effect=_fake_stage), \
             patch("core.runner.cleanup_temp"), \
             patch("core.pipeline.run_pipeline", side_effect=lambda **kwargs: [{"fsa_dir": kwargs["fsa_dir"]}]):
            entries = run_pipeline_job_collect(
                fsa_dir=None,
                base_outdir=Path("/tmp/out"),
                out_folder_name="ASSAY_REPORTS",
                scope="all",
                needle="",
                files=files,
            )

        self.assertEqual(len(staged_dirs), 2)
        self.assertEqual(len(staged_dirs[0][1]), CHUNK_SIZE)
        self.assertEqual(len(staged_dirs[1][1]), 2)
        self.assertEqual([e["fsa_dir"] for e in entries], [chunk_dir for chunk_dir, _ in staged_dirs])

    def test_run_pipeline_job_collect_can_stage_all_files_once(self):
        from core.runner import CHUNK_SIZE, run_pipeline_job_collect

        files = [Path(f"/tmp/sample_{i}.fsa") for i in range(CHUNK_SIZE + 2)]
        staged_dirs = []

        def _fake_stage(chunk):
            tmp = Path(f"/tmp/staged_{len(staged_dirs)}")
            staged_dirs.append((tmp, list(chunk)))
            return tmp

        with patch("core.runner.stage_files", side_effect=_fake_stage), \
             patch("core.runner.cleanup_temp"), \
             patch("core.pipeline.run_pipeline", side_effect=lambda **kwargs: [{"fsa_dir": kwargs["fsa_dir"]}]):
            entries = run_pipeline_job_collect(
                fsa_dir=None,
                base_outdir=Path("/tmp/out"),
                out_folder_name="ASSAY_REPORTS",
                scope="all",
                needle="",
                files=files,
                chunk_files=False,
            )

        self.assertEqual(len(staged_dirs), 1)
        self.assertEqual(staged_dirs[0][1], files)
        self.assertEqual([e["fsa_dir"] for e in entries], [staged_dirs[0][0]])

    def test_run_pipeline_job_collect_rejects_blank_custom_filter(self):
        from core.runner import run_pipeline_job_collect

        with self.assertRaises(ValueError):
            run_pipeline_job_collect(
                fsa_dir=Path("/tmp/in"),
                base_outdir=Path("/tmp/out"),
                out_folder_name="ASSAY_REPORTS",
                scope="custom",
                needle="",
                files=[Path("/tmp/sample_1.fsa")],
            )

    def test_run_pipeline_job_collect_rejects_custom_filter_with_no_matching_explicit_files(self):
        from core.runner import run_pipeline_job_collect

        with self.assertRaises(ValueError):
            run_pipeline_job_collect(
                fsa_dir=None,
                base_outdir=Path("/tmp/out"),
                out_folder_name="ASSAY_REPORTS",
                scope="custom",
                needle="target",
                files=[Path("/tmp/sample_1.fsa")],
            )

    def test_run_pipeline_job_collect_rejects_custom_filter_with_no_matching_folder_files(self):
        from core.runner import run_pipeline_job_collect

        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "sample_1.fsa").write_text("", encoding="utf-8")

            with self.assertRaises(ValueError):
                run_pipeline_job_collect(
                    fsa_dir=base,
                    base_outdir=Path("/tmp/out"),
                    out_folder_name="ASSAY_REPORTS",
                    scope="custom",
                    needle="target",
                )

    def test_build_filtered_input_uses_shared_staging_helper(self):
        from core.runner import build_filtered_input

        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            matched = base / "sample_a.fsa"
            other = base / "sample_b.txt"
            matched.write_text("", encoding="utf-8")
            other.write_text("", encoding="utf-8")

            staged = Path("/tmp/staged-filtered")
            with patch("core.runner.stage_files", return_value=staged) as mock_stage:
                result = build_filtered_input(base, "sample")

        self.assertEqual(result, staged)
        mock_stage.assert_called_once()
        self.assertEqual([p.name for p in mock_stage.call_args.args[0]], ["sample_a.fsa"])

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

    def test_flt3_pipeline_run_orchestrates_selection_reports_and_finalize(self):
        from core.analyses.flt3 import pipeline as flt3_pipeline

        fsa_dir = Path("/tmp/flt3-input")
        outdir = Path("/tmp/flt3-output")
        selected_entries = [
            {"fsa": type("DummyFsa", (), {"file_name": "sample_d835.fsa"})(), "assay": "FLT3-D835", "selection_key": "d835"},
            {"fsa": type("DummyFsa", (), {"file_name": "sample_itd.fsa"})(), "assay": "FLT3-ITD", "selection_key": "itd"},
        ]
        files = [Path("/tmp/a.fsa"), Path("/tmp/b.fsa"), Path("/tmp/c.fsa")]
        classified_meta = [
            {"selection_key": "d835", "assay": "FLT3-D835"},
            {"selection_key": "d835", "assay": "FLT3-D835"},
            {"selection_key": "itd", "assay": "FLT3-ITD"},
        ]

        with patch.object(flt3_pipeline, "normalize_pipeline_paths", return_value=(fsa_dir, outdir)), \
             patch.object(flt3_pipeline, "_scan_files", return_value=files), \
             patch.object(flt3_pipeline, "classify_fsa", side_effect=classified_meta), \
             patch.object(flt3_pipeline, "_select_best_entry", side_effect=selected_entries) as mock_select, \
             patch.object(flt3_pipeline, "_calculate_ratios") as mock_ratios, \
             patch.object(flt3_pipeline, "generate_flt3_peak_report") as mock_peak, \
             patch.object(flt3_pipeline, "generate_flt3_bp_validation_report") as mock_bp, \
             patch.object(flt3_pipeline, "update_flt3_qc_trends") as mock_trends, \
             patch.object(flt3_pipeline, "finalize_pipeline_run", return_value=["finalized"]) as mock_finalize:
            result = flt3_pipeline.run_pipeline(
                fsa_dir,
                base_outdir=Path("/tmp/base"),
                assay_folder_name="REPORTS",
                return_entries=True,
                make_dit_reports=False,
                mode="all",
            )

        self.assertEqual(result, ["finalized"])
        self.assertEqual(mock_select.call_count, 2)
        mock_ratios.assert_called_once_with(selected_entries)
        mock_peak.assert_called_once_with(selected_entries, outdir)
        mock_bp.assert_called_once_with(selected_entries, outdir)
        mock_trends.assert_called_once_with(outdir / flt3_pipeline.FLT3_QC_TRENDS_FILENAME, selected_entries)
        mock_finalize.assert_called_once_with(
            selected_entries,
            outdir,
            return_entries=True,
            make_dit_reports=False,
            mode="all",
        )

    def test_flt3_summary_tables_render_validation_strings_without_real_data(self):
        from core.html_reports import _build_flt3_summary_table
        import pandas as pd

        d835_peaks = pd.DataFrame(
            [
                {"basepairs": 80.0, "peaks": 1000.0, "area": 8000.0, "label": "WT"},
                {"basepairs": 129.0, "peaks": 400.0, "area": 2400.0, "label": "MUT"},
                {"basepairs": 150.0, "peaks": 140.0, "area": 1100.0, "label": "unspecific"},
            ]
        )
        itd_peaks = pd.DataFrame(
            [
                {
                    "peak_id": "pk_wt",
                    "basepairs": 329.8,
                    "peaks": 1200.0,
                    "area": 5000.0,
                    "label": "WT",
                    "area_DATA1": 3200.0,
                    "area_DATA2": 1800.0,
                },
                {
                    "peak_id": "pk_mut",
                    "basepairs": 360.2,
                    "peaks": 650.0,
                    "area": 2100.0,
                    "label": "ITD",
                    "area_DATA1": 1600.0,
                    "area_DATA2": 500.0,
                },
            ]
        )

        d835_html = _build_flt3_summary_table(
            {
                "assay": "FLT3-D835",
                "ratio": 0.30,
                "ratio_numerator_area": 2400.0,
                "ratio_denominator_area": 8000.0,
                "primary_peak_channel": "DATA3",
                "peaks_by_channel": {"DATA3": d835_peaks},
            }
        )
        itd_html = _build_flt3_summary_table(
            {
                "assay": "FLT3-ITD",
                "ratio": 0.32,
                "ratio_numerator_area": 2100.0,
                "ratio_denominator_area": 5000.0,
                "ratio_mode": "manual",
                "primary_peak_channel": "DATA1",
                "peaks_by_channel": {"DATA1": itd_peaks},
                "selected_wt_peak_id": "pk_wt",
                "selected_wt_channel": "DATA1",
                "selected_mutant_peak_ids": ["pk_mut"],
                "selected_mutant_channels": ["DATA2"],
            }
        )

        # D835 without ratio_mode="manual" now returns empty (no placeholder table)
        self.assertEqual(d835_html, "")
        # ITD with ratio_mode="manual" still renders the full validation table
        self.assertIn("Validering", itd_html)
        self.assertIn("Mut/(Mut+WT)", itd_html)
        self.assertIn("WT-topp", itd_html)

if __name__ == "__main__":
    unittest.main()
