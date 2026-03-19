import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from core.analyses.flt3.classification import classify_fsa
from core.analyses.flt3.config import ASSAY_REFERENCE_LABEL, ASSAY_REFERENCE_RANGES
from core.analyses.flt3.pipeline import (
    FLT3_QC_TRENDS_FILENAME,
    _build_control_qc_row,
    _calculate_ratios,
    _build_flt3_qc_trend_frames,
    _interpret_entry,
    _peak_area_half_width_bp,
    _resolve_peak_area,
    _scan_files,
    _select_best_entry,
    generate_flt3_peak_report,
    run_pipeline,
    update_flt3_qc_trends,
)
from core.html_reports import _build_flt3_summary_table, _flt3_report_blocks


def _full_meta(**overrides):
    base = {
        "injection_time": 3,
        "injection_voltage": 2000,
        "well_id": "A09",
        "run_name": "Run_3730DNA",
        "run_date": "2026-03-16",
        "run_time": "14:19:33",
        "injection_protocol": "D_3sek_2500_POP7_36cm",
    }
    base.update(overrides)
    return base


class TestFlt3PipelineHardening(unittest.TestCase):
    def test_control_qc_row_marks_negative_control_with_mutant_signal_as_fail(self):
        peaks = pd.DataFrame(
            [
                {"basepairs": 80.0, "peaks": 900.0, "area": 8000.0, "label": "WT"},
                {"basepairs": 129.0, "peaks": 420.0, "area": 2400.0, "label": "MUT"},
            ]
        )
        entry = {
            "fsa": type("DummyFsa", (), {"file_name": "NEG_control.fsa"})(),
            "group": "negative_control",
            "assay": "FLT3-D835",
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
            "ratio": 0.3,
            "well_id": "A01",
            "selected_injection": "3s",
            "injection_time": 3,
            "selection_reason": "",
        }

        row = _build_control_qc_row(entry)

        self.assertEqual(row["Status"], "FAIL")
        self.assertIn("Unexpected mutant peaks found", row["Details"])
        self.assertEqual(row["Expectation"], "Ingen mutant/ITD-topper forventet")

    def test_run_pipeline_does_not_write_legacy_flt3_injection_reports(self):
        peaks = pd.DataFrame(
            [
                {"basepairs": 80.0, "peaks": 1000.0, "area": 9000.0, "label": "WT"},
                {"basepairs": 129.0, "peaks": 600.0, "area": 4200.0, "label": "MUT"},
            ]
        )
        entry = {
            "fsa": type("DummyFsa", (), {"file_name": "POS_control.fsa"})(),
            "group": "positive_control",
            "assay": "FLT3-D835",
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
            "ratio": 0.525,
            "ratio_numerator_area": 4200.0,
            "ratio_denominator_area": 8000.0,
            "ladder_qc_status": "ok",
            "ladder_fit_note": "All expected ladder steps were fitted.",
            "well_id": "B03",
            "selected_injection": "3s",
            "injection_time": 3,
            "selection_reason": "Preferred 3s injection selected",
            "analysis_type": "standard",
            "protocol_injection_time": 3,
            "source_run_dir": "0623",
            "sizing_method": "spline",
        }

        with TemporaryDirectory() as tmp, \
             patch("core.analyses.flt3.pipeline.normalize_pipeline_paths", return_value=(Path("/tmp/flt3-in"), Path(tmp))), \
             patch("core.analyses.flt3.pipeline._scan_files", return_value=[Path("/tmp/control.fsa")]), \
             patch("core.analyses.flt3.pipeline.classify_fsa", return_value={"selection_key": "control"}), \
             patch("core.analyses.flt3.pipeline._select_best_entry", return_value=entry), \
             patch("core.analyses.flt3.pipeline.generate_flt3_peak_report"), \
             patch("core.analyses.flt3.pipeline.generate_flt3_bp_validation_report"), \
             patch("core.analyses.flt3.pipeline.finalize_pipeline_run", return_value=[entry]):
            outdir = Path(tmp)
            run_pipeline(Path("/tmp/flt3-in"), return_entries=True, make_dit_reports=False)

            self.assertFalse((outdir / "QC_FLT3_Injections.csv").exists())
            self.assertFalse((outdir / "QC_FLT3_Injections.html").exists())
            self.assertTrue((outdir / FLT3_QC_TRENDS_FILENAME).exists())

    def test_recursive_scan_finds_nested_fsa_files(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            nested = root / "0622"
            nested.mkdir()
            (nested / "26OUM04232_ITD__130326_A01.fsa").write_text("")
            (root / "water_blank.fsa").write_text("")
            (nested / "V__130326_A11.fsa").write_text("")

            files = _scan_files(root)

        self.assertEqual([p.name for p in files], ["26OUM04232_ITD__130326_A01.fsa"])

    def test_bare_flt3_filename_defaults_to_d835(self):
        with patch(
            "core.analyses.flt3.classification.get_injection_metadata",
            return_value=_full_meta(well_id="A09", injection_time=3),
        ):
            result = classify_fsa(Path("/tmp/0623/26OUM04232__130326_A09_H9H1DI0C.fsa"))

        self.assertIsNotNone(result)
        self.assertEqual(result["assay"], "FLT3-D835")
        self.assertEqual(result["well_id"], "A09")
        self.assertEqual(result["source_run_dir"], "0623")
        self.assertIn("A09", result["selection_key"])
        self.assertEqual(result["protocol_injection_time"], 3)

    def test_selection_key_keeps_replicate_wells_separate(self):
        with patch(
            "core.analyses.flt3.classification.get_injection_metadata",
            side_effect=[
                _full_meta(well_id="A01", injection_time=1, injection_protocol="D_1sek"),
                _full_meta(well_id="A02", injection_time=1, injection_protocol="D_1sek"),
            ],
        ):
            first = classify_fsa(Path("/tmp/0622/26OUM04232_ITD__130326_A01.fsa"))
            second = classify_fsa(Path("/tmp/0622/26OUM04232_ITD__130326_A02.fsa"))

        self.assertNotEqual(first["selection_key"], second["selection_key"])

    def test_select_best_entry_prefers_1s_for_itd_and_3s_for_d835(self):
        def fake_build_entry(path, meta):
            return {
                "fsa": type("DummyFsa", (), {"file_name": path.name})(),
                "assay": meta["assay"],
                "analysis_type": meta["analysis_type"],
                "injection_time": meta["injection_time"],
                "selected_injection": f"{meta['injection_time']}s",
                "selected_injection_time": meta["injection_time"],
                "preferred_injection_time": meta["injection_time"],
                "selection_reason": "",
                "source_run_dir": meta.get("source_run_dir", ""),
                "well_id": meta.get("well_id"),
                "parallel": meta.get("parallel"),
                "selection_key": meta.get("selection_key"),
                "group": meta.get("group", "sample"),
                "ladder_qc_status": "ok",
                "peak_qc_pass": True,
                "peak_qc_status": "ok",
                "alternate_injections": [],
                "alternate_injections_summary": "",
                "ratio": 0.0,
                "ratio_numerator_area": 0.0,
                "ratio_denominator_area": 0.0,
                "mutant_fraction": 0.0,
            }

        itd_candidates = [
            (Path("sample_1s.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 1, "source_run_dir": "0622"}),
            (Path("sample_3s.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 3, "source_run_dir": "0623"}),
        ]
        d835_candidates = [
            (Path("sample_1s.fsa"), {"assay": "FLT3-D835", "analysis_type": "standard", "injection_time": 1, "source_run_dir": "0622"}),
            (Path("sample_3s.fsa"), {"assay": "FLT3-D835", "analysis_type": "standard", "injection_time": 3, "source_run_dir": "0623"}),
        ]

        with patch("core.analyses.flt3.pipeline._build_entry_from_candidate", side_effect=fake_build_entry):
            itd_entry = _select_best_entry(itd_candidates)
            d835_entry = _select_best_entry(d835_candidates)

        self.assertEqual(itd_entry["selected_injection"], "1s")
        self.assertEqual(d835_entry["selected_injection"], "3s")

    def test_select_best_entry_falls_back_when_preferred_fails(self):
        def fake_build_entry(path, meta):
            if meta["injection_time"] == 3:
                return {
                    "fsa": type("DummyFsa", (), {"file_name": path.name})(),
                    "assay": meta["assay"],
                    "analysis_type": meta["analysis_type"],
                    "injection_time": meta["injection_time"],
                    "selected_injection": "3s",
                    "selected_injection_time": 3,
                    "preferred_injection_time": 3,
                    "selection_reason": "",
                    "source_run_dir": meta.get("source_run_dir", ""),
                    "well_id": meta.get("well_id"),
                    "parallel": meta.get("parallel"),
                    "selection_key": meta.get("selection_key"),
                    "group": "sample",
                    "ladder_qc_status": "ladder_qc_failed",
                    "peak_qc_pass": False,
                    "peak_qc_status": "no_relevant_peaks",
                    "alternate_injections": [],
                    "alternate_injections_summary": "",
                    "ratio": 0.0,
                    "ratio_numerator_area": 0.0,
                    "ratio_denominator_area": 0.0,
                    "mutant_fraction": 0.0,
                }
            return {
                "fsa": type("DummyFsa", (), {"file_name": path.name})(),
                "assay": meta["assay"],
                "analysis_type": meta["analysis_type"],
                "injection_time": 1,
                "selected_injection": "1s",
                "selected_injection_time": 1,
                "preferred_injection_time": 3,
                "selection_reason": "",
                "source_run_dir": meta.get("source_run_dir", ""),
                "well_id": meta.get("well_id"),
                "parallel": meta.get("parallel"),
                "selection_key": meta.get("selection_key"),
                "group": "sample",
                "ladder_qc_status": "ok",
                "peak_qc_pass": True,
                "peak_qc_status": "ok",
                "alternate_injections": [],
                "alternate_injections_summary": "",
                "ratio": 0.0,
                "ratio_numerator_area": 0.0,
                "ratio_denominator_area": 0.0,
                "mutant_fraction": 0.0,
            }

        candidates = [
            (Path("sample_1s.fsa"), {"assay": "FLT3-D835", "analysis_type": "standard", "injection_time": 1, "source_run_dir": "0622"}),
            (Path("sample_3s.fsa"), {"assay": "FLT3-D835", "analysis_type": "standard", "injection_time": 3, "source_run_dir": "0623"}),
        ]

        with patch("core.analyses.flt3.pipeline._build_entry_from_candidate", side_effect=fake_build_entry):
            entry = _select_best_entry(candidates)

        self.assertEqual(entry["selected_injection"], "1s")
        self.assertIn("Preferred 3s failed", entry["selection_reason"])

    def test_d835_summary_table_shows_real_ratio_and_selection_metadata(self):
        peaks = pd.DataFrame(
            [
                {"peak_id": "pk_wt", "basepairs": 80.0, "peaks": 1000.0, "area": 8000.0, "label": "WT"},
                {"peak_id": "pk_mut", "basepairs": 129.0, "peaks": 400.0, "area": 2400.0, "label": "MUT"},
                {"peak_id": "pk_digest", "basepairs": 150.0, "peaks": 120.0, "area": 900.0, "label": "unspecific"},
            ]
        )
        entry = {
            "assay": "FLT3-D835",
            "ratio": 0.3,
            "ratio_numerator_area": 2400.0,
            "ratio_denominator_area": 8000.0,
            "ratio_mode": "manual",
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
            "selected_wt_peak_id": "pk_wt",
            "selected_wt_channel": "DATA3",
            "selected_mutant_peak_ids": ["pk_mut"],
            "selected_mutant_channels": ["DATA3"],
            "analysis_type": "standard",
            "protocol_injection_time": 3,
            "selected_injection": "3s",
            "source_run_dir": "0623",
            "selection_reason": "Preferred 3s injection selected",
            "sizing_method": "spline",
        }

        html = _build_flt3_summary_table(entry)

        self.assertIn("0.3000", html)
        self.assertNotIn("Injeksjonsvalg:", html)
        self.assertNotIn("Digest-status", html)
        self.assertIn("150.0 bp", html)

    def test_itd_reference_window_is_300_to_1000(self):
        self.assertEqual(ASSAY_REFERENCE_RANGES["FLT3-ITD"], [(300.0, 1000.0)])
        self.assertIn("300-1000 bp", ASSAY_REFERENCE_LABEL["FLT3-ITD"])

    def test_d835_reference_window_and_label_are_report_friendly(self):
        self.assertEqual(ASSAY_REFERENCE_RANGES["FLT3-D835"], [(50.0, 250.0)])
        self.assertIn("50-250 bp", ASSAY_REFERENCE_LABEL["FLT3-D835"])
        self.assertIn("Mutert >129 bp", ASSAY_REFERENCE_LABEL["FLT3-D835"])

    def test_flt3_report_blocks_show_ratio_before_d835_before_other_itd(self):
        ratio_entry = {"assay": "FLT3-ITD", "analysis_type": "ratio_quant"}
        d835_entry = {"assay": "FLT3-D835", "analysis_type": "standard"}
        itd_entry = {"assay": "FLT3-ITD", "analysis_type": "undiluted"}

        blocks = _flt3_report_blocks(
            {
                "FLT3-ITD": [itd_entry, ratio_entry],
                "FLT3-D835": [d835_entry],
            }
        )

        self.assertEqual(
            [(assay_key, title, len(entries)) for assay_key, title, entries in blocks],
            [
                ("FLT3-ITD", "FLT3-ITD-ratio", 1),
                ("FLT3-D835", "FLT3-D835", 1),
                ("FLT3-ITD", "FLT3-ITD", 1),
            ],
        )

    def test_d835_area_windows_use_narrower_label_specific_widths(self):
        self.assertEqual(_peak_area_half_width_bp("FLT3-D835", "WT", 80.0), 1.2)
        self.assertEqual(_peak_area_half_width_bp("FLT3-D835", "MUT", 129.0), 0.5)
        self.assertEqual(_peak_area_half_width_bp("FLT3-D835", "unspecific", 150.0), 0.8)
        self.assertEqual(_peak_area_half_width_bp("FLT3-ITD", "WT", 330.0), 2.0)
        self.assertEqual(_peak_area_half_width_bp("FLT3-ITD", "ITD", 350.0), 1.0)

    def test_itd_peak_area_prefers_strongest_single_channel(self):
        self.assertEqual(
            _resolve_peak_area("FLT3-ITD", combined_area=596576.0, channel_areas={"DATA1": 384822.0, "DATA2": 160385.0}),
            384822.0,
        )
        self.assertEqual(
            _resolve_peak_area("FLT3-D835", combined_area=1612.0, channel_areas={"DATA3": 797.0}),
            1612.0,
        )

    def test_manual_itd_ratio_uses_selected_blue_wt_and_green_mutant(self):
        peaks = pd.DataFrame(
            [
                {
                    "peak_id": "pk_wt",
                    "basepairs": 330.0,
                    "peaks": 9000.0,
                    "area": 14000.0,
                    "area_DATA1": 10000.0,
                    "area_DATA2": 4000.0,
                    "label": "WT",
                    "keep": True,
                },
                {
                    "peak_id": "pk_mut",
                    "basepairs": 346.0,
                    "peaks": 3800.0,
                    "area": 7000.0,
                    "area_DATA1": 1000.0,
                    "area_DATA2": 5000.0,
                    "label": "ITD",
                    "keep": True,
                },
            ]
        )
        entry = {
            "fsa": type("DummyFsa", (), {"file_name": "sample_itd.fsa"})(),
            "assay": "FLT3-ITD",
            "analysis_type": "standard",
            "primary_peak_channel": "DATA1",
            "peaks_by_channel": {"DATA1": peaks},
            "manual_ratio_selection": {
                "enabled": True,
                "version": 1,
                "wt": {"peak_id": "pk_wt", "channel": "DATA1"},
                "mutants": [{"peak_id": "pk_mut", "channel": "DATA2"}],
            },
        }

        _calculate_ratios([entry])

        self.assertEqual(entry["ratio_mode"], "manual")
        self.assertTrue(entry["manual_ratio_selection_valid"])
        self.assertEqual(entry["selected_wt_peak_id"], "pk_wt")
        self.assertEqual(entry["selected_mutant_peak_ids"], ["pk_mut"])
        self.assertAlmostEqual(entry["ratio_denominator_area"], 10000.0, places=4)
        self.assertAlmostEqual(entry["ratio_numerator_area"], 5000.0, places=4)
        self.assertAlmostEqual(entry["ratio"], 0.5, places=4)
        self.assertAlmostEqual(entry["mutant_fraction"], 1 / 3, places=4)
        self.assertEqual(_interpret_entry(entry), "Positiv FLT3-ITD")

        with TemporaryDirectory() as tmp:
            generate_flt3_peak_report([entry], Path(tmp))
            report = pd.read_csv(Path(tmp) / "Final_Detailed_Peak_Report.csv")

        self.assertEqual(report.iloc[0]["RatioMode"], "manual")
        self.assertEqual(report.iloc[0]["SelectedWT_PeakID"], "pk_wt")
        self.assertEqual(report.iloc[0]["SelectedMutant_PeakIDs"], "pk_mut")
        self.assertAlmostEqual(float(report.iloc[0]["RatioNumeratorArea"]), 5000.0, places=4)
        self.assertAlmostEqual(float(report.iloc[0]["RatioDenominatorArea"]), 10000.0, places=4)

    def test_invalid_manual_itd_selection_does_not_fall_back_to_auto(self):
        peaks = pd.DataFrame(
            [
                {
                    "peak_id": "pk_wt",
                    "basepairs": 330.0,
                    "peaks": 8200.0,
                    "area": 12000.0,
                    "area_DATA1": 12000.0,
                    "area_DATA2": 3000.0,
                    "label": "WT",
                    "keep": True,
                },
                {
                    "peak_id": "pk_mut",
                    "basepairs": 350.0,
                    "peaks": 2100.0,
                    "area": 4500.0,
                    "area_DATA1": 500.0,
                    "area_DATA2": 4500.0,
                    "label": "ITD",
                    "keep": True,
                },
            ]
        )
        invalid_manual_entry = {
            "fsa": type("DummyFsa", (), {"file_name": "sample_itd_invalid.fsa"})(),
            "assay": "FLT3-ITD",
            "analysis_type": "standard",
            "primary_peak_channel": "DATA1",
            "peaks_by_channel": {"DATA1": peaks},
            "manual_ratio_selection": {
                "enabled": True,
                "version": 1,
                "wt": {"peak_id": "missing_wt", "channel": "DATA1"},
                "mutants": [{"peak_id": "missing_mut", "channel": "DATA2"}],
            },
        }

        _calculate_ratios([invalid_manual_entry])

        self.assertEqual(invalid_manual_entry["ratio_mode"], "manual_required")
        self.assertFalse(invalid_manual_entry["manual_ratio_selection_valid"])
        self.assertEqual(invalid_manual_entry["ratio"], 0.0)
        self.assertEqual(invalid_manual_entry["ratio_numerator_area"], 0.0)
        self.assertEqual(invalid_manual_entry["ratio_denominator_area"], 0.0)
        self.assertIn("Ingen gyldige manuelle mutantpeaks", invalid_manual_entry["manual_ratio_selection_reason"])
        self.assertEqual(_interpret_entry(invalid_manual_entry), "Ingen FLT3-ITD pavist")

    def test_manual_d835_ratio_uses_inferred_wt_and_selected_mutant(self):
        peaks = pd.DataFrame(
            [
                {
                    "peak_id": "pk_wt",
                    "basepairs": 80.1,
                    "peaks": 2100.0,
                    "area": 9000.0,
                    "label": "WT",
                    "keep": True,
                },
                {
                    "peak_id": "pk_mut",
                    "basepairs": 129.0,
                    "peaks": 620.0,
                    "area": 2700.0,
                    "label": "MUT",
                    "keep": True,
                },
                {
                    "peak_id": "pk_digest",
                    "basepairs": 150.0,
                    "peaks": 180.0,
                    "area": 950.0,
                    "label": "unspecific",
                    "keep": True,
                },
            ]
        )
        entry = {
            "fsa": type("DummyFsa", (), {"file_name": "sample_d835.fsa"})(),
            "assay": "FLT3-D835",
            "analysis_type": "standard",
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
            "manual_ratio_selection": {
                "enabled": True,
                "version": 2,
                "mutants": [{"peak_id": "pk_mut", "channel": "DATA3"}],
            },
        }

        _calculate_ratios([entry])

        self.assertEqual(entry["ratio_mode"], "manual")
        self.assertTrue(entry["manual_ratio_selection_valid"])
        self.assertEqual(entry["selected_wt_peak_id"], "pk_wt")
        self.assertEqual(entry["selected_mutant_peak_ids"], ["pk_mut"])
        self.assertAlmostEqual(entry["ratio_denominator_area"], 9000.0, places=4)
        self.assertAlmostEqual(entry["ratio_numerator_area"], 2700.0, places=4)
        self.assertAlmostEqual(entry["ratio"], 0.3, places=4)
        self.assertEqual(_interpret_entry(entry), "Positiv FLT3-D835")

    def test_control_qc_row_uses_manual_selection_areas(self):
        peaks = pd.DataFrame(
            [
                {
                    "peak_id": "pk_wt",
                    "basepairs": 330.0,
                    "peaks": 9100.0,
                    "area": 16000.0,
                    "area_DATA1": 11000.0,
                    "area_DATA2": 5000.0,
                    "label": "WT",
                    "keep": True,
                },
                {
                    "peak_id": "pk_mut",
                    "basepairs": 349.5,
                    "peaks": 3900.0,
                    "area": 8000.0,
                    "area_DATA1": 1200.0,
                    "area_DATA2": 6200.0,
                    "label": "ITD",
                    "keep": True,
                },
            ]
        )
        entry = {
            "fsa": type("DummyFsa", (), {"file_name": "control_itd.fsa"})(),
            "group": "positive_control",
            "assay": "FLT3-ITD",
            "analysis_type": "standard",
            "primary_peak_channel": "DATA1",
            "peaks_by_channel": {"DATA1": peaks},
            "manual_ratio_selection": {
                "enabled": True,
                "version": 1,
                "wt": {"peak_id": "pk_wt", "channel": "DATA1"},
                "mutants": [{"peak_id": "pk_mut", "channel": "DATA2"}],
            },
        }

        _calculate_ratios([entry])
        row = _build_control_qc_row(entry)

        self.assertEqual(row["Status"], "PASS")
        self.assertAlmostEqual(row["WT_Area"], 11000.0, places=4)
        self.assertAlmostEqual(row["Mutant_Area"], 6200.0, places=4)
        self.assertAlmostEqual(row["Ratio"], 6200.0 / 11000.0, places=4)

    def test_small_standard_itd_shoulders_do_not_trigger_positive_interpretation(self):
        peaks = pd.DataFrame(
            [
                {"basepairs": 328.0, "peaks": 10000.0, "area": 100000.0, "label": "WT"},
                {"basepairs": 336.5, "peaks": 180.0, "area": 2500.0, "label": "ITD"},
                {"basepairs": 337.4, "peaks": 160.0, "area": 1800.0, "label": "ITD"},
            ]
        )
        entry = {
            "assay": "FLT3-ITD",
            "analysis_type": "standard",
            "primary_peak_channel": "DATA1",
            "peaks_by_channel": {"DATA1": peaks},
        }

        _calculate_ratios([entry])

        self.assertEqual(entry["ratio"], 0.0)
        self.assertEqual(_interpret_entry(entry), "Ingen FLT3-ITD pavist")

    def test_update_flt3_qc_trends_writes_and_dedupes_controls(self):
        with TemporaryDirectory() as tmp:
            excel_path = Path(tmp) / FLT3_QC_TRENDS_FILENAME
            peaks_initial = pd.DataFrame(
                [
                    {"basepairs": 80.0, "peaks": 1000.0, "area": 8000.0, "label": "WT", "keep": True},
                    {"basepairs": 129.0, "peaks": 400.0, "area": 2400.0, "label": "MUT", "keep": True},
                ]
            )
            peaks_updated = pd.DataFrame(
                [
                    {"basepairs": 80.0, "peaks": 1200.0, "area": 9000.0, "label": "WT", "keep": True},
                    {"basepairs": 129.0, "peaks": 500.0, "area": 3000.0, "label": "MUT", "keep": True},
                    {"basepairs": 150.0, "peaks": 120.0, "area": 950.0, "label": "unspecific", "keep": True},
                ]
            )
            base_entry = {
                "fsa": type("DummyFsa", (), {"file_name": "IVS-P001_D8365_kutting__310725_F05.fsa"})(),
                "group": "positive_control",
                "assay": "FLT3-D835",
                "analysis_type": "TKD_digested",
                "dit": "",
                "specimen_id": "IVS-P001",
                "well_id": "F05",
                "run_date": "2026-03-16",
                "run_time": "14:19:33",
                "run_name": "Run_3730DNA",
                "source_run_dir": "0623",
                "injection_protocol": "D_3sek_2500_POP7_36cm",
                "injection_time": 3,
                "selected_injection": "3s",
                "preferred_injection_time": 3,
                "protocol_injection_time": 3,
                "selection_reason": "Preferred 3s injection selected",
                "alternate_injections_summary": "",
                "sizing_method": "spline",
                "ladder": "ROX400HD",
                "ladder_qc_status": "ok",
                "ladder_r2": 0.9987,
                "peak_qc_status": "ok",
                "primary_peak_channel": "DATA3",
                "ratio_numerator_area": 2400.0,
                "ratio_denominator_area": 8000.0,
                "ratio": 0.3,
                "mutant_fraction": 0.2308,
                "peaks_by_channel": {"DATA3": peaks_initial},
            }

            update_flt3_qc_trends(excel_path, [base_entry])

            updated_entry = dict(base_entry)
            updated_entry["selection_reason"] = "Preferred 3s injection selected after rerun"
            updated_entry["ratio_numerator_area"] = 3000.0
            updated_entry["ratio_denominator_area"] = 9000.0
            updated_entry["ratio"] = 0.3333
            updated_entry["mutant_fraction"] = 0.25
            updated_entry["peaks_by_channel"] = {"DATA3": peaks_updated}

            update_flt3_qc_trends(excel_path, [updated_entry])

            runs = pd.read_excel(excel_path, sheet_name="Control_Runs", engine="openpyxl")
            peak_rows = pd.read_excel(excel_path, sheet_name="Control_Peaks", engine="openpyxl")

        self.assertEqual(len(runs), 1)
        self.assertEqual(runs.iloc[0]["File"], "IVS-P001_D8365_kutting__310725_F05.fsa")
        self.assertEqual(runs.iloc[0]["SelectionReason"], "Preferred 3s injection selected after rerun")
        self.assertAlmostEqual(float(runs.iloc[0]["Ratio"]), 0.3333, places=4)
        self.assertEqual(len(peak_rows), 3)
        self.assertEqual(sorted(peak_rows["PeakRank"].tolist()), [1, 2, 3])

    def test_build_flt3_qc_trend_frames_filters_non_controls(self):
        control_entry = {
            "fsa": type("DummyFsa", (), {"file_name": "NTC_RATIO__310725_E04.fsa"})(),
            "group": "negative_control",
            "assay": "FLT3-ITD",
            "analysis_type": "ratio_quant",
            "dit": "",
            "specimen_id": "NTC",
            "well_id": "E04",
            "run_date": "2026-03-16",
            "run_time": "14:19:33",
            "run_name": "Run_3730DNA",
            "source_run_dir": "0623",
            "injection_protocol": "D_1sek",
            "injection_time": 1,
            "selected_injection": "1s",
            "preferred_injection_time": 1,
            "protocol_injection_time": 1,
            "selection_reason": "Preferred 1s injection selected",
            "alternate_injections_summary": "",
            "sizing_method": "spline",
            "ladder": "ROX500",
            "ladder_qc_status": "ok",
            "ladder_r2": 0.9991,
            "peak_qc_status": "ok",
            "primary_peak_channel": "DATA1",
            "ratio_numerator_area": 0.0,
            "ratio_denominator_area": 0.0,
            "ratio": 0.0,
            "mutant_fraction": 0.0,
            "peaks_by_channel": {"DATA1": pd.DataFrame(columns=["basepairs", "peaks", "area", "label", "keep"])},
        }
        sample_entry = dict(control_entry)
        sample_entry["fsa"] = type("DummyFsa", (), {"file_name": "25OUM04232_ITD__130326_A01.fsa"})()
        sample_entry["group"] = "sample"
        sample_entry["specimen_id"] = "25OUM04232"

        df_runs, df_peaks = _build_flt3_qc_trend_frames([control_entry, sample_entry])

        self.assertEqual(len(df_runs), 1)
        self.assertEqual(df_runs.iloc[0]["ControlGroup"], "negative_control")
        self.assertTrue(df_peaks.empty)

    def test_run_pipeline_updates_flt3_qc_trends(self):
        fsa_dir = Path("/tmp/flt3-input")
        assay_dir = Path("/tmp/flt3-output")
        selected_entry = {
            "fsa": type("DummyFsa", (), {"file_name": "IVS-P001_D8365_kutting__310725_F05.fsa"})(),
            "assay": "FLT3-D835",
            "selection_key": "d835",
            "group": "positive_control",
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": pd.DataFrame(columns=["basepairs", "peaks", "area", "label", "keep"])},
            "ratio_numerator_area": 0.0,
            "ratio_denominator_area": 0.0,
            "ratio": 0.0,
            "mutant_fraction": 0.0,
        }

        with patch("core.analyses.flt3.pipeline.normalize_pipeline_paths", return_value=(fsa_dir, assay_dir)), \
             patch("core.analyses.flt3.pipeline._scan_files", return_value=[Path("/tmp/a.fsa")]), \
             patch("core.analyses.flt3.pipeline.classify_fsa", return_value={"selection_key": "d835"}), \
             patch("core.analyses.flt3.pipeline._select_best_entry", return_value=selected_entry), \
             patch("core.analyses.flt3.pipeline._calculate_ratios"), \
             patch("core.analyses.flt3.pipeline.generate_flt3_peak_report"), \
             patch("core.analyses.flt3.pipeline.generate_flt3_bp_validation_report"), \
             patch("core.analyses.flt3.pipeline.update_flt3_qc_trends") as mock_trends, \
             patch("core.analyses.flt3.pipeline.finalize_pipeline_run", return_value=["done"]):
            result = run_pipeline(fsa_dir, return_entries=True, make_dit_reports=False)

        self.assertEqual(result, ["done"])
        mock_trends.assert_called_once_with(assay_dir / FLT3_QC_TRENDS_FILENAME, [selected_entry])


if __name__ == "__main__":
    unittest.main()
