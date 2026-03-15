import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from core.analyses.flt3.classification import classify_fsa
from core.analyses.flt3.config import ASSAY_REFERENCE_LABEL, ASSAY_REFERENCE_RANGES
from core.analyses.flt3.pipeline import (
    _calculate_ratios,
    _interpret_entry,
    _peak_area_half_width_bp,
    _resolve_peak_area,
    _scan_files,
    _select_best_entry,
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
                {"basepairs": 80.0, "peaks": 1000.0, "area": 8000.0, "label": "WT"},
                {"basepairs": 129.0, "peaks": 400.0, "area": 2400.0, "label": "MUT"},
                {"basepairs": 150.0, "peaks": 120.0, "area": 900.0, "label": "unspecific"},
            ]
        )
        entry = {
            "assay": "FLT3-D835",
            "ratio": 0.3,
            "ratio_numerator_area": 2400.0,
            "ratio_denominator_area": 8000.0,
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
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


if __name__ == "__main__":
    unittest.main()
