import unittest
import copy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from config import APP_SETTINGS
from core.analyses.flt3.classification import classify_fsa
from core.analyses.flt3.config import ASSAY_REFERENCE_LABEL, ASSAY_REFERENCE_RANGES, ROX_LADDER as FLT3_ROX_LADDER
from core.analyses.flt3.pipeline import (
    FLT3_QC_TRENDS_FILENAME,
    FLT3_NPM1_QC_TRACKER_FILENAME,
    _analyse_fsa_candidate,
    _attempt_flt3_template_fit,
    _build_control_qc_row,
    _build_flt3_npm1_tracker_frames,
    _calculate_ratios,
    _build_flt3_qc_trend_frames,
    _candidate_flt3_template_keys,
    _flt3_fit_is_geometrically_invalid,
    _flt3_high_end_anchors_are_plausible,
    _flt3_template_key_allowed_for_fsa,
    _late_trace_peak_candidates,
    _flt3_mapping_shape_penalty,
    _rank_flt3_template_keys_for_fsa,
    _resolved_flt3_template_key,
    _select_flt3_high_end_anchor_combo,
    _template_mapping_payload_for_anchors,
    _flt3_template_key,
    _interpret_entry,
    _peak_area_half_width_bp,
    _resolve_peak_area,
    _scan_files,
    _select_best_entry,
    generate_flt3_peak_report,
    run_pipeline,
    update_flt3_npm1_qc_tracker_workbook,
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


def _dummy_template_fsa(
    peak_times: list[int],
    peak_intensities: list[float] | None = None,
    *,
    ladder_steps: list[float] | None = None,
    best_size_standard: list[float] | None = None,
):
    intensities = peak_intensities or [100.0] * len(peak_times)
    size = max(5000, max(peak_times) + 50)
    trace = np.zeros(size, dtype=float)
    for time, intensity in zip(peak_times, intensities):
        trace[int(time)] = float(intensity)
    steps = np.asarray(
        ladder_steps
        or [35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0, 400.0, 450.0, 490.0, 500.0],
        dtype=float,
    )
    fsa = type("DummyFsa", (), {})()
    fsa.size_standard = trace
    fsa.expected_ladder_steps = steps.copy()
    fsa.ladder_steps = steps.copy()
    if best_size_standard is None:
        fsa.best_size_standard = np.asarray(peak_times[: len(steps)], dtype=float)
    else:
        fsa.best_size_standard = np.asarray(best_size_standard, dtype=float)
    return fsa


class TestFlt3PipelineHardening(unittest.TestCase):
    def setUp(self):
        self._settings_backup = copy.deepcopy(APP_SETTINGS)

    def tearDown(self):
        APP_SETTINGS.clear()
        APP_SETTINGS.update(self._settings_backup)

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
             patch("core.analyses.flt3.pipeline.update_flt3_npm1_qc_tracker_workbook"), \
             patch("core.analyses.flt3.pipeline.finalize_pipeline_run", return_value=[entry]), \
             patch.dict("config.APP_SETTINGS", {"analyses": {"flt3": {"batch": {"tracking_excel_path": ""}}}}, clear=False):
            outdir = Path(tmp)
            run_pipeline(Path("/tmp/flt3-in"), return_entries=True, make_dit_reports=False)

            self.assertFalse((outdir / "QC_FLT3_Injections.csv").exists())
            self.assertFalse((outdir / "QC_FLT3_Injections.html").exists())
            self.assertTrue((outdir / FLT3_QC_TRENDS_FILENAME).exists())

    def test_run_pipeline_uses_configured_flt3_tracking_workbook_path(self):
        entry = {
            "fsa": type("DummyFsa", (), {"file_name": "POS_control.fsa"})(),
            "group": "positive_control",
            "assay": "FLT3-D835",
        }
        custom_excel = Path("/tmp/shared/flt3/custom.xlsx")

        with patch("core.analyses.flt3.pipeline.normalize_pipeline_paths", return_value=(Path("/tmp/flt3-in"), Path("/tmp/flt3-out"))), \
             patch("core.analyses.flt3.pipeline._scan_files", return_value=[Path("/tmp/control.fsa")]), \
             patch("core.analyses.flt3.pipeline.classify_fsa", return_value={"selection_key": "control"}), \
             patch("core.analyses.flt3.pipeline._select_best_entry", return_value=entry), \
             patch("core.analyses.flt3.pipeline._calculate_ratios"), \
             patch("core.analyses.flt3.pipeline.generate_flt3_peak_report"), \
             patch("core.analyses.flt3.pipeline.generate_flt3_bp_validation_report"), \
             patch("core.analyses.flt3.pipeline.update_flt3_qc_trends") as mock_trends, \
             patch("core.analyses.flt3.pipeline.update_flt3_npm1_qc_tracker_workbook") as mock_tracker, \
             patch("core.analyses.flt3.pipeline.finalize_pipeline_run", return_value=[entry]), \
             patch.dict("config.APP_SETTINGS", {"analyses": {"flt3": {"batch": {"tracking_excel_path": str(custom_excel)}}}}, clear=False):
            run_pipeline(Path("/tmp/flt3-in"), return_entries=True, make_dit_reports=False)

        mock_trends.assert_called_once_with(custom_excel, [entry])
        mock_tracker.assert_called_once_with(custom_excel.parent / FLT3_NPM1_QC_TRACKER_FILENAME, [entry])

    def test_recursive_scan_finds_nested_fsa_files(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            nested = root / "0622"
            nested.mkdir()
            (nested / "26OUM04232_ITD__130326_A01.fsa").write_text("ok", encoding="utf-8")
            (root / "water_blank.fsa").write_text("ok", encoding="utf-8")
            (nested / "V__130326_A11.fsa").write_text("ok", encoding="utf-8")

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

    def test_analyse_fsa_candidate_uses_configured_flt3_rox_ladder(self):
        dummy_fsa = type("DummyFsa", (), {"_flt3_sizing_method": "preconfigured"})()

        with patch("core.analyses.flt3.pipeline.analyse_fsa_rox", return_value=dummy_fsa) as mock_analyse, \
             patch("core.analyses.flt3.pipeline._attempt_flt3_template_fit", return_value=None):
            result = _analyse_fsa_candidate(Path("/tmp/sample_itd.fsa"), "DATA1", "FLT3-ITD", "standard")

        self.assertIs(result, dummy_fsa)
        mock_analyse.assert_called_once()
        call_args = mock_analyse.call_args
        self.assertEqual(call_args.args[:2], (Path("/tmp/sample_itd.fsa"), "DATA1"))
        self.assertEqual(call_args.kwargs["ladder_name"], FLT3_ROX_LADDER)

    def test_analyse_fsa_candidate_prefers_template_rescue_when_available(self):
        auto_fsa = type("DummyFsa", (), {"file_name": "auto.fsa"})()
        rescued_fsa = type("DummyFsa", (), {"file_name": "rescued.fsa"})()

        with patch("core.analyses.flt3.pipeline.analyse_fsa_rox", return_value=auto_fsa), \
             patch("core.analyses.flt3.pipeline._attempt_flt3_template_fit", return_value=rescued_fsa) as mock_rescue, \
             patch("core.analyses.flt3.pipeline._infer_sizing_method", return_value="spline"):
            result = _analyse_fsa_candidate(Path("/tmp/sample_itd.fsa"), "DATA1", "FLT3-ITD", "10x_diluted")

        self.assertIs(result, rescued_fsa)
        mock_rescue.assert_called_once_with(auto_fsa, "FLT3-ITD", "10x_diluted")

    def test_flt3_template_key_uses_group_specific_profiles(self):
        self.assertEqual(_flt3_template_key("FLT3-ITD", "25x_diluted"), ("FLT3-ITD", "25x_diluted"))
        self.assertEqual(_flt3_template_key("FLT3-D835", "standard"), ("FLT3-D835", "standard"))
        self.assertEqual(_flt3_template_key("FLT3-ITD", "standard"), ("FLT3-ITD", "standard"))
        self.assertEqual(_flt3_template_key("FLT3-ITD", "unknown"), ("FLT3-ITD", "10x_diluted"))

    def test_candidate_template_keys_include_compact_d835_family(self):
        self.assertEqual(
            _candidate_flt3_template_keys("FLT3-D835", "TKD_digested"),
            [("FLT3-D835", "standard"), ("FLT3-D835", "standard_compact")],
        )

    def test_resolved_template_key_chooses_compact_family_for_4100_high_end_cluster(self):
        fsa = _dummy_template_fsa(
            [1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854, 4077, 4123],
            peak_intensities=[200.0] * 16 + [1439.0, 1438.0],
            best_size_standard=[1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854],
        )
        candidate_df = pd.DataFrame(
            [
                {"time": 3287.0, "intensity": 1228.0, "source": "auto"},
                {"time": 3581.0, "intensity": 1270.0, "source": "auto"},
                {"time": 3854.0, "intensity": 1337.0, "source": "auto"},
            ]
        )

        with patch("core.analyses.flt3.pipeline.get_ladder_candidates", return_value=candidate_df):
            self.assertEqual(
                _resolved_flt3_template_key(fsa, "FLT3-ITD", "10x_diluted"),
                ("FLT3-ITD", "25x_diluted_compact"),
            )

    def test_resolved_template_key_keeps_late_d835_family_when_500_is_near_4590(self):
        fsa = _dummy_template_fsa(
            [1578, 1652, 1799, 1944, 2174, 2233, 2293, 2538, 2839, 3163, 3408, 3469, 3786, 4285, 4536, 4588],
            peak_intensities=[300.0] * 13 + [427.0, 457.0, 459.0],
        )
        candidate_df = pd.DataFrame(
            [
                {"time": 4285.0, "intensity": 427.0, "source": "trace"},
                {"time": 4536.0, "intensity": 457.0, "source": "trace"},
                {"time": 4588.0, "intensity": 459.0, "source": "trace"},
            ]
        )

        with patch("core.analyses.flt3.pipeline.get_ladder_candidates", return_value=candidate_df):
            self.assertEqual(
                _resolved_flt3_template_key(fsa, "FLT3-D835", "standard"),
                ("FLT3-D835", "standard"),
            )

    def test_ranked_template_keys_move_late_family_ahead_of_false_early_auto_fit(self):
        fsa = _dummy_template_fsa(
            [1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854, 4077, 4123],
            peak_intensities=[150.0] * 16 + [1439.0, 1438.0],
            best_size_standard=[1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854],
        )
        candidate_df = pd.DataFrame(
            [
                {"time": 3287.0, "intensity": 1228.0, "source": "auto"},
                {"time": 3581.0, "intensity": 1270.0, "source": "auto"},
                {"time": 3854.0, "intensity": 1337.0, "source": "auto"},
            ]
        )

        with patch("core.analyses.flt3.pipeline.get_ladder_candidates", return_value=candidate_df):
            ranked = _rank_flt3_template_keys_for_fsa(fsa, "FLT3-ITD", "10x_diluted")

        self.assertEqual(ranked[0], ("FLT3-ITD", "25x_diluted_compact"))
        self.assertNotEqual(ranked[0], ("FLT3-ITD", "10x_diluted"))

    def test_flt3_high_end_anchor_guard_rejects_implausible_early_auto_fit(self):
        self.assertFalse(_flt3_high_end_anchors_are_plausible(3480.0, 3802.0, 4335.0, 4385.0))
        self.assertTrue(_flt3_high_end_anchors_are_plausible(4334.0, 4386.0, 4325.0, 4376.5))

    def test_late_trace_peak_candidates_only_keeps_positive_late_peaks(self):
        trace = pd.Series(
            [0.0] * 3800
            + [12.0, 0.0, -5.0]
            + [0.0] * 200
            + [150.0, 0.0, -3.0]
            + [0.0] * 50
            + [180.0, 0.0, -2.0]
            + [0.0] * 50
        ).to_numpy()

        rows = _late_trace_peak_candidates(trace, min_time=3800.0, min_intensity=10.0)

        self.assertEqual(rows["time"].round().astype(int).tolist(), [3800, 4003, 4056])
        self.assertTrue((rows["intensity"] > 0).all())

    def test_select_high_end_anchor_combo_prefers_late_gap_consistent_family(self):
        rows = pd.DataFrame(
            [
                {"time": 3432.0, "intensity": 167.0, "source": "trace"},
                {"time": 3749.0, "intensity": 129.0, "source": "trace"},
                {"time": 4041.0, "intensity": 141.0, "source": "trace"},
                {"time": 4279.0, "intensity": 159.0, "source": "trace"},
                {"time": 4329.0, "intensity": 188.0, "source": "trace"},
            ]
        )

        combo = _select_flt3_high_end_anchor_combo(
            rows,
            template_450=4162.0,
            template_490=4410.0,
            template_500=4461.0,
        )

        self.assertEqual(combo, (4041.0, 4279.0, 4329.0))

    def test_attempt_template_fit_can_replace_worse_manual_adjustment(self):
        fsa = _dummy_template_fsa(
            [1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854],
            best_size_standard=[1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854],
        )
        setattr(fsa, "ladder_fit_strategy", "manual_adjustment")

        def fake_apply_mapping(trial, payload):
            setattr(trial, "rescued_fit", True)
            setattr(trial, "best_size_standard", np.asarray([1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3854, 4077, 4123], dtype=float))
            return trial

        def fake_qc(obj):
            if getattr(obj, "rescued_fit", False):
                return {"r2": 0.999975, "max_abs_error_bp": 1.41, "mean_abs_error_bp": 0.62}
            return {"r2": 0.992000, "max_abs_error_bp": 12.50, "mean_abs_error_bp": 3.10}

        payload = {
            "mapping": {},
            "mapping_times": {13: 3854.0, 14: 4077.0, 15: 4123.0},
            "manual_candidates": [],
            "template_key": ("FLT3-ITD", "25x_diluted_compact"),
        }

        with patch("core.analyses.flt3.pipeline.get_ladder_candidates", return_value=pd.DataFrame(columns=["time", "intensity", "source"])), \
             patch("core.analyses.flt3.pipeline._rank_flt3_template_keys_for_fsa", return_value=[("FLT3-ITD", "25x_diluted_compact")]), \
             patch("core.analyses.flt3.pipeline._rank_flt3_high_end_anchor_combos", return_value=[(3854.0, 4077.0, 4123.0)]), \
             patch("core.analyses.flt3.pipeline._template_mapping_payload", return_value=payload), \
             patch("core.analyses.flt3.pipeline._template_mapping_payload_for_anchors", return_value=payload), \
             patch("core.analyses.flt3.pipeline.apply_manual_ladder_mapping", side_effect=fake_apply_mapping), \
             patch("core.analyses.flt3.pipeline.compute_ladder_qc_metrics", side_effect=fake_qc):
            rescued = _attempt_flt3_template_fit(fsa, "FLT3-ITD", "10x_diluted")

        self.assertIsNotNone(rescued)
        self.assertEqual(getattr(rescued, "ladder_fit_strategy", ""), "flt3_template_rescue")

    def test_template_rescue_allows_three_bp_max_residual_without_review_flag(self):
        fsa = _dummy_template_fsa(
            [1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854],
            best_size_standard=[1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3287, 3581, 3854],
        )

        def fake_apply_mapping(trial, payload):
            setattr(trial, "rescued_fit", True)
            setattr(trial, "best_size_standard", np.asarray([1490, 1520, 1589, 1731, 1865, 2083, 2138, 2193, 2421, 2699, 3004, 3165, 3229, 3854, 4077, 4123], dtype=float))
            return trial

        def fake_qc(obj):
            if getattr(obj, "rescued_fit", False):
                return {"r2": 0.999950, "max_abs_error_bp": 2.95, "mean_abs_error_bp": 0.85}
            return {"r2": 0.992000, "max_abs_error_bp": 12.50, "mean_abs_error_bp": 3.10}

        payload = {
            "mapping": {},
            "mapping_times": {13: 3854.0, 14: 4077.0, 15: 4123.0},
            "manual_candidates": [],
            "template_key": ("FLT3-ITD", "25x_diluted_compact"),
        }

        with patch("core.analyses.flt3.pipeline.get_ladder_candidates", return_value=pd.DataFrame(columns=["time", "intensity", "source"])), \
             patch("core.analyses.flt3.pipeline._rank_flt3_template_keys_for_fsa", return_value=[("FLT3-ITD", "25x_diluted_compact")]), \
             patch("core.analyses.flt3.pipeline._rank_flt3_high_end_anchor_combos", return_value=[(3854.0, 4077.0, 4123.0)]), \
             patch("core.analyses.flt3.pipeline._template_mapping_payload", return_value=payload), \
             patch("core.analyses.flt3.pipeline._template_mapping_payload_for_anchors", return_value=payload), \
             patch("core.analyses.flt3.pipeline.apply_manual_ladder_mapping", side_effect=fake_apply_mapping), \
             patch("core.analyses.flt3.pipeline.compute_ladder_qc_metrics", side_effect=fake_qc):
            rescued = _attempt_flt3_template_fit(fsa, "FLT3-ITD", "10x_diluted")

        self.assertIsNotNone(rescued)
        self.assertFalse(bool(getattr(rescued, "ladder_review_required", False)))

    def test_mapping_shape_penalty_prefers_manual_like_spacing(self):
        template_times = np.asarray(
            [1580.0, 1652.0, 1804.0, 1948.0, 2181.0, 2240.0, 2299.0, 2544.0, 2842.0, 3174.0, 3420.0, 3484.0, 3807.0, 4104.0, 4346.0, 4396.0],
            dtype=float,
        )
        manual_like = {idx: float(value) for idx, value in enumerate(template_times)}
        distorted = dict(manual_like)
        distorted[10] = 3355.0
        distorted[11] = 3490.0
        distorted[13] = 3950.0
        distorted[14] = 4215.0
        distorted[15] = 4295.0

        self.assertLess(
            _flt3_mapping_shape_penalty(manual_like, template_times),
            _flt3_mapping_shape_penalty(distorted, template_times),
        )

    def test_fit_is_geometrically_invalid_for_perfect_but_wrong_high_end(self):
        template_times = np.asarray(
            [1510.0, 1578.0, 1719.0, 1854.0, 2070.0, 2125.0, 2180.0, 2408.0, 2686.0, 2990.0, 3215.0, 3273.0, 3568.0, 3841.0, 4064.0, 4110.0],
            dtype=float,
        )
        expected_steps = np.asarray([35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0, 400.0, 450.0, 490.0, 500.0], dtype=float)
        wrong_mapping = {
            0: 1501.0,
            1: 1559.0,
            2: 1714.0,
            3: 1873.0,
            4: 2024.0,
            5: 2268.0,
            6: 2330.0,
            7: 2392.0,
            8: 2650.0,
            9: 2963.0,
            10: 3311.0,
            11: 3498.0,
            12: 3212.0,
            13: 3269.0,
            14: 3561.0,
            15: 3832.0,
        }
        trace = np.zeros(5000, dtype=float)
        self.assertTrue(
            _flt3_fit_is_geometrically_invalid(wrong_mapping, expected_steps, template_times, trace)
        )

    def test_compact_template_is_blocked_for_known_noncompact_march_runs(self):
        fsa = type("DummyFsa", (), {})()
        fsa.file_name = "26OUM04273_ITD_X10__130326_C04_H9H1DI0C.fsa"

        self.assertFalse(
            _flt3_template_key_allowed_for_fsa(fsa, ("FLT3-ITD", "25x_diluted_compact"))
        )
        self.assertTrue(
            _flt3_template_key_allowed_for_fsa(fsa, ("FLT3-ITD", "10x_diluted"))
        )

    def test_template_payload_prefers_gap_consistent_35_peak_over_blob(self):
        template_times = np.asarray(
            [1575.0, 1647.0, 1798.5, 1942.0, 2174.5, 2233.5, 2292.5, 2537.5, 2834.5, 3166.0, 3411.5, 3475.5, 3797.5, 4093.5, 4335.0, 4385.0],
            dtype=float,
        )
        expected_steps = np.asarray(
            [35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0, 400.0, 450.0, 490.0, 500.0],
            dtype=float,
        )
        true_times = [1554, 1624, 1773, 1915, 2143, 2201, 2260, 2501, 2795, 3121, 3364, 3427, 3745, 4040, 4281, 4331]
        trace = np.zeros(5000, dtype=float)
        for time in true_times:
            trace[time] = 450.0
        trace[1522] = 1400.0
        fsa = _dummy_template_fsa(true_times, best_size_standard=true_times)
        fsa.size_standard = trace
        candidates = pd.DataFrame([{"time": 1522.0, "intensity": 1400.0, "source": "auto"}])

        payload = _template_mapping_payload_for_anchors(
            fsa,
            ("FLT3-ITD", "10x_diluted"),
            template_times,
            expected_steps,
            candidates,
            4040.0,
            4281.0,
            4331.0,
        )

        self.assertIsNotNone(payload)
        self.assertEqual(payload["mapping_times"][0], 1554.0)

    def test_template_payload_prefers_340_350_pair_over_late_blob(self):
        template_times = np.asarray(
            [1618.0, 1686.0, 1842.5, 1986.5, 2221.5, 2281.5, 2340.0, 2588.0, 2883.5, 3220.5, 3465.5, 3532.0, 3859.5, 4162.0, 4410.0, 4461.0],
            dtype=float,
        )
        expected_steps = np.asarray(
            [35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0, 400.0, 450.0, 490.0, 500.0],
            dtype=float,
        )
        true_times = [1610, 1681, 1836, 1982, 2219, 2278, 2338, 2587, 2887, 3225, 3473, 3539, 3866, 4165, 4408, 4458]
        trace = np.zeros(5000, dtype=float)
        for time in true_times:
            trace[time] = 750.0
        trace[3598] = 3000.0
        fsa = _dummy_template_fsa(true_times, best_size_standard=true_times)
        fsa.size_standard = trace
        candidates = pd.DataFrame(
            [
                {"time": 3598.0, "intensity": 3000.0, "source": "auto"},
                {"time": 3486.0, "intensity": 2600.0, "source": "auto"},
            ]
        )

        payload = _template_mapping_payload_for_anchors(
            fsa,
            ("FLT3-ITD", "25x_diluted"),
            template_times,
            expected_steps,
            candidates,
            4165.0,
            4408.0,
            4458.0,
        )

        self.assertIsNotNone(payload)
        self.assertEqual(payload["mapping_times"][10], 3473.0)
        self.assertEqual(payload["mapping_times"][11], 3539.0)

    def test_template_payload_prefers_gap_consistent_close_triplet(self):
        template_times = np.asarray(
            [1575.0, 1647.0, 1798.5, 1942.0, 2174.5, 2233.5, 2292.5, 2537.5, 2834.5, 3166.0, 3411.5, 3475.5, 3797.5, 4093.5, 4335.0, 4385.0],
            dtype=float,
        )
        expected_steps = np.asarray(
            [35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0, 400.0, 450.0, 490.0, 500.0],
            dtype=float,
        )
        candidate_times = [1491, 1510, 1546, 1615, 1764, 1904, 2130, 2188, 2246, 2485, 2775, 3097, 3336, 3398, 3711, 4001, 4236, 4286]
        trace = np.zeros(5000, dtype=float)
        for time in candidate_times:
            trace[time] = 350.0
        trace[1491] = 20000.0
        fsa = _dummy_template_fsa(candidate_times, best_size_standard=candidate_times[:16])
        fsa.size_standard = trace
        candidates = pd.DataFrame(
            [{"time": float(time), "intensity": float(trace[time]), "source": "auto"} for time in candidate_times]
        )

        payload = _template_mapping_payload_for_anchors(
            fsa,
            ("FLT3-ITD", "10x_diluted"),
            template_times,
            expected_steps,
            candidates,
            4001.0,
            4236.0,
            4286.0,
        )

        self.assertIsNotNone(payload)
        self.assertEqual(payload["mapping_times"][4], 2130.0)
        self.assertEqual(payload["mapping_times"][5], 2188.0)
        self.assertEqual(payload["mapping_times"][6], 2246.0)

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

    def test_select_best_entry_stays_on_required_injection_even_when_it_fails(self):
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

        self.assertEqual(entry["selected_injection"], "3s")
        self.assertIn("No candidate passed QC; kept 3s", entry["selection_reason"])

    def test_select_best_entry_chooses_best_preferred_candidate_by_ladder_r2(self):
        def fake_build_entry(path, meta):
            return {
                "fsa": type("DummyFsa", (), {"file_name": path.name})(),
                "assay": meta["assay"],
                "analysis_type": meta["analysis_type"],
                "injection_time": meta["injection_time"],
                "selected_injection": f"{meta['injection_time']}s",
                "selected_injection_time": int(meta["injection_time"]),
                "preferred_injection_time": 1,
                "selection_reason": "",
                "source_run_dir": meta.get("source_run_dir", ""),
                "well_id": meta.get("well_id"),
                "parallel": meta.get("parallel"),
                "selection_key": meta.get("selection_key"),
                "group": "sample",
                "ladder_qc_status": "ok",
                "ladder_r2": meta["ladder_r2"],
                "n_ladder_steps": meta.get("n_ladder_steps", 16),
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
            (Path("sample_1s_low.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 1, "ladder_r2": 0.991, "source_run_dir": "0622"}),
            (Path("sample_1s_high.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 1, "ladder_r2": 0.999, "source_run_dir": "0623"}),
            (Path("sample_3s.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 3, "ladder_r2": 1.0, "source_run_dir": "0624"}),
        ]

        with patch("core.analyses.flt3.pipeline._build_entry_from_candidate", side_effect=fake_build_entry):
            entry = _select_best_entry(candidates)

        self.assertEqual(entry["selected_injection"], "1s")
        self.assertEqual(entry["fsa"].file_name, "sample_1s_high.fsa")

    def test_select_best_entry_treats_manual_adjustment_as_acceptable(self):
        def fake_build_entry(path, meta):
            return {
                "fsa": type("DummyFsa", (), {"file_name": path.name})(),
                "assay": meta["assay"],
                "analysis_type": meta["analysis_type"],
                "injection_time": meta["injection_time"],
                "selected_injection": f"{meta['injection_time']}s",
                "selected_injection_time": int(meta["injection_time"]),
                "preferred_injection_time": 1,
                "selection_reason": "",
                "source_run_dir": meta.get("source_run_dir", ""),
                "well_id": meta.get("well_id"),
                "parallel": meta.get("parallel"),
                "selection_key": meta.get("selection_key"),
                "group": "sample",
                "ladder_qc_status": meta["ladder_qc_status"],
                "ladder_r2": meta.get("ladder_r2", 0.998),
                "n_ladder_steps": meta.get("n_ladder_steps", 16),
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
            (Path("sample_manual_1s.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 1, "ladder_qc_status": "manual_adjustment", "ladder_r2": 0.95}),
            (Path("sample_failed_1s.fsa"), {"assay": "FLT3-ITD", "analysis_type": "standard", "injection_time": 1, "ladder_qc_status": "ladder_qc_failed", "ladder_r2": 0.99}),
        ]

        with patch("core.analyses.flt3.pipeline._build_entry_from_candidate", side_effect=fake_build_entry):
            entry = _select_best_entry(candidates)

        self.assertEqual(entry["fsa"].file_name, "sample_manual_1s.fsa")
        self.assertEqual(entry["ladder_qc_status"], "manual_adjustment")

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

    def test_build_flt3_npm1_tracker_frames_emits_control_and_ladder_markers(self):
        peaks = pd.DataFrame(
            [
                {"basepairs": 80.2, "peaks": 1200.0, "area": 12000.0, "label": "WT", "keep": True},
                {"basepairs": 128.8, "peaks": 650.0, "area": 6400.0, "label": "MUT", "keep": True},
            ]
        )
        sample_bp = np.array([74.8, 79.4, 80.2, 98.9, 100.2, 128.8, 149.3, 199.2], dtype=float)
        fsa = type("DummyFsa", (), {"file_name": "IVS-P001_D8365_kutting__310725_F05.fsa"})()
        fsa.ladder_steps = np.asarray([75.0, 100.0, 150.0, 200.0], dtype=float)
        fsa.best_size_standard = np.asarray([10.0, 20.0, 30.0, 40.0], dtype=float)
        fsa.size_standard = np.asarray([0.0] * 50, dtype=float)
        fsa.size_standard[10] = 120.0
        fsa.size_standard[20] = 180.0
        fsa.size_standard[30] = 160.0
        fsa.size_standard[40] = 140.0
        fsa.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": np.asarray([10.0, 20.0, 30.0, 40.0], dtype=float),
                "basepairs": np.asarray([74.8, 100.2, 149.3, 199.2], dtype=float),
            }
        )
        entry = {
            "fsa": fsa,
            "group": "positive_control",
            "assay": "FLT3-D835",
            "specimen_id": "IVS-P001",
            "run_date": "2026-03-16",
            "source_run_dir": "0623",
            "run_name": "Run_3730DNA",
            "analysis_type": "TKD_digested",
            "selected_injection": "3s",
            "selected_injection_time": 3,
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
            "ratio": 0.5333,
            "ladder_qc_status": "ok",
            "peak_qc_status": "ok",
        }

        patient_df, control_df, marker_df = _build_flt3_npm1_tracker_frames([entry])

        self.assertTrue(patient_df.empty)
        self.assertEqual(list(control_df["Control"]), ["PK"])
        self.assertEqual(list(control_df["SampleKind"]), ["control"])
        marker_names = set(marker_df["MarkerName"].tolist())
        self.assertIn("IVSP001_D835_128_129", marker_names)
        self.assertIn("D835_Ladder_139", marker_names)
        d835_mut = marker_df.loc[marker_df["MarkerName"] == "IVSP001_D835_128_129"].iloc[0]
        self.assertTrue(bool(d835_mut["OK"]))
        self.assertAlmostEqual(float(d835_mut["FoundBP"]), 128.8, places=2)

    def test_update_flt3_npm1_qc_tracker_workbook_writes_expected_sheets(self):
        peaks = pd.DataFrame(
            [
                {"basepairs": 299.1, "peaks": 1400.0, "area": 24000.0, "label": "WT", "keep": True},
                {"basepairs": 304.2, "peaks": 900.0, "area": 12500.0, "label": "MUT", "keep": True},
            ]
        )
        fsa = type("DummyFsa", (), {"file_name": "IVS-0000_NPM1___310725_G07.fsa"})()
        fsa.ladder_steps = np.asarray([139.0, 250.0, 350.0], dtype=float)
        fsa.best_size_standard = np.asarray([100.0, 200.0, 300.0], dtype=float)
        fsa.size_standard = np.asarray([0.0] * 400, dtype=float)
        fsa.size_standard[100] = 180.0
        fsa.size_standard[200] = 220.0
        fsa.size_standard[300] = 200.0
        fsa.sample_data_with_basepairs = pd.DataFrame(
            {
                "time": np.asarray([100.0, 200.0, 300.0], dtype=float),
                "basepairs": np.asarray([139.4, 249.2, 349.4], dtype=float),
            }
        )
        entry = {
            "fsa": fsa,
            "group": "reactive_control",
            "assay": "NPM1",
            "specimen_id": "IVS-0000",
            "run_date": "2026-03-16",
            "source_run_dir": "0623",
            "run_name": "Run_3730DNA",
            "analysis_type": "standard",
            "selected_injection": "3s",
            "primary_peak_channel": "DATA3",
            "peaks_by_channel": {"DATA3": peaks},
            "ratio": 0.12,
            "ladder_qc_status": "ok",
            "peak_qc_status": "ok",
        }

        with TemporaryDirectory() as tmp:
            excel_path = Path(tmp) / FLT3_NPM1_QC_TRACKER_FILENAME
            update_flt3_npm1_qc_tracker_workbook(excel_path, [entry])
            sheet_names = pd.ExcelFile(excel_path, engine="openpyxl").sheet_names
            control_df = pd.read_excel(excel_path, sheet_name="Control_Runs", engine="openpyxl")
            patient_df = pd.read_excel(excel_path, sheet_name="Patient_Runs", engine="openpyxl")
            marker_df = pd.read_excel(excel_path, sheet_name="PK_Peaks", engine="openpyxl")

        self.assertIn("Dashboard", sheet_names)
        self.assertEqual(set(sheet_names) & {"Patient_Runs", "Control_Runs", "PK_Peaks"}, {"Patient_Runs", "Control_Runs", "PK_Peaks"})
        self.assertTrue(patient_df.empty)
        self.assertEqual(control_df.iloc[0]["Control"], "RK")
        self.assertEqual(control_df.iloc[0]["Assay"], "NPM1")
        self.assertIn("IVS0000_NPM1_299", set(marker_df["MarkerName"].tolist()))
        self.assertIn("NPM1_Ladder_350", set(marker_df["MarkerName"].tolist()))
        self.assertEqual(float(marker_df.loc[marker_df["MarkerName"] == "IVS0000_NPM1_299", "Area"].iloc[0]), 24000.0)

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
             patch("core.analyses.flt3.pipeline.update_flt3_npm1_qc_tracker_workbook") as mock_tracker, \
             patch("core.analyses.flt3.pipeline.finalize_pipeline_run", return_value=["done"]), \
             patch.dict("config.APP_SETTINGS", {"analyses": {"flt3": {"batch": {"tracking_excel_path": ""}}}}, clear=False):
            result = run_pipeline(fsa_dir, return_entries=True, make_dit_reports=False)

        self.assertEqual(result, ["done"])
        mock_trends.assert_called_once_with(assay_dir / FLT3_QC_TRENDS_FILENAME, [selected_entry])
        mock_tracker.assert_called_once_with(assay_dir / FLT3_NPM1_QC_TRACKER_FILENAME, [selected_entry])


if __name__ == "__main__":
    unittest.main()
