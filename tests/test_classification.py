import unittest
from pathlib import Path
from core.classification import detect_assay, strip_stage_prefix, classify_fsa
from config import APP_SETTINGS

class TestClassification(unittest.TestCase):
    def setUp(self):
        self._active_analysis = APP_SETTINGS.get("active_analysis", "clonality")
        APP_SETTINGS["active_analysis"] = "clonality"

    def tearDown(self):
        APP_SETTINGS["active_analysis"] = self._active_analysis

    def test_detect_assay_sl(self):
        self.assertEqual(detect_assay("sample_sl.fsa"), "SL")
        self.assertEqual(detect_assay("sample_SIZELADDER.fsa"), "SL")
        self.assertEqual(detect_assay("sample-size-ladder.fsa"), "SL")

    def test_detect_assay_tcrb_mix(self):
        self.assertEqual(detect_assay("TCRB mix A.fsa"), "TCRbA")
        self.assertEqual(detect_assay("TRB_Mix_B.fsa"), "TCRbB")
        self.assertEqual(detect_assay("tcrb-mix-c.fsa"), "TCRbC")

    def test_detect_assay_tcrg_mix(self):
        self.assertEqual(detect_assay("TCRG mix A.fsa"), "TCRgA")
        self.assertEqual(detect_assay("TRG_Mix_B.fsa"), "TCRgB")

    def test_detect_assay_igh(self):
        self.assertEqual(detect_assay("sample_FR1.fsa"), "FR1")
        self.assertEqual(detect_assay("sample_FR2.fsa"), "FR2")
        self.assertEqual(detect_assay("sample_FR3.fsa"), "FR3")

    def test_detect_assay_dhjh(self):
        self.assertEqual(detect_assay("sample_DHJH_D.fsa"), "DHJH_D")
        self.assertEqual(detect_assay("sample_DHJHmixE.fsa"), "DHJH_E")

    def test_detect_assay_liz(self):
        self.assertEqual(detect_assay("sample_IGK.fsa"), "IGK")
        self.assertEqual(detect_assay("sample_KDE.fsa"), "KDE")

    def test_strip_stage_prefix(self):
        filename = "12345_abcdef12_sample_FR1.fsa"
        self.assertEqual(strip_stage_prefix(filename), "sample_FR1.fsa")

    def test_classify_fsa_pk_prefix(self):
        p = Path("PK1_sample_FR1.fsa")
        res = classify_fsa(p)
        self.assertIsNotNone(res)
        assay, group, ladder, trace_channels, peak_channels, primary_ch, bp_min, bp_max = res
        self.assertEqual(assay, "FR1")
        self.assertEqual(group, "positive")
        self.assertEqual(ladder, "ROX")

    def test_classify_fsa_nk_prefix(self):
        p = Path("NK_sample_IGK.fsa")
        res = classify_fsa(p)
        self.assertIsNotNone(res)
        self.assertEqual(res[0], "IGK")
        self.assertEqual(res[1], "negative")
        self.assertEqual(res[2], "LIZ")

    def test_detect_assay_dispatch_for_flt3(self):
        APP_SETTINGS["active_analysis"] = "flt3"
        self.assertEqual(detect_assay("patient_itd_ratio_p1.fsa"), "FLT3-ITD")

    def test_classify_fsa_dispatch_for_flt3(self):
        APP_SETTINGS["active_analysis"] = "flt3"
        from unittest.mock import patch
        from core.analyses.flt3 import classification as flt3_classification

        with patch.object(flt3_classification, "get_injection_metadata", return_value={"injection_time": 3, "injection_voltage": 15}):
            res = classify_fsa(Path("ivs-p001_itd_p1_ufort.fsa"))
        self.assertIsNotNone(res)
        self.assertEqual(res["group"], "positive_control")
        self.assertEqual(res["parallel"], "p1")
        self.assertEqual(res["injection_time"], 3)

    def test_flt3_npm1_filename_dispatch(self):
        APP_SETTINGS["active_analysis"] = "flt3"
        self.assertEqual(detect_assay("26OUM00042_npm1_p2.fsa"), "NPM1")

    def test_flt3_parallel_regex_avoids_false_matches(self):
        APP_SETTINGS["active_analysis"] = "flt3"
        from unittest.mock import patch
        from core.analyses.flt3 import classification as flt3_classification

        with patch.object(flt3_classification, "get_injection_metadata", return_value={"injection_time": 3, "injection_voltage": 15}):
            res = classify_fsa(Path("26OUM00042_npm1_run2026.fsa"))
        self.assertIsNotNone(res)
        self.assertIsNone(res["parallel"])

    def test_flt3_protocol_injection_uses_preferred_run_type(self):
        APP_SETTINGS["active_analysis"] = "flt3"
        from unittest.mock import patch
        from core.analyses.flt3 import classification as flt3_classification

        with patch.object(flt3_classification, "get_injection_metadata", return_value={"injection_time": 3, "injection_voltage": 15}):
            res = classify_fsa(Path("25OUM11314_p1_ITD_1-10__310725_A02.fsa"))
        self.assertEqual(res["injection_time"], 3)
        self.assertEqual(res["protocol_injection_time"], 1)

if __name__ == "__main__":
    unittest.main()
