import unittest
from pathlib import Path
from core.classification import detect_assay, strip_stage_prefix, classify_fsa

class TestClassification(unittest.TestCase):
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

if __name__ == "__main__":
    unittest.main()
