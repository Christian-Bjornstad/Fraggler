import unittest

from core.html_reports import _render_file_summary_table


class _DummyFsa:
    def __init__(self, file_name: str):
        self.file_name = file_name


class TestHtmlReportLadderStatus(unittest.TestCase):
    def test_summary_table_renders_compact_ladder_status_badges(self):
        entries = [
            {
                "fsa": _DummyFsa("01_ok.fsa"),
                "assay": "SL",
                "ladder": "ROX400HD",
                "bp_min": 100.0,
                "bp_max": 500.0,
                "ladder_qc_status": "ok",
                "ladder_r2": 0.99995,
                "ladder_fit_note": "All expected ladder steps were fitted.",
            },
            {
                "fsa": _DummyFsa("02_partial.fsa"),
                "assay": "SL",
                "ladder": "ROX400HD",
                "bp_min": 100.0,
                "bp_max": 500.0,
                "ladder_qc_status": "review_required",
                "ladder_r2": 0.99993,
                "ladder_fit_note": "High-end rescue used the stable top 15/21 ladder steps because the lower ROX region was unreliable.",
            },
            {
                "fsa": _DummyFsa("03_manual.fsa"),
                "assay": "SL",
                "ladder": "ROX400HD",
                "bp_min": 100.0,
                "bp_max": 500.0,
                "ladder_qc_status": "manual_adjustment",
                "ladder_r2": 0.99999,
                "ladder_fit_note": "Manual ladder adjustment applied from saved sidecar.",
            },
        ]

        html_lines: list[str] = []
        _render_file_summary_table(entries, html_lines)
        html = "\n".join(html_lines)

        self.assertIn(">OK</span>", html)
        self.assertIn(">Warning</span>", html)
        self.assertIn(">Manual</span>", html)
        self.assertNotIn("review_required", html)
        self.assertNotIn("manual_adjustment", html)
        self.assertIn("lower ROX region was unreliable", html)


if __name__ == "__main__":
    unittest.main()
