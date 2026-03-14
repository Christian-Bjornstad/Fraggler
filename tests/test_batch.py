import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

from core.batch import generate_jobs, group_files_by_patient


class TestBatch(unittest.TestCase):
    def test_group_files_by_patient_separates_qc(self):
        files = [
            Path("25OUM10166_sample_FR1.fsa"),
            Path("25OUM10166_sample_FR2.fsa"),
            Path("PK_sample_FR1.fsa"),
        ]
        grouped = group_files_by_patient(files, r"\d{2}OUM\d{5}")
        self.assertEqual(len(grouped["25OUM10166"]), 2)
        self.assertEqual([p.name for p in grouped["QC"]], ["PK_sample_FR1.fsa"])

    def test_generate_jobs_aggregates_patient_and_qc_jobs(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            folder_a = base / "run_a"
            folder_b = base / "run_b"
            folder_a.mkdir()
            folder_b.mkdir()

            for path in [
                folder_a / "25OUM10166_sample_FR1.fsa",
                folder_b / "25OUM10166_sample_FR2.fsa",
                folder_b / "NK_sample_IGK.fsa",
            ]:
                path.write_text("")

            jobs = generate_jobs([base], aggregate_patients=True)
            names = {(job["name"], job["type"]) for job in jobs}
            self.assertIn(("25OUM10166", "pipeline"), names)
            self.assertIn(("QC", "qc"), names)


if __name__ == "__main__":
    unittest.main()
