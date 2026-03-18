import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch

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

    def test_run_batch_jobs_returns_failed_summary(self):
        from core.batch import run_batch_jobs

        jobs = [
            {"name": "ok", "type": "pipeline", "path": None, "files": [Path("a.fsa")]},
            {"name": "bad", "type": "pipeline", "path": None, "files": [Path("b.fsa")]},
        ]
        updates = []

        def _pipeline_job(**kwargs):
            if kwargs["files"] == [Path("b.fsa")]:
                raise RuntimeError("boom")

        with patch("core.batch.run_pipeline_job", side_effect=_pipeline_job):
            result = run_batch_jobs(
                jobs=jobs,
                output_base=Path("/tmp/out"),
                out_folder_tmpl="ASSAY_REPORTS",
                outfile_html_tmpl="QC_REPORT_{name}.html",
                excel_name_tmpl="Fraggler_QC_Trends.xlsx",
                pipeline_scope="all",
                assay_filter="",
                aggregate_dit_reports=False,
                continue_on_error=True,
                update_callback=lambda *args: updates.append(args),
            )

        self.assertEqual(result["failed_jobs"], ["bad"])
        self.assertEqual(result["completed_jobs"], ["ok"])
        self.assertIn((1, 2, "ok", "success"), updates)
        self.assertTrue(any(update[2] == "bad" and str(update[3]).startswith("error:") for update in updates))

    def test_run_batch_jobs_marks_unknown_job_type_as_failed(self):
        from core.batch import run_batch_jobs

        updates = []
        result = run_batch_jobs(
            jobs=[{"name": "weird", "type": "unknown", "path": None, "files": []}],
            output_base=Path("/tmp/out"),
            out_folder_tmpl="ASSAY_REPORTS",
            outfile_html_tmpl="QC_REPORT_{name}.html",
            excel_name_tmpl="Fraggler_QC_Trends.xlsx",
            pipeline_scope="all",
            assay_filter="",
            aggregate_dit_reports=False,
            continue_on_error=True,
            update_callback=lambda *args: updates.append(args),
        )

        self.assertEqual(result["completed_jobs"], [])
        self.assertEqual(result["failed_jobs"], ["weird"])
        self.assertTrue(any(update[2] == "weird" and str(update[3]).startswith("error:") for update in updates))


if __name__ == "__main__":
    unittest.main()
