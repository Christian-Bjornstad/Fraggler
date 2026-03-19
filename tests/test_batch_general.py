import copy
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from config import APP_SETTINGS
from core.batch import generate_jobs, run_batch_jobs
from core.runner import run_pipeline_job


class TestBatchGeneral(unittest.TestCase):
    def setUp(self):
        self._settings_backup = copy.deepcopy(APP_SETTINGS)
        APP_SETTINGS.setdefault("qc", {})
        APP_SETTINGS.setdefault("analyses", {})
        APP_SETTINGS["active_analysis"] = "general"
        APP_SETTINGS["analyses"]["general"] = {
            "batch": {
                "base_input_dir": "/tmp/input",
                "output_base": "/tmp/output",
                "aggregate_by_patient": True,
                "patient_id_regex": r"\d{2}OUM\d{5}",
                "aggregate_dit_reports": True,
            },
            "pipeline": {
                "mode": "all",
                "assay_filter_substring": "",
            },
        }

    def tearDown(self):
        APP_SETTINGS.clear()
        APP_SETTINGS.update(self._settings_backup)

    def test_generate_jobs_general_avoids_patient_aggregation(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp) / "input"
            a = root / "a"
            b = root / "b"
            a.mkdir(parents=True)
            b.mkdir(parents=True)
            (a / "one.fsa").write_text("", encoding="utf-8")
            (b / "two.fsa").write_text("", encoding="utf-8")

            jobs = generate_jobs([root], aggregate_patients=True)

        self.assertEqual([j["type"] for j in jobs], ["pipeline", "pipeline"])
        self.assertTrue(all(j["files"] == [] for j in jobs))

    def test_run_batch_jobs_disables_dit_aggregation_for_general(self):
        jobs = [{"name": "a", "type": "pipeline", "path": Path("/tmp/input/a"), "files": []}]

        with patch("core.batch.run_pipeline_job") as mock_run, \
             patch("core.batch.run_pipeline_job_collect") as mock_collect, \
             patch("core.batch.run_qc_job"), \
             patch("core.batch.run_dit_job"):
            run_batch_jobs(
                jobs=jobs,
                output_base=Path("/tmp/output"),
                out_folder_tmpl="ASSAY_REPORTS",
                outfile_html_tmpl="QC_REPORT_{name}.html",
                excel_name_tmpl="Fraggler_QC_TRENDS_{name}.xlsx",
                pipeline_scope="all",
                assay_filter="",
                aggregate_dit_reports=True,
                continue_on_error=True,
                update_callback=None,
            )

        mock_run.assert_called_once()
        mock_collect.assert_not_called()

    def test_run_pipeline_job_general_stages_explicit_files_once(self):
        files = [Path("/tmp/input/a.fsa"), Path("/tmp/input/b.fsa")]
        staged_dir = Path("/tmp/general-staged")

        with patch("core.runner.stage_files", return_value=staged_dir) as mock_stage, \
             patch("core.runner.cleanup_temp") as mock_cleanup, \
             patch("core.pipeline.run_pipeline") as mock_run:
            run_pipeline_job(
                fsa_dir=None,
                base_outdir=Path("/tmp/output"),
                out_folder_name="GENERAL_REPORTS",
                scope="all",
                needle="",
                files=files,
            )

        mock_stage.assert_called_once_with(files)
        mock_run.assert_called_once_with(
            fsa_dir=staged_dir,
            base_outdir=Path("/tmp/output"),
            assay_folder_name="GENERAL_REPORTS",
            mode="all",
        )
        mock_cleanup.assert_called_once_with(staged_dir)
