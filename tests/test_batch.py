import unittest
import os
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch

from core.batch import generate_jobs, group_files_by_patient
from core.utils import is_water_file


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

    def test_generate_jobs_excludes_v_and_vann_water_files(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            folder = base / "run_a"
            folder.mkdir()

            for path in [
                folder / "25OUM10166_sample_FR1.fsa",
                folder / "V.fsa",
                folder / "Vann_kontroll.fsa",
            ]:
                path.write_text("")

            jobs = generate_jobs([base], aggregate_patients=True)

        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0]["name"], "25OUM10166")
        self.assertEqual([p.name for p in jobs[0]["files"]], ["25OUM10166_sample_FR1.fsa"])

    def test_generate_jobs_non_aggregated_excludes_water_only_folders(self):
        with TemporaryDirectory() as tmp:
            base = Path(tmp)
            patient = base / "patient_run"
            water = base / "water_run"
            patient.mkdir()
            water.mkdir()

            (patient / "25OUM10166_sample_FR1.fsa").write_text("")
            (water / "Vann.fsa").write_text("")

            jobs = generate_jobs([base], aggregate_patients=False)

        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0]["name"], "patient_run")
        self.assertEqual([p.name for p in jobs[0]["files"]], ["25OUM10166_sample_FR1.fsa"])

    def test_group_files_by_patient_skips_water_files(self):
        files = [
            Path("25OUM10166_sample_FR1.fsa"),
            Path("V.fsa"),
            Path("Vann_blank.fsa"),
            Path("PK_sample_FR1.fsa"),
        ]

        grouped = group_files_by_patient(files, r"\d{2}OUM\d{5}")

        self.assertEqual([p.name for p in grouped["25OUM10166"]], ["25OUM10166_sample_FR1.fsa"])
        self.assertEqual([p.name for p in grouped["QC"]], ["PK_sample_FR1.fsa"])
        self.assertEqual(set(grouped.keys()), {"25OUM10166", "QC"})

    def test_is_water_file_matches_v_and_vann_patterns(self):
        self.assertTrue(is_water_file("V.fsa"))
        self.assertTrue(is_water_file("V_001.fsa"))
        self.assertTrue(is_water_file("Vann_kontroll.fsa"))
        self.assertFalse(is_water_file("25OUM10166_sample_FR1.fsa"))

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

    def test_run_batch_jobs_uses_single_collect_pass_for_flt3_aggregated_jobs(self):
        from core.batch import run_batch_jobs

        jobs = [
            {"name": "patient", "type": "pipeline", "path": None, "files": [Path("a.fsa"), Path("b.fsa")]},
        ]

        with patch("core.batch.run_pipeline_job_collect", return_value=[{"ok": True}]) as mock_collect, \
             patch("config.APP_SETTINGS", {"active_analysis": "flt3", "qc": {}}):
            result = run_batch_jobs(
                jobs=jobs,
                output_base=Path("/tmp/out"),
                out_folder_tmpl="ASSAY_REPORTS",
                outfile_html_tmpl="QC_REPORT_{name}.html",
                excel_name_tmpl="Fraggler_QC_Trends.xlsx",
                pipeline_scope="all",
                assay_filter="",
                aggregate_dit_reports=True,
                continue_on_error=True,
                update_callback=None,
            )

        self.assertEqual(result["failed_jobs"], [])
        self.assertEqual(result["completed_jobs"], ["patient"])
        self.assertFalse(mock_collect.call_args.kwargs["chunk_files"])


class TestQtBatchTableScrolling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from PyQt6.QtWidgets import QApplication

        cls._app = QApplication.instance() or QApplication([])

    def test_detected_jobs_table_wheel_scrolls(self):
        from PyQt6.QtCore import QPoint, QPointF, Qt
        from PyQt6.QtGui import QWheelEvent
        from PyQt6.QtWidgets import QApplication
        from gui_qt.tabs.tab_batch import TabBatch

        widget = TabBatch()
        widget.resize(1200, 800)
        widget._detected_jobs = [
            {"name": f"job_{i:03d}", "type": "pipeline", "path": None, "files": [Path("a.fsa")]}
            for i in range(120)
        ]
        widget._job_states = {job["name"]: "pending" for job in widget._detected_jobs}
        widget._rebuild_table()
        widget.show()
        QApplication.processEvents()

        bar = widget.table.verticalScrollBar()
        self.assertGreater(bar.maximum(), 0)
        start = bar.value()

        event = QWheelEvent(
            QPointF(20.0, 20.0),
            QPointF(20.0, 20.0),
            QPoint(0, 0),
            QPoint(0, -120),
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
            Qt.ScrollPhase.ScrollUpdate,
            False,
        )
        QApplication.sendEvent(widget.table.viewport(), event)
        QApplication.processEvents()

        self.assertTrue(event.isAccepted())
        self.assertGreater(bar.value(), start)

        widget.close()


if __name__ == "__main__":
    unittest.main()
