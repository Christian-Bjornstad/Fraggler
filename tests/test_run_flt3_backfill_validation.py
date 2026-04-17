from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts.run_flt3_backfill_validation import (
    FLT3_NPM1_QC_TRACKER_FILENAME,
    _archive_metric_paths,
    build_tracking_workbook,
)


def test_archive_metric_paths_dedupes_and_skips_missing_files():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        existing = root / "A" / "sample1.fsa"
        existing.parent.mkdir(parents=True, exist_ok=True)
        existing.write_text("ok", encoding="utf-8")
        metrics_csv = root / "flt3_ladder_metrics.csv"
        metrics_csv.write_text(
            "path,status\n"
            f"{existing},ok\n"
            f"{existing},review_required\n"
            f"{root / 'missing.fsa'},analysis_error\n",
            encoding="utf-8",
        )

        paths = _archive_metric_paths(metrics_csv)

    assert paths == [existing]


def test_build_tracking_workbook_replays_archive_files_into_tracker():
    with TemporaryDirectory() as tmp:
        run_dir = Path(tmp)
        file_a = run_dir / "folderA" / "a.fsa"
        file_b = run_dir / "folderA" / "b.fsa"
        file_c = run_dir / "folderB" / "c.fsa"
        for path in (file_a, file_b, file_c):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("ok", encoding="utf-8")

        metrics_csv = run_dir / "flt3_ladder_metrics.csv"
        metrics_csv.write_text(
            "path,status\n"
            f"{file_a},ok\n"
            f"{file_b},review_required\n"
            f"{file_c},ok\n",
            encoding="utf-8",
        )

        captured_calls = []

        def _fake_collect(**kwargs):
            captured_calls.append(kwargs)
            tracker_path = Path(kwargs["tracking_excel_path"])
            tracker_path.parent.mkdir(parents=True, exist_ok=True)
            tracker_path.write_text("tracker", encoding="utf-8")
            return []

        with patch("scripts.run_flt3_backfill_validation.run_pipeline_job_collect", side_effect=_fake_collect):
            tracker_path = build_tracking_workbook(run_dir)

        assert tracker_path == run_dir / FLT3_NPM1_QC_TRACKER_FILENAME
        assert tracker_path.exists()
        assert len(captured_calls) == 2
        assert {call["fsa_dir"] for call in captured_calls} == {file_a.parent, file_c.parent}
        assert all(call["tracking_excel_path"] == tracker_path for call in captured_calls)
        assert all(call["update_tracking_workbook"] is True for call in captured_calls)
        assert all(call["chunk_files"] is False for call in captured_calls)
