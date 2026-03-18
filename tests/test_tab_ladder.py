import os
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication

from gui_qt.tabs.tab_ladder import TabLadder


_APP = QApplication.instance() or QApplication([])


class TestTabLadder(unittest.TestCase):
    def test_source_file_list_uses_compact_height(self):
        widget = TabLadder()
        self.assertEqual(widget.file_list.minimumHeight(), 220)

    def test_file_selection_schedules_background_metadata_load(self):
        widget = TabLadder()
        started = []

        with patch.object(widget.threadpool, "start", side_effect=lambda worker: started.append(worker)), \
             patch("gui_qt.tabs.tab_ladder.detect_fsa_for_ladder") as mock_detect, \
             patch("gui_qt.tabs.tab_ladder.load_adjustable_fsa") as mock_load:
            widget._update_current_file(Path("/tmp/example.fsa"))

        self.assertEqual(len(started), 2)
        self.assertTrue(widget._metadata_loading)
        mock_detect.assert_not_called()
        mock_load.assert_not_called()

    def test_set_analysis_does_not_trigger_metadata_reload_for_current_file(self):
        widget = TabLadder()
        widget._current_file = Path("/tmp/example.fsa")

        with patch.object(widget, "_refresh_current_metadata") as mock_refresh:
            widget.set_analysis("flt3")

        mock_refresh.assert_not_called()
        self.assertIsNone(widget._current_meta)
        self.assertIsNone(widget._current_fsa)


if __name__ == "__main__":
    unittest.main()
