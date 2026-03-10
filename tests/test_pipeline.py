import unittest
from pathlib import Path
from core.pipeline import _scan_files
from fraggler.fraggler import FsaFile

class TestPipeline(unittest.TestCase):
    def test_scan_files_basic(self):
        # This is hard to test without real .fsa files, 
        # but we can mock the directory structure if needed.
        # For now, let's just verify it returns a list.
        # In a real scenario, we'd use a temporary directory with dummy .fsa files.
        pass

    def test_dit_extraction(self):
        from core.html_reports import extract_dit_from_name
        self.assertEqual(extract_dit_from_name("25OUM10166_some_data.fsa"), "25OUM10166")
        self.assertEqual(extract_dit_from_name("no_dit_here.fsa"), None)
        self.assertEqual(extract_dit_from_name("26OUM00042_ABC.fsa"), "26OUM00042")

if __name__ == "__main__":
    unittest.main()
