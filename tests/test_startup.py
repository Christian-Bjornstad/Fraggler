import contextlib
import importlib
import io
import time
import unittest


class TestStartup(unittest.TestCase):
    def test_fraggler_import_has_no_startup_banner(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            module = importlib.import_module("fraggler.fraggler")
            importlib.reload(module)
        self.assertNotIn("Starting fraggler, importing libraries", buf.getvalue())

    def test_qt_entry_defaults_to_no_legacy_panel_server(self):
        import qt_app

        self.assertFalse(qt_app.LEGACY_PANEL_ENABLED)
        self.assertEqual(qt_app.LEGACY_PANEL_PORT, 5078)

    def test_core_pipeline_import_smoke(self):
        start = time.perf_counter()
        module = importlib.import_module("core.pipeline")
        elapsed = time.perf_counter() - start
        self.assertTrue(hasattr(module, "run_pipeline"))
        self.assertLess(elapsed, 5.0)


if __name__ == "__main__":
    unittest.main()
