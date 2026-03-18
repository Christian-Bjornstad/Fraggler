import importlib
import importlib.util
import os
import sys
import unittest
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch


class TestBuildQt(unittest.TestCase):
    def test_build_pyinstaller_args_include_runtime_hook_hidden_imports_and_datas(self):
        import build_qt

        with patch.object(build_qt.sys, "platform", "linux"), \
             patch.object(build_qt, "_collect_linux_binaries", return_value=["--add-binary=/tmp/lib.so:."]):
            args = build_qt._build_pyinstaller_args()

        self.assertIn("qt_app.py", args)
        self.assertIn(f"--runtime-hook={build_qt.HOOK_DIR / 'runtime_desktop.py'}", args)
        self.assertIn("--hidden-import=core.analyses.clonality.pipeline", args)
        self.assertIn("--hidden-import=core.analyses.flt3.pipeline", args)
        self.assertIn(build_qt._format_data_arg("assets", "assets"), args)
        self.assertIn(build_qt._format_data_arg("app.py", "."), args)
        self.assertIn("--add-binary=/tmp/lib.so:.", args)

    def test_zip_path_preserves_root_name_for_directory_contents(self):
        import build_qt

        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = root / "Fraggler"
            (src / "_internal").mkdir(parents=True)
            (src / "Fraggler").write_text("launcher", encoding="utf-8")
            (src / "_internal" / "qt.conf").write_text("[Paths]\nPrefix = .\n", encoding="utf-8")
            zip_path = root / "bundle.zip"

            build_qt._zip_path(src, zip_path, root_name="Fraggler_Linux")

            with zipfile.ZipFile(zip_path) as zf:
                names = sorted(zf.namelist())

        self.assertIn("Fraggler_Linux/Fraggler", names)
        self.assertIn("Fraggler_Linux/_internal/qt.conf", names)

    def test_post_build_linux_stages_release_bundle_with_readme_and_launcher(self):
        import build_qt

        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            dist_dir = root / "dist"
            release_dir = dist_dir / "releases"
            source_dir = dist_dir / build_qt.APP_NAME
            source_dir.mkdir(parents=True)
            launcher = source_dir / build_qt.APP_NAME
            launcher.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            (source_dir / "_internal").mkdir()
            legacy_guide = root / "LINUX_GUIDE.md"
            legacy_guide.write_text("legacy guide", encoding="utf-8")

            with patch.object(build_qt, "DIST_DIR", dist_dir), \
                 patch.object(build_qt, "RELEASE_DIR", release_dir), \
                 patch.object(build_qt, "LEGACY_LINUX_GUIDE", legacy_guide):
                build_qt._post_build_linux()

            staged_dir = dist_dir / "Fraggler_Linux"
            staged_zip = release_dir / "Fraggler_Linux_offline.zip"

            self.assertTrue((staged_dir / "README.txt").exists())
            self.assertTrue((staged_dir / "LINUX_GUIDE.md").exists())
            self.assertTrue((staged_dir / build_qt.APP_NAME).exists())
            self.assertTrue(os.access(staged_dir / build_qt.APP_NAME, os.X_OK))
            self.assertTrue(staged_zip.exists())

            with zipfile.ZipFile(staged_zip) as zf:
                names = set(zf.namelist())

        self.assertIn("Fraggler_Linux/README.txt", names)
        self.assertIn("Fraggler_Linux/LINUX_GUIDE.md", names)
        self.assertIn(f"Fraggler_Linux/{build_qt.APP_NAME}", names)

    def test_runtime_hook_disables_legacy_panel_by_default(self):
        module_name = "runtime_desktop_test"
        module_path = Path("/Users/christian/Desktop/OUS/packaging/hooks/runtime_desktop.py")
        old_panel = os.environ.pop("FRAGGLER_ENABLE_LEGACY_PANEL", None)
        old_qpa = os.environ.pop("QT_QPA_PLATFORM", None)

        try:
            with patch.object(sys, "platform", "darwin"):
                spec = importlib.util.spec_from_file_location(module_name, module_path)
                module = importlib.util.module_from_spec(spec)
                assert spec.loader is not None
                spec.loader.exec_module(module)
            self.assertIsNotNone(module)
            self.assertEqual(os.environ.get("FRAGGLER_ENABLE_LEGACY_PANEL"), "0")
            self.assertNotIn("QT_QPA_PLATFORM", os.environ)
        finally:
            sys.modules.pop(module_name, None)
            if old_panel is None:
                os.environ.pop("FRAGGLER_ENABLE_LEGACY_PANEL", None)
            else:
                os.environ["FRAGGLER_ENABLE_LEGACY_PANEL"] = old_panel
            if old_qpa is None:
                os.environ.pop("QT_QPA_PLATFORM", None)
            else:
                os.environ["QT_QPA_PLATFORM"] = old_qpa


if __name__ == "__main__":
    unittest.main()
