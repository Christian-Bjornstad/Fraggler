import unittest
from tempfile import TemporaryDirectory
from pathlib import Path

import yaml

from config import load_settings, save_settings


class TestConfig(unittest.TestCase):
    def test_load_settings_migrates_legacy_default_output(self):
        with TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / "settings.yaml"
            cfg_path.write_text(yaml.safe_dump({"default_output": "/tmp/out"}), encoding="utf-8")

            settings = load_settings(cfg_path, env={})
            self.assertEqual(settings["general"]["default_output"], "/tmp/out")
            self.assertEqual(settings["batch"]["output_base"], "/tmp/out")
            self.assertEqual(settings["pipeline"]["output_base"], "/tmp/out")

    def test_load_settings_applies_nested_env_override(self):
        settings = load_settings(
            Path("/does/not/exist.yaml"),
            env={
                "FRAGGLER_BATCH_OUTPUT_BASE": "/env/out",
                "FRAGGLER_QC_MIN_R2_OK": "0.991",
                "FRAGGLER_GENERAL_AUTHOR": "Lab",
            },
        )
        self.assertEqual(settings["batch"]["output_base"], "/env/out")
        self.assertEqual(settings["general"]["author"], "Lab")
        self.assertAlmostEqual(settings["qc"]["min_r2_ok"], 0.991)

    def test_save_settings_preserves_legacy_key(self):
        with TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / "settings.yaml"
            settings = load_settings(cfg_path, env={})
            settings["general"]["default_output"] = "/saved/out"
            save_settings(settings, cfg_path)
            payload = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["default_output"], "/saved/out")


if __name__ == "__main__":
    unittest.main()
