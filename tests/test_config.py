import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch

import yaml

import config
from config import get_analysis_settings, load_settings, save_settings


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

    def test_load_settings_migrates_legacy_batch_and_pipeline_into_analysis_profiles(self):
        with TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / "settings.yaml"
            cfg_path.write_text(
                yaml.safe_dump(
                    {
                        "batch": {
                            "base_input_dir": "/legacy/in",
                            "output_base": "/legacy/out",
                            "aggregate_by_patient": False,
                            "patient_id_regex": "ABC",
                            "aggregate_dit_reports": False,
                        },
                        "pipeline": {
                            "mode": "custom",
                            "assay_filter_substring": "FLT3",
                        },
                    }
                ),
                encoding="utf-8",
            )

            settings = load_settings(cfg_path, env={})

        for analysis_id in ("clonality", "flt3"):
            profile = settings["analyses"][analysis_id]
            self.assertEqual(profile["batch"]["base_input_dir"], "/legacy/in")
            self.assertEqual(profile["batch"]["output_base"], "/legacy/out")
            self.assertFalse(profile["batch"]["aggregate_by_patient"])
            self.assertEqual(profile["batch"]["patient_id_regex"], "ABC")
            self.assertFalse(profile["batch"]["aggregate_dit_reports"])
            self.assertEqual(profile["pipeline"]["mode"], "custom")
            self.assertEqual(profile["pipeline"]["assay_filter_substring"], "FLT3")

    def test_get_analysis_settings_returns_analysis_specific_profile(self):
        settings = load_settings(Path("/does/not/exist.yaml"), env={})
        settings["analyses"]["flt3"]["batch"]["base_input_dir"] = "/flt3/in"
        settings["analyses"]["clonality"]["batch"]["base_input_dir"] = "/clonality/in"

        self.assertEqual(
            get_analysis_settings("flt3", settings)["batch"]["base_input_dir"],
            "/flt3/in",
        )
        self.assertEqual(
            get_analysis_settings("clonality", settings)["batch"]["base_input_dir"],
            "/clonality/in",
        )

    def test_load_settings_warns_and_tracks_error_on_invalid_yaml(self):
        with TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / "settings.yaml"
            cfg_path.write_text("{invalid", encoding="utf-8")

            with self.assertWarns(RuntimeWarning):
                settings = load_settings(cfg_path, env={})

        self.assertEqual(settings["active_analysis"], "clonality")
        self.assertIsNotNone(config.LAST_SETTINGS_LOAD_ERROR)

    def test_save_settings_warns_and_tracks_error_on_write_failure(self):
        with TemporaryDirectory() as tmp:
            cfg_path = Path(tmp) / "settings.yaml"
            settings = load_settings(Path("/does/not/exist.yaml"), env={})

            with patch("builtins.open", side_effect=OSError("disk full")):
                with self.assertWarns(RuntimeWarning):
                    save_settings(settings, cfg_path)

        self.assertIsNotNone(config.LAST_SETTINGS_SAVE_ERROR)

    def test_invalid_active_analysis_falls_back_to_default(self):
        settings = load_settings(
            Path("/does/not/exist.yaml"),
            env={"FRAGGLER_ACTIVE_ANALYSIS": "does-not-exist"},
        )
        self.assertEqual(settings["active_analysis"], "clonality")


if __name__ == "__main__":
    unittest.main()
