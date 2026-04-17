from pathlib import Path


REPO_ROOT = Path("/Users/christian/Desktop/OUS")


def test_third_party_notice_mentions_upstream_and_license():
    notice = (REPO_ROOT / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")

    assert "willros/fraggler" in notice
    assert "MIT" in notice
    assert "fraggler/" in notice
    assert "core/analysis.py" in notice


def test_upstream_mit_license_file_matches_expected_header():
    license_text = (REPO_ROOT / "LICENSES" / "fraggler_MIT.txt").read_text(encoding="utf-8")

    assert license_text.startswith("MIT License")
    assert "Clinical Genomic Umea" in license_text
    assert "Permission is hereby granted" in license_text


def test_readme_points_to_notice_and_license_files():
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")

    assert "THIRD_PARTY_NOTICES.md" in readme
    assert "LICENSES/fraggler_MIT.txt" in readme
    assert "willros/fraggler" in readme
