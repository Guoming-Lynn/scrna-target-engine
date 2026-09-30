"""legacy stays out of tests and out of the published wheel."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_legacy_is_not_collected_or_packaged() -> None:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'testpaths = ["tests"]' in text
    assert 'norecursedirs = ["legacy"]' in text
    assert 'only-include = ["scripts"]' in text
    assert "legacy" not in text.split("only-include", 1)[1]
