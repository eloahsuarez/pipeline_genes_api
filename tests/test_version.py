import re
from pathlib import Path

from version import APP_VERSION


def test_app_version_uses_semantic_versioning():
    assert re.fullmatch(r"\d+\.\d+\.\d+", APP_VERSION)


def test_readme_declares_the_same_app_version():
    readme = Path(__file__).parents[1].joinpath("README.md").read_text(encoding="utf-8")
    match = re.search(r"Versão atual: \*\*(\d+\.\d+\.\d+)\*\*", readme)

    assert match is not None
    assert match.group(1) == APP_VERSION
