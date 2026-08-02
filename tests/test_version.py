import re

from version import APP_VERSION


def test_app_version_uses_semantic_versioning():
    assert re.fullmatch(r"\d+\.\d+\.\d+", APP_VERSION)
