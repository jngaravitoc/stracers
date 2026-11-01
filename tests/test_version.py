import pathlib
import re
from importlib.metadata import version

import stracers

CHANGELOG = pathlib.Path(__file__).parents[1] / "CHANGELOG.md"


def test_metadata_matches_package():
    """pyproject.toml reads the version from stracers.__version__."""
    assert version("stracers") == stracers.__version__


def test_changelog_has_current_version():
    """The newest dated release in CHANGELOG.md is the package version."""
    releases = re.findall(r"^## \[(\d+\.\d+\.\d+)\] - \d{4}-\d{2}-\d{2}$", CHANGELOG.read_text(), re.M)
    assert releases[0] == stracers.__version__
