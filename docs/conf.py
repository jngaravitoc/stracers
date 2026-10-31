import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

import stracers  # noqa: E402

project = "stracers"
author = "Nico Garavito-Camargo and Chervin Laporte"
copyright = "2026, " + author
release = stracers.__version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

html_theme = "sphinx_rtd_theme"
exclude_patterns = ["_build"]
