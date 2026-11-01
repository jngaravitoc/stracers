# Sphinx configuration for stracers, following nba's docs (sphinx-astropy
# defaults, pydata-sphinx-theme, automodapi API pages). Build with `make html`.
import csv
import pathlib
import shutil
import sys

try:
    from sphinx_astropy.conf.v1 import *  # noqa: F401,F403
except ImportError:
    print("ERROR: building the stracers documentation requires sphinx-astropy: "
          "pip install -e '.[docs]'")
    sys.exit(1)

import stracers
from stracers import mw

docs_root = pathlib.Path(__file__).parent.resolve()

project = "stracers"
author = "Nicolas Garavito-Camargo and Chervin Laporte"
copyright = "2026, Nicolas Garavito-Camargo and Chervin Laporte"
version = release = stracers.__version__

highlight_language = "python3"
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
templates_path = ["_templates"]
default_role = "obj"

# API pages
numpydoc_show_class_members = False
numpydoc_xref_param_type = True
automodapi_toctreedirnm = "api"
automodsumm_inherited_members = False

intersphinx_mapping.update({  # noqa: F405
    "h5py": ("https://docs.h5py.org/en/stable/", None),
})

# Notebooks and Markdown pages with myst-nb. The tutorial needs no simulation
# data, so it is run when the documentation is built.
extensions += ["myst_nb"]  # noqa: F405
source_suffix = {".rst": "restructuredtext", ".md": "myst-nb", ".ipynb": "myst-nb"}
nb_execution_mode = "auto"
nb_execution_timeout = 300
myst_enable_extensions = ["dollarmath"]

# The tutorials live in tutorials/ at the top of the repository; copy the ones
# shown in the docs next to tutorials.rst at build time (ignored by git).
TUTORIALS = ["hernquist_tracers.ipynb"]
(docs_root / "tutorials").mkdir(exist_ok=True)
for name in TUTORIALS:
    shutil.copy(docs_root.parent / "tutorials" / name, docs_root / "tutorials" / name)

# Tables of the Milky Way profiles shown in mw/index.rst, written from the
# package data at build time (ignored by git).
def _fmt(entry, name):
    value = entry[name]
    if value is None:
        return ""
    lo, hi = entry.get(name + "_err_minus"), entry.get(name + "_err_plus")
    if lo is None:
        return "{:g}".format(value)
    if lo == hi:
        return "{:g} ± {:g}".format(value, lo)
    return ":math:`{:g}^{{+{:g}}}_{{-{:g}}}`".format(value, hi, lo)


def _range(entry):
    if entry["r_min"] is None:
        return "<{:g}".format(entry["r_max"])
    return "{:g}–{:g}".format(entry["r_min"], entry["r_max"])


with open(docs_root / "mw" / "spl.csv", "w", newline="") as f_spl, \
        open(docs_root / "mw" / "bpl.csv", "w", newline="") as f_bpl:
    spl, bpl = csv.writer(f_spl), csv.writer(f_bpl)
    for e in mw.load_table():
        if e["model"] == "spl":
            spl.writerow([e["key"], _fmt(e, "slope"), _range(e), e["tracer"], e["reference"]])
        else:
            bpl.writerow([e["key"], _fmt(e, "r_break"), _fmt(e, "inner_slope"),
                          _fmt(e, "outer_slope"), _range(e), e["tracer"], e["reference"]])

# HTML
html_theme = "pydata_sphinx_theme"
html_title = f"stracers v{version}"
html_baseurl = "https://jngaravitoc.github.io/tracers_triaxial/"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "github_url": "https://github.com/jngaravitoc/tracers_triaxial",
    "navbar_align": "left",
    "show_toc_level": 2,
}
html_sidebars = {"install": [], "getting_started": [], "changes": []}
