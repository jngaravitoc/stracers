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
TUTORIALS = ["hernquist_tracers.ipynb", "tracer_kinematics.ipynb", "simulation_tracers.ipynb"]
# the simulation tutorial needs local simulation data: it is committed with its
# outputs and shown as it is
nb_execution_excludepatterns = ["tutorials/simulation_tracers.ipynb"]
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


def _link(entry):
    # anonymous hyperlink: the same reference text appears in several rows
    return "`{} <{}>`__".format(entry["reference"], mw.reference_url(entry["key"]))


with open(docs_root / "mw" / "spl.csv", "w", newline="") as f_spl, \
        open(docs_root / "mw" / "bpl.csv", "w", newline="") as f_bpl, \
        open(docs_root / "mw" / "tpl.csv", "w", newline="") as f_tpl:
    spl, bpl, tpl = csv.writer(f_spl), csv.writer(f_bpl), csv.writer(f_tpl)
    for e in mw.load_table():
        if e["model"] == "spl":
            spl.writerow([e["key"], _fmt(e, "slope"), _range(e), e["tracer"], _link(e)])
        elif e["model"] == "bpl":
            bpl.writerow([e["key"], _fmt(e, "r_break"), _fmt(e, "inner_slope"),
                          _fmt(e, "outer_slope"), _range(e), e["tracer"], _link(e)])
        else:
            tpl.writerow([e["key"], _fmt(e, "r_break"), _fmt(e, "r_break2"), _fmt(e, "inner_slope"),
                          _fmt(e, "middle_slope"), _fmt(e, "outer_slope"), _range(e), e["tracer"],
                          _link(e)])

# List of the references of the compilation, with their arXiv (or DOI) links
_refs = {}
for e in mw.load_table():
    ident = "arXiv:" + e["arxiv"] if e["arxiv"] else "doi:" + e["doi"]
    _refs.setdefault(e["reference"], (ident, mw.reference_url(e["key"])))
with open(docs_root / "mw" / "references.inc", "w") as f:
    for ref, (ident, url) in sorted(_refs.items(), key=lambda kv: kv[0].lower()):
        f.write("- {}, `{} <{}>`__\n".format(ref, ident, url))

# Figures made with the plot directive (mw/index.rst), in LaTeX style
from stracers.plotting import LATEX_STYLE  # noqa: E402

plot_rcparams = dict(LATEX_STYLE, **{"savefig.bbox": "tight"})  # keep the legend
plot_formats = [("png", 110)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_include_source = False

# HTML
html_theme = "pydata_sphinx_theme"
html_title = f"stracers v{version}"
html_baseurl = "https://jngaravitoc.github.io/stracers/"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "github_url": "https://github.com/jngaravitoc/stracers",
    "navbar_align": "left",
    "show_toc_level": 2,
}
html_sidebars = {"install": [], "getting_started": [], "changes": []}
