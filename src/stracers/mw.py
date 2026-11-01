"""
Milky Way stellar halo density profiles from the literature.

The profiles are the simple and broken power laws compiled in Table 6 of
Medina et al. 2024, MNRAS, 531, 4762 (doi:10.1093/mnras/stae1137), and the
doubly broken power law of Li et al. 2025 (arXiv:2512.01350), stored in
``stracers/data/mw_stellar_halo_profiles.csv``. They are spherical number
density profiles of halo stars, fitted only over the radial range listed for
each one (``r_min`` to ``r_max``).

Examples
--------
>>> from stracers import mw
>>> mw.available()[:3]
['wetterer1996', 'vivas2006', 'bell2008']
>>> nu = mw.mw_profile("medina2024")
"""

import csv
import functools
from importlib import resources

import numpy as np

from .profiles import broken_power_law, power_law, triple_power_law

__all__ = ["load_table", "available", "get_entry", "reference_url", "mw_profile", "profile_from_params"]

_DATA_FILE = "mw_stellar_halo_profiles.csv"
_TEXT_COLUMNS = ("key", "model", "tracer", "reference", "note", "arxiv", "doi")
_PARAMETERS = {
    "spl": ("slope",),
    "bpl": ("r_break", "inner_slope", "outer_slope"),
    "tpl": ("r_break", "r_break2", "inner_slope", "middle_slope", "outer_slope"),
}


def _parse(row):
    entry = {}
    for name, value in row.items():
        if name in _TEXT_COLUMNS:
            entry[name] = value
        else:
            entry[name] = float(value) if value != "" else None
    return entry


@functools.lru_cache(maxsize=1)
def _table():
    """The parsed data file, read once (callers get copies of the rows)."""
    text = resources.files("stracers").joinpath("data", _DATA_FILE).read_text()
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    return tuple(_parse(row) for row in csv.DictReader(lines))


def load_table():
    """
    Return the compiled profiles as a list of dicts, one per row of the data
    file. Missing values are ``None``.
    """
    return [dict(entry) for entry in _table()]


def available(model=None):
    """
    Keys of the available profiles, optionally only ``model`` = "spl"
    (simple power law), "bpl" (broken power law) or "tpl" (power law with
    two breaks).
    """
    return [e["key"] for e in load_table() if model is None or e["model"] == model]


def get_entry(key):
    """Return the table row (a dict) for ``key``."""
    for entry in _table():
        if entry["key"] == key:
            return dict(entry)
    raise ValueError("Unknown profile {!r}; see stracers.mw.available()".format(key))


def reference_url(key):
    """
    Link to the paper of a profile: its arXiv abstract page, or its DOI for
    the papers that are not on arXiv.
    """
    entry = get_entry(key)
    if entry["arxiv"]:
        return "https://arxiv.org/abs/" + entry["arxiv"]
    return "https://doi.org/" + entry["doi"]


def _draw(value, err_minus, err_plus, rng, positive=False):
    """
    Draw a parameter quoted as value (-err_minus, +err_plus), read as the
    median and the 16th/84th percentiles: half-normals of width err_minus
    below and err_plus above the value, each with probability 1/2. With
    ``positive``, non-positive draws are rejected and drawn again (for radii).
    """
    if err_minus is None or err_plus is None:
        return value
    while True:
        if rng.uniform() < 0.5:
            draw = value - abs(rng.normal(0.0, err_minus))
        else:
            draw = value + abs(rng.normal(0.0, err_plus))
        if draw > 0 or not positive:
            return draw


def mw_profile(key, r_core=1.0, smoothness=0.1, rng=None):
    """
    Tracer profile nu(r) for a literature Milky Way stellar halo profile.

    Parameters
    ----------
    key : str
        Profile key, e.g. "medina2024" (see :func:`available`).
    r_core : float
        Core radius [kpc] that flattens the profile at the centre, where none
        of these profiles were fitted and where power laws steeper than -3
        have a divergent mass. With ``r_core=0`` such profiles put most of the
        stellar mass on the few innermost particles. For slopes steeper than
        -3 the mass inside a few ``r_core`` depends strongly on its value, so
        choose it deliberately.
    smoothness : float
        Width (in dex) of the transition at the breaks of broken power laws.
        0 gives the sharp break of the published fits, which puts a narrow
        feature into f(E) that the energy binning cannot resolve (visible as
        a dip or excess near the break, large for steep outer slopes). The
        default, 0.1, removes it. Ignored for simple power laws.
    rng : numpy.random.Generator, optional
        If given, each parameter is drawn independently from its quoted
        uncertainties, read as median and 16th/84th percentiles (two
        half-normals joined at the median). Parameters without uncertainties
        are kept fixed, and non-positive break radii are drawn again (this
        matters for wide uncertainties, e.g. 26.5 +/- 8.9 kpc), as are pairs
        of break radii out of order.

    Returns
    -------
    nu : callable
        nu(r), with attribute ``params`` holding the parameters used.
    """
    entry = get_entry(key)
    while True:
        params = {}
        for name in _PARAMETERS[entry["model"]]:
            value = entry[name]
            if rng is not None:
                value = _draw(value, entry[name + "_err_minus"], entry[name + "_err_plus"], rng,
                              positive=name.startswith("r_break"))
            params[name] = value
        if params.get("r_break", 0.0) < params.get("r_break2", np.inf):
            break

    return profile_from_params(key, dict(params, r_core=r_core, smoothness=smoothness))


def profile_from_params(key, params):
    """
    Rebuild the profile of fit ``key`` from its parameters, e.g. the
    ``params`` attribute of a profile from :func:`mw_profile` (including
    parameters drawn with ``rng``), as stored with a saved tracer population.

    Parameters
    ----------
    key : str
        Profile key (see :func:`available`).
    params : dict
        ``slope`` (simple power law), ``r_break``, ``inner_slope`` and
        ``outer_slope`` (broken power law), or ``r_break``, ``r_break2``,
        ``inner_slope``, ``middle_slope`` and ``outer_slope`` (power law with
        two breaks), plus ``r_core`` and ``smoothness``.

    Returns
    -------
    nu : callable
        nu(r), with attribute ``params``.
    """
    model = get_entry(key)["model"]
    params = {k: float(v) for k, v in params.items()}
    r_core, smoothness = params.get("r_core", 0.0), params.get("smoothness", 0.0)
    if model == "spl":

        def nu(r):
            return power_law(r, params["slope"], r_core=r_core)

    elif model == "tpl":

        def nu(r):
            return triple_power_law(
                r, params["r_break"], params["r_break2"], params["inner_slope"],
                params["middle_slope"], params["outer_slope"], r_core=r_core, smoothness=smoothness,
            )

    else:

        def nu(r):
            return broken_power_law(
                r, params["r_break"], params["inner_slope"], params["outer_slope"],
                r_core=r_core, smoothness=smoothness,
            )

    nu.params = dict(params, r_core=r_core, smoothness=smoothness)
    nu.__name__ = "mw_profile({!r})".format(key)
    return nu
