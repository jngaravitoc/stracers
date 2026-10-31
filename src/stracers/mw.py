"""
Milky Way stellar halo density profiles from the literature.

The profiles are the simple and broken power laws compiled in Table 6 of
Medina et al. 2024, MNRAS, 531, 4762 (doi:10.1093/mnras/stae1137), stored in
``stracers/data/mw_stellar_halo_profiles.csv``. They are spherical number
density profiles of halo stars, fitted only over the radial range listed for
each one (``r_min`` to ``r_max``).

Example
-------
>>> from stracers import mw
>>> mw.available()[:3]
['wetterer1996', 'vivas2006', 'bell2008']
>>> nu = mw.mw_profile("medina2024")
"""

import csv
from importlib import resources

import numpy as np

from .profiles import broken_power_law, power_law

__all__ = ["load_table", "available", "get_entry", "mw_profile"]

_DATA_FILE = "mw_stellar_halo_profiles.csv"
_TEXT_COLUMNS = ("key", "model", "tracer", "reference", "note")
_PARAMETERS = {"spl": ("slope",), "bpl": ("r_break", "inner_slope", "outer_slope")}


def _parse(row):
    entry = {}
    for name, value in row.items():
        if name in _TEXT_COLUMNS:
            entry[name] = value
        else:
            entry[name] = float(value) if value != "" else None
    return entry


def load_table():
    """
    Return the compiled profiles as a list of dicts, one per row of the data
    file. Missing values are ``None``.
    """
    text = resources.files("stracers").joinpath("data", _DATA_FILE).read_text()
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    return [_parse(row) for row in csv.DictReader(lines)]


def available(model=None):
    """
    Keys of the available profiles, optionally only ``model`` = "spl"
    (simple power law) or "bpl" (broken power law).
    """
    return [e["key"] for e in load_table() if model is None or e["model"] == model]


def get_entry(key):
    """Return the table row (a dict) for ``key``."""
    for entry in load_table():
        if entry["key"] == key:
            return entry
    raise ValueError("Unknown profile {!r}; see stracers.mw.available()".format(key))


def _draw(value, err_minus, err_plus, rng):
    """
    Draw a parameter quoted as value (-err_minus, +err_plus), read as the
    median and the 16th/84th percentiles: half-normals of width err_minus
    below and err_plus above the value, each with probability 1/2.
    """
    if err_minus is None or err_plus is None:
        return value
    if rng.uniform() < 0.5:
        return value - abs(rng.normal(0.0, err_minus))
    return value + abs(rng.normal(0.0, err_plus))


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
        Width (in dex) of the transition at the break of broken power laws.
        0 gives the sharp break of the published fits, which puts a narrow
        feature into f(E) that the energy binning cannot resolve (visible as
        a dip or excess near the break, large for steep outer slopes). The
        default, 0.1, removes it. Ignored for simple power laws.
    rng : numpy.random.Generator, optional
        If given, each parameter is drawn independently from its quoted
        uncertainties, read as median and 16th/84th percentiles (two
        half-normals joined at the median). Parameters without uncertainties
        are kept fixed.

    Returns
    -------
    nu : callable
        nu(r), with attribute ``params`` holding the parameters used.
    """
    entry = get_entry(key)
    params = {}
    for name in _PARAMETERS[entry["model"]]:
        value = entry[name]
        if rng is not None:
            value = _draw(value, entry[name + "_err_minus"], entry[name + "_err_plus"], rng)
        params[name] = value

    if entry["model"] == "spl":

        def nu(r):
            return power_law(r, params["slope"], r_core=r_core)

    else:

        def nu(r):
            return broken_power_law(
                r, params["r_break"], params["inner_slope"], params["outer_slope"],
                r_core=r_core, smoothness=smoothness,
            )

    nu.params = dict(params, r_core=r_core, smoothness=smoothness)
    nu.__name__ = "mw_profile({!r})".format(key)
    return nu
