"""
stracers: stellar-tracer weights for dark matter particles in N-body halos,
following Laporte et al. (2013).
"""

from .profiles import (
    broken_power_law,
    double_power_law,
    einasto,
    get_profile,
    hernquist,
    mass_fraction_within,
    nfw,
    plummer,
    power_law,
    truncate,
)
from .weights import WeightsResult, compute_weights
from .io import load_weights, match_weights, save_weights
from . import halos, mw
from . import nstars
from .nstars import number_of_stars, number_of_stars_table

__version__ = "0.2.1"

__all__ = [
    "compute_weights",
    "WeightsResult",
    "plummer",
    "hernquist",
    "nfw",
    "einasto",
    "power_law",
    "broken_power_law",
    "double_power_law",
    "get_profile",
    "truncate",
    "mass_fraction_within",
    "save_weights",
    "load_weights",
    "match_weights",
    "mw",
    "halos",
    "nstars",
    "number_of_stars",
    "number_of_stars_table",
]
