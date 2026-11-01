"""
Mock populations of stellar tracers.

:func:`paint_tracer` paints a tracer population, e.g. RR Lyrae, BHB stars or
K giants, onto the dark matter particles of a halo: it chooses a literature
profile for that tracer from :mod:`stracers.mw`, computes the weights for the
whole stellar halo, and keeps the particles in a radial range, with the
expected number of stars of that type on each one.

The radial range is a selection applied *after* the weighting, like an
observational cut: truncating the profile instead would change the
equilibrium distribution function and the kinematics of the stars. Radii are
in kpc, like the profiles of :mod:`stracers.mw`.

Example
-------
>>> pop = paint_tracer(pos, vel, pot, mass, "RRLyr", r_range=(20, 100), n_stars=3000)  # doctest: +SKIP
>>> pop.key, pop.total                                                                 # doctest: +SKIP
('medina2024', 3000.0)
"""

import re
import warnings
from dataclasses import dataclass

import numpy as np

from . import kinematics, mw
from .weights import WeightsResult, compute_weights

__all__ = ["available_tracers", "tracer_fits", "choose_fit", "paint_tracer", "TracerPopulation"]

# canonical name -> (aliases, tracers in the ``tracer`` column of the table)
_TRACERS = {
    "RR Lyrae": (("rrlyr", "rrl", "rrlyrae", "rrlyra"), ("RRLs", "RRab")),
    "RRab": (("rrab",), ("RRab",)),
    "BHB": (("bhb", "bluehorizontalbranch"), ("BHB stars", "BHB and BS stars")),
    "K giants": (("kgiant",), ("K giants",)),
    "M giants": (("mgiant",), ("M giants",)),
    "Giants": (("giant",), ("Giant stars", "K giants", "M giants")),
    "MSTO": (("msto", "mainsequenceturnoff"), ("MSTO stars", "near MSTO stars")),
    "F stars": (("fstar", "f"), ("F stars",)),
    "A-type stars": (("atype", "atypestar", "astar"), ("A-type stars",)),
}


def _normalise(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _resolve(tracer):
    """Canonical name of a tracer, from any of its spellings."""
    key = _normalise(tracer)
    for name, (aliases, _) in _TRACERS.items():
        if key in aliases or (key.endswith("s") and key[:-1] in aliases):
            return name
    raise ValueError(
        "Unknown tracer {!r}; choose from {}".format(tracer, ", ".join(_TRACERS))
    )


def available_tracers():
    """
    The tracer names and the keys of their literature fits, as a dict
    ``{name: [keys]}``. Names are matched case- and punctuation-insensitively,
    so e.g. ``"K-giant"``, ``"k giants"`` and ``"KGiants"`` are the same.
    """
    return {name: tracer_fits(name) for name in _TRACERS}


def tracer_fits(tracer):
    """Keys of the fits in :mod:`stracers.mw` traced by ``tracer``."""
    populations = _TRACERS[_resolve(tracer)][1]
    return [e["key"] for e in mw.load_table() if e["tracer"] in populations]


def _year(entry):
    match = re.search(r"(\d{4})", entry["reference"])
    return int(match.group(1)) if match else 0


def _rank(entry):
    """Newest first, then broken over simple power laws, then the widest range."""
    r_lo = entry["r_min"] if entry["r_min"] else 1e-3
    return (_year(entry), entry["model"] == "bpl", np.log(entry["r_max"] / r_lo))


def choose_fit(tracer, r_range):
    """
    Choose the literature fit for ``tracer`` in the radial range ``r_range``.

    The fit is the most recent whose fitted range covers ``r_range`` (fits
    quoted only with an outer radius, "<X kpc", count as starting at 0);
    ties go to broken over simple power laws, then to the widest fitted range.
    If no fit covers the range, the one covering the largest part of it (in
    log r) is chosen, with a warning: part of the range then uses an
    extrapolated profile.

    Parameters
    ----------
    tracer : str
        Tracer name, e.g. ``"RRLyr"``, ``"BHB"`` or ``"K-giant"`` (see
        :func:`available_tracers`).
    r_range : (float, float)
        Radial range [kpc].

    Returns
    -------
    str
        Key of the fit in :mod:`stracers.mw`.
    """
    r1, r2 = map(float, r_range)
    if not 0 < r1 < r2:
        raise ValueError("r_range must be (r1, r2) with 0 < r1 < r2")
    entries = [mw.get_entry(k) for k in tracer_fits(tracer)]
    covering = [e for e in entries if (e["r_min"] or 0.0) <= r1 and e["r_max"] >= r2]
    if covering:
        return max(covering, key=_rank)["key"]

    def overlap(e):
        lo, hi = max(r1, e["r_min"] or 0.0), min(r2, e["r_max"])
        return max(0.0, np.log(hi / lo)) / np.log(r2 / r1) if hi > lo else 0.0

    best = max(entries, key=lambda e: (round(overlap(e), 6),) + _rank(e))
    warnings.warn(
        "No {} fit covers {:g}-{:g} kpc; using {} (fitted over {}-{:g} kpc, {:.0%} of the "
        "range): the rest of the range uses the extrapolated profile.".format(
            _resolve(tracer), r1, r2, best["key"],
            "{:g}".format(best["r_min"]) if best["r_min"] else "0", best["r_max"], overlap(best),
        ),
        RuntimeWarning,
    )
    return best["key"]


@dataclass
class TracerPopulation:
    """
    A mock tracer population, returned by :func:`paint_tracer`.

    Attributes
    ----------
    tracer : str
        Canonical tracer name.
    key : str
        Key of the literature fit used (:mod:`stracers.mw`).
    reference : str
        Reference of that fit.
    url : str
        Link to its paper (:func:`stracers.mw.reference_url`).
    fit_range : tuple
        Radial range of the fit [kpc] (``None`` for a missing inner radius).
    r_range : tuple
        Radial range of the population [kpc].
    profile : callable
        The tracer profile used for the weights.
    in_range : ndarray of bool
        Particles inside ``r_range``.
    n_stars : ndarray
        Expected number of stars of the tracer on each particle (zero outside
        ``r_range``); fractions of the tracer population inside ``rmax`` if no
        normalisation was given.
    result : WeightsResult
        The weights of the whole stellar halo, for diagnostics.
    pos, vel : ndarray
        Positions (relative to the centre) and velocities of the particles.
    """

    tracer: str
    key: str
    reference: str
    url: str
    fit_range: tuple
    r_range: tuple
    profile: object
    in_range: np.ndarray
    n_stars: np.ndarray
    result: WeightsResult
    pos: np.ndarray
    vel: np.ndarray

    @property
    def total(self):
        """Expected number of stars of the tracer in ``r_range``."""
        return float(np.sum(self.n_stars))

    def velocity_profile(self, bins=10, **kwargs):
        """
        Velocity profile of the tracer in ``r_range``
        (:func:`stracers.kinematics.velocity_profile`, weighted by
        ``n_stars``). An integer ``bins`` spans ``r_range``.
        """
        sel = self.in_range
        if np.ndim(bins) == 0:
            kwargs.setdefault("rmin", self.r_range[0])
            kwargs.setdefault("rmax", self.r_range[1])
        return kinematics.velocity_profile(self.pos[sel], self.vel[sel], self.n_stars[sel],
                                           bins=bins, **kwargs)

    def draw_stars(self, rng=None):
        """
        Draw a mock catalogue: a Poisson number of stars on each particle,
        with mean ``n_stars``. Returns the particle index of each star, so
        ``pos[idx]`` and ``vel[idx]`` are the stars' positions and velocities.
        """
        rng = np.random.default_rng(rng)
        counts = rng.poisson(self.n_stars)
        return np.repeat(np.arange(len(self.n_stars)), counts)


def paint_tracer(pos, vel, pot, mass, tracer, r_range, *, key=None, n_stars=None, n_ref=None,
                 ref_range=None, profile_kwargs=None, center=None, **weights_kwargs):
    """
    Paint a mock population of a stellar tracer onto a halo.

    Parameters
    ----------
    pos, vel, pot, mass :
        Positions [kpc], velocities, potential and masses of the dark matter
        particles, as in :func:`stracers.compute_weights`.
    tracer : str
        Tracer name, e.g. ``"RRLyr"``, ``"BHB"``, ``"K-giant"`` (see
        :func:`available_tracers`).
    r_range : (float, float)
        Radial range of the population [kpc].
    key : str, optional
        Literature fit to use; by default :func:`choose_fit` picks the most
        recent one covering ``r_range``.
    n_stars : float, optional
        Number of stars of the tracer in ``r_range``.
    n_ref : float, optional
        Normalise instead to ``n_ref`` stars in ``ref_range`` [kpc], e.g. an
        observed count.
    profile_kwargs : dict, optional
        Passed to :func:`stracers.mw.mw_profile` (``r_core``, ``smoothness``,
        ``rng``).
    center : (3,) array, optional
        Centre of the halo.
    **weights_kwargs :
        Passed to :func:`stracers.compute_weights`. ``rmax`` defaults to
        ``max(300, r_range[1])`` kpc.

    Returns
    -------
    TracerPopulation
    """
    if n_stars is not None and n_ref is not None:
        raise ValueError("Give either n_stars or n_ref/ref_range, not both")
    if n_ref is not None and ref_range is None:
        raise ValueError("n_ref needs ref_range")
    name = _resolve(tracer)
    r1, r2 = map(float, r_range)
    if key is None:
        key = choose_fit(tracer, (r1, r2))
    elif key not in tracer_fits(tracer):
        raise ValueError("{!r} is not a {} fit; choose from {}".format(key, name, tracer_fits(tracer)))

    rmax = weights_kwargs.pop("rmax", max(300.0, r2))
    if r2 > rmax or (ref_range is not None and ref_range[1] > rmax):
        raise ValueError("r_range and ref_range must be inside rmax = {:g}".format(rmax))

    profile = mw.mw_profile(key, **(profile_kwargs or {}))
    result = compute_weights(pos, vel, pot, mass, profile, stellar_mass=1.0, rmax=rmax,
                             center=center, **weights_kwargs)

    pos = np.asarray(pos, dtype=float)
    if center is not None:
        pos = pos - np.asarray(center, dtype=float)
    r = np.linalg.norm(pos, axis=1)
    frac = result.weights * np.broadcast_to(np.asarray(mass, dtype=float), r.shape)
    in_range = (r >= r1) & (r <= r2)

    scale = 1.0
    if n_stars is not None:
        scale = n_stars / np.sum(frac[in_range])
    elif n_ref is not None:
        in_ref = (r >= ref_range[0]) & (r <= ref_range[1])
        scale = n_ref / np.sum(frac[in_ref])

    entry = mw.get_entry(key)
    return TracerPopulation(
        tracer=name,
        key=key,
        reference=entry["reference"],
        url=mw.reference_url(key),
        fit_range=(entry["r_min"], entry["r_max"]),
        r_range=(r1, r2),
        profile=profile,
        in_range=in_range,
        n_stars=np.where(in_range, frac * scale, 0.0),
        result=result,
        pos=pos,
        vel=np.asarray(vel, dtype=float),
    )
