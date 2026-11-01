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

import json
import re
import warnings
from dataclasses import dataclass, replace

import numpy as np

from . import kinematics, mw
from .io import match_weights
from .weights import compute_weights

__all__ = ["available_tracers", "tracer_fits", "choose_fit", "paint_tracer", "TracerPopulation",
           "load_population"]

_FORMAT = "stracers.TracerPopulation"

# canonical name -> (aliases, tracers in the ``tracer`` column of the table)
_TRACERS = {
    "RR Lyrae": (("rrlyr", "rrl", "rrlyrae", "rrlyra"), ("RRLs", "RRab")),
    "RRab": (("rrab",), ("RRab",)),
    "BHB": (("bhb", "bhbstar", "bluehorizontalbranch", "bluehorizontalbranchstar"),
            ("BHB stars", "BHB and BS stars")),
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


def _check_range(r_range, name="r_range"):
    """(r1, r2) as floats, with 0 < r1 < r2."""
    r1, r2 = map(float, r_range)
    if not 0 < r1 < r2:
        raise ValueError("{} must be (r1, r2) with 0 < r1 < r2, got {}".format(name, tuple(r_range)))
    return r1, r2


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
    r1, r2 = _check_range(r_range)
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
    result : WeightsResult or None
        The weights of the whole stellar halo, for diagnostics (``None`` for a
        population loaded from a file).
    pos, vel : ndarray
        Positions (relative to ``center``) and velocities of the particles.
    n_stars_all : ndarray
        Stars on every particle of the whole population, before the radial
        selection. They are what is carried to other snapshots.
    scale : float
        Normalisation: number of stars per unit fraction of the population.
    ids : ndarray or None
        Particle IDs, needed to follow the population to other snapshots.
    center : ndarray or None
        Halo centre used for the radial selection.
    snapshot : str or None
        Label of the snapshot.
    matched_fraction : float
        Fraction of the stars of the population found in this snapshot (1 for
        the snapshot where it was painted).
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
    result: object
    pos: np.ndarray
    vel: np.ndarray
    n_stars_all: np.ndarray = None
    scale: float = 1.0
    ids: np.ndarray = None
    center: np.ndarray = None
    snapshot: str = None
    matched_fraction: float = 1.0

    def _require_particles(self):
        if self.pos is None:
            raise ValueError("This population has no particles yet: place it on a snapshot "
                             "with at_snapshot(ids, pos, vel)")

    @property
    def total(self):
        """Expected number of stars of the tracer in ``r_range``."""
        self._require_particles()
        return float(np.sum(self.n_stars))

    def at_snapshot(self, ids, pos, vel, center=None, r_range=None, snapshot=None):
        """
        The same population in another snapshot of the simulation.

        The stars of each particle are carried by its ID (particles missing
        from the snapshot lose theirs), and the radial selection is applied
        again where the particles are now, about ``center``. The stars are
        conserved, not renormalised: the number in the range changes as they
        move.

        Parameters
        ----------
        ids : (N,) array
            Particle IDs of the snapshot.
        pos, vel : (N, 3) arrays
            Positions [kpc] and velocities of the particles.
        center : (3,) array, optional
            Halo centre in this snapshot.
        r_range : (float, float), optional
            Radial range of the selection [kpc]; defaults to that of the
            population.
        snapshot : str, optional
            Label of the snapshot.

        Returns
        -------
        TracerPopulation
        """
        if self.ids is None:
            raise ValueError("Pass ids to paint_tracer to follow the population to other snapshots")
        r1, r2 = _check_range(self.r_range if r_range is None else r_range)
        pos = np.asarray(pos, dtype=float)
        ids = _check_ids(ids, len(pos))
        if center is not None:
            center = np.asarray(center, dtype=float)
            pos = pos - center
        n_all = match_weights(self.ids, self.n_stars_all, ids)
        total_all = np.sum(self.n_stars_all)
        matched = float(np.sum(n_all) / total_all) if total_all > 0 else 1.0
        if matched < 0.99:
            warnings.warn("Only {:.1%} of the stars of the population are on particles of this "
                          "snapshot".format(matched), RuntimeWarning)
        r = np.linalg.norm(pos, axis=1)
        in_range = (r >= r1) & (r <= r2)
        return replace(
            self, r_range=(r1, r2), in_range=in_range, n_stars=np.where(in_range, n_all, 0.0),
            pos=pos, vel=np.asarray(vel, dtype=float), n_stars_all=n_all, ids=ids,
            center=center, snapshot=snapshot, matched_fraction=matched,
        )

    def catalogue(self, kind="stars", rng=None, columns=None, vcenter=None):
        """
        Mock catalogue in memory, as a dict of columns
        (:func:`stracers.catalogue.catalogue_columns`). Use
        :func:`stracers.write_catalogue` for large catalogues.
        """
        from .catalogue import catalogue_columns

        return catalogue_columns(self, kind=kind, rng=rng, columns=columns, vcenter=vcenter)

    def save(self, path):
        """
        Save the population to an HDF5 file (needs h5py): the particle IDs, the
        stars of each particle and the metadata. Load it with
        :func:`load_population` and place it on any snapshot with
        :meth:`at_snapshot`, without computing the weights again.
        """
        import h5py

        if self.ids is None:
            raise ValueError("Pass ids to paint_tracer to save the population")
        with h5py.File(path, "w") as f:
            f.create_dataset("ids", data=self.ids)
            f.create_dataset("n_stars_all", data=self.n_stars_all)
            f.attrs.update(
                format=_FORMAT, tracer=self.tracer, key=self.key, reference=self.reference,
                url=self.url, r_range=self.r_range, scale=self.scale,
                fit_range=[np.nan if v is None else v for v in self.fit_range],
                profile_params=json.dumps(getattr(self.profile, "params", {})),
                snapshot="" if self.snapshot is None else str(self.snapshot),
                stracers_version=_version(),
            )

    def velocity_profile(self, bins=10, **kwargs):
        """
        Velocity profile of the tracer in ``r_range``
        (:func:`stracers.kinematics.velocity_profile`, weighted by
        ``n_stars``). An integer ``bins`` spans ``r_range``.
        """
        self._require_particles()
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
        self._require_particles()
        rng = np.random.default_rng(rng)
        counts = rng.poisson(self.n_stars)
        return np.repeat(np.arange(len(self.n_stars)), counts)


def _check_ids(ids, n, name="ids"):
    """IDs as an array, one per particle and without duplicates."""
    ids = np.asarray(ids)
    if ids.shape != (n,):
        raise ValueError("{} must have one value per particle".format(name))
    if np.unique(ids).size != ids.size:
        raise ValueError("{} contains duplicated IDs: the stars of a particle would be "
                         "copied onto each of them".format(name))
    return ids


def _version():
    from . import __version__

    return __version__


def load_population(path):
    """
    Load a population saved with :meth:`TracerPopulation.save`.

    The population has no particles until it is placed on a snapshot with
    :meth:`TracerPopulation.at_snapshot`; the weights are not computed again.
    """
    import h5py

    with h5py.File(path, "r") as f:
        attrs = dict(f.attrs)
        if attrs.get("format") != _FORMAT:
            raise ValueError("{} is not a stracers tracer population file".format(path))
        ids, n_all = f["ids"][()], f["n_stars_all"][()]
    fit_range = tuple(None if np.isnan(v) else float(v) for v in attrs["fit_range"])
    return TracerPopulation(
        tracer=str(attrs["tracer"]), key=str(attrs["key"]), reference=str(attrs["reference"]),
        url=str(attrs["url"]), fit_range=fit_range, r_range=tuple(float(v) for v in attrs["r_range"]),
        profile=mw.profile_from_params(str(attrs["key"]), json.loads(attrs["profile_params"])),
        in_range=None, n_stars=None, result=None, pos=None, vel=None,
        n_stars_all=n_all, scale=float(attrs["scale"]), ids=ids,
        snapshot=str(attrs["snapshot"]) or None,
    )


def _tracer_fraction(frac, sel, name, r_range):
    """Fraction of the tracer population in a selection; must be positive to normalise."""
    total = np.sum(frac[sel])
    if not (np.isfinite(total) and total > 0):
        raise ValueError(
            "No tracer stars in {} = {:g}-{:g} kpc ({} particles, all with zero weight or none "
            "selected), so it cannot be used to normalise the population".format(
                name, r_range[0], r_range[1], int(np.sum(sel)))
        )
    return total


def paint_tracer(pos, vel, pot, mass, tracer, r_range, *, key=None, n_stars=None, n_ref=None,
                 ref_range=None, profile_kwargs=None, center=None, ids=None, snapshot=None,
                 **weights_kwargs):
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
    ids : (N,) array, optional
        Particle IDs. Needed to follow the population to other snapshots
        (:meth:`TracerPopulation.at_snapshot`) and to save it.
    snapshot : str, optional
        Label of the snapshot, kept as metadata.
    **weights_kwargs :
        Passed to :func:`stracers.compute_weights`. ``rmax`` defaults to
        ``max(300, r_range[1])`` kpc, and ``beyond_rmax`` to ``"energy"``:
        particles beyond ``rmax`` also carry stars, according to their energy,
        so that the population is stationary and stars on orbits reaching
        beyond ``rmax`` are not lost when it is followed in time.

    Returns
    -------
    TracerPopulation
    """
    if n_stars is not None and n_ref is not None:
        raise ValueError("Give either n_stars or n_ref/ref_range, not both")
    if n_ref is not None and ref_range is None:
        raise ValueError("n_ref needs ref_range")
    name = _resolve(tracer)
    # validated here, before a key given by the user bypasses choose_fit
    r1, r2 = _check_range(r_range)
    if ref_range is not None:
        ref_range = _check_range(ref_range, "ref_range")
    if key is None:
        key = choose_fit(tracer, (r1, r2))
    elif key not in tracer_fits(tracer):
        raise ValueError("{!r} is not a {} fit; choose from {}".format(key, name, tracer_fits(tracer)))

    rmax = weights_kwargs.pop("rmax", max(300.0, r2))
    # weights as a function of energy for every particle, also beyond rmax,
    # so that the population is stationary and can be followed in time
    weights_kwargs.setdefault("beyond_rmax", "energy")
    if r2 > rmax or (ref_range is not None and ref_range[1] > rmax):
        raise ValueError("r_range and ref_range must be inside rmax = {:g}".format(rmax))

    profile = mw.mw_profile(key, **(profile_kwargs or {}))
    result = compute_weights(pos, vel, pot, mass, profile, stellar_mass=1.0, rmax=rmax,
                             center=center, **weights_kwargs)

    pos = np.asarray(pos, dtype=float)
    if center is not None:
        center = np.asarray(center, dtype=float)
        pos = pos - center
    if ids is not None:
        ids = _check_ids(ids, len(pos))
    r = np.linalg.norm(pos, axis=1)
    frac = result.weights * np.broadcast_to(np.asarray(mass, dtype=float), r.shape)
    in_range = (r >= r1) & (r <= r2)

    scale = 1.0
    if n_stars is not None:
        scale = n_stars / _tracer_fraction(frac, in_range, "r_range", (r1, r2))
    elif n_ref is not None:
        in_ref = (r >= ref_range[0]) & (r <= ref_range[1])
        scale = n_ref / _tracer_fraction(frac, in_ref, "ref_range", ref_range)

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
        n_stars_all=frac * scale,
        scale=float(scale),
        ids=ids,
        center=center,
        snapshot=snapshot,
    )
