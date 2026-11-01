r"""
Number of stars of a stellar halo in a range of galactocentric radius.

For a number density profile n(r_q), with ellipsoidal radius
:math:`r_q = \sqrt{R^2 + z^2/q^2}` (axisymmetric, flattening q), the number of
stars with galactocentric radius between r_min and r_max is

.. math::
    N = \int_{r_{min}}^{r_{max}} r^2\, dr \int d\Omega\;
        n\!\left(r \sqrt{\sin^2\theta + \cos^2\theta / q^2}\right),

as used, for example, in Medina et al. (2024, Sect. 5.2.1 and Fig. 7) to
predict the number of RR Lyrae beyond a given radius. The profile can be any
callable or the key of a Milky Way profile in :mod:`stracers.mw`.
"""

import warnings

import numpy as np
from scipy.integrate import quad

from . import mw

__all__ = ["number_of_stars", "number_of_stars_table", "FULL_SKY_DEG2"]

FULL_SKY_DEG2 = 4.0 * np.pi * (180.0 / np.pi) ** 2  # 41252.96 deg^2

_MU, _MU_W = np.polynomial.legendre.leggauss(64)
_MU = 0.5 * (_MU + 1.0)  # nodes on [0, 1] for mu = cos(theta)
_MU_W = 0.5 * _MU_W


def _shell(nu, r, q):
    """Integral of nu over the sphere of radius r, 4 pi <nu>."""
    if q == 1.0:
        return 4.0 * np.pi * nu(r)
    s = np.sqrt(1.0 - _MU**2 + _MU**2 / q**2)
    return 4.0 * np.pi * np.sum(_MU_W * nu(r * s))


def _integral(nu, r_min, r_max, q, rtol=1e-2):
    """
    Integral of nu over r_min < r < r_max. For r_max = inf the integral is
    taken to 10^8 r_min and rejected if radii beyond 10^6 r_min contribute
    more than ``rtol`` of it.
    """

    def integrand(lnr):
        r = np.exp(lnr)
        return r**3 * _shell(nu, r, q)

    if r_min > 0:
        lo = np.log(r_min)
        if np.isfinite(r_max):
            return quad(integrand, lo, np.log(r_max), limit=500)[0]
        part_1 = quad(integrand, lo, lo + np.log(1e6), limit=500)[0]
        part_2 = quad(integrand, lo + np.log(1e6), lo + np.log(1e8), limit=500)[0]
        if not np.isfinite(part_1 + part_2) or part_2 > rtol * (part_1 + part_2):
            raise ValueError(
                "The number of stars out to infinity does not converge (outer "
                "slope -3 or shallower?); give a finite r_max"
            )
        return part_1 + part_2

    # r_min = 0: integrate the centre in linear r, then the rest in ln r
    r_split = min(1e-3, r_max)
    inner = quad(lambda r: r**2 * _shell(nu, r, q), 0.0, r_split, limit=500)[0]
    if r_max <= r_split:
        return inner
    return inner + _integral(nu, r_split, r_max, q, rtol)


def number_of_stars(
    nu,
    r_min,
    r_max=np.inf,
    *,
    q=1.0,
    n0=None,
    r0=None,
    n_ref=None,
    ref_range=None,
    area_deg2=None,
    r_core=1.0,
    smoothness=0.1,
    rng=None,
):
    """
    Number of stars with galactocentric radius between ``r_min`` and ``r_max``.

    Parameters
    ----------
    nu : callable or str
        Number density profile as a function of the (ellipsoidal) radius
        [kpc], e.g. from :mod:`stracers.profiles`, or the key of a Milky Way
        profile in :mod:`stracers.mw` (e.g. "medina2024"; see
        :func:`stracers.mw.available`). A callable's own normalisation is
        used unless ``n0`` or ``n_ref`` is given; for a key one of them is
        required.
    r_min : float or array_like
        Inner radius [kpc]. An array gives the cumulative counts for each
        value (e.g. N(> r) curves as in Medina et al. 2024, Fig. 7).
    r_max : float
        Outer radius [kpc], i.e. the outer edge of the halo. It may be
        ``np.inf`` if the counts converge. It matters a lot for shallow
        profiles.
    q : float
        Flattening of the halo; ``nu`` is evaluated at
        r_q = sqrt(R^2 + z^2 / q^2). q < 1 is oblate.
    n0, r0 : float, optional
        Normalise to a number density ``n0`` [stars / kpc^3] at radius
        ``r0`` [kpc] in the plane (z = 0), e.g. the local density of RR Lyrae
        at the solar radius.
    n_ref : float, optional
        Normalise instead so that there are ``n_ref`` stars (over the whole
        sky) in ``ref_range`` = (r_min, r_max), e.g. ``n_ref=3000,
        ref_range=(100, 300)``.
    area_deg2 : float, optional
        Survey area [deg^2]. The full-sky number is scaled by
        ``area_deg2 / 41253``, i.e. assuming the stars are spread uniformly
        over the sky. This ignores the direction dependence of a flattened
        halo and the survey footprint.
    r_core, smoothness, rng :
        Only used when ``nu`` is a key; passed to
        :func:`stracers.mw.mw_profile`. ``rng`` draws the profile parameters
        within their quoted uncertainties.

    Returns
    -------
    float or ndarray
        Number of stars (in the survey area if ``area_deg2`` is given).

    Notes
    -----
    The literature profiles are fits over a limited radial range (columns
    ``r_min``/``r_max`` of the table, see :func:`stracers.mw.get_entry`);
    counts outside that range are extrapolations.
    """
    if isinstance(nu, str):
        if n0 is None and n_ref is None:
            raise ValueError(
                "The Milky Way profiles have an arbitrary normalisation "
                "(nu = 1 at the break or at 1 kpc); give n0/r0 or n_ref/ref_range"
            )
        nu = mw.mw_profile(nu, r_core=r_core, smoothness=smoothness, rng=rng)
    if n0 is not None and n_ref is not None:
        raise ValueError("Give either n0/r0 or n_ref/ref_range, not both")
    if q <= 0:
        raise ValueError("q must be positive")

    scale = 1.0
    if n0 is not None:
        if r0 is None:
            raise ValueError("n0 needs r0")
        scale = n0 / float(nu(np.asarray(r0, dtype=float)))
    elif n_ref is not None:
        if ref_range is None:
            raise ValueError("n_ref needs ref_range")
        scale = n_ref / _integral(nu, ref_range[0], ref_range[1], q)

    if area_deg2 is not None:
        scale *= area_deg2 / FULL_SKY_DEG2

    r_min_arr = np.atleast_1d(np.asarray(r_min, dtype=float))
    if np.any(r_min_arr < 0) or np.any(r_min_arr > r_max):
        raise ValueError("Need 0 <= r_min <= r_max")
    counts = np.array([scale * _integral(nu, rm, r_max, q) for rm in r_min_arr])
    return counts[0] if np.ndim(r_min) == 0 else counts


def number_of_stars_table(r_min, r_max=np.inf, *, keys=None, **kwargs):
    """
    :func:`number_of_stars` for every Milky Way profile in the table.

    Parameters
    ----------
    r_min, r_max :
        As in :func:`number_of_stars`.
    keys : list of str, optional
        Profiles to use; defaults to all of :func:`stracers.mw.available`.
    **kwargs :
        Passed to :func:`number_of_stars`; a normalisation (``n0``/``r0`` or
        ``n_ref``/``ref_range``) is required. The same normalisation is
        applied to every profile.

    Returns
    -------
    dict
        ``{key: N}``. Profiles whose counts diverge for the requested range
        (e.g. outer slope -3 or shallower with ``r_max=np.inf``) get ``nan``,
        with a warning listing them.
    """
    if keys is None:
        keys = mw.available()
    out, diverging = {}, []
    for key in keys:
        try:
            out[key] = number_of_stars(key, r_min, r_max, **kwargs)
        except ValueError as err:
            if "does not converge" not in str(err):
                raise
            out[key] = np.full(np.shape(r_min), np.nan) if np.ndim(r_min) else np.nan
            diverging.append(key)
    if diverging:
        warnings.warn(
            "Counts do not converge for {}; set to nan. Use a finite r_max.".format(
                ", ".join(diverging)
            ),
            RuntimeWarning,
        )
    return out
