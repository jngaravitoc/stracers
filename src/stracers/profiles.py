"""
Tracer (stellar halo) density profiles.

The Plummer, Hernquist, NFW and Einasto profiles are normalised to unit
total mass, i.e. they return nu(r) = rho(r) / M_tracer. The power-law
profiles are relative densities (their mass may diverge); the normalisation
of a tracer profile does not matter for :func:`stracers.compute_weights`.
Lengths are in kpc.

Any callable ``nu(r)`` can be used as a tracer profile in
:func:`stracers.compute_weights`; the functions here are provided for
convenience and can be bound to their parameters with :func:`get_profile`.
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

__all__ = [
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
]


def plummer(r, a):
    """Plummer profile with scale length ``a``."""
    r = np.asarray(r, dtype=float)
    return 3.0 / (4.0 * np.pi * a**3) * (1.0 + (r / a) ** 2) ** (-2.5)


def hernquist(r, a):
    """Hernquist (1990) profile with scale length ``a``."""
    r = np.asarray(r, dtype=float)
    return a / (2.0 * np.pi * r * (r + a) ** 3)


def nfw(r, rs, rvir):
    """
    NFW profile with scale radius ``rs``, normalised to unit mass inside
    ``rvir`` (the NFW mass diverges, so a truncation radius is required).
    """
    r = np.asarray(r, dtype=float)
    c = rvir / rs
    f_c = np.log(1.0 + c) - c / (1.0 + c)
    x = r / rs
    return 1.0 / (4.0 * np.pi * rs**3 * f_c * x * (1.0 + x) ** 2)


def einasto(r, n, r_eff):
    """
    Einasto profile with index ``n`` and effective radius ``r_eff``.

    Uses the d_n approximation of Merritt et al. (2006),
    https://arxiv.org/abs/astro-ph/0509417, valid for n > 0.5.
    """
    if n <= 0.5:
        raise ValueError("Einasto profile is only implemented for n > 0.5")
    r = np.asarray(r, dtype=float)
    d_n = 3.0 * n - 1.0 / 3.0 + 0.0079 / n
    # total mass for rho_e = 1
    m_tot = 4.0 * np.pi * r_eff**3 * n * np.exp(d_n) * d_n ** (-3.0 * n) * gamma(3.0 * n)
    return np.exp(-d_n * ((r / r_eff) ** (1.0 / n) - 1.0)) / m_tot


def power_law(r, slope, r_core=0.0):
    """
    Power law nu(r) = (r_s / 1 kpc)^slope with r_s = sqrt(r^2 + r_core^2).

    ``r_core`` > 0 flattens the profile inside ``r_core`` so that its mass
    does not diverge at the centre.
    """
    r = np.asarray(r, dtype=float)
    return np.sqrt(r**2 + r_core**2) ** slope


def broken_power_law(r, r_break, inner_slope, outer_slope, r_core=0.0, smoothness=0.0):
    """
    Broken power law with slopes ``inner_slope`` inside ``r_break`` and
    ``outer_slope`` outside, normalised to nu(r_break) = 1.

    With x = r_s / r_break and r_s = sqrt(r^2 + r_core^2):

    - ``smoothness`` = 0: nu = x^inner_slope for x < 1 and x^outer_slope for
      x >= 1 (continuous, with a kink at the break);
    - ``smoothness`` = s > 0: the :func:`double_power_law` with
      alpha = 1/s, which turns over smoothly over ~s dex around the break and
      avoids a spike in d2nu/dpsi2.

    ``r_core`` > 0 flattens the profile inside ``r_core``.
    """
    if smoothness > 0:
        return double_power_law(
            r, r_break, inner_slope, outer_slope, alpha=1.0 / smoothness, r_core=r_core
        )
    r = np.asarray(r, dtype=float)
    x = np.sqrt(r**2 + r_core**2) / r_break
    return np.where(x < 1.0, x**inner_slope, x**outer_slope)


def double_power_law(r, r_break, inner_slope, outer_slope, alpha=1.0, r_core=0.0):
    """
    Smooth double power law (Zhao 1996 "alpha-beta-gamma" form),

        nu = x^inner_slope [(1 + x^alpha) / 2]^((outer_slope - inner_slope) / alpha),

    with x = r_s / r_break and r_s = sqrt(r^2 + r_core^2), normalised to
    nu(r_break) = 1. The logarithmic slope goes from ``inner_slope`` at
    r << r_break to ``outer_slope`` at r >> r_break; ``alpha`` sets how sharp
    the transition is (larger is sharper). Hernquist (-1, -4, alpha=1), NFW
    (-1, -3, alpha=1) and Plummer (0, -5, alpha=2) are special cases.
    """
    r = np.asarray(r, dtype=float)
    x = np.sqrt(r**2 + r_core**2) / r_break
    return x**inner_slope * (0.5 * (1.0 + x**alpha)) ** ((outer_slope - inner_slope) / alpha)


_PROFILES = {
    "plummer": plummer,
    "hernquist": hernquist,
    "nfw": nfw,
    "einasto": einasto,
    "power_law": power_law,
    "broken_power_law": broken_power_law,
    "double_power_law": double_power_law,
    # short names
    "spl": power_law,
    "bpl": broken_power_law,
    "dpl": double_power_law,
}


def get_profile(name, **params):
    """
    Return a callable ``nu(r)`` for a named profile.

    Examples
    --------
    >>> nu = get_profile("hernquist", a=2.0)
    """
    try:
        func = _PROFILES[name.lower()]
    except KeyError:
        raise ValueError(
            "Unknown profile {!r}; choose from {}".format(name, sorted(_PROFILES))
        ) from None

    def nu(r):
        return func(r, **params)

    nu.__name__ = "{}({})".format(
        name.lower(), ", ".join("{}={}".format(k, v) for k, v in params.items())
    )
    return nu


def truncate(nu, r_t, sharpness=2.0):
    """
    Taper a profile smoothly beyond ``r_t``:
    nu_t(r) = nu(r) exp(-(r / r_t)^sharpness).

    The result is a new profile that can be passed to
    :func:`stracers.compute_weights`. A larger ``sharpness`` gives a sharper
    edge, but also makes a negative distribution function near the edge more
    likely. Tapering also gives a finite total mass to profiles that lack one
    (e.g. NFW).
    """

    def nu_t(r):
        r = np.asarray(r, dtype=float)
        return nu(r) * np.exp(-((r / r_t) ** sharpness))

    nu_t.__name__ = "truncate({}, r_t={}, sharpness={})".format(
        getattr(nu, "__name__", "nu"), r_t, sharpness
    )
    return nu_t


def mass_fraction_within(nu, rmax, rtol=1e-3):
    """
    Fraction of the total mass of the profile ``nu(r)`` that lies inside
    ``rmax``. ``nu`` does not need to be normalised.

    Raises ``ValueError`` if the total mass does not converge, or converges
    so slowly that it depends on radii > 10^6 rmax. For power laws this
    happens when the outer slope is -3 or shallower (e.g. NFW, or many of
    the simple power-law fits to the Milky Way halo).
    """

    def shell(lnr):
        r = np.exp(lnr)
        return 4.0 * np.pi * r**3 * nu(r)

    inner = quad(lambda r: 4.0 * np.pi * r**2 * nu(r), 0.0, rmax, limit=500)[0]
    lnr0 = np.log(rmax)
    tail_1 = quad(shell, lnr0, lnr0 + np.log(1e6), limit=500)[0]
    tail_2 = tail_1 + quad(shell, lnr0 + np.log(1e6), lnr0 + np.log(1e8), limit=500)[0]
    total = inner + tail_2
    if not np.isfinite(total) or total <= 0 or (tail_2 - tail_1) > rtol * total:
        raise ValueError(
            "The total mass of the tracer profile does not converge (outer "
            "slope -3 or shallower?). Give the stellar mass inside rmax "
            "(within='rmax'), taper the profile with stracers.truncate(nu, r_t), "
            "or truncate the stellar halo with r_trunc."
        )
    return inner / total
