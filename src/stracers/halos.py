"""
Toy equilibrium halos for tests and examples.

:func:`sample_hernquist` draws an isotropic Hernquist (1990) halo from its
analytic distribution function, with the analytic potential. It is what the
test suite and the documentation use in place of a simulation snapshot.
Units: kpc, km/s, Msun.
"""

import numpy as np

from .weights import G_GADGET

__all__ = ["hernquist_df", "hernquist_potential", "sample_hernquist"]


def hernquist_df(E, M, a, G=G_GADGET):
    """Isotropic Hernquist (1990) distribution function f(E), for E < 0."""
    vg = np.sqrt(G * M / a)
    q = np.sqrt(np.clip(-np.asarray(E, dtype=float) * a / (G * M), 0, 1 - 1e-12))
    pre = M / (8.0 * np.sqrt(2.0) * np.pi**3 * a**3 * vg**3)
    return pre / (1 - q**2) ** 2.5 * (
        3 * np.arcsin(q) + q * np.sqrt(1 - q**2) * (1 - 2 * q**2) * (8 * q**4 - 8 * q**2 - 3)
    )


def hernquist_potential(r, M, a, G=G_GADGET):
    """Hernquist potential -G M / (r + a) [(km/s)^2]."""
    return -G * M / (np.asarray(r, dtype=float) + a)


def sample_hernquist(n, M, a, seed=None, G=G_GADGET):
    """
    Sample an isotropic Hernquist halo in equilibrium.

    Parameters
    ----------
    n : int
        Number of particles.
    M : float
        Halo mass [Msun].
    a : float
        Scale radius [kpc].
    seed : int or numpy.random.Generator, optional
        Random seed or generator, for reproducible halos.
    G : float
        Gravitational constant [kpc (km/s)^2 / Msun].

    Returns
    -------
    pos, vel : ndarray, shape (n, 3)
        Positions [kpc] and velocities [km/s].
    pot : ndarray, shape (n,)
        Potential of each particle [(km/s)^2].

    Examples
    --------
    >>> from stracers.halos import sample_hernquist
    >>> pos, vel, pot = sample_hernquist(1000, 1e10, 10.0, seed=0)
    >>> pos.shape, vel.shape, pot.shape
    ((1000, 3), (1000, 3), (1000,))
    """
    rng = np.random.default_rng(seed)
    u = rng.uniform(0, 1, n)
    r = a * np.sqrt(u) / (1 - np.sqrt(u))

    phi = hernquist_potential(r, M, a, G)
    v_esc = np.sqrt(-2 * phi)
    f_max = hernquist_df(phi, M, a, G)
    v = np.empty(n)
    todo = np.arange(n)
    while len(todo):
        # proposal p(v) ~ v^2 on [0, v_esc]; accept with f(E) / f(phi)
        v_try = v_esc[todo] * rng.uniform(0, 1, len(todo)) ** (1 / 3)
        E = phi[todo] + 0.5 * v_try**2
        ok = rng.uniform(0, 1, len(todo)) < hernquist_df(E, M, a, G) / f_max[todo]
        v[todo[ok]] = v_try[ok]
        todo = todo[~ok]

    def directions(k):
        z = rng.uniform(-1, 1, k)
        t = rng.uniform(0, 2 * np.pi, k)
        s = np.sqrt(1 - z**2)
        return np.column_stack([s * np.cos(t), s * np.sin(t), z])

    pos = r[:, None] * directions(n)
    vel = v[:, None] * directions(n)
    return pos, vel, phi
