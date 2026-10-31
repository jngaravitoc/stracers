import numpy as np
import pytest

G = 4.30071e-6  # kpc (km/s)^2 / Msun


def hernquist_df(E, M, a):
    """Isotropic Hernquist (1990) distribution function f(E), E < 0."""
    vg = np.sqrt(G * M / a)
    q = np.sqrt(np.clip(-E * a / (G * M), 0, 1 - 1e-12))
    pre = M / (8.0 * np.sqrt(2.0) * np.pi**3 * a**3 * vg**3)
    return pre / (1 - q**2) ** 2.5 * (
        3 * np.arcsin(q) + q * np.sqrt(1 - q**2) * (1 - 2 * q**2) * (8 * q**4 - 8 * q**2 - 3)
    )


def hernquist_potential(r, M, a):
    return -G * M / (r + a)


def sample_hernquist(n, M, a, rng):
    """Positions and velocities of an isotropic Hernquist halo in equilibrium."""
    u = rng.uniform(0, 1, n)
    r = a * np.sqrt(u) / (1 - np.sqrt(u))

    phi = hernquist_potential(r, M, a)
    v_esc = np.sqrt(-2 * phi)
    f_max = hernquist_df(phi, M, a)
    v = np.empty(n)
    todo = np.arange(n)
    while len(todo):
        # proposal p(v) ~ v^2 on [0, v_esc]; accept with f(E) / f(phi)
        v_try = v_esc[todo] * rng.uniform(0, 1, len(todo)) ** (1 / 3)
        E = phi[todo] + 0.5 * v_try**2
        ok = rng.uniform(0, 1, len(todo)) < hernquist_df(E, M, a) / f_max[todo]
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


@pytest.fixture(scope="session")
def hernquist_halo():
    """Isotropic Hernquist halo: M = 1e10 Msun, a = 10 kpc, 2e5 particles."""
    rng = np.random.default_rng(42)
    n, M, a = 200_000, 1e10, 10.0
    pos, vel, pot = sample_hernquist(n, M, a, rng)
    return {
        "pos": pos,
        "vel": vel,
        "pot": pot,
        "mass": M / n,
        "ids": rng.permutation(n) + 1,
        "M": M,
        "a": a,
    }


def weighted_density(r, mass, edges):
    """Spherical density profile from (weighted) particle masses."""
    m_shell, _ = np.histogram(r, bins=edges, weights=mass)
    vol = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    centers = np.sqrt(edges[1:] * edges[:-1])
    return centers, m_shell / vol
