import numpy as np
import pytest

from stracers.halos import hernquist_df, hernquist_potential, sample_hernquist  # noqa: F401
from stracers.weights import G_GADGET as G  # noqa: F401


@pytest.fixture(scope="session")
def hernquist_halo():
    """Isotropic Hernquist halo: M = 1e10 Msun, a = 10 kpc, 2e5 particles."""
    rng = np.random.default_rng(42)
    n, M, a = 200_000, 1e10, 10.0
    pos, vel, pot = sample_hernquist(n, M, a, seed=rng)
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
