"""Cross-check of the weights with the method of npaint (Errani & Penarrubia 2020).

npaint (https://github.com/rerrani/nbopy, arXiv:1906.01642) paints stars on an
N-body model with probabilities f_*(E) / f_DM(E), where both distribution
functions are Eddington inversions of analytic densities in the analytic
potential. stracers instead estimates f_DM = N(E) / g(E) from the particles and
uses their own potentials. The two must agree on an equilibrium halo.

npaint is GPL-3 and is not copied here. This file re-implements its method
independently, for its default models: a Hernquist halo (G = M = a = 1) and the
tracer nu(r) = r^-gS (rsS^aS + r^aS)^((gS - bS) / aS) with aS = 2, bS = 5,
gS = 0.1, rsS = 0.2. Run once on the same 2e5-particle halo, npaint and this
reference gave the same weights to 5e-4 (median per particle).
"""

import numpy as np
import pytest
from scipy.integrate import quad

from stracers import compute_weights
from stracers.halos import hernquist_df, sample_hernquist
from stracers.kinematics import velocity_profile

A_S, B_S, G_S, R_S = 2.0, 5.0, 0.1, 0.2     # npaint's default tracer (aS, bS, gS, rsS)
RMAX = 100.0                                  # npaint's Rmax


def nu_tracer(r, gamma=G_S):
    """npaint's tracer density (not normalised)."""
    return r ** (-gamma) * (R_S**A_S + r**A_S) ** ((gamma - B_S) / A_S)


def d2nu_dpsi2(psi, gamma=G_S):
    """Exact d2nu/dpsi2 of the tracer in the Hernquist potential psi = 1 / (1 + r)."""
    r = 1.0 / psi - 1.0
    s = R_S**A_S + r**A_S
    # first and second derivatives of ln(nu) with respect to r
    L1 = -gamma / r + (gamma - B_S) * r ** (A_S - 1) / s
    L2 = gamma / r**2 + (gamma - B_S) * (
        (A_S - 1) * r ** (A_S - 2) * s - A_S * r ** (2 * A_S - 2)
    ) / s**2
    nu = nu_tracer(r, gamma)
    d1, d2 = nu * L1, nu * (L2 + L1**2)
    p1, p2 = -1.0 / (1 + r) ** 2, 2.0 / (1 + r) ** 3
    return (d2 * p1 - d1 * p2) / p1**3


def eddington(eps, d2):
    """f(eps) = 1/(sqrt(8) pi^2) int_0^eps d2(psi) dpsi / sqrt(eps - psi), with psi = eps (1 - u^2)."""
    integrand = lambda u: d2(eps * (1 - u * u)) * 2 * np.sqrt(eps)
    return quad(integrand, 0.0, 1.0, limit=400, epsrel=1e-7)[0] / (np.sqrt(8) * np.pi**2)


def reference_weights(E, n_grid=400):
    """npaint's probabilities f_*(E) / f_DM(E), not normalised; zero for E >= 0."""
    eps = np.geomspace(1e-6, 1 - 1e-6, n_grid)
    f_star = np.array([eddington(e, d2nu_dpsi2) for e in eps])
    ratio = f_star / hernquist_df(-eps, 1.0, 1.0, G=1.0)
    w = np.interp(-E, eps, ratio, left=0.0)
    w[E >= 0] = 0.0
    return w


@pytest.fixture(scope="module")
def painted():
    """npaint's halo (G = M = a = 1, 2e5 particles), painted by stracers and by the reference."""
    n = 200_000
    pos, vel, pot = sample_hernquist(n, 1.0, 1.0, seed=7, G=1.0)
    E = pot + 0.5 * np.sum(vel**2, axis=1)
    r = np.linalg.norm(pos, axis=1)
    sel = (r < RMAX) & (E < 0)

    res = compute_weights(pos, vel, pot, 1.0 / n, nu_tracer, rmax=RMAX, G=1.0)
    w_ref = np.where(sel, reference_weights(E), 0.0)
    w_str = np.where(sel, res.weights, 0.0)
    return {
        "pos": pos, "vel": vel, "E": E, "r": r, "res": res,
        "w_ref": w_ref / w_ref.sum(), "w_str": w_str / w_str.sum(),
    }


def test_reference_eddington_recovers_hernquist_df():
    """The reference inversion applied to the Hernquist density gives its analytic DF."""
    def d2_hernquist(psi):
        r = 1.0 / psi - 1.0
        nu = 1.0 / (2 * np.pi * r * (1 + r) ** 3)
        L1 = -1.0 / r - 3.0 / (1 + r)
        L2 = 1.0 / r**2 + 3.0 / (1 + r) ** 2
        d1, d2 = nu * L1, nu * (L2 + L1**2)
        p1, p2 = -1.0 / (1 + r) ** 2, 2.0 / (1 + r) ** 3
        return (d2 * p1 - d1 * p2) / p1**3

    for eps in (0.05, 0.3, 0.7, 0.95):
        f = eddington(eps, d2_hernquist)
        assert f == pytest.approx(hernquist_df(-eps, 1.0, 1.0, G=1.0), rel=1e-6)


def test_npaint_tracer_has_positive_df_but_a_core_does_not():
    """In a cusp (psi_0 - psi ~ r), f_* >= 0 needs gS > 0 (npaint's note on gS)."""
    eps = np.linspace(0.05, 0.999, 30)
    f_cusp = [eddington(e, d2nu_dpsi2) for e in eps]
    f_core = [eddington(e, lambda p: d2nu_dpsi2(p, gamma=0.0)) for e in eps]
    assert min(f_cusp) > 0
    assert min(f_core) < 0
    assert f_core[-1] < 0                     # negative at the most bound energies


def test_no_negative_df_bins(painted):
    assert painted["res"].n_negative_f == 0


def test_shell_masses_match_npaint(painted):
    edges = np.geomspace(0.02, 5.0, 9)
    m_str, _ = np.histogram(painted["r"], edges, weights=painted["w_str"])
    m_ref, _ = np.histogram(painted["r"], edges, weights=painted["w_ref"])
    np.testing.assert_allclose(m_str / m_ref, 1.0, atol=0.03)


def test_energy_distribution_matches_npaint(painted):
    """Stellar mass in deciles of the particle energies (most weakly bound decile looser)."""
    E = painted["E"]
    bound = painted["w_ref"] > 0
    edges = np.quantile(E[bound], np.linspace(0, 1, 11))
    edges[-1] = 0.0
    m_str, _ = np.histogram(E, edges, weights=painted["w_str"])
    m_ref, _ = np.histogram(E, edges, weights=painted["w_ref"])
    ratio = m_str / m_ref
    np.testing.assert_allclose(ratio[:-1], 1.0, atol=0.05)    # most bound first
    assert ratio[-1] == pytest.approx(1.0, abs=0.08)          # truncation at rmax


def test_velocity_dispersions_match_npaint(painted):
    bins = np.geomspace(0.05, 5.0, 7)
    p_str = velocity_profile(painted["pos"], painted["vel"], painted["w_str"], bins=bins)
    p_ref = velocity_profile(painted["pos"], painted["vel"], painted["w_ref"], bins=bins)
    np.testing.assert_allclose(p_str.sigma_r, p_ref.sigma_r, rtol=0.01)
    np.testing.assert_allclose(p_str.sigma_3d, p_ref.sigma_3d, rtol=0.01)
