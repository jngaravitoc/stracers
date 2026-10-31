import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from conftest import G, hernquist_df, weighted_density
from stracers import compute_weights, get_profile, truncate
from stracers.weights import density_of_states


def run(halo, profile, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return compute_weights(
            halo["pos"], halo["vel"], halo["pot"], halo["mass"], profile, **kwargs
        )


def particle_energy(halo):
    return 0.5 * np.sum(halo["vel"] ** 2, axis=1) + halo["pot"]


def test_eddington_recovers_hernquist_df(hernquist_halo):
    """For the halo's own profile, f(E) must match the analytic DF."""
    halo = hernquist_halo
    res = run(halo, get_profile("hernquist", a=halo["a"]))
    counts = np.histogram(particle_energy(halo)[res.mask], res.E_edges)[0]
    ratio = res.f_E / hernquist_df(res.E_centers, halo["M"], halo["a"])
    sel = counts > 200
    ratio /= np.median(ratio[sel])
    # energies resolved by the radial bins vs. those relying on the
    # power-law extrapolation beyond rmax
    inside = sel & (-res.E_centers > res.psi.min() + 0.05 * np.ptp(res.psi))
    assert np.max(np.abs(ratio[inside] - 1)) < 0.03
    assert np.max(np.abs(ratio[sel & ~inside] - 1)) < 0.15


def test_self_consistent_weights_are_constant(hernquist_halo):
    halo = hernquist_halo
    res = run(halo, get_profile("hernquist", a=halo["a"]))
    counts = np.histogram(particle_energy(halo)[res.mask], res.E_edges)[0]
    w = res.w_bins[counts > 500]
    assert np.std(w) / np.mean(w) < 0.05


def test_density_of_states_matches_analytic(hernquist_halo):
    """g(E) of Laporte+13 eq. (3) for the Hernquist potential, truncated at
    the outermost radial bin like the particle sample."""
    halo = hernquist_halo
    M, a = halo["M"], halo["a"]
    res = run(halo, get_profile("hernquist", a=2.0))
    r = res.r_bins
    r_out = r[-1] * np.sqrt(r[-1] / r[-2])

    def g_true(e):
        r_e = min(-G * M / e - a, r_out)
        integrand = lambda x: x**2 * np.sqrt(2 * max(e + G * M / (x + a), 0.0))
        return (4 * np.pi) ** 2 * quad(integrand, 0, r_e, limit=200)[0]

    counts = np.histogram(particle_energy(halo)[res.mask], res.E_edges)[0]
    sel = counts > 200
    ratio = res.g_E[sel] / np.array([g_true(e) for e in res.E_centers[sel]])
    assert np.max(np.abs(ratio - 1)) < 0.03

    # the inner sphere r < r_1 matters for the most bound energies
    dr = r * np.log(r[1] / r[0])
    E_deep = res.E_centers[1:4]
    g_no_inner = density_of_states(E_deep, r, dr, -res.psi)
    g_inner = density_of_states(E_deep, r, dr, -res.psi, r_in=r[0] - 0.5 * dr[0])
    truth = np.array([g_true(e) for e in E_deep])
    assert np.all(np.abs(g_inner / truth - 1) < np.abs(g_no_inner / truth - 1))


def test_energy_uses_particle_potential(hernquist_halo):
    """Laporte+13: <Phi>(r) sets nu(psi), but E uses each particle's own
    potential. Deepening the potential of some particles must make them more
    bound and so change their weights."""
    halo = hernquist_halo
    nu = get_profile("hernquist", a=2.0)
    ref = run(halo, nu)
    pot = halo["pot"].copy()
    r = np.linalg.norm(halo["pos"], axis=1)
    sub = np.where((r > 5) & (r < 6))[0][:50]
    pot[sub] *= 1.2
    new = compute_weights(halo["pos"], halo["vel"], pot, halo["mass"], nu)
    assert np.all(new.weights[sub] != ref.weights[sub])
    assert np.mean(new.weights[sub]) > np.mean(ref.weights[sub])


@pytest.mark.parametrize(
    "name, params",
    [
        ("hernquist", {"a": 2.0}),
        ("plummer", {"a": 3.0}),
        ("einasto", {"n": 3.0, "r_eff": 15.0}),
        ("spl", {"slope": -3.5, "r_core": 2.0}),
        ("dpl", {"r_break": 25.0, "inner_slope": -2.5, "outer_slope": -4.5, "alpha": 2.0}),
    ],
)
def test_recovers_tracer_profile(hernquist_halo, name, params):
    halo = hernquist_halo
    nu = get_profile(name, **params)
    res = run(halo, nu, stellar_mass=1e8)

    r = np.linalg.norm(halo["pos"], axis=1)
    edges = np.logspace(np.log10(0.5), np.log10(100.0), 12)
    rc, rho = weighted_density(r, res.weights * halo["mass"], edges)
    ratio = rho / nu(rc)
    ratio /= np.median(ratio)
    assert np.max(np.abs(ratio - 1)) < 0.1


def test_gradient_mode_runs(hernquist_halo):
    halo = hernquist_halo
    nu = get_profile("hernquist", a=2.0)
    res = run(halo, nu, derivative="gradient", dlogr=0.05)
    r = np.linalg.norm(halo["pos"], axis=1)
    rc, rho = weighted_density(r, res.weights * halo["mass"], np.logspace(0, 2, 8))
    ratio = rho / nu(rc)
    assert np.max(np.abs(ratio / np.median(ratio) - 1)) < 0.2


def test_normalisation_and_mask(hernquist_halo):
    halo = hernquist_halo
    rmax = 100.0
    res = run(halo, get_profile("plummer", a=3.0), stellar_mass=5e7, rmax=rmax)
    r = np.linalg.norm(halo["pos"], axis=1)

    assert res.weights.shape == r.shape
    assert np.sum(res.weights * halo["mass"]) == pytest.approx(5e7)
    assert np.all(res.weights[r >= rmax] == 0)
    assert np.array_equal(res.mask, r < rmax)
    assert np.all(res.weights >= 0)


def test_center_shift(hernquist_halo):
    halo = hernquist_halo
    nu = get_profile("hernquist", a=2.0)
    offset = np.array([50.0, -20.0, 10.0])
    ref = run(halo, nu)
    shifted = compute_weights(
        halo["pos"] + offset, halo["vel"], halo["pot"], halo["mass"], nu, center=offset
    )
    np.testing.assert_allclose(shifted.weights, ref.weights, rtol=1e-8)


def test_truncation_correction_is_a_constant_shift(hernquist_halo):
    halo = hernquist_halo
    nu = get_profile("hernquist", a=2.0)
    ref = run(halo, nu, boundary_term=False)
    corr = run(halo, nu, boundary_term=False, truncation_correction=True)
    np.testing.assert_allclose(corr.weights, ref.weights, rtol=1e-6)


def test_bad_derivative_option(hernquist_halo):
    halo = hernquist_halo
    with pytest.raises(ValueError):
        run(halo, get_profile("hernquist", a=2.0), derivative="magic")


def test_stellar_fraction_inside_rmax(hernquist_halo):
    halo = hernquist_halo
    r = np.linalg.norm(halo["pos"], axis=1)
    res = run(halo, get_profile("hernquist", a=2.0), stellar_fraction=0.01, rmax=100.0)
    m_dm = halo["mass"] * np.sum((r > 0) & (r < 100.0))
    assert np.sum(res.weights * halo["mass"]) == pytest.approx(0.01 * m_dm)
    assert res.stellar_mass == pytest.approx(0.01 * m_dm)


def test_total_stellar_mass(hernquist_halo):
    halo = hernquist_halo
    a, rmax = 2.0, 100.0
    res = run(halo, get_profile("hernquist", a=a), stellar_mass=1e8, within="total", rmax=rmax)
    inside = 1e8 * (rmax / (rmax + a)) ** 2  # Hernquist M(<r) / M
    assert np.sum(res.weights * halo["mass"]) == pytest.approx(inside, rel=1e-6)


def test_total_stellar_fraction_uses_all_particles(hernquist_halo):
    halo = hernquist_halo
    a, rmax = 2.0, 100.0
    res = run(halo, get_profile("hernquist", a=a), stellar_fraction=0.01, within="total", rmax=rmax)
    m_dm_total = halo["mass"] * len(halo["pos"])
    inside = 0.01 * m_dm_total * (rmax / (rmax + a)) ** 2
    assert res.stellar_mass == pytest.approx(inside, rel=1e-6)


def test_stellar_mass_option_errors(hernquist_halo):
    halo = hernquist_halo
    nu = get_profile("hernquist", a=2.0)
    with pytest.raises(ValueError, match="not both"):
        run(halo, nu, stellar_mass=1e8, stellar_fraction=0.01)
    with pytest.raises(ValueError, match="within"):
        run(halo, nu, within="virial")
    with pytest.raises(ValueError, match="does not converge"):
        run(halo, get_profile("nfw", rs=20.0, rvir=200.0), within="total")


def test_energy_truncation(hernquist_halo):
    """r_trunc: no stellar mass beyond r_t, nu(r) followed well inside it."""
    halo = hernquist_halo
    nu = get_profile("hernquist", a=5.0)
    r_t = 50.0
    res = run(halo, nu, stellar_mass=1e8, r_trunc=r_t)
    r = np.linalg.norm(halo["pos"], axis=1)

    assert np.sum(res.weights[r > 1.01 * r_t]) == 0.0
    assert np.sum(res.weights * halo["mass"]) == pytest.approx(1e8)

    rc, rho = weighted_density(r, res.weights * halo["mass"], np.logspace(np.log10(0.5), 1, 8))
    ratio = rho / nu(rc)
    assert np.max(np.abs(ratio / np.median(ratio) - 1)) < 0.1


def test_energy_truncation_total_mass_is_inside(hernquist_halo):
    halo = hernquist_halo
    res = run(halo, get_profile("hernquist", a=5.0), stellar_mass=1e8,
              within="total", r_trunc=50.0, rmax=100.0)
    assert res.stellar_mass == pytest.approx(1e8)


def test_energy_truncation_range(hernquist_halo):
    halo = hernquist_halo
    nu = get_profile("hernquist", a=5.0)
    with pytest.raises(ValueError, match="larger than rmax"):
        run(halo, nu, r_trunc=200.0, rmax=100.0)
    with pytest.raises(ValueError, match="radial range"):
        run(halo, nu, r_trunc=1e-4)


@pytest.mark.parametrize("sharpness", [2.0, 8.0])
def test_tapered_profile(hernquist_halo, sharpness):
    """truncate(): the tapered profile is recovered, including its edge."""
    halo = hernquist_halo
    r_t = 50.0
    nu_t = truncate(get_profile("hernquist", a=5.0), r_t, sharpness)
    res = run(halo, nu_t, stellar_mass=1.0)
    r = np.linalg.norm(halo["pos"], axis=1)

    rc, rho = weighted_density(r, res.weights * halo["mass"], np.logspace(np.log10(0.5), np.log10(r_t), 10))
    ratio = rho / nu_t(rc)
    assert np.max(np.abs(ratio / np.median(ratio) - 1)) < 0.1

    expected = quad(lambda x: 4 * np.pi * x**2 * nu_t(x), r_t, 300.0)[0] / quad(
        lambda x: 4 * np.pi * x**2 * nu_t(x), 0, 300.0
    )[0]
    assert np.sum(res.weights[r > r_t] * halo["mass"]) == pytest.approx(expected, rel=0.15)


def test_total_mass_power_law_profiles(hernquist_halo):
    """within='total' for user power laws: analytic BPL enclosed mass, divergent
    SPL rejected, and accepted again once tapered."""
    halo = hernquist_halo
    rb, gi, go, rmax = 18.1, -2.05, -4.47, 100.0
    m_in = rb**3 / (3 + gi)
    frac = (m_in + rb**3 * (1 - (rmax / rb) ** (3 + go)) / -(3 + go)) / (m_in + rb**3 / -(3 + go))
    nu = get_profile("bpl", r_break=rb, inner_slope=gi, outer_slope=go, smoothness=0.0)
    res = run(halo, nu, stellar_mass=1e9, within="total", rmax=rmax)
    assert res.stellar_mass == pytest.approx(1e9 * frac, rel=1e-5)

    spl = get_profile("spl", slope=-2.5, r_core=1.0)
    with pytest.raises(ValueError, match="does not converge"):
        run(halo, spl, stellar_mass=1e9, within="total")
    res = run(halo, truncate(spl, 150.0), stellar_mass=1e9, within="total", rmax=rmax)
    assert 0 < res.stellar_mass < 1e9
