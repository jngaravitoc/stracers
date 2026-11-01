import warnings

import numpy as np
import pytest
from scipy.integrate import quad, trapezoid

from conftest import G, hernquist_df, weighted_density
from stracers import compute_weights, get_profile, truncate
from stracers.weights import density_of_states, eddington


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
    # the outermost bins rely on the power-law extension of nu(psi) beyond the
    # grid; the most weakly bound one is ~17% off
    assert np.max(np.abs(ratio[sel & ~inside] - 1)) < 0.2


def test_self_consistent_weights_are_constant(hernquist_halo):
    halo = hernquist_halo
    res = run(halo, get_profile("hernquist", a=halo["a"]))
    counts = np.histogram(particle_energy(halo)[res.mask], res.E_edges)[0]
    w = res.w_bins[counts > 500]
    assert np.std(w) / np.mean(w) < 0.05


def test_density_of_states_matches_analytic(hernquist_halo):
    """g(E) of Laporte+13 eq. (3) for the Hernquist potential, integrated out
    to the outermost selected particle, i.e. over the same volume as the
    particles counted in N(E) (not a boundary derived from the fine grid)."""
    halo = hernquist_halo
    M, a = halo["M"], halo["a"]
    res = run(halo, get_profile("hernquist", a=2.0))
    r_part = np.linalg.norm(halo["pos"], axis=1)
    r_out = r_part[res.mask].max()
    assert res.r_range[1] == pytest.approx(r_out)

    def g_true(e):
        r_e = min(-G * M / e - a, r_out)
        integrand = lambda x: x**2 * np.sqrt(2 * max(e + G * M / (x + a), 0.0))
        return (4 * np.pi) ** 2 * quad(integrand, 0, r_e, limit=200)[0]

    counts = np.histogram(particle_energy(halo)[res.mask], res.E_edges)[0]
    sel = counts > 200
    ratio = res.g_E[sel] / np.array([g_true(e) for e in res.E_centers[sel]])
    assert np.max(np.abs(ratio - 1)) < 0.03

    # the inner sphere r < r_range[0] matters for the most bound energies
    r = res.r_bins
    dlnr = np.log(r[1] / r[0])
    dr = r * (np.exp(dlnr / 2) - np.exp(-dlnr / 2))
    E_deep = res.E_centers[1:4]
    g_no_inner = density_of_states(E_deep, r, dr, -res.psi)
    g_inner = density_of_states(E_deep, r, dr, -res.psi, r_in=res.r_range[0])
    truth = np.array([g_true(e) for e in E_deep])
    assert np.all(np.abs(g_inner / truth - 1) < np.abs(g_no_inner / truth - 1))


def test_weakly_bound_weights_not_biased(hernquist_halo):
    """Regression: the fine grid used to stop at the last coarse bin centre,
    so g(E) missed the outer half of the last shell and the self-consistent
    weights of the most weakly bound particles were ~30% too low."""
    halo = hernquist_halo
    res = run(halo, get_profile("hernquist", a=halo["a"]))
    counts = np.histogram(particle_energy(halo)[res.mask], res.E_edges)[0]
    w = res.w_bins[counts > 200]
    w = w / np.median(w)
    assert np.max(np.abs(w[-5:] - 1)) < 0.08


def test_number_counts_and_density_of_states_same_volume(hernquist_halo):
    """With sparse outer bins dropped, N(E) only counts the particles inside
    the volume covered by g(E)."""
    halo = hernquist_halo
    res = run(halo, get_profile("hernquist", a=2.0), min_per_bin=5000)
    r = np.linalg.norm(halo["pos"], axis=1)
    r_hi = res.r_range[1]
    assert r_hi < r[res.mask].max()  # some outer bins were dropped
    counted = np.sum(res.N_E * np.diff(res.E_edges))
    assert counted == pytest.approx(halo["mass"] * np.sum(res.mask & (r <= r_hi)))
    assert np.all(np.isfinite(res.weights)) and np.all(res.weights >= 0)


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


def test_profile_vanishing_at_rmax_no_warning(hernquist_halo):
    """A taper that underflows to zero at rmax must not warn about 0/0."""
    halo = hernquist_halo
    nu_t = truncate(get_profile("hernquist", a=2.0), r_t=50.0, sharpness=4)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        warnings.filterwarnings("ignore", message=".*f\\(E\\) < 0.*")
        res = compute_weights(halo["pos"], halo["vel"], halo["pot"], halo["mass"], nu_t)
    assert np.all(np.isfinite(res.weights))


def test_gradient_mode_excludes_dropped_intermediate_shell(hernquist_halo):
    """derivative="gradient": g(E) only covers the bins kept, so particles in
    a dropped intermediate shell must not be counted in N(E). The spline grid
    is continuous across the gap, so there they are counted."""
    halo = hernquist_halo
    r = np.linalg.norm(halo["pos"], axis=1)
    keep = ~((r > 20.0) & (r < 32.0)) | (np.arange(len(r)) % 2000 == 0)  # empty a shell
    sub = {k: (v[keep] if np.ndim(v) else v) for k, v in halo.items() if k != "ids"}
    r_sub = r[keep]
    nu = get_profile("hernquist", a=2.0)

    res = run(sub, nu, derivative="gradient")
    lo, hi = res.bin_edges[:, 0], res.bin_edges[:, 1]
    gaps = lo[1:] > hi[:-1] * (1 + 1e-12)
    assert np.any(gaps), "the test needs a dropped intermediate bin"
    in_bins = np.zeros(len(r_sub), dtype=bool)
    for a, b in res.bin_edges:
        in_bins |= (r_sub >= a) & (r_sub <= b)
    covered = res.mask & ((r_sub < res.r_range[0]) | in_bins)
    counted = np.sum(res.N_E * np.diff(res.E_edges))
    assert counted == pytest.approx(sub["mass"] * covered.sum())
    assert covered.sum() < res.mask.sum()  # the shell's particles were left out

    res = run(sub, nu, derivative="spline")
    counted = np.sum(res.N_E * np.diff(res.E_edges))
    assert counted == pytest.approx(sub["mass"] * np.sum(res.mask & (r_sub <= res.r_range[1])))


def test_outer_slope_jump_is_accounted():
    """The power-law extension with a slope different from the interior
    derivative puts a jump in dnu/dpsi at psi_min; its contribution must be
    included so that the DF does not depend on the slope chosen when nu(psi)
    is an exact power law, and is consistent otherwise."""
    k = 3.0
    psi = np.linspace(1.0, 10.0, 4000)[::-1]           # decreasing with radius
    nu = psi**k
    dnu = k * psi ** (k - 1)
    d2 = k * (k - 1) * psi ** (k - 2)
    eps = np.array([0.5, 2.0, 5.0, 9.0])
    f_local = eddington(eps, psi, d2, nu, dnu)
    f_exact = eddington(eps, psi, d2, nu, dnu, outer_slope=k)
    np.testing.assert_allclose(f_exact, f_local, rtol=1e-12)
    # the analytic DF of nu = psi^k: k (k-1) B(k-1, 1/2) eps^(k-3/2) / (sqrt(8) pi^2)
    from scipy.special import beta
    f_true = k * (k - 1) * beta(k - 1, 0.5) * eps ** (k - 1.5) / (np.sqrt(8) * np.pi**2)
    np.testing.assert_allclose(f_exact, f_true, rtol=1e-3)

    # with a different extension slope, nu is continuous at psi_min and the
    # jump term keeps the result equal to the DF of that piecewise profile:
    # interior psi^k joined to (psi/psi_min)^k2 below psi_min
    k2 = 4.0
    f_jump = eddington(eps, psi, d2, nu, dnu, outer_slope=k2)
    pm, num = psi[-1], nu[-1]
    grid = np.linspace(1e-6, 10.0, 400001)
    nu_pw = np.where(grid >= pm, grid**k, num * (grid / pm) ** k2)
    d1_pw = np.gradient(nu_pw, grid)
    # Eddington in the integrated-by-parts form avoids the delta function:
    # f = d/deps [ int_0^eps dnu/dpsi / sqrt(eps - psi) dpsi ] / (sqrt(8) pi^2)
    def inner(e):
        g = grid[grid < e]
        return trapezoid(np.interp(g, grid, d1_pw) / np.sqrt(e - g + 1e-300), g)
    h = 1e-3
    f_ref = np.array([(inner(e + h) - inner(e - h)) / (2 * h) for e in eps]) / (np.sqrt(8) * np.pi**2)
    above = eps > pm
    np.testing.assert_allclose(f_jump[above], f_ref[above], rtol=2e-2)
