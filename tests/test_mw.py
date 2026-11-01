import warnings

import numpy as np
import pytest

from conftest import weighted_density
from stracers import broken_power_law, compute_weights, mw, power_law, triple_power_law


def test_table_contents():
    table = mw.load_table()
    assert len(table) == 32
    keys = [e["key"] for e in table]
    assert len(set(keys)) == len(keys)
    assert len(mw.available("spl")) == 16
    assert len(mw.available("bpl")) == 15
    assert mw.available("tpl") == ["li2026"]
    for e in table:
        if e["model"] == "spl":
            assert e["slope"] is not None
        elif e["model"] == "tpl":
            assert None not in (e["r_break"], e["r_break2"], e["inner_slope"], e["middle_slope"],
                                e["outer_slope"])
        else:
            assert None not in (e["r_break"], e["inner_slope"], e["outer_slope"])
        assert e["r_max"] is not None and e["reference"]


def test_medina2024_entry():
    e = mw.get_entry("medina2024")
    assert e["r_break"] == 18.1
    assert (e["r_break_err_minus"], e["r_break_err_plus"]) == (1.1, 2.1)
    assert (e["inner_slope"], e["outer_slope"]) == (-2.05, -4.47)
    assert (e["r_min"], e["r_max"]) == (7.0, 145.0)
    with pytest.raises(ValueError):
        mw.get_entry("nonexistent2099")


def test_li2026_entry():
    """Li et al. 2026, Table 1, whole halo."""
    e = mw.get_entry("li2026")
    assert (e["r_break"], e["r_break2"]) == (16.0, 76.3)
    assert (e["r_break_err_minus"], e["r_break_err_plus"]) == (0.2, 0.25)
    assert (e["r_break2_err_minus"], e["r_break2_err_plus"]) == (1.94, 2.37)
    assert (e["inner_slope"], e["middle_slope"], e["outer_slope"]) == (-1.5, -3.45, -5.2)
    assert (e["r_min"], e["r_max"], e["tracer"]) == (8.0, 200.0, "K giants")
    assert mw.reference_url("li2026") == "https://arxiv.org/abs/2512.01350"
    nu = mw.mw_profile("li2026", r_core=0.0, smoothness=0.0)
    r = np.array([10.0, 12.0, 30.0, 40.0, 100.0, 120.0])
    assert np.allclose(np.diff(np.log(nu(r)))[::2] / np.diff(np.log(r))[::2], [-1.5, -3.45, -5.2])
    # drawn parameters keep the breaks in order
    rng = np.random.default_rng(3)
    for _ in range(200):
        p = mw.mw_profile("li2026", rng=rng).params
        assert 0 < p["r_break"] < p["r_break2"]


def test_power_law_slopes():
    r = np.array([10.0, 20.0])
    assert np.diff(np.log(power_law(r, -3.5)))[0] / np.log(2) == pytest.approx(-3.5)
    # a core flattens the centre
    assert power_law(0.0, -3.5, r_core=2.0) == pytest.approx(2.0**-3.5)


@pytest.mark.parametrize("smoothness", [0.0, 0.1])
def test_broken_power_law(smoothness):
    def slope(r):
        return np.gradient(np.log(broken_power_law(r, 20.0, -2.0, -4.5, smoothness=smoothness)), np.log(r))

    assert slope(np.array([0.1, 0.11, 0.12]))[1] == pytest.approx(-2.0, abs=1e-3)
    assert slope(np.array([2000.0, 2100.0, 2200.0]))[1] == pytest.approx(-4.5, abs=1e-3)
    assert broken_power_law(20.0, 20.0, -2.0, -4.5) == 1.0


@pytest.mark.parametrize("smoothness", [0.0, 0.1])
def test_triple_power_law(smoothness):
    def slope(r):
        nu = triple_power_law(r, 10.0, 100.0, -1.5, -3.5, -5.0, smoothness=smoothness)
        return np.gradient(np.log(nu), np.log(r))[1]

    assert slope(np.array([0.1, 0.11, 0.12])) == pytest.approx(-1.5, abs=1e-3)
    assert slope(np.array([30.0, 31.0, 32.0])) == pytest.approx(-3.5, abs=1e-2)
    assert slope(np.array([1e4, 1.1e4, 1.2e4])) == pytest.approx(-5.0, abs=1e-3)
    assert triple_power_law(10.0, 10.0, 100.0, -1.5, -3.5, -5.0, smoothness=smoothness) == pytest.approx(1.0)
    # continuous at the second break
    r = 100.0 * np.array([1 - 1e-9, 1 + 1e-9])
    nu = triple_power_law(r, 10.0, 100.0, -1.5, -3.5, -5.0, smoothness=smoothness)
    assert nu[0] == pytest.approx(nu[1], rel=1e-6)
    with pytest.raises(ValueError):
        triple_power_law(1.0, 100.0, 10.0, -1.5, -3.5, -5.0)


def test_mw_profile_parameters_and_sampling():
    nu = mw.mw_profile("medina2024")
    assert nu.params["r_break"] == 18.1

    rng = np.random.default_rng(1)
    rb = np.array([mw.mw_profile("medina2024", rng=rng).params["r_break"] for _ in range(4000)])
    # median and 16/84th percentiles reproduce 18.1 (-1.1, +2.1)
    np.testing.assert_allclose(np.percentile(rb, [16, 50, 84]), [17.0, 18.1, 20.2], atol=0.15)
    # parameters without uncertainties stay fixed
    w = mw.mw_profile("watkins2009", rng=np.random.default_rng(0)).params
    assert (w["r_break"], w["inner_slope"], w["outer_slope"]) == (23.0, -2.4, -4.5)
    # reproducible
    a = mw.mw_profile("feng2024", rng=np.random.default_rng(5)).params
    b = mw.mw_profile("feng2024", rng=np.random.default_rng(5)).params
    assert a == b


@pytest.mark.parametrize("key", ["medina2024", "feng2024", "fukushima2019", "li2026"])
def test_mw_profile_recovered(hernquist_halo, key):
    halo = hernquist_halo
    nu = mw.mw_profile(key)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        res = compute_weights(halo["pos"], halo["vel"], halo["pot"], halo["mass"], nu)
    e = mw.get_entry(key)
    r = np.linalg.norm(halo["pos"], axis=1)
    edges = np.logspace(np.log10(max(e["r_min"], 1.0)), np.log10(min(e["r_max"], 150.0)), 8)
    rc, rho = weighted_density(r, res.weights * halo["mass"], edges)
    ratio = rho / nu(rc)
    assert np.max(np.abs(ratio / np.median(ratio) - 1)) < 0.1


def test_steep_outer_slope_gradient_mode_is_finite(hernquist_halo):
    """Regression: outer slope -15 used to overflow in outer_eddington."""
    halo = hernquist_halo
    nu = mw.mw_profile("fukushima2019", smoothness=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        res = compute_weights(halo["pos"], halo["vel"], halo["pot"], halo["mass"], nu,
                              derivative="gradient")
    assert np.all(np.isfinite(res.weights))


def test_every_profile_has_a_link():
    """Each fit has an arXiv identifier, or a DOI for the papers not on arXiv."""
    for e in mw.load_table():
        assert bool(e["arxiv"]) != bool(e["doi"]), e["key"]
    assert mw.reference_url("medina2024") == "https://arxiv.org/abs/2402.14055"
    assert mw.reference_url("vivas2006") == "https://arxiv.org/abs/astro-ph/0604359"
    assert mw.reference_url("wetterer1996") == "https://doi.org/10.1086/118076"
    # the two fits of one paper share its link
    assert mw.reference_url("xue2015_spl") == mw.reference_url("xue2015_bpl")


def test_plot_mw_profiles():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from stracers.plotting import plot_mw_profiles

    tracers = ["BHB stars", "BHB and BS stars"]
    keys = [e["key"] for e in mw.load_table() if e["tracer"] in tracers]
    ax = plot_mw_profiles(tracers=tracers, r_extrap=200.0)
    solid = [l for l in ax.get_lines() if l.get_linestyle() == "-"]
    dashed = [l for l in ax.get_lines() if l.get_linestyle() == "--"]
    assert len(solid) == len(keys)
    # fits that end before 200 kpc are extrapolated with a dashed line up to it
    n_extrap = sum(mw.get_entry(k)["r_max"] < 200.0 for k in keys)
    assert len(dashed) == n_extrap
    assert all(np.isclose(l.get_xdata()[-1], 200.0) for l in dashed)
    # solid lines span the fitted range, and profiles are 1 at r_norm = 30 kpc
    fuk = solid[keys.index("fukushima2019")]
    assert np.isclose(fuk.get_xdata()[0], 35.0) and np.isclose(fuk.get_xdata()[-1], 360.0)
    das = solid[keys.index("das2016")]                   # fitted over 10-70 kpc
    assert np.isclose(np.interp(np.log(30.0), np.log(das.get_xdata()), das.get_ydata()), 1.0, rtol=1e-3)
    # all the text is in LaTeX math mode
    texts = [t.get_text() for t in ax.get_legend().get_texts()]
    texts += [ax.get_xlabel(), ax.get_ylabel()]
    assert all(t.startswith("$") and t.endswith("$") for t in texts)
    assert not any("&" in t for t in texts)
    plt.close(ax.figure)


def test_drawn_break_radii_are_positive():
    """Regression: faccioli2014 (r_break = 26.5 +/- 8.9 kpc) could draw
    non-positive break radii, giving NaN profiles."""
    rng = np.random.default_rng(0)
    rb = np.array([mw.mw_profile("faccioli2014", rng=rng).params["r_break"] for _ in range(20000)])
    assert np.all(rb > 0)
    assert abs(np.median(rb) - 26.5) < 0.5
    nu = mw.mw_profile("faccioli2014", rng=np.random.default_rng(int(np.argmin(rb))))
    assert np.all(np.isfinite(nu(np.geomspace(0.1, 300, 50))))
