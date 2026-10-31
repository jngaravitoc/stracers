import warnings

import numpy as np
import pytest

from conftest import weighted_density
from stracers import broken_power_law, compute_weights, mw, power_law


def test_table_contents():
    table = mw.load_table()
    assert len(table) == 31
    keys = [e["key"] for e in table]
    assert len(set(keys)) == len(keys)
    assert len(mw.available("spl")) == 16
    assert len(mw.available("bpl")) == 15
    for e in table:
        if e["model"] == "spl":
            assert e["slope"] is not None
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


@pytest.mark.parametrize("key", ["medina2024", "feng2024", "fukushima2019"])
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
