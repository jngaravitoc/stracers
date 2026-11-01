import numpy as np
import pytest
from scipy.integrate import quad

from stracers import get_profile, mw, number_of_stars, number_of_stars_table
from stracers.nstars import FULL_SKY_DEG2

DIVERGENT = {
    "wetterer1996", "vivas2006", "bell2008", "depropris2010",
    "deason2011_spl", "sesar2013", "iorio2018", "thomas2018",
}


def spl_counts(n, n0, r0, r_min):
    """Analytic N(> r_min) for n0 (r / r0)^n, n < -3."""
    return 4 * np.pi * n0 * r0 ** (-n) * r_min ** (3 + n) / -(3 + n)


def test_simple_power_law_analytic():
    nu = get_profile("spl", slope=-4.47)
    got = number_of_stars(nu, 100.0, n0=5.0, r0=8.2)
    assert got == pytest.approx(spl_counts(-4.47, 5.0, 8.2, 100.0), rel=1e-8)


def test_flattened_power_law():
    n, q = -4.47, 0.6
    flat = quad(lambda mu: (1 - mu**2 + mu**2 / q**2) ** (n / 2), 0, 1)[0]
    got = number_of_stars(get_profile("spl", slope=n), 100.0, q=q, n0=5.0, r0=8.2)
    assert got == pytest.approx(spl_counts(n, 5.0, 8.2, 100.0) * flat, rel=1e-8)


def test_finite_range_and_array():
    nu = get_profile("spl", slope=-4.0)
    n = number_of_stars(nu, [50.0, 100.0], 200.0, n0=1.0, r0=1.0)
    exact = 4 * np.pi * (np.array([50.0, 100.0]) ** -1 - 200.0**-1)
    np.testing.assert_allclose(n, exact, rtol=1e-8)


def test_medina2024_numbers():
    """Medina et al. 2024, Sect. 5.2.1: ~3000 RRLs beyond 100 kpc, 0.073 per
    deg^2, ~25 in 350 deg^2 and ~350 beyond 230 kpc (for a 300 kpc edge)."""
    kw = dict(n_ref=3000.0, ref_range=(100.0, 300.0))
    n100, n230 = number_of_stars("medina2024", [100.0, 230.0], 300.0, **kw)
    assert n100 == pytest.approx(3000.0)
    assert n100 / FULL_SKY_DEG2 == pytest.approx(0.073, abs=5e-4)
    assert number_of_stars("medina2024", 100.0, 300.0, area_deg2=350.0, **kw) == pytest.approx(25, abs=1)
    assert n230 == pytest.approx(350.0, rel=0.1)
    # flattening only rescales a power law: the ratio is unchanged
    flat = number_of_stars("medina2024", 230.0, 300.0, q=0.6, **kw)
    assert flat == pytest.approx(n230, rel=1e-3)


def test_key_matches_callable():
    kw = dict(n_ref=1.0, ref_range=(50.0, 200.0))
    by_key = number_of_stars("stringer2021", 100.0, 200.0, **kw)
    by_fn = number_of_stars(mw.mw_profile("stringer2021"), 100.0, 200.0, **kw)
    assert by_key == pytest.approx(by_fn)
    # rng draws the parameters within their uncertainties
    a = number_of_stars("stringer2021", 100.0, 200.0, rng=np.random.default_rng(1), **kw)
    b = number_of_stars("stringer2021", 100.0, 200.0, rng=np.random.default_rng(2), **kw)
    assert a != b


def test_errors():
    with pytest.raises(ValueError, match="arbitrary normalisation"):
        number_of_stars("medina2024", 100.0)
    with pytest.raises(ValueError, match="Unknown profile"):
        number_of_stars("nonexistent2099", 100.0, n_ref=1.0, ref_range=(10.0, 20.0))
    nu = get_profile("spl", slope=-2.5)
    with pytest.raises(ValueError, match="does not converge"):
        number_of_stars(nu, 100.0, n0=1.0, r0=8.0)
    with pytest.raises(ValueError, match="needs r0"):
        number_of_stars(nu, 100.0, 200.0, n0=1.0)
    with pytest.raises(ValueError, match="not both"):
        number_of_stars(nu, 100.0, 200.0, n0=1.0, r0=8.0, n_ref=1.0, ref_range=(1.0, 2.0))


def test_table():
    kw = dict(n_ref=3000.0, ref_range=(100.0, 300.0))
    finite = number_of_stars_table(230.0, 300.0, **kw)
    assert set(finite) == set(mw.available())
    assert all(np.isfinite(v) and v > 0 for v in finite.values())

    with pytest.warns(RuntimeWarning, match="do not converge"):
        infinite = number_of_stars_table(230.0, np.inf, **kw)
    assert {k for k, v in infinite.items() if np.isnan(v)} == DIVERGENT

    # same normalisation: steeper simple power laws leave fewer stars beyond 230 kpc
    spl = sorted(mw.available("spl"), key=lambda k: mw.get_entry(k)["slope"])
    counts = [finite[k] for k in spl]
    assert all(np.diff(counts) >= -1e-9)
