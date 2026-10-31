import numpy as np
import pytest
from scipy.integrate import quad

from stracers import einasto, get_profile, hernquist, mass_fraction_within, nfw, plummer, truncate


def enclosed(nu, rmax):
    return quad(lambda r: 4 * np.pi * r**2 * nu(r), 0, rmax, limit=200)[0]


@pytest.mark.parametrize(
    "nu",
    [
        lambda r: plummer(r, a=2.0),
        lambda r: hernquist(r, a=5.0),
        lambda r: einasto(r, n=4.0, r_eff=20.0),
    ],
)
def test_unit_mass(nu):
    assert enclosed(nu, np.inf) == pytest.approx(1.0, rel=1e-4)


def test_nfw_unit_mass_inside_rvir():
    assert enclosed(lambda r: nfw(r, rs=20.0, rvir=200.0), 200.0) == pytest.approx(1.0, rel=1e-6)


def test_get_profile():
    r = np.logspace(-1, 2, 10)
    np.testing.assert_allclose(get_profile("Hernquist", a=3.0)(r), hernquist(r, 3.0))
    with pytest.raises(ValueError):
        get_profile("king", rc=1.0)


def test_einasto_rejects_small_n():
    with pytest.raises(ValueError):
        einasto(1.0, n=0.4, r_eff=1.0)


def test_mass_fraction_within():
    assert mass_fraction_within(lambda r: hernquist(r, a=2.0), 100.0) == pytest.approx(
        (100.0 / 102.0) ** 2, rel=1e-6
    )
    # normalisation of nu does not matter
    assert mass_fraction_within(lambda r: 7.0 * plummer(r, a=3.0), 10.0) == pytest.approx(
        10.0**3 / (10.0**2 + 3.0**2) ** 1.5, rel=1e-6
    )
    with pytest.raises(ValueError):
        mass_fraction_within(lambda r: nfw(r, rs=20.0, rvir=200.0), 100.0)


def test_truncate():
    nu = lambda r: hernquist(r, a=2.0)
    nu_t = truncate(nu, 10.0, sharpness=2.0)
    r = np.array([1.0, 10.0, 30.0])
    np.testing.assert_allclose(nu_t(r), nu(r) * np.exp(-((r / 10.0) ** 2)))
    # tapering gives NFW a finite total mass
    assert 0 < mass_fraction_within(truncate(lambda r: nfw(r, 20.0, 200.0), 100.0), 100.0) < 1


def test_double_power_law_special_cases():
    from stracers import double_power_law

    r = np.logspace(-2, 3, 30)

    def proportional(a, b):
        ratio = a / b
        np.testing.assert_allclose(ratio, ratio[0], rtol=1e-10)

    proportional(double_power_law(r, 3.0, -1.0, -4.0, alpha=1.0), hernquist(r, 3.0))
    proportional(double_power_law(r, 20.0, -1.0, -3.0, alpha=1.0), nfw(r, 20.0, 200.0))
    proportional(double_power_law(r, 2.0, 0.0, -5.0, alpha=2.0), plummer(r, 2.0))
    assert double_power_law(7.0, 7.0, -2.0, -4.0, alpha=3.0) == pytest.approx(1.0)


def test_double_power_law_slopes_and_sharpness():
    from stracers import broken_power_law, double_power_law

    def slope(f, r):
        rr = r * np.array([0.99, 1.0, 1.01])
        return np.gradient(np.log(f(rr)), np.log(rr))[1]

    f = lambda r: double_power_law(r, 20.0, -2.5, -4.5, alpha=2.0)
    assert slope(f, 1e-3) == pytest.approx(-2.5, abs=1e-3)
    assert slope(f, 1e5) == pytest.approx(-4.5, abs=1e-3)
    # large alpha approaches the sharp broken power law (offset ~ 2^(-delta/alpha))
    r = np.array([5.0, 10.0, 40.0, 80.0])
    np.testing.assert_allclose(
        double_power_law(r, 20.0, -2.5, -4.5, alpha=200.0),
        broken_power_law(r, 20.0, -2.5, -4.5),
        rtol=1e-2,
    )


def test_short_names():
    r = np.array([3.0, 30.0])
    np.testing.assert_allclose(get_profile("spl", slope=-3.5, r_core=1.0)(r), (r**2 + 1) ** -1.75)
    np.testing.assert_allclose(
        get_profile("bpl", r_break=20.0, inner_slope=-2.0, outer_slope=-4.0)(r),
        [(3 / 20) ** -2.0, (30 / 20) ** -4.0],
    )
    assert get_profile("dpl", r_break=20.0, inner_slope=-2.0, outer_slope=-4.0, alpha=2.0)(20.0) == 1.0
