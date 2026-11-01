import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from conftest import G
from stracers import compute_weights, get_profile
from stracers.kinematics import (
    anisotropy,
    bulk_velocity,
    spherical_velocities,
    velocity_profile,
    weighted_dispersion,
    weighted_mean,
)


def test_spherical_velocities_simple_cases():
    pos = np.array([[2.0, 0, 0], [0, 3.0, 0], [1.0, 0, 0], [0, 0, 5.0]])
    vel = np.array([[7.0, 0, 0], [-4.0, 0, 0], [0, 0, 2.0], [1.0, 0, 0]])
    v_r, v_theta, v_phi = spherical_velocities(pos, vel)
    np.testing.assert_allclose(v_r, [7.0, 0, 0, 0])
    np.testing.assert_allclose(v_phi, [0, 4.0, 0, 0])      # x-motion at +y is clockwise
    np.testing.assert_allclose(v_theta, [0, 0, -2.0, 0])   # +z at the equator is towards theta=0
    # on the z axis phi is undefined: the motion is left out of v_theta and v_phi
    np.testing.assert_allclose([v_theta[3], v_phi[3]], [0, 0])


def test_spherical_velocities_preserve_speed():
    rng = np.random.default_rng(0)
    pos, vel = rng.normal(size=(1000, 3)), rng.normal(size=(1000, 3)) * 100
    v_r, v_theta, v_phi = spherical_velocities(pos, vel)
    np.testing.assert_allclose(v_r**2 + v_theta**2 + v_phi**2, np.sum(vel**2, axis=1))


def test_weighted_statistics():
    rng = np.random.default_rng(1)
    x = rng.normal(3.0, 2.0, 500)
    assert weighted_mean(x) == pytest.approx(np.mean(x))
    assert weighted_dispersion(x, np.full(500, 7.0)) == pytest.approx(np.std(x))
    # integer weights are equivalent to repeating the values
    w = rng.integers(0, 4, 500)
    rep = np.repeat(x, w)
    assert weighted_mean(x, w) == pytest.approx(np.mean(rep))
    assert weighted_dispersion(x, w) == pytest.approx(np.std(rep))
    vel = rng.normal(size=(500, 3)) + [10.0, -5.0, 0.0]
    np.testing.assert_allclose(bulk_velocity(vel, w), np.mean(np.repeat(vel, w, axis=0), axis=0))


def test_anisotropy():
    assert anisotropy(100.0, 100.0, 100.0) == pytest.approx(0.0)
    assert anisotropy(100.0, 0.0, 0.0) == pytest.approx(1.0)
    assert anisotropy(100.0, 200.0, 200.0) == pytest.approx(-3.0)


def jeans_sigma_r(nu, r, M, a):
    """Isotropic Jeans dispersion of a tracer nu in a Hernquist potential."""
    return np.sqrt([quad(lambda x: nu(x) * G * M / (x + a) ** 2, ri, np.inf, limit=200)[0] / nu(ri)
                    for ri in r])


@pytest.mark.parametrize("tracer", [None, ("hernquist", {"a": 2.0}), ("plummer", {"a": 3.0})])
def test_profile_matches_jeans(hernquist_halo, tracer):
    """The weighted dispersions of tracers in the isotropic Hernquist halo
    follow the isotropic Jeans equation for their density profile."""
    halo = hernquist_halo
    if tracer is None:
        nu, weights = get_profile("hernquist", a=halo["a"]), None
    else:
        nu = get_profile(tracer[0], **tracer[1])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            weights = compute_weights(halo["pos"], halo["vel"], halo["pot"], halo["mass"], nu).weights
    prof = velocity_profile(halo["pos"], halo["vel"], weights, bins=np.geomspace(1.0, 100.0, 7))
    sigma = jeans_sigma_r(nu, prof.r, halo["M"], halo["a"])
    for s in (prof.sigma_r, prof.sigma_theta, prof.sigma_phi):
        assert np.max(np.abs(s / sigma - 1)) < 0.05
    assert np.max(np.abs(prof.beta)) < 0.08
    assert np.max(np.abs(prof.mean_vr / sigma)) < 0.02
    assert np.all(prof.n_eff > 1000)


def test_profile_recovers_rotation(hernquist_halo):
    halo = hernquist_halo
    pos, vel = halo["pos"], halo["vel"].copy()
    R = np.hypot(pos[:, 0], pos[:, 1])
    v_rot = 50.0
    vel[:, 0] += -v_rot * pos[:, 1] / R
    vel[:, 1] += v_rot * pos[:, 0] / R
    prof = velocity_profile(pos, vel, bins=np.geomspace(2.0, 100.0, 5))
    # the added velocity is along phi-hat, so <v_phi> = v_rot (the halo has none)
    np.testing.assert_allclose(prof.mean_vphi, v_rot, rtol=0.01)
    np.testing.assert_allclose(prof.mean_vr, 0.0, atol=1.0)


def test_profile_center_and_vcenter(hernquist_halo):
    halo = hernquist_halo
    w = np.linspace(0, 1, len(halo["pos"]))
    ref = velocity_profile(halo["pos"], halo["vel"], w, bins=6, rmin=1.0, rmax=100.0)
    shifted = velocity_profile(halo["pos"] + [30.0, -10.0, 5.0], halo["vel"] + [200.0, 0.0, -50.0], w,
                               bins=6, rmin=1.0, rmax=100.0,
                               center=[30.0, -10.0, 5.0], vcenter=[200.0, 0.0, -50.0])
    for name in ("sigma_r", "sigma_theta", "sigma_phi", "mean_vr", "mass", "n_eff"):
        np.testing.assert_allclose(getattr(shifted, name), getattr(ref, name), rtol=1e-8, atol=1e-8)


def test_profile_empty_bins_and_errors(hernquist_halo):
    halo = hernquist_halo
    r_out = np.linalg.norm(halo["pos"], axis=1).max()
    prof = velocity_profile(halo["pos"], halo["vel"], bins=r_out * np.array([2.0, 3.0, 4.0]))
    assert np.all(prof.n == 0) and np.all(np.isnan(prof.sigma_r))
    with pytest.raises(ValueError):
        velocity_profile(halo["pos"], halo["vel"], np.ones(3))
    with pytest.raises(ValueError):
        velocity_profile(halo["pos"], halo["vel"], -np.ones(len(halo["pos"])))
