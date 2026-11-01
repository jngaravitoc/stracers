"""
Tests on a 1M-particle sample of the dark matter halo of an isolated Milky Way
simulation (Gadget-4, 10^8 particles, snapshots at t = 0 to 5 time units of
0.98 Gyr), made with tests/make_sim_sample.py.

The sample is not in the repository: the tests are skipped unless it exists
at $STRACERS_SIM_SAMPLE (default ~/.cache/stracers/mw_iso_b1_1M.hdf5).
"""

import os
import pathlib
import warnings

import numpy as np
import pytest
from scipy.integrate import quad

import stracers

SAMPLE = pathlib.Path(os.environ.get(
    "STRACERS_SIM_SAMPLE", pathlib.Path.home() / ".cache" / "stracers" / "mw_iso_b1_1M.hdf5"))
pytestmark = pytest.mark.skipif(not SAMPLE.exists(), reason="simulation sample not available")

EDGES = np.geomspace(20.0, 100.0, 6)


def load(snapshot):
    h5py = pytest.importorskip("h5py")
    with h5py.File(SAMPLE, "r") as f:
        g = f[snapshot]
        return dict(pos=g["pos"][()].astype(float), vel=g["vel"][()].astype(float),
                    pot=g["pot"][()].astype(float), ids=g["ids"][()], t=float(g.attrs["time"]),
                    mass=float(g.attrs["mass"]), center=g.attrs["center"], vcenter=g.attrs["vcenter"])


@pytest.fixture(scope="module")
def snaps():
    return {s: load(s) for s in ("000", "044", "110", "250")}


def paint(snap, tracer, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return stracers.paint_tracer(snap["pos"], snap["vel"], snap["pot"], snap["mass"], tracer,
                                     (20, 100), n_stars=3000, ids=snap["ids"], center=snap["center"],
                                     **kwargs)


def observe(pop, snap, label=None):
    return pop.at_snapshot(snap["ids"], snap["pos"], snap["vel"], center=snap["center"], snapshot=label)


def shells(pop):
    return np.histogram(np.linalg.norm(pop.pos, axis=1), EDGES, weights=pop.n_stars)[0]


@pytest.mark.parametrize("tracer, tol", [("RRLyr", 0.03), ("BHB", 0.06)])
def test_profile_recovered_in_simulation(snaps, tracer, tol):
    """The stars painted on the simulated (radially anisotropic) halo follow
    the literature profile of the tracer: within 1% (RR Lyrae, Medina+2024)
    and 4% (BHB, Thomas+2018) in shells over 20-100 kpc."""
    pop = paint(snaps["000"], tracer)
    shell = np.array([quad(lambda x: 4 * np.pi * x**2 * pop.profile(x), lo, hi)[0]
                      for lo, hi in zip(EDGES[:-1], EDGES[1:])])
    np.testing.assert_allclose(shells(pop) / (shell / shell.sum() * 3000), 1.0, atol=tol)


@pytest.mark.parametrize("tracer", ["RRLyr", "BHB"])
def test_population_is_stationary(snaps, tracer):
    """A population painted on a relaxed snapshot (t = 2.2) of the isolated
    halo keeps its number of stars, shell counts and dispersions until t = 5:
    an f(E) population in a static potential is a steady state."""
    pop = paint(snaps["110"], tracer)
    later = observe(pop, snaps["250"])
    assert later.matched_fraction == pytest.approx(1.0)
    assert later.total == pytest.approx(3000, rel=0.01)
    np.testing.assert_allclose(shells(later) / shells(pop), 1.0, atol=0.03)
    s0 = pop.velocity_profile(bins=EDGES, vcenter=snaps["110"]["vcenter"]).sigma_r
    s1 = later.velocity_profile(bins=EDGES, vcenter=snaps["250"]["vcenter"]).sigma_r
    np.testing.assert_allclose(s1 / s0, 1.0, atol=0.02)


def test_population_from_initial_conditions(snaps):
    """Painted on the initial conditions, the population relaxes with the halo
    during the first ~1 Gyr (a few per cent fewer stars in the range), then
    stays constant."""
    pop = paint(snaps["000"], "RRLyr")
    n = {s: observe(pop, snaps[s]).total for s in ("044", "110", "250")}
    assert all(0.9 < v / 3000 < 1.0 for v in n.values())
    assert n["250"] == pytest.approx(n["110"], rel=0.005)
    assert n["250"] == pytest.approx(n["044"], rel=0.01)


def test_stars_beyond_rmax_are_followed(snaps):
    """With weights only inside rmax, stars whose orbits reach beyond rmax are
    lost over time (the shallow BHB profile loses ~6% in its outer shell);
    with weights by energy beyond rmax (the paint_tracer default) they are
    not."""
    zero = paint(snaps["110"], "BHB", beyond_rmax="zero")
    energy = paint(snaps["110"], "BHB")
    loss_zero = 1 - shells(observe(zero, snaps["250"]))[-1] / shells(zero)[-1]
    loss_energy = 1 - shells(observe(energy, snaps["250"]))[-1] / shells(energy)[-1]
    assert loss_zero > 0.04
    assert abs(loss_energy) < 0.03


def test_saved_population_and_catalogue(snaps, tmp_path):
    """Paint once on the initial conditions, save, load and observe a later
    snapshot; write its mock catalogue and read it back."""
    pop = paint(snaps["000"], "RRLyr", snapshot="000")
    pop.save(tmp_path / "rrl.h5")
    loaded = stracers.load_population(tmp_path / "rrl.h5")
    a, b = observe(pop, snaps["250"]), observe(loaded, snaps["250"], label="250")
    np.testing.assert_allclose(b.n_stars, a.n_stars)

    n_h5 = stracers.write_catalogue(b, tmp_path / "rrl_250.h5", rng=0, vcenter=snaps["250"]["vcenter"])
    n_csv = stracers.write_catalogue(b, tmp_path / "rrl_250.csv", rng=0, vcenter=snaps["250"]["vcenter"])
    cols, meta = stracers.read_catalogue(tmp_path / "rrl_250.h5")
    cols_csv, _ = stracers.read_catalogue(tmp_path / "rrl_250.csv")
    assert n_h5 == n_csv == len(cols["x"])
    assert abs(n_h5 - b.total) < 5 * np.sqrt(b.total)
    np.testing.assert_array_equal(cols["id"], cols_csv["id"])
    assert np.all((cols["r"] >= 20) & (cols["r"] <= 100))
    assert meta["key"] == "medina2024"
    assert meta["snapshot"] == "250"
