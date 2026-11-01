import warnings

import numpy as np
import pytest

from stracers import paint_tracer, read_catalogue, write_catalogue
from stracers.catalogue import COLUMNS
from stracers.kinematics import spherical_velocities


@pytest.fixture(scope="module")
def population(hernquist_halo):
    halo = hernquist_halo
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return paint_tracer(halo["pos"], halo["vel"], halo["pot"], halo["mass"], "RRLyr", (20, 100),
                            n_stars=5000, ids=halo["ids"], snapshot="000")


@pytest.mark.parametrize("ext", [".h5", ".csv"])
def test_stars_catalogue_roundtrip(population, tmp_path, ext):
    path = tmp_path / ("stars" + ext)
    n = write_catalogue(population, path, rng=1, vcenter=[1.0, 2.0, 3.0])
    cols, meta = read_catalogue(path)
    assert list(cols) == [c for c in COLUMNS if c != "n_stars"]
    assert len(cols["x"]) == n == meta["n_rows"]
    assert abs(n - 5000) < 5 * np.sqrt(5000)                     # Poisson around the expected total
    # each star is on a particle of the population, with its phase-space position
    idx = np.searchsorted(population.ids, cols["id"], sorter=np.argsort(population.ids))
    p = np.argsort(population.ids)[idx]
    assert np.all(population.in_range[p])
    tol = dict(rtol=1e-5, atol=1e-4)
    np.testing.assert_allclose(cols["x"], population.pos[p, 0], **tol)
    vel = population.vel[p] - [1.0, 2.0, 3.0]
    np.testing.assert_allclose(cols["vz"], vel[:, 2], **tol)
    v_r, v_theta, v_phi = spherical_velocities(population.pos[p], vel)
    np.testing.assert_allclose(cols["v_phi"], v_phi, **tol)
    np.testing.assert_allclose(cols["r"], np.linalg.norm(population.pos[p], axis=1), **tol)
    for key in ("tracer", "key", "reference", "url", "kind", "seed", "snapshot"):
        assert meta[key] is not None, key
    assert meta["kind"] == "stars" and meta["seed"] == 1 and meta["vcenter"] == [1.0, 2.0, 3.0]
    assert meta["expected_stars"] == pytest.approx(5000)


def test_chunk_size_does_not_change_catalogue(population, tmp_path):
    a = write_catalogue(population, tmp_path / "a.h5", rng=7, chunk_size=10**9)
    b = write_catalogue(population, tmp_path / "b.h5", rng=7, chunk_size=97)
    ca, _ = read_catalogue(tmp_path / "a.h5")
    cb, _ = read_catalogue(tmp_path / "b.h5")
    assert a == b
    for c in ca:
        np.testing.assert_array_equal(ca[c], cb[c])
    mem = population.catalogue(rng=7)
    np.testing.assert_array_equal(mem["id"], ca["id"])


def test_particles_catalogue_and_columns(population, tmp_path):
    n = write_catalogue(population, tmp_path / "p.csv", kind="particles", columns=["id", "r", "n_stars"])
    cols, meta = read_catalogue(tmp_path / "p.csv")
    assert n == population.in_range.sum() and list(cols) == ["id", "r", "n_stars"]
    assert cols["n_stars"].sum() == pytest.approx(5000, rel=1e-6)
    assert meta["kind"] == "particles" and meta["units"]["r"] == "kpc"


def test_catalogue_errors(population, tmp_path):
    path = tmp_path / "c.h5"
    write_catalogue(population, path, rng=0)
    with pytest.raises(FileExistsError):
        write_catalogue(population, path, rng=0)
    write_catalogue(population, path, rng=0, overwrite=True)
    with pytest.raises(ValueError, match="format"):
        write_catalogue(population, tmp_path / "c.fits")
    with pytest.raises(ValueError, match="columns"):
        write_catalogue(population, tmp_path / "d.h5", columns=["x", "colour"])
    with pytest.raises(ValueError, match="columns"):
        write_catalogue(population, tmp_path / "e.h5", kind="stars", columns=["n_stars"])
    with pytest.raises(ValueError, match="kind"):
        write_catalogue(population, tmp_path / "f.h5", kind="galaxies")
