import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from stracers import compute_weights, mw, paint_tracer
from stracers.kinematics import velocity_profile
from stracers.tracers import available_tracers, choose_fit, tracer_fits


def paint(halo, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return paint_tracer(halo["pos"], halo["vel"], halo["pot"], halo["mass"], *args, **kwargs)


def test_tracer_names():
    rrl = tracer_fits("RRLyr")
    assert rrl == tracer_fits("rr lyrae") == tracer_fits("RRLs")
    assert {mw.get_entry(k)["tracer"] for k in rrl} == {"RRLs", "RRab"} and len(rrl) == 15
    assert tracer_fits("K-giant") == tracer_fits("k giants") == tracer_fits("KGiants")
    assert set(tracer_fits("K-giant")) == {"xue2015_spl", "xue2015_bpl", "lopezcorredoira2024_k"}
    assert len(tracer_fits("BHB")) == 7
    assert set(available_tracers()) >= {"RR Lyrae", "BHB", "K giants", "M giants", "MSTO"}
    with pytest.raises(ValueError, match="Unknown tracer"):
        tracer_fits("Cepheids")


def test_choose_fit():
    # most recent fit covering the range; ties go to broken power laws
    assert choose_fit("RRLyr", (20, 100)) == "medina2024"   # feng2024 is also 2024 but an SPL
    assert choose_fit("BHB", (20, 100)) == "thomas2018"     # the only BHB fit covering 20-100 kpc
    assert choose_fit("RRLyr", (150, 250)) == "feng2024"    # medina2024 stops at 145 kpc
    # no K-giant fit reaches 100 kpc: the largest overlap is chosen, with a warning
    with pytest.warns(RuntimeWarning, match="No K giants fit covers"):
        assert choose_fit("K-giant", (20, 100)) == "xue2015_bpl"
    with pytest.raises(ValueError):
        choose_fit("BHB", (100, 20))


def test_paint_tracer_normalisation(hernquist_halo):
    halo = hernquist_halo
    pop = paint(halo, "RRLyr", (20, 100), n_stars=3000)
    r = np.linalg.norm(halo["pos"], axis=1)
    assert pop.key == "medina2024" and pop.tracer == "RR Lyrae"
    assert pop.url == "https://arxiv.org/abs/2402.14055"
    assert pop.total == pytest.approx(3000)
    assert np.array_equal(pop.in_range, (r >= 20) & (r <= 100))
    assert np.all(pop.n_stars[~pop.in_range] == 0)

    # n_ref: the scale (stars per unit fraction of the population) is set by the
    # stars in ref_range, outside r_range here, so they are not in n_stars
    ref = paint(halo, "RRLyr", (20, 100), n_ref=500, ref_range=(100, 200))
    frac = ref.result.weights * halo["mass"]
    scale = ref.n_stars[ref.in_range] / frac[ref.in_range]
    np.testing.assert_allclose(scale, scale[0])
    assert scale[0] * np.sum(frac[(r >= 100) & (r <= 200)]) == pytest.approx(500)

    unnorm = paint(halo, "RRLyr", (20, 100))
    assert unnorm.total == pytest.approx(np.sum(unnorm.result.weights[unnorm.in_range] * halo["mass"]))


def test_paint_tracer_uses_whole_halo_weights(hernquist_halo):
    """The weights are those of the whole stellar halo, and the range is only a
    selection: the stars follow the profile across the range."""
    halo = hernquist_halo
    pop = paint(halo, "BHB", (20, 100), n_stars=1000)
    nu = mw.mw_profile(pop.key)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        w = compute_weights(halo["pos"], halo["vel"], halo["pot"], halo["mass"], nu).weights
    ratio = pop.n_stars[pop.in_range] / w[pop.in_range]
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-10)

    r = np.linalg.norm(halo["pos"], axis=1)
    edges = np.geomspace(20, 100, 5)
    got = np.histogram(r, edges, weights=pop.n_stars)[0]
    shell = np.array([quad(lambda x: 4 * np.pi * x**2 * nu(x), lo, hi)[0]
                      for lo, hi in zip(edges[:-1], edges[1:])])
    np.testing.assert_allclose(got / (shell * 1000 / shell.sum()), 1.0, atol=0.05)


def test_population_kinematics_and_catalogue(hernquist_halo):
    halo = hernquist_halo
    pop = paint(halo, "RRLyr", (20, 100), n_stars=3000)
    prof = pop.velocity_profile(bins=4)
    sel = pop.in_range
    ref = velocity_profile(halo["pos"][sel], halo["vel"][sel], pop.n_stars[sel], bins=4,
                           rmin=20, rmax=100)
    np.testing.assert_allclose(prof.sigma_r, ref.sigma_r)
    assert prof.edges[0] == pytest.approx(20) and prof.edges[-1] == pytest.approx(100)

    idx = pop.draw_stars(rng=1)
    assert np.array_equal(idx, pop.draw_stars(rng=1))
    assert abs(len(idx) - 3000) < 5 * np.sqrt(3000)
    assert np.all(pop.in_range[idx])


def test_paint_tracer_options_and_errors(hernquist_halo):
    halo = hernquist_halo
    pop = paint(halo, "RRLyr", (20, 100), key="feng2024", n_stars=10)
    assert pop.key == "feng2024"
    with pytest.raises(ValueError, match="is not a"):
        paint(halo, "BHB", (20, 100), key="medina2024")
    with pytest.raises(ValueError, match="not both"):
        paint(halo, "BHB", (20, 100), n_stars=10, n_ref=10, ref_range=(10, 20))
    with pytest.raises(ValueError, match="inside rmax"):
        paint(halo, "BHB", (20, 100), rmax=50.0)
    with pytest.raises(ValueError, match="needs ref_range"):
        paint(halo, "BHB", (20, 100), n_ref=10)
