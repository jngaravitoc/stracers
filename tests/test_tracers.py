import pathlib
import re
import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from stracers import compute_weights, load_population, mw, paint_tracer
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
    assert set(tracer_fits("K-giant")) == {"xue2015_spl", "xue2015_bpl", "lopezcorredoira2024_k",
                                            "li2025"}
    assert len(tracer_fits("BHB")) == 7
    assert set(available_tracers()) >= {"RR Lyrae", "BHB", "K giants", "M giants", "MSTO"}
    with pytest.raises(ValueError, match="Unknown tracer"):
        tracer_fits("Cepheids")


def test_choose_fit():
    # most recent fit covering the range; ties go to broken power laws
    assert choose_fit("RRLyr", (20, 100)) == "medina2024"   # feng2024 is also 2024 but an SPL
    assert choose_fit("BHB", (20, 100)) == "thomas2018"     # the only BHB fit covering 20-100 kpc
    assert choose_fit("RRLyr", (150, 250)) == "feng2024"    # medina2024 stops at 145 kpc
    assert choose_fit("K-giant", (20, 100)) == "li2025"     # two breaks, 8-200 kpc
    # no M-giant fit reaches 100 kpc: the largest overlap is chosen, with a warning
    with pytest.warns(RuntimeWarning, match="No M giants fit covers"):
        assert choose_fit("M giants", (20, 100)) == "lopezcorredoira2024_m"
    with pytest.warns(RuntimeWarning, match="No K giants fit covers"):
        assert choose_fit("K-giant", (20, 300)) == "li2025"
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


def test_documented_spellings_resolve():
    """Every spelling listed in the "Tracer names" table of the user guide is
    accepted (regression: "bhb stars" was documented but rejected)."""
    page = (pathlib.Path(__file__).parents[1] / "docs" / "tracers" / "index.rst").read_text()
    table = page.split("Tracer names")[1].split("Choosing the fit")[0]
    spellings = re.findall(r'``"([^"]+)"``', table)
    assert "bhb stars" in spellings and len(spellings) >= 15
    for name in spellings:
        assert tracer_fits(name), name


def test_range_validated_with_explicit_key(hernquist_halo):
    """A key given by the user must not bypass the radial-range checks."""
    halo = hernquist_halo
    with pytest.raises(ValueError, match="0 < r1 < r2"):
        paint(halo, "BHB", (100, 20), key="thomas2018", n_stars=10)
    with pytest.raises(ValueError, match="0 < r1 < r2"):
        paint(halo, "BHB", (0, 20), key="thomas2018")
    with pytest.raises(ValueError, match="ref_range"):
        paint(halo, "BHB", (20, 100), key="thomas2018", n_ref=10, ref_range=(50, 30))


def test_normalisation_needs_tracer_stars(hernquist_halo):
    """Normalising to a range without weighted particles raises instead of
    giving NaN or infinite star counts."""
    halo = hernquist_halo
    r = np.linalg.norm(halo["pos"], axis=1)
    keep = (r < 40) | (r > 50)                               # an empty shell
    sub = {k: (v[keep] if np.ndim(v) else v) for k, v in halo.items() if k != "ids"}
    with pytest.raises(ValueError, match="No tracer stars in r_range"):
        paint(sub, "RRLyr", (41, 49), n_stars=100)
    with pytest.raises(ValueError, match="No tracer stars in ref_range"):
        paint(sub, "RRLyr", (20, 100), n_ref=100, ref_range=(41, 49))
    # without normalisation an empty range is allowed: it just has no stars
    pop = paint(sub, "RRLyr", (41, 49))
    assert pop.total == 0 and len(pop.draw_stars(rng=0)) == 0


def test_at_snapshot_follows_ids(hernquist_halo):
    """The stars follow their particles by ID; they are conserved, not
    renormalised, and the radial selection is applied again."""
    halo = hernquist_halo
    ids = halo["ids"]
    pop = paint(halo, "RRLyr", (20, 100), n_stars=3000, ids=ids, snapshot="000")
    assert pop.matched_fraction == 1.0 and pop.snapshot == "000"

    # a "later snapshot": the same particles in another order, moved
    rng = np.random.default_rng(3)
    order = rng.permutation(len(ids))
    shift = np.array([5.0, -3.0, 1.0])
    pos2 = halo["pos"][order] * 1.1 + shift            # all radii grow by 10%
    later = pop.at_snapshot(ids[order], pos2, halo["vel"][order], center=shift, snapshot="010")

    np.testing.assert_allclose(later.n_stars_all, pop.n_stars_all[order])
    assert np.sum(later.n_stars_all) == pytest.approx(np.sum(pop.n_stars_all))   # conserved
    r2 = np.linalg.norm(pos2 - shift, axis=1)
    assert np.array_equal(later.in_range, (r2 >= 20) & (r2 <= 100))
    assert later.total != pytest.approx(3000)          # not renormalised
    assert later.snapshot == "010" and later.matched_fraction == pytest.approx(1.0)

    # missing particles lose their stars, with a warning
    keep = np.arange(len(ids)) % 4 != 0
    with pytest.warns(RuntimeWarning, match="of the stars"):
        part = pop.at_snapshot(ids[keep], halo["pos"][keep], halo["vel"][keep])
    assert part.matched_fraction == pytest.approx(np.sum(pop.n_stars_all[keep]) / np.sum(pop.n_stars_all))

    nope = paint(halo, "RRLyr", (20, 100), n_stars=3000)
    with pytest.raises(ValueError, match="Pass ids"):
        nope.at_snapshot(ids, halo["pos"], halo["vel"])


def test_save_and_load_population(hernquist_halo, tmp_path):
    halo = hernquist_halo
    pop = paint(halo, "BHB", (20, 100), n_stars=1000, ids=halo["ids"], snapshot="000",
                profile_kwargs={"rng": np.random.default_rng(5)})
    pop.save(tmp_path / "bhb.h5")
    loaded = load_population(tmp_path / "bhb.h5")
    for name in ("tracer", "key", "reference", "url", "r_range", "fit_range", "scale", "snapshot"):
        assert getattr(loaded, name) == getattr(pop, name), name
    assert loaded.profile.params == pytest.approx(pop.profile.params)   # drawn parameters kept
    r = np.geomspace(1, 200, 9)
    np.testing.assert_allclose(loaded.profile(r), pop.profile(r))
    with pytest.raises(ValueError, match="at_snapshot"):
        loaded.total

    a = pop.at_snapshot(halo["ids"], halo["pos"], halo["vel"])
    b = loaded.at_snapshot(halo["ids"], halo["pos"], halo["vel"])
    np.testing.assert_allclose(b.n_stars, a.n_stars)
    assert b.total == pytest.approx(1000)

    with pytest.raises(ValueError, match="Pass ids"):
        paint(halo, "BHB", (20, 100)).save(tmp_path / "x.h5")


def test_duplicate_ids_rejected(hernquist_halo):
    """Regression: duplicated IDs in a later snapshot copied a particle's stars
    onto each of them (matched_fraction > 1)."""
    halo = hernquist_halo
    pop = paint(halo, "RRLyr", (20, 100), n_stars=1000, ids=halo["ids"])
    dup = halo["ids"].copy()
    dup[1] = dup[0]
    with pytest.raises(ValueError, match="duplicated"):
        pop.at_snapshot(dup, halo["pos"], halo["vel"])
    with pytest.raises(ValueError, match="duplicated"):
        paint(halo, "RRLyr", (20, 100), ids=dup)
