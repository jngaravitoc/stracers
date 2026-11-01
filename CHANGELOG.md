# Changelog

All notable changes to this project will be documented in this file. Versions follow
[semantic versioning](https://semver.org) (MAJOR.MINOR.PATCH); before 1.0.0, a new minor version may
change the API. Each release is tagged `vX.Y.Z` in git.

## [Unreleased]
### Added
- `CONTRIBUTING.md`: development setup, workflow, tests, documentation and release steps
- `tests/test_npaint.py`: the weights compared with the method of npaint (Errani & Peñarrubia
  2020), an independent re-implementation with analytic distribution functions for its default
  Hernquist halo and tracer; shell masses agree within 2%, velocity dispersions within 0.3%

### Fixed
- The weights page attributed to An & Evans (2006) a minimum inner slope for isotropic tracers
  (half the cusp slope of the potential); that relation is about anisotropy. In a cuspy halo, a
  non-negative isotropic DF near the centre only needs a non-zero inner slope; a core gives f < 0

## [0.3.0] - 2026-11-01
### Added
- `CITATION.cff` with Laporte et al. (2013) and Garavito-Camargo et al. (2019), and a "Citing"
  section in the README
- GitHub workflow running the tests on Python 3.9 with NumPy < 2 and Python 3.12 with NumPy >= 2
- MIT license in the package metadata
- `arxiv` and `doi` columns in the Milky Way profile data file (all 28 papers; DOIs only for the two
  not on arXiv) and `stracers.mw.reference_url`
- `stracers.plotting.plot_mw_profiles`: the Milky Way profiles over their fitted ranges (solid) and
  extrapolated to 200 kpc (dashed), normalised at 30 kpc, with all text in LaTeX math mode;
  `stracers.plotting.LATEX_STYLE` for Computer Modern fonts
- Milky Way page of the documentation: links to the papers in the tables, a list of references
  with arXiv links, and figures for RR Lyrae, BHB stars, and giants and other tracers
- `stracers.kinematics`: velocities and velocity dispersions of the tracers with the weights:
  `spherical_velocities`, `weighted_mean`, `weighted_dispersion`, `bulk_velocity`, `anisotropy` and
  `velocity_profile` (radial profiles of mean velocities, dispersions, anisotropy and effective
  number of particles), tested against the Jeans equation; a user-guide page and the tutorial
  `tutorials/tracer_kinematics.ipynb`
- `stracers.tracers`: mock tracer populations. `paint_tracer(pos, vel, pot, mass, "BHB",
  r_range=(20, 100), n_stars=...)` chooses the most recent literature fit for the tracer covering the
  range (`choose_fit`), computes the weights of the whole stellar halo and selects the range, giving
  the number of stars of the tracer on each particle; `TracerPopulation.velocity_profile` and
  `draw_stars` (mock catalogues). Tracer names such as "RRLyr", "BHB" and "K-giant" are matched
  flexibly (`available_tracers`, `tracer_fits`)
- Following tracer populations through a simulation: `paint_tracer(..., ids=...)` and
  `TracerPopulation.at_snapshot`, which carries the stars of each particle by ID to another snapshot
  and selects the radial range again (stars are conserved, not renormalised); `TracerPopulation.save`
  and `load_population` to paint once and observe snapshots separately; `mw.profile_from_params`
- `compute_weights(..., beyond_rmax="energy")`: particles beyond `rmax` get the weight of their
  energy, so that the stellar halo is stationary; the default in `paint_tracer`
- Mock catalogues, `stracers.catalogue`: `write_catalogue` (HDF5 or CSV, one row per mock star or per
  particle, written in chunks for catalogues of ~1e8 stars, reproducible for any chunk size) and
  `read_catalogue`, with metadata; `TracerPopulation.catalogue` in memory
- Tests on a 1M-particle sample of an isolated Milky Way simulation (`tests/make_sim_sample.py`,
  `tests/test_simulation.py`, skipped without the data) and the tutorial
  `tutorials/simulation_tracers.ipynb`
- Units page in the documentation: consistent units, what assumes kpc, cosmological snapshots, the
  zero point of the potential, and two checks of the inputs

### Changed
- The Milky Way profile table is parsed once and cached: drawing profiles with `mw_profile` in a
  loop is ~40 times faster

### Fixed
- Unbound particles (E >= 0) could carry stars when their energy bin straddled E = 0; they now
  never do, and the energy grid has an edge at E = 0
- `r_trunc` zeroed the whole energy bin containing phi(r_trunc), removing particles below the cut
  on coarse grids; the grid now has an exact edge at phi(r_trunc) (`energy_grid`), and particles
  exactly at the cut belong to the bin below it (`energy_bins`), so the cut is exact
- `double_power_law` (and smoothed broken power laws) overflowed for sharp transitions (large
  `alpha`), dropping to zero far from the break; it is now evaluated in log space (and stays
  finite at r = 0 for a flat inner slope, e.g. Plummer)
- The smoothing spline of the potential failed with 3 or 4 radial bins; its degree is now lowered
  for few bins
- Drawing Milky Way profile parameters within their uncertainties could give non-positive break
  radii (e.g. faccioli2014, 26.5 +/- 8.9 kpc) and NaN profiles; they are now drawn again

## [0.2.1] - 2026-10-31
### Fixed
- With `derivative="spline"` (the default) the fine radial grid ran between the first and last bin
  centres, so the density of states missed the outer half of the last radial bin while N(E) counted
  the particles in it. g(E) was ~20% too low and the weights of the most weakly bound particles
  ~30% too low (self-consistent test halo). The fine grid now covers the radial range of the bins,
  and N(E) only counts particles inside the volume covered by g(E)
- The inner sphere of g(E) uses the inner edge of the first radial bin
- With `derivative="gradient"`, particles in a dropped intermediate radial bin are left out of
  N(E), as their shell is not in g(E)
- The power-law extension of the Eddington integral beyond the grid uses a secant slope over the
  outermost radial bin instead of the derivative at the last grid point, and the contribution of
  the resulting jump in dnu/dpsi at the end of the grid is included

### Changed
- `radial_potential` also returns the edges of the bins kept; `smooth_potential` takes a radial
  range `r_range`; `eddington` takes an optional `outer_slope`; `WeightsResult` has `r_range` and
  `bin_edges` fields

## [0.2.0] - 2026-10-31
The untested `Stracers` package was rewritten as the `stracers` library.

### Breaking changes
- Package renamed `Stracers` to `stracers`, with a new API; `setup.py` replaced by `pyproject.toml`
- The old package, the kinematics code (`tracers_properties.py`) and the notebooks moved, unchanged,
  to `legacy/`; `examples/` renamed `tutorials/`

### Added
- `compute_weights`: Laporte et al. (2013) weights from plain arrays, with a shared energy grid for
  N(E), g(E) and f(E), a spline-smoothed spherically averaged potential, exact integration of the
  Eddington kernel, a power-law extension of the Eddington integral beyond `rmax`, and g(E)
  integrated from r = 0
- Stellar mass as a mass or a fraction of the DM mass, inside `rmax` or in total
  (`stellar_mass`, `stellar_fraction`, `within`); energy truncation of the stellar halo (`r_trunc`)
- Tracer profiles: Plummer, Hernquist, NFW, Einasto, simple, broken and double power laws,
  `truncate` and `mass_fraction_within`
- `stracers.mw`: the 31 Milky Way stellar halo profiles compiled by Medina et al. (2024, Table 6),
  with sampling within their uncertainties
- `stracers.io`: save and load weights, and match them across snapshots by particle ID
- `stracers.nstars`: `number_of_stars` and `number_of_stars_table`, the number of stars in a range of
  galactocentric radius for any profile or Milky Way profile key, with flattening, normalisation by a
  local density or a reference count, and survey-area scaling
- `stracers.halos`: `sample_hernquist`, an isotropic Hernquist halo sampled from its analytic
  distribution function, used by the tests and the documentation
- Tests against an analytic Hernquist halo
- Documentation on the structure of the nba documentation (sphinx-astropy, pydata-sphinx-theme,
  automodapi), with doctested examples and the tutorial run at build time, published with GitHub
  Pages
- The version is defined once, in `stracers.__version__`

### Fixed (relative to the old code)
- g(E) had `dr` inside the square root
- f(E) was reversed twice relative to the energy bins
- Weights were cast with an off-by-one error and normalised to the DM mass instead of the stellar
  mass
- Empty radial bins were given a potential of zero

## [0.1.0] - 2019-02-11
The original `Stracers` package (`setup.py` version 0.1), with the weights of Laporte et al. (2013)
and the kinematics of the weighted particles. It is kept, unchanged, in `legacy/`.
