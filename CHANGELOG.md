# Changelog

All notable changes to this project will be documented in this file. Versions follow
[semantic versioning](https://semver.org) (MAJOR.MINOR.PATCH); before 1.0.0, a new minor version may
change the API. Each release is tagged `vX.Y.Z` in git.

## [Unreleased]

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
