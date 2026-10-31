# stracers

Assign stellar-tracer weights to the dark matter particles of an N-body halo,
following the particle-tagging method of Laporte et al. (2013).

Each DM particle gets a weight that depends only on its energy,

    w(E) = f_tracer(E) / f_DM(E) = f_tracer(E) g(E) / N(E)

- f_tracer is the Eddington distribution function of a tracer population with
  density ν(r) in the spherically averaged halo potential.
- g(E) is the density of states.
- N(E) is the differential energy distribution of the DM particles.

Weighting the DM particles by w makes them trace a stellar halo with profile
ν(r). Because weights belong to particle IDs, the tracers can be followed
through later snapshots.

## Install

```bash
pip install -e ".[io,plot]"      # or ".[test]" to run the tests
```

## Usage

The library works on plain numpy arrays, so any snapshot format can be used.

```python
import numpy as np
import stracers

# pos [kpc], vel [km/s], pot [(km/s)^2], mass [Msun], ids: read from your snapshot
nu = stracers.get_profile("hernquist", a=2.0)     # or any callable nu(r)
res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8, rmax=300)

stellar_mass_per_particle = res.weights * mass      # sums to 1e8

# or: 1% of the DM mass inside rmax, or a total (to infinity) stellar mass
res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_fraction=0.01)
res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e9, within="total")

# truncated stellar halo: energy truncation at 50 kpc, or a smoothly tapered profile
res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8, r_trunc=50.0)
res = stracers.compute_weights(pos, vel, pot, mass, stracers.truncate(nu, 50.0), stellar_mass=1e8)

stracers.save_weights("weights.h5", ids, res.weights, stellar_mass=1e8)
ids_w, w, attrs = stracers.load_weights("weights.h5")
w_later = stracers.match_weights(ids_w, w, ids_later_snapshot)
```

`res` also contains the diagnostics: ψ(r), ν, dν/dψ, d²ν/dψ², N(E), g(E), f(E)
and the binned weights. `stracers.plotting.plot_diagnostics(res)` plots them.

See [examples/hernquist_tracers.ipynb](examples/hernquist_tracers.ipynb).

Available profiles: `plummer(a)`, `hernquist(a)`, `nfw(rs, rvir)`, and
`einasto(n, r_eff)` (normalised to unit mass), and power laws with your own
slopes and break radius: `power_law` / `"spl"`, `broken_power_law` / `"bpl"`
and the smooth double power law `double_power_law` / `"dpl"`, all with an
optional core radius. For example,
`get_profile("dpl", r_break=20, inner_slope=-2.5, outer_slope=-4.5, alpha=2)`.
`truncate(nu, r_t)` tapers any profile beyond `r_t`.

`stracers.mw` provides the 31 Milky Way stellar halo profiles (simple and
broken power laws) compiled in Medina et al. (2024, MNRAS 531, 4762, Table 6),
stored in `src/stracers/data/mw_stellar_halo_profiles.csv`:

```python
from stracers import mw
nu = mw.mw_profile("medina2024")                      # cored (1 kpc), smoothed break
nu = mw.mw_profile("medina2024", rng=np.random.default_rng(0))   # parameters drawn within errors
```

## Units

kpc, km/s, Msun, potential in (km/s)², G = 4.30071e-6 kpc (km/s)² / Msun
(Gadget conventions). The potential should be zero at infinity.

## Notes on the method

- The potential is averaged in logarithmic radial shells (`dlogr`). Shells with
  fewer than `min_per_bin` particles are dropped. By default
  (`derivative="spline"`) a smoothing spline is fitted to Φ(ln r) before ν is
  differentiated with respect to ψ. `derivative="gradient"` reproduces the
  original notebook, but its derivatives are much noisier.
- The Eddington integral is evaluated with exact integration of the
  1/√(ε−ψ) kernel on each ψ segment.
- With `boundary_term=True` (the default), the part of the integral beyond
  `rmax` is included by extrapolating ν(ψ) as a power law. Without it the
  tracer profile is underestimated near `rmax`.
- Bins where f(E) < 0 get zero weight, and the code warns. This is physical
  when the tracer is inconsistent with an isotropic DF in the halo potential.
  For example, a cored tracer (Plummer) in a ρ ∝ r⁻¹ cusp: An & Evans (2006)
  showed that an isotropic tracer needs an inner slope of at least half the
  potential's cusp slope.
- The method assumes equal-mass particles and a halo close to spherical in
  energy space (the potential is spherically averaged, as in Laporte et al.
  2013).

## Documentation

The full documentation is in [docs/](docs/). It covers the method (Laporte+13),
how the code works step by step, a usage guide and the API reference. To build
it:

```bash
pip install -e ".[docs]"
make -C docs html        # open docs/_build/html/index.html
```

## Tests

```bash
pytest
```

The tests sample an isotropic Hernquist halo from its analytic distribution
function. They check that the Eddington inversion recovers the analytic DF,
that self-consistent weights are constant, and that the Hernquist, Plummer and
Einasto tracer profiles are recovered to within 10%.

## Legacy code

The original package and notebooks, including `tracers_properties.py` with the
kinematics analysis, are kept unchanged in [legacy/](legacy/) for reference.
