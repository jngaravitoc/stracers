# stracers

Paint stellar traicers in dark matter particles of an N-body halo,
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

## Documentation

The documentation is in [docs/](docs/) and is published at
<https://jngaravitoc.github.io/stracers/>. It has installation
instructions, a getting-started guide, a tutorial, and a user guide with one
page per module (the method, how the code works, and the API). 


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

See [tutorials/hernquist_tracers.ipynb](tutorials/hernquist_tracers.ipynb).

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

## Kinematics of the tracers

The stars move with the particles that carry them, so their velocities and velocity dispersions are
weighted averages over the particles:

```python
prof = stracers.velocity_profile(pos, vel, res.weights, bins=np.geomspace(1, 100, 11))
prof.sigma_r, prof.sigma_theta, prof.sigma_phi, prof.beta, prof.mean_vphi, prof.n_eff
v_r, v_theta, v_phi = stracers.kinematics.spherical_velocities(pos, vel)
```

See [tutorials/tracer_kinematics.ipynb](tutorials/tracer_kinematics.ipynb).

## Mock tracer populations

Paint a population of a given tracer in a radial range, using the most recent literature profile for
it that covers the range:

```python
pop = stracers.paint_tracer(pos, vel, pot, mass, "RRLyr", r_range=(20, 100), n_stars=3000)
pop.key, pop.n_stars                    # fit used (medina2024), stars of the tracer per particle
prof = pop.velocity_profile(bins=8)     # kinematics of the population
idx = pop.draw_stars(rng=0)             # mock catalogue: particle index of each star
```

Tracers: "RR Lyrae" ("RRLyr"), "BHB", "K giants" ("K-giant"), "M giants", "Giants", "MSTO",
"F stars", "A-type stars".

Paint once and observe later snapshots, and write mock catalogues (HDF5 or CSV, up to ~1e8 stars):

```python
pop = stracers.paint_tracer(pos, vel, pot, mass, "BHB", (20, 100), n_stars=3000, ids=ids)
pop.save("bhb.h5")                                         # paint once ...
later = stracers.load_population("bhb.h5").at_snapshot(ids_t, pos_t, vel_t, center=center_t)
stracers.write_catalogue(later, "bhb_t.h5", rng=0)         # ... observe any snapshot
cols, meta = stracers.read_catalogue("bhb_t.h5")
```

## Counting stars

`stracers.number_of_stars` gives the number of stars in a range of
galactocentric radius for any profile, or any key in `stracers.mw`, with
optional flattening and survey area. For example, this reproduces Medina et
al. (2024, Sect. 5.2.1):

```python
stracers.number_of_stars("medina2024", [100, 230], r_max=300,
                         n_ref=3000, ref_range=(100, 300))      # -> [3000, 356]
stracers.number_of_stars_table(230, r_max=300, n_ref=3000, ref_range=(100, 300))
```

## Units

The weights and kinematics work in any **consistent** units: the potential must
be in the units of the velocities squared, and zero at infinity, and all radii
in one length unit. The docs use kpc, km/s, Msun and (km/s)²; the Milky Way
profiles, star counts and some defaults (`rmax=300`) assume kpc. See the
[Units](https://jngaravitoc.github.io/stracers/units.html) page for
cosmological snapshots, N-body units, and two quick checks of the inputs.
