# stracers

[![Tests](https://github.com/jngaravitoc/stracers/actions/workflows/tests.yml/badge.svg)](https://github.com/jngaravitoc/stracers/actions/workflows/tests.yml)
[![Documentation](https://github.com/jngaravitoc/stracers/actions/workflows/docs.yml/badge.svg)](https://jngaravitoc.github.io/stracers/)

Paint stellar tracers in dark matter particles of an N-body halo,
following the particle-tagging method of Laporte et al. (2013) that generalizes the method in Bullock & Johnston (2005). 

In short, the code weights the DM particles by w makes them trace a stellar halo with profile
ν(r). Because weights belong to particle IDs, the tracers can be followed
through later snapshots.

## Install

`stracers` needs Python 3.9 or later. It is not on PyPI yet; install it from
GitHub:

```bash
git clone https://github.com/jngaravitoc/stracers.git
cd stracers
python -m pip install .           # or: python -m pip install -e .  (editable, to work on the code)
```

The only required dependencies are NumPy (>= 1.20) and SciPy (>= 1.7), which
pip installs with the package. The code works on plain NumPy arrays, so it
does not depend on any snapshot reader. Optional features need extra
packages, installed with extras:

| Extra   | Installs                       | Needed for                                                                 |
|---------|--------------------------------|----------------------------------------------------------------------------|
| `io`    | h5py                           | saving and loading weights and tracer populations, HDF5 mock catalogues    |
| `plot`  | matplotlib                     | `stracers.plotting` (diagnostic plots, Milky Way profiles)                 |
| `test`  | pytest, h5py                   | running the tests                                                          |
| `docs`  | Sphinx and its extensions      | building the documentation                                                 |

For example:

```bash
python -m pip install ".[io,plot]"               # everything for using the library
python -m pip install -e ".[io,plot,test,docs]"  # development
```

CSV catalogues and everything else work without the extras.

## Documentation

The documentation is in [docs/](docs/) and is published at
<https://jngaravitoc.github.io/stracers/>. It has installation
instructions, a getting-started guide, a tutorial, and a user guide with one
page per module (the method, how the code works, and the API). 


## Usage

See [tutorials/hernquist_tracers.ipynb](tutorials/hernquist_tracers.ipynb).


Available profiles: `plummer(a)`, `hernquist(a)`, `nfw(rs, rvir)`, and
`einasto(n, r_eff)` (normalised to unit mass), and power laws with your own
slopes and break radius: `power_law` / `"spl"`, `broken_power_law` / `"bpl"`
the smooth double power law `double_power_law` / `"dpl"`, and the power law with
two breaks `triple_power_law` / `"tpl"`, all with an
optional core radius. For example,
`get_profile("dpl", r_break=20, inner_slope=-2.5, outer_slope=-4.5, alpha=2)`.
`truncate(nu, r_t)` tapers any profile beyond `r_t`.

`stracers.mw` provides 32 Milky Way stellar halo profiles, stored in
`src/stracers/data/mw_stellar_halo_profiles.csv`: the 31 simple and broken power
laws compiled in Medina et al. (2024, MNRAS 531, 4762, Table 6), and the power law
with two breaks of Li et al. (2026, DESI DR2 K giants; `"li2026"`).

Tracers: "RR Lyrae" ("RRLyr"), "BHB", "K giants" ("K-giant"), "M giants", "Giants", "MSTO",
"F stars", "A-type stars".


## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, tests, documentation and release
steps.

## Citing

If you use `stracers`, please cite Laporte et al. (2013, MNRAS, 435, 901;
[arXiv:1301.5319](https://arxiv.org/abs/1301.5319)), who introduced the
particle-tagging method, and Garavito-Camargo et al. (2019, ApJ, 884, 51;
[arXiv:1902.05089](https://arxiv.org/abs/1902.05089)). The metadata are in
[CITATION.cff](CITATION.cff), and GitHub's "Cite this repository" button gives
them as BibTeX or APA. If you use the Milky Way profiles, please also cite
Medina et al. (2024) and the papers of the fits you use.
