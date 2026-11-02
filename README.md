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

The tutorials, in [tutorials/](tutorials/), show the library at work. They are also rendered,
with their outputs, in the documentation:

| Tutorial | What it shows |
|---|---|
| [Stellar tracers in a DM halo](https://jngaravitoc.github.io/stracers/tutorials/hernquist_tracers.html) ([notebook](tutorials/hernquist_tracers.ipynb)) | Weights for a Hernquist stellar halo painted on an isotropic Hernquist DM halo, the recovered density profile, and following the tracers to another snapshot by particle ID. Start here. |
| [Kinematics of the stellar tracers](https://jngaravitoc.github.io/stracers/tutorials/tracer_kinematics.html) ([notebook](tutorials/tracer_kinematics.ipynb)) | Weighted mean velocities, velocity dispersions and anisotropy of two stellar haloes, checked against the Jeans equation; rotation, and tracer populations in a radial range. |
| [Tracer populations in a simulation](https://jngaravitoc.github.io/stracers/tutorials/simulation_tracers.html) ([notebook](tutorials/simulation_tracers.ipynb)) | RR Lyrae and BHB populations painted on an isolated Milky Way simulation and followed for ~5 Gyr, stars beyond `rmax`, saving a population, and mock catalogues. Uses non-public data, so it is shown as it was run. |

### Tracer density profiles

The stellar tracer profile ν(r) can be any function of radius. These are built in, as functions or
by name with `get_profile(name, **params)`:

| Name | Function | Parameters | Profile |
|---|---|---|---|
| `"plummer"` | `plummer` | `a` | Plummer, normalised to unit mass |
| `"hernquist"` | `hernquist` | `a` | Hernquist (1990), normalised to unit mass |
| `"nfw"` | `nfw` | `rs`, `rvir` | NFW, normalised to unit mass inside `rvir` |
| `"einasto"` | `einasto` | `n`, `r_eff` | Einasto, normalised to unit mass |
| `"spl"` | `power_law` | `slope`, `r_core` | simple power law |
| `"bpl"` | `broken_power_law` | `r_break`, `inner_slope`, `outer_slope`, `r_core`, `smoothness` | power law with one break, sharp or smoothed |
| `"dpl"` | `double_power_law` | `r_break`, `inner_slope`, `outer_slope`, `alpha`, `r_core` | smooth double power law (α-β-γ) |
| `"tpl"` | `triple_power_law` | `r_break`, `r_break2`, `inner_slope`, `middle_slope`, `outer_slope`, `r_core`, `smoothness` | power law with two breaks, sharp or smoothed |

The power laws are relative densities; `r_core` flattens them at the centre. For example,
`get_profile("dpl", r_break=20, inner_slope=-2.5, outer_slope=-4.5, alpha=2)`. `truncate(nu, r_t)`
tapers any profile beyond `r_t`, and `mw.mw_profile(key)` gives the Milky Way fits below.

### Milky Way stellar halo profiles

`stracers.mw` provides 32 Milky Way stellar halo profiles, stored in
`src/stracers/data/mw_stellar_halo_profiles.csv`: the 31 simple and broken power
laws compiled in Medina et al. (2024, MNRAS 531, 4762, Table 6), and the power law
with two breaks of Li et al. (2026, DESI DR2 K giants; `"li2026"`).

Slopes γ are those of the number density, d ln n / d ln r, and break radii r<sub>b</sub> are in kpc.
SPL: simple power law; BPL: broken power law; TPL: power law with two breaks. Each fit is valid
only over its range.

| Tracer | Key | Model | Parameters | Range [kpc] | Reference |
|---|---|---|---|---|---|
| RRLs | `wetterer1996` | SPL | γ = −3.03 ± 0.08 | 1–80 | [Wetterer & McGraw 1996](https://doi.org/10.1086/118076) |
| RRLs | `vivas2006` | SPL | γ = −2.8 | 4–60 | [Vivas & Zinn 2006](https://arxiv.org/abs/astro-ph/0604359) |
| RRLs | `keller2008` | BPL | r<sub>b</sub> = 45; γ = −2.48 ± 0.09, −5 ± 0.3 | 10–80 | [Keller et al. 2008](https://arxiv.org/abs/0711.2337) |
| RRLs | `watkins2009` | BPL | r<sub>b</sub> = 23; γ = −2.4, −4.5 | 5–100 | [Watkins et al. 2009](https://arxiv.org/abs/0906.0498) |
| RRLs | `akhter2012` | BPL | r<sub>b</sub> = 45; γ = −2.78 ± 0.02, −5 ± 0.2 | 10–80 | [Akhter et al. 2012](https://doi.org/10.1088/0004-637X/756/1/23) |
| RRab | `sesar2013` | SPL | γ = −2.42 ± 0.13 | 5–30 | [Sesar et al. 2013](https://arxiv.org/abs/1305.2160) |
| RRLs | `faccioli2014` | BPL | r<sub>b</sub> = 26.5 ± 8.9; γ = −2.7 ± 0.6, −3.6 ± 0.4 | 9–49 | [Faccioli et al. 2014](https://arxiv.org/abs/1404.6318) |
| RRLs | `zinn2014` | BPL | r<sub>b</sub> = 24; γ = −2.8 ± 0.5, −5.4 ± 0.5 | 5–60 | [Zinn et al. 2014](https://arxiv.org/abs/1312.1602) |
| RRab | `cohen2017` | SPL | γ = −4 | 50–90 | [Cohen et al. 2017](https://arxiv.org/abs/1710.01276) |
| RRLs | `iorio2018` | SPL | γ = −2.96 ± 0.05 | < 28 | [Iorio et al. 2018](https://arxiv.org/abs/1707.03833) |
| RRLs | `medina2018` | SPL | γ = −4.17 (−0.2/+0.18) | 17–145 | [Medina et al. 2018](https://arxiv.org/abs/1802.01581) |
| RRLs | `stringer2021` | BPL | r<sub>b</sub> = 32.1 (−0.9/+1.1); γ = −2.54 ± 0.09, −5.42 (−0.14/+0.13) | 9–100 | [Stringer et al. 2021](https://arxiv.org/abs/2011.13930) |
| RRLs | `chen2023` | BPL | r<sub>b</sub> = 31.11 (−5.88/+7.61); γ = −2.46 (−0.2/+0.18), −3.99 (−1.33/+0.75) | < 80 | [Chen et al. 2023](https://arxiv.org/abs/2307.15385) |
| RRLs | `feng2024` | SPL | γ = −4.09 ± 0.1 | 20–300 | [Feng et al. 2024](https://arxiv.org/abs/2402.05418) |
| RRLs | `medina2024` | BPL | r<sub>b</sub> = 18.1 (−1.1/+2.1); γ = −2.05 (−0.15/+0.13), −4.47 (−0.18/+0.11) | 7–145 | [Medina et al. 2024](https://arxiv.org/abs/2402.14055) |
| BHB stars | `depropris2010` | SPL | γ = −2.5 ± 0.2 | 10–90 | [De Propris, Harrison, & Mares 2010](https://arxiv.org/abs/1007.0013) |
| BHB and BS stars | `deason2011_bpl` | BPL | r<sub>b</sub> = 27 ± 1; γ = −2.3 ± 0.1, −4.6 (−0.1/+0.2) | 10–45 | [Deason, Belokurov, & Evans 2011](https://arxiv.org/abs/1104.3220) |
| BHB stars | `deason2011_spl` | SPL | γ = −2.7 ± 0.5 | < 40 | [Deason, Belokurov, & Evans 2011](https://arxiv.org/abs/1104.3220) |
| BHB stars | `das2016` | BPL | r<sub>b</sub> = 29.87 (−3.55/+2.8); γ = −3.61 (−0.16/+0.15), −4.75 (−0.28/+0.3) | 10–70 | [Das, Williams, & Binney 2016](https://arxiv.org/abs/1608.07297) |
| BHB stars | `thomas2018` | BPL | r<sub>b</sub> = 41.4 (−2.4/+2.5); γ = −4.24 ± 0.08, −3.21 ± 0.07 | 15–220 | [Thomas et al. 2018](https://arxiv.org/abs/1809.08245) |
| BHB stars | `fukushima2019` | BPL | r<sub>b</sub> = 160 (−19/+18); γ = −2.9 ± 0.3, −15 (−3.7/+4.5) | 35–360 | [Fukushima et al. 2019](https://arxiv.org/abs/1904.04966) |
| BHB stars | `yu2024` | SPL | γ = −4.28 (−0.13/+0.12) | 20–70 | [Yu et al. 2024](https://arxiv.org/abs/2402.00104) |
| K giants | `xue2015_bpl` | BPL | r<sub>b</sub> = 18 ± 1; γ = −2.1 ± 0.3, −3.8 ± 0.1 | 10–80 | [Xue et al. 2015](https://arxiv.org/abs/1506.06144) |
| K giants | `xue2015_spl` | SPL | γ = −3.4 ± 0.1 | 10–80 | [Xue et al. 2015](https://arxiv.org/abs/1506.06144) |
| K giants | `lopezcorredoira2024_k` | SPL | γ = −4.6 ± 0.4 | 25–90 | [Lopez-Corredoira et al. 2024](https://arxiv.org/abs/2402.05157) |
| K giants | `li2026` | TPL | r<sub>b</sub> = 16 (−0.2/+0.25), 76.3 (−1.94/+2.37); γ = −1.5 ± 0.07, −3.45 ± 0.02, −5.2 ± 0.14 | 8–200 | [Li et al. 2026](https://arxiv.org/abs/2512.01350) |
| M giants | `lopezcorredoira2024_m` | SPL | γ = −4.5 ± 0.2 | 25–90 | [Lopez-Corredoira et al. 2024](https://arxiv.org/abs/2402.05157) |
| Giant stars | `slater2016` | SPL | γ = −3.5 ± 0.2 | 30–90 | [Slater et al. 2016](https://arxiv.org/abs/1610.02034) |
| MSTO stars | `bell2008` | SPL | γ = −3 | 5–40 | [Bell et al. 2008](https://arxiv.org/abs/0706.0004) |
| near MSTO stars | `sesar2011` | BPL | r<sub>b</sub> = 28; γ = −2.6 ± 0.04, −3.8 ± 0.1 | 12–40 | [Sesar, Juric, & Ivezic 2011](https://arxiv.org/abs/1011.4487) |
| F stars | `piladiez2015` | BPL | r<sub>b</sub> = 20; γ = −2.5 ± 0.4, −4.9 ± 0.4 | 10–60 | [Pila-Diez et al. 2015](https://arxiv.org/abs/1502.02460) |
| A-type stars | `deason2014` | SPL | γ = −6 | 50–100 | [Deason et al. 2014](https://arxiv.org/abs/1403.7205) |

`deason2014` is an upper limit (γ < −6), `vivas2006` an approximate value, and `thomas2018` has an
outer slope shallower than its inner one. `li2026` is a triaxial fit (axis ratios 1 : 0.85 : 0.74);
its radius is the ellipsoidal radius. Ranges "< X" are quoted only with an outer radius.

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
