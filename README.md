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

Built-in profiles, by name in `get_profile(name, **params)`:

- `"plummer"`: Plummer
- `"hernquist"`: Hernquist
- `"nfw"`: NFW
- `"einasto"`: Einasto
- `"spl"`: simple power law
- `"bpl"`: broken power law
- `"dpl"`: double power law
- `"tpl"`: power law with two breaks

Any function ν(r) can also be used. The parameters of each profile are in the
[documentation](https://jngaravitoc.github.io/stracers/profiles/index.html).

### Milky Way stellar halo profiles

`stracers.mw` provides 32 Milky Way stellar halo profiles, stored in
`src/stracers/data/mw_stellar_halo_profiles.csv`: the 31 simple and broken power
laws compiled in Medina et al. (2024, MNRAS 531, 4762, Table 6), and the power law
with two breaks of Li et al. (2026, DESI DR2 K giants; `"li2026"`).

SPL: simple power law; BPL: broken power law; TPL: power law with two breaks. The parameters and
radial range of each fit are on the
[Milky Way page](https://jngaravitoc.github.io/stracers/mw/index.html) of the documentation.

| Tracer | Key | Model | Reference |
|---|---|---|---|
| RRLs | `wetterer1996` | SPL | [Wetterer & McGraw 1996](https://doi.org/10.1086/118076) |
| RRLs | `vivas2006` | SPL | [Vivas & Zinn 2006](https://arxiv.org/abs/astro-ph/0604359) |
| RRLs | `keller2008` | BPL | [Keller et al. 2008](https://arxiv.org/abs/0711.2337) |
| RRLs | `watkins2009` | BPL | [Watkins et al. 2009](https://arxiv.org/abs/0906.0498) |
| RRLs | `akhter2012` | BPL | [Akhter et al. 2012](https://doi.org/10.1088/0004-637X/756/1/23) |
| RRab | `sesar2013` | SPL | [Sesar et al. 2013](https://arxiv.org/abs/1305.2160) |
| RRLs | `faccioli2014` | BPL | [Faccioli et al. 2014](https://arxiv.org/abs/1404.6318) |
| RRLs | `zinn2014` | BPL | [Zinn et al. 2014](https://arxiv.org/abs/1312.1602) |
| RRab | `cohen2017` | SPL | [Cohen et al. 2017](https://arxiv.org/abs/1710.01276) |
| RRLs | `iorio2018` | SPL | [Iorio et al. 2018](https://arxiv.org/abs/1707.03833) |
| RRLs | `medina2018` | SPL | [Medina et al. 2018](https://arxiv.org/abs/1802.01581) |
| RRLs | `stringer2021` | BPL | [Stringer et al. 2021](https://arxiv.org/abs/2011.13930) |
| RRLs | `chen2023` | BPL | [Chen et al. 2023](https://arxiv.org/abs/2307.15385) |
| RRLs | `feng2024` | SPL | [Feng et al. 2024](https://arxiv.org/abs/2402.05418) |
| RRLs | `medina2024` | BPL | [Medina et al. 2024](https://arxiv.org/abs/2402.14055) |
| BHB stars | `depropris2010` | SPL | [De Propris, Harrison, & Mares 2010](https://arxiv.org/abs/1007.0013) |
| BHB and BS stars | `deason2011_bpl` | BPL | [Deason, Belokurov, & Evans 2011](https://arxiv.org/abs/1104.3220) |
| BHB stars | `deason2011_spl` | SPL | [Deason, Belokurov, & Evans 2011](https://arxiv.org/abs/1104.3220) |
| BHB stars | `das2016` | BPL | [Das, Williams, & Binney 2016](https://arxiv.org/abs/1608.07297) |
| BHB stars | `thomas2018` | BPL | [Thomas et al. 2018](https://arxiv.org/abs/1809.08245) |
| BHB stars | `fukushima2019` | BPL | [Fukushima et al. 2019](https://arxiv.org/abs/1904.04966) |
| BHB stars | `yu2024` | SPL | [Yu et al. 2024](https://arxiv.org/abs/2402.00104) |
| K giants | `xue2015_bpl` | BPL | [Xue et al. 2015](https://arxiv.org/abs/1506.06144) |
| K giants | `xue2015_spl` | SPL | [Xue et al. 2015](https://arxiv.org/abs/1506.06144) |
| K giants | `lopezcorredoira2024_k` | SPL | [Lopez-Corredoira et al. 2024](https://arxiv.org/abs/2402.05157) |
| K giants | `li2026` | TPL | [Li et al. 2026](https://arxiv.org/abs/2512.01350) |
| M giants | `lopezcorredoira2024_m` | SPL | [Lopez-Corredoira et al. 2024](https://arxiv.org/abs/2402.05157) |
| Giant stars | `slater2016` | SPL | [Slater et al. 2016](https://arxiv.org/abs/1610.02034) |
| MSTO stars | `bell2008` | SPL | [Bell et al. 2008](https://arxiv.org/abs/0706.0004) |
| near MSTO stars | `sesar2011` | BPL | [Sesar, Juric, & Ivezic 2011](https://arxiv.org/abs/1011.4487) |
| F stars | `piladiez2015` | BPL | [Pila-Diez et al. 2015](https://arxiv.org/abs/1502.02460) |
| A-type stars | `deason2014` | SPL | [Deason et al. 2014](https://arxiv.org/abs/1403.7205) |

`li2026` is a triaxial fit (axis ratios 1 : 0.85 : 0.74); its radius is the ellipsoidal radius.

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
