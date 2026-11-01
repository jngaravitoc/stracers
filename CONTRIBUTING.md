# Contributing to stracers

Bug reports, questions and feature requests are welcome as
[GitHub issues](https://github.com/jngaravitoc/stracers/issues). This page describes how to work on the
code, run the tests and documentation, and make a release.

## Development setup

```bash
git clone https://github.com/jngaravitoc/stracers.git
cd stracers
python -m pip install -e ".[test,io,plot,docs]"
```

`stracers` needs Python 3.9 or later. Its only required dependencies are NumPy and SciPy; h5py
(`io`), matplotlib (`plot`) and the documentation tools (`docs`) are optional.

## Workflow

1. Create a branch from an up-to-date `master` for each change.
2. Add tests for the change. For a bug fix, add a regression test and check that it fails without
   the fix.
3. Add an entry to `CHANGELOG.md` under `## [Unreleased]`, in the `### Added`, `### Changed` or
   `### Fixed` section.
4. Open a pull request against `master`. GitHub runs the tests (Python 3.9 with NumPy < 2, and
   Python 3.12 with NumPy >= 2) and builds the documentation, and Copilot reviews the changes.
   Address the review comments in the same branch.
5. Merge once the checks pass. Only release commits go to `master` directly.

## Tests

```bash
pytest
```

The tests use toy haloes sampled from analytic distribution functions (`stracers.halos`), so they need
no data. The tests in `tests/test_simulation.py` and the tutorial `tutorials/simulation_tracers.ipynb`
use a 10^6-particle sample of an isolated Milky Way simulation, which is not public; they are skipped
without it. To make the sample from the snapshots:

```bash
python tests/make_sim_sample.py [snapshot_dir] [output_file]
```

It is written to `$STRACERS_SIM_SAMPLE` (default `~/.cache/stracers/mw_iso_b1_1M.hdf5`), reading the
snapshots in chunks (about 3 minutes and 600 MB of memory for six snapshots of 10^8 particles).

## Documentation

The documentation (`docs/`) follows the structure of the nba documentation: Sphinx with
sphinx-astropy, pydata-sphinx-theme and automodapi, and myst-nb for the tutorials. The examples in the
pages are doctests, and the build must have no warnings:

```bash
pytest --doctest-plus --doctest-rst docs/          # the examples in the pages
make -C docs html SPHINXOPTS="-W --keep-going"     # docs/_build/html/index.html
```

The toy-halo tutorials are run when the documentation is built. The simulation tutorial needs local
data: run it yourself (`jupyter nbconvert --to notebook --execute --inplace`) and commit it with its
outputs; it is not executed by the build. The documentation is published at
<https://jngaravitoc.github.io/stracers/> by the Documentation workflow on every push to `master`.

## Conventions

- **Docstrings** follow the NumPy format; new public functions go in the API section of their
  user-guide page.
- **Units**: the code works in any consistent units, and the documentation uses kpc, km/s, Msun and
  (km/s)^2 (see the Units page). Anything that assumes particular units must say so.
- **Literature data** (`src/stracers/data/mw_stellar_halo_profiles.csv`): check every value and
  reference against the paper, and give its arXiv identifier (or its DOI if it is not on arXiv).

## Releases

The version is defined once, in `src/stracers/__init__.py` (`__version__`), and follows
[semantic versioning](https://semver.org): before 1.0.0, a new minor version may change the API. To
make a release, once everything for it is merged and the checks pass on `master`:

1. In `CHANGELOG.md`, move the entries under `## [Unreleased]` to a new `## [X.Y.Z] - YYYY-MM-DD`
   section, and leave an empty `## [Unreleased]` above it.
2. Set `__version__ = "X.Y.Z"` in `src/stracers/__init__.py`, and `version` and `date-released` in
   `CITATION.cff`.
3. Run `pytest`: `tests/test_version.py` checks that the installed version, `__version__` and the
   changelog agree (reinstall with `pip install -e .` first).
4. Commit on `master` ("Release X.Y.Z"), tag and push:

   ```bash
   git tag -a vX.Y.Z -m "stracers X.Y.Z"
   git push origin master vX.Y.Z
   ```

5. Check that the Tests and Documentation workflows pass for the release commit and that the
   documentation shows the new version.

## Citing

If you use `stracers`, please cite the papers in [CITATION.cff](CITATION.cff); see the README.
