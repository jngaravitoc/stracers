.. module:: stracers.io

************************************************
Saving weights and catalogues (`stracers.io`)
************************************************

Introduction
============

Weights are computed once, usually on the initial conditions of a simulation,
and belong to particle IDs. `stracers.io` saves them with their IDs and maps
them onto the particles of any other snapshot, so that the stars can be
followed as the simulation evolves. Saving and loading need h5py.

Saving and loading
==================

`~stracers.save_weights` writes the IDs and weights to an HDF5 file, with any
extra keyword arguments as file attributes. Store the parameters used, so
that every file can be traced back to how it was made. `~stracers.load_weights`
reads them back::

    >>> import os, tempfile
    >>> import numpy as np
    >>> import stracers
    >>> ids = np.arange(1, 6)
    >>> weights = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    >>> fname = os.path.join(tempfile.mkdtemp(), "weights.h5")
    >>> stracers.save_weights(fname, ids, weights, stellar_mass=1e8, profile="hernquist", a=2.0)
    >>> ids_w, w, attrs = stracers.load_weights(fname)
    >>> print(ids_w, w, attrs["profile"], attrs["a"])
    [1 2 3 4 5] [0.  0.5 1.  1.5 2. ] hernquist 2.0

Matching weights to another snapshot
====================================

`~stracers.match_weights` returns the weights in the order of the particles of
another snapshot. Particles that have no weight, for example because they
were outside ``rmax`` when the weights were computed and were not saved, get
``fill_value`` (0 by default)::

    >>> ids_later = np.array([5, 3, 7, 1])     # shuffled, with a new particle 7
    >>> stracers.match_weights(ids_w, w, ids_later)
    array([2., 1., 0., 0.])

Mock catalogues
===============

`~stracers.write_catalogue` writes a mock catalogue of a
`~stracers.TracerPopulation` (from `~stracers.paint_tracer`, possibly placed
on a later snapshot) to an HDF5 (``.h5``, ``.hdf5``) or CSV (``.csv``) file, and
`~stracers.read_catalogue` reads it back. There are two kinds of catalogue:

- ``kind="stars"`` (default): one row per mock star, drawn as a Poisson number
  of stars on each particle. The stars of a particle share its position and
  velocity.
- ``kind="particles"``: one row per particle in the radial range, with its
  expected number of stars ``n_stars``. It is smaller, and exact for
  statistics.

The columns are ``id``, ``x``, ``y``, ``z`` [kpc], ``vx``, ``vy``, ``vz`` [km/s],
``r``, ``v_r``, ``v_theta``, ``v_phi`` (and ``n_stars``), relative to the halo
centre of the population, with the velocity ``vcenter`` subtracted::

    >>> import warnings
    >>> from stracers.halos import sample_hernquist
    >>> pos, vel, pot = sample_hernquist(100_000, 1e10, 10.0, seed=1)
    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter("ignore")
    ...     pop = stracers.paint_tracer(pos, vel, pot, 1e5, "RRLyr", (20, 100), n_stars=3000,
    ...                                 ids=np.arange(100_000))
    >>> folder = tempfile.mkdtemp()
    >>> stracers.write_catalogue(pop, os.path.join(folder, "rrl_stars.h5"), rng=1)
    2988
    >>> cols, meta = stracers.read_catalogue(os.path.join(folder, "rrl_stars.h5"))
    >>> list(cols)
    ['id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'r', 'v_r', 'v_theta', 'v_phi']
    >>> meta["kind"], meta["key"], meta["seed"]
    ('stars', 'medina2024', 1)

A particle catalogue with a subset of the columns, as CSV::

    >>> stracers.write_catalogue(pop, os.path.join(folder, "rrl_particles.csv"),
    ...                          kind="particles", columns=["id", "r", "v_r", "n_stars"])
    38085
    >>> cols, meta = stracers.read_catalogue(os.path.join(folder, "rrl_particles.csv"))
    >>> print(f"{cols['n_stars'].sum():.0f}")
    3000

The metadata stores the tracer, fit, reference and link, the ranges, the kind,
the number of rows and expected stars, the seed, the snapshot label, the
centre and ``vcenter``, the units, the version of ``stracers`` and the date.
They are HDF5 attributes, or the first line of a CSV file.

Large catalogues
----------------

- **Memory.** Catalogues are written in chunks of at most ``chunk_size`` rows
  (default 10\ :sup:`7`), even when a single particle carries more stars than
  that, so memory stays at about one chunk whatever the size of the catalogue.
- **IDs.** Particle IDs keep their integer type in both formats (also 64-bit
  IDs above 2\ :sup:`53`); the column types are stored in the metadata.
- **Size.** With all the columns a star takes ~45 bytes in float32, so 10\ :sup:`8`
  stars take ~4.5 GB. Use ``columns=`` to keep only what is needed,
  ``kind="particles"`` for statistics, and HDF5 (CSV warns above 10\ :sup:`7`
  rows). HDF5 ``compression`` (``"gzip"`` by default, ``"lzf"`` or `None`)
  trades speed for size; it works best when particles carry many stars, which
  repeat the particle's columns.
- **Example.** 10\ :sup:`8` RR Lyrae stars drawn from a 10\ :sup:`6`-particle
  sample of a Milky Way halo (~1300 stars per particle in 20–100 kpc) took
  28 s and 4.4 GB without compression, and 38 s and 0.03 GB with gzip, with a
  peak memory of 2 GB in both cases.
- **Reproducibility.** With an integer ``rng`` the catalogue is the same for
  any ``chunk_size``.

API
===

.. automodapi:: stracers.io
   :no-heading:
   :no-main-docstr:
   :headings: =-

Catalogues
----------

.. automodapi:: stracers.catalogue
   :no-heading:
   :no-main-docstr:
   :headings: ^"
