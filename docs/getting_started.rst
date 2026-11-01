***************
Getting started
***************

A typical use reads the dark matter particles of a halo from a snapshot,
chooses a stellar density profile, and computes the weight of each particle:
the stellar mass it carries is its weight times its mass.

A halo
======

`stracers.compute_weights` takes plain arrays: positions [kpc], velocities
[km/s], the potential of each particle [(km/s)²], zero at infinity, and the
particle mass [M☉]. Any snapshot reader can provide them, in any consistent
set of units (see :doc:`units`). Here an isotropic
Hernquist halo of 10\ :sup:`10` M☉ and scale radius 10 kpc, sampled with
`stracers.halos.sample_hernquist`, stands in for a snapshot::

    >>> import numpy as np
    >>> import stracers
    >>> from stracers.halos import sample_hernquist
    >>> n, M = 100_000, 1e10
    >>> pos, vel, pot = sample_hernquist(n, M, 10.0, seed=1)
    >>> mass = M / n

Weights
=======

The stellar halo is a Hernquist profile with a scale radius of 2 kpc and a
stellar mass of 10\ :sup:`8` M☉ inside ``rmax``::

    >>> nu = stracers.get_profile("hernquist", a=2.0)
    >>> res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8, rmax=250.0)
    >>> m_star = res.weights * mass
    >>> print(f"{m_star.sum():.3g}")
    1e+08

``res`` is a `~stracers.WeightsResult`: besides ``res.weights`` it holds every
intermediate quantity of the calculation, for diagnostics.

Checking the stellar profile
============================

The weighted particles follow the chosen profile. The stellar mass in shells
between 1 and 50 kpc is within a few percent of the mass of the profile in
the same shells, computed from its enclosed mass fraction
(`~stracers.mass_fraction_within`) and scaled to 10\ :sup:`8` M☉ inside
``rmax``::

    >>> r = np.linalg.norm(pos, axis=1)
    >>> edges = np.logspace(0, np.log10(50), 6)
    >>> m_shell = np.histogram(r, edges, weights=m_star)[0]
    >>> enclosed = np.array([stracers.mass_fraction_within(nu, x) for x in edges])
    >>> expected = 1e8 * np.diff(enclosed) / stracers.mass_fraction_within(nu, 250.0)
    >>> print(np.round(m_shell / expected, 2))
    [0.98 1.02 1.01 1.02 1.  ]

Following the stars
===================

The weights belong to particle IDs. Save them once, and map them onto the
particles of any later snapshot::

    >>> import os, tempfile
    >>> ids = np.arange(n)
    >>> fname = os.path.join(tempfile.mkdtemp(), "weights.h5")
    >>> stracers.save_weights(fname, ids, res.weights, stellar_mass=1e8)
    >>> ids_w, w, attrs = stracers.load_weights(fname)
    >>> ids_later = np.random.default_rng(0).permutation(ids)   # IDs of a later snapshot
    >>> w_later = stracers.match_weights(ids_w, w, ids_later)
    >>> print(f"{np.sum(w_later * mass):.3g}")
    1e+08

Next steps
==========

- :doc:`weights/index` explains the method, how the weights are computed, and
  how to choose the parameters.
- :doc:`profiles/index` and :doc:`mw/index` list the stellar profiles,
  including the observed Milky Way stellar halo.
- :doc:`nstars/index` counts the stars in a range of distances.
- The :doc:`tutorials` show the whole calculation with plots.
