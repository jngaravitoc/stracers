.. module:: stracers.io

***********************************************
Saving and following weights (`stracers.io`)
***********************************************

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

API
===

.. automodapi:: stracers.io
   :no-heading:
   :no-main-docstr:
   :headings: =-
