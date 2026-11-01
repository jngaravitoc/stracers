.. module:: stracers.halos

***********************************
Toy halos (`stracers.halos`)
***********************************

Introduction
============

`stracers.halos` samples an isotropic Hernquist (1990) halo in equilibrium,
from its analytic distribution function and with its analytic potential. The
test suite, the examples in this documentation and the :doc:`../tutorials` use
it in place of a simulation snapshot. Because the distribution function is
known, it is also the reference that the weights are validated against.

Sampling a halo
===============

`~stracers.halos.sample_hernquist` returns positions [kpc], velocities [km/s]
and the potential of each particle [(km/s)²]::

    >>> import numpy as np
    >>> from stracers.halos import sample_hernquist, hernquist_potential
    >>> pos, vel, pot = sample_hernquist(10_000, 1e10, 10.0, seed=0)
    >>> pos.shape, vel.shape, pot.shape
    ((10000, 3), (10000, 3), (10000,))
    >>> r = np.linalg.norm(pos, axis=1)
    >>> bool(np.allclose(pot, hernquist_potential(r, 1e10, 10.0)))
    True

Half of the mass of a Hernquist halo is inside
:math:`(1 + \sqrt{2})\, a`::

    >>> print(f"{np.mean(r < (1 + np.sqrt(2)) * 10.0):.2f}")
    0.50

All particles are bound, :math:`E = v^2/2 + \Phi < 0`::

    >>> bool(np.all(0.5 * np.sum(vel**2, axis=1) + pot < 0))
    True

References
==========

- Hernquist, L. 1990, ApJ, 356, 359

API
===

.. automodapi:: stracers.halos
   :no-heading:
   :no-main-docstr:
   :headings: =-
