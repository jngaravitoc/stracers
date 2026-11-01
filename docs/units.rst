*****
Units
*****

``stracers`` does not convert units: it works on the arrays it is given. Most
of it does not depend on the units chosen, as long as they are **consistent**.
This page lists what that means, which parts assume particular units, and how
to check the inputs from a simulation.

The rule
========

=====================  ===================================================  ===========================
Input                  Must be                                              Default units in the docs
=====================  ===================================================  ===========================
``pos``, radii         one length unit, also used for profile radii,        kpc
                       ``rmax``, ``r_trunc``, ``r_core``, ``bins``, …
``vel``                one velocity unit                                    km/s
``pot``                potential per unit mass, in **velocity² units**      (km/s)²
                       (the units of ``vel`` squared), zero at infinity
``mass``               one mass unit, also used for ``stellar_mass``        M☉
=====================  ===================================================  ===========================

`~stracers.compute_weights` only uses :math:`E = \tfrac{1}{2}v^2 + \Phi` and
radii, so the weights are the same in any consistent set of units. For
example, a code using N-body units with :math:`G = 1` can pass its outputs
directly, with ``rmax`` and the profile radii in its length unit. The
kinematics (`stracers.kinematics`) are also unit-agnostic: velocities and
dispersions come out in the units of ``vel``.

What assumes particular units
=============================

- **Milky Way profiles** (`stracers.mw`): break radii, fitted ranges and the
  default ``r_core=1`` are in kpc. The positions must be in kpc to use them.
- **Star counts** (`stracers.nstars`): radii in kpc, densities in
  stars / kpc\ :sup:`3`, survey areas in deg².
- **Defaults in kpc**: ``rmax=300`` in `~stracers.compute_weights` and the radii
  in the examples.
- **The gravitational constant**: ``stracers.weights.G_GADGET`` =
  4.30071 × 10\ :sup:`-6` kpc (km/s)² / M☉. It is used only by the optional
  ``truncation_correction`` (pass ``G`` in your units) and by the toy halo of
  `stracers.halos`.

Two common mistakes
===================

**Mixed units.** Converting the positions to kpc, but not the potential, or
the velocities but not the potential, gives wrong weights without any error.
Convert all the quantities, or none.

**Cosmological snapshots.** Codes such as Gadget write comoving quantities
with factors of the Hubble parameter :math:`h` and the scale factor
:math:`a`. In the standard Gadget conventions, the physical position is
:math:`a\, x / h` (with :math:`x` in kpc/h), the peculiar velocity is
:math:`\sqrt{a}\, u` (with :math:`u` the velocity in the file), and the mass is
:math:`10^{10} m / h` M☉. The potential has its own convention: check the
documentation of the code and convert it to physical (km/s)² as well.

The zero point of the potential
===============================

The extension of the Eddington integral beyond ``rmax`` (``boundary_term``)
assumes :math:`\Phi \to 0` at infinity, as for an isolated halo. Potentials
from periodic boxes, or written with an offset, need to be shifted. In a
cosmological simulation the potential of the halo particles alone (e.g.
recomputed with a tree code) is the cleanest choice.

Checking the inputs
===================

Two quick checks catch most unit and zero-point problems. Here they are on
the toy halo of `stracers.halos`, which is in kpc, km/s, M☉ and (km/s)²::

    >>> import numpy as np
    >>> from stracers.halos import sample_hernquist
    >>> from stracers.weights import G_GADGET as G      # kpc (km/s)^2 / Msun
    >>> pos, vel, pot = sample_hernquist(100_000, 1e10, 10.0, seed=1)
    >>> mass = np.full(100_000, 1e5)
    >>> r = np.linalg.norm(pos, axis=1)
    >>> rmax = 300.0

**1. Most particles inside rmax are bound.** In a halo in equilibrium,
nearly all the particles inside the virial radius have :math:`E < 0`. A small
bound fraction means the velocities are too large for the potential (wrong
units), or the zero point is too high::

    >>> E = 0.5 * np.sum(vel**2, axis=1) + pot
    >>> print(f"{np.mean(E[r < rmax] < 0):.2f}")
    1.00

**2. The potential at the outskirts is about** :math:`-G M(<r)/r`. With
:math:`G` in the units of the inputs, the ratio is close to 1 (a little above 1
when there is mass beyond ``rmax``). A ratio far from 1 means that the
potential, positions or masses are not in the units assumed for :math:`G`, or
that the zero point is shifted::

    >>> shell = (r > 0.8 * rmax) & (r < rmax)
    >>> m_in = mass[r < 0.8 * rmax].sum()
    >>> print(f"{np.median(pot[shell]) / (-G * m_in / np.median(r[shell])):.2f}")
    1.05

For the same halo with deliberately wrong inputs the checks give:

=======================================  ==============  ======================
Input                                    Bound fraction  :math:`\Phi / (-GM/r)`
=======================================  ==============  ======================
correct                                  1.00            1.05
potential in (m/s)², velocities in km/s  1.00            10\ :sup:`6`
velocities in m/s                        0.00            1.05
positions in Mpc                         1.00            0.00
zero point raised by 500 (km/s)²         0.72            −2.3
zero point lowered by 2000 (km/s)²       1.00            14.5
=======================================  ==============  ======================

Each check misses some problems, so use both.
