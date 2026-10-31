The method
==========

This page describes the weighting scheme of Laporte et al. (2013, Sect. 2;
the source is in ``doc/method.tex``). It generalises the scheme of Bullock &
Johnston (2005) to DM haloes that are not spherical.

Why the weights depend only on energy
-------------------------------------

Dark matter, like stars, is collisionless. Its distribution function (DF) obeys
the collisionless Boltzmann equation. If the halo is in a steady state, Jeans'
theorem says the DF can be written as a function of the isolating integrals of
motion. In any potential, including a triaxial one, a DF of the form
:math:`f = f(H)` is allowed, where :math:`H = \tfrac{1}{2}v^2 + \Phi` is the
energy per unit mass.

Real haloes are never exactly in a steady state, because accretion keeps
reshaping them. Between mergers, however, their virialised regions are in
*quasi-equilibrium* (Natarajan et al. 1997). So, as long as we stay inside the
virialised region, we can give the stars a DF of the form
:math:`f_* = f_*(E)`. Stars painted on this way are in equilibrium in the halo
and keep their profile when the simulation is evolved.

The weight function
-------------------

Each particle with energy :math:`E` carries both DM and stars, through the
weight

.. math::
   :label: weight

   w(E) = \frac{N_*(E)}{N(E)} = \frac{f_*(E)\, g(E)}{N(E)},

where

- :math:`N(E)\,dE` is the *differential energy distribution* of the DM: the
  mass of DM particles with energy between :math:`E` and :math:`E + dE`. It is
  measured directly from the simulation.
- :math:`g(E)` is the *density of states*: the phase-space volume per unit
  energy.
- :math:`f_*(E)` is the stellar DF, and :math:`N_*(E) = f_*(E) g(E)` is the
  differential energy distribution of the stars.

Since :math:`N(E) = f_{DM}(E) g(E)`, the weight is just the ratio of the two
DFs, :math:`w = f_*/f_{DM}`. The stellar mass carried by particle :math:`i` is

.. math::

   m_{*,i} = w(E_i)\, m_{DM}.

The stellar density profile
---------------------------

The stellar tracer density :math:`\nu(r) = \rho_*(r)/M_*` can be any profile
that has a non-negative DF in the halo potential (see `Limits`_).
Laporte et al. (2013) use a Hernquist (1990) profile,

.. math::

   \nu(r) = \frac{a}{2\pi\, r\, (r + a)^3},

where the scale radius is related to the 3D half-mass radius by
:math:`a = r_e / (\sqrt{2} + 1)`, and the projected half-mass radius is
:math:`R_e \simeq 1.33\, r_e`. ``stracers`` also provides Plummer, NFW and
Einasto profiles and simple, broken and double power laws
(:mod:`stracers.profiles`),
the observed Milky Way stellar halo profiles compiled by Medina et al. (2024)
(:mod:`stracers.mw`), and accepts any function
:math:`\nu(r)`.

Eddington's formula
-------------------

For a spherical, isotropic system, the DF of a tracer with density
:math:`\nu` in a potential :math:`\Phi` follows from Eddington's formula,

.. math::
   :label: eddington

   f_*(\mathcal{E}) = \frac{1}{\sqrt{8}\,\pi^2}\left[
   \int_0^{\mathcal{E}} \frac{d^2\nu}{d\Psi^2}
   \frac{d\Psi}{\sqrt{\mathcal{E} - \Psi}}
   + \frac{1}{\sqrt{\mathcal{E}}}
   \left.\frac{d\nu}{d\Psi}\right|_{\Psi = 0}\right],

where :math:`\Psi = -\Phi + \Phi_0` is the relative potential and
:math:`\mathcal{E} = -E + \Phi_0 = \Psi - v^2/2` is the relative energy.
Following Laporte et al. we take :math:`\Phi_0 = 0`, so the potential must go
to zero at infinity. This holds for the potentials written out by codes such
as Gadget. Particles with :math:`\mathcal{E} \le 0` are unbound and have
:math:`f_* = 0`.

The tracer density enters only through :math:`\nu(\Psi)`. This requires a
one-to-one relation between :math:`r` and :math:`\Psi`, i.e. a potential that
decreases monotonically with radius.

Density of states
-----------------

For a spherical potential,

.. math::
   :label: gE

   g(E) = (4\pi)^2 \int_0^{r_E} r^2 \sqrt{2\,(E - \Phi(r))}\, dr,

where :math:`r_E` is the radius at which :math:`\Phi(r_E) = E`.

Triaxial haloes
---------------

Equations :eq:`eddington` and :eq:`gE` hold only for spherical systems. In a
triaxial halo :math:`\Phi(\mathbf{r})` takes many values at a given radius, so
there is no single :math:`\nu(\Psi)`. Laporte et al. handle this in two steps:

1. To relate :math:`\nu` and :math:`\Psi` (and to compute :math:`g(E)`), use
   the **spherically averaged** potential :math:`\Phi(r) = \langle
   \Phi(\mathbf{r}) \rangle`.
2. To define each particle's energy :math:`E = \tfrac{1}{2}v^2 +
   \Phi(\mathbf{r})`, keep the **actual potential of the particle**.

Because the weight depends on the true energy, the stellar profile has
equidensity surfaces that follow the triaxial potential, not spheres.

Truncated stellar haloes
------------------------

A stellar halo can be given an outer edge :math:`r_t` in two ways.

**Energy truncation.** Set :math:`f_*(\mathcal{E}) = 0` for
:math:`\mathcal{E} < \Psi(r_t)`, i.e. for :math:`E > \Phi(r_t)`. In a
spherical potential such stars could reach :math:`r_t`, so removing them
means no star ever goes beyond :math:`r_t`, and the truncated halo is still in
equilibrium. The price is that the density is no longer exactly :math:`\nu`:
it falls to zero at :math:`r_t` and lies below :math:`\nu` close to it,
because the stars removed also spent part of their orbits inside
:math:`r_t`. This is the same idea as lowered (King-type) models.

**Tapered profile.** Replace :math:`\nu(r)` by a profile that falls smoothly
beyond :math:`r_t`, for example :math:`\nu(r)\, e^{-(r/r_t)^p}`, and apply
Eddington's formula as usual. The stars follow the tapered profile, including
the small amount of mass it has beyond :math:`r_t`. Very sharp tapers make a
negative :math:`f_*` near the edge more likely.

Limits
------

- **Equilibrium.** Use only the virialised region of the halo (``rmax``
  :math:`\approx r_{vir}`). Haloes that are merging or strongly perturbed are
  not in equilibrium.
- **Isotropy.** :math:`f_* = f_*(E)` describes an isotropic tracer. Not every
  :math:`\nu` has a non-negative isotropic DF in a given potential. An & Evans
  (2006) showed that the inner slope of an isotropic tracer,
  :math:`\gamma_* = -d\ln\nu/d\ln r`, must be at least half the slope of the
  potential's density cusp. A cored tracer (e.g. Plummer) in an NFW or
  Hernquist cusp (:math:`\gamma = 1`) therefore gives :math:`f_* < 0` at the
  most bound energies. ``stracers`` sets those weights to zero and warns.
- **Equal-mass particles.** :math:`N(E)` is built from particle masses, and the
  weights assume all DM particles have the same mass.
- **Resolution.** The tracer profile can only be resolved where there are
  enough DM particles. Laporte et al. note that compact stellar profiles get
  noisy below ~1 kpc in cosmological haloes.
