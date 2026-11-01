.. module:: stracers.weights

*************************************
Tracer weights (`stracers.weights`)
*************************************

Introduction
============

`stracers.weights` computes the weights of the particle-tagging scheme of
Laporte et al. (2013, Sect. 2; the source is in ``doc/method.tex``), which
generalises that of Bullock & Johnston (2005) to haloes that are not
spherical. `~stracers.compute_weights` takes the dark matter particles of a
halo and a stellar density profile, and returns a `~stracers.WeightsResult`
with the weight of each particle and every intermediate quantity:

====================  ==============================================================
Input                 Description
====================  ==============================================================
``pos``, ``vel``      (N, 3) positions [kpc] and velocities [km/s]
``pot``               (N,) potential of each particle [(km/s)²], zero at infinity
``mass``              particle mass [M☉] (scalar, or (N,) array of equal masses)
``tracer_profile``    a function ν(r); its normalisation does not matter
====================  ==============================================================

The examples on this page use an isotropic Hernquist halo of 10\ :sup:`10` M☉
and scale radius 10 kpc, sampled with `stracers.halos.sample_hernquist`::

    >>> import numpy as np
    >>> import stracers
    >>> from stracers.halos import sample_hernquist
    >>> n, M = 100_000, 1e10
    >>> pos, vel, pot = sample_hernquist(n, M, 10.0, seed=1)
    >>> mass = M / n
    >>> r = np.linalg.norm(pos, axis=1)
    >>> nu = stracers.get_profile("hernquist", a=2.0)

The method
==========

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

Each particle with energy :math:`E` carries both dark matter and stars, through
the weight

.. math::
   :label: weight

   w(E) = \frac{N_*(E)}{N(E)} = \frac{f_*(E)\, g(E)}{N(E)},

where

- :math:`N(E)\,dE` is the *differential energy distribution* of the dark
  matter: the mass of particles with energy between :math:`E` and
  :math:`E + dE`. It is measured directly from the simulation.
- :math:`g(E)` is the *density of states*: the phase-space volume per unit
  energy.
- :math:`f_*(E)` is the stellar DF, and :math:`N_*(E) = f_*(E) g(E)` is the
  differential energy distribution of the stars.

Since :math:`N(E) = f_{DM}(E) g(E)`, the weight is the ratio of the two DFs,
:math:`w = f_*/f_{DM}`. The stellar mass carried by particle :math:`i` is
:math:`m_{*,i} = w(E_i)\, m_{DM}`.

The stellar profile
-------------------

The stellar tracer density :math:`\nu(r) = \rho_*(r)/M_*` can be any profile
that has a non-negative DF in the halo potential (see :ref:`weights-limits`).
Laporte et al. (2013) use a Hernquist (1990) profile,

.. math::

   \nu(r) = \frac{a}{2\pi\, r\, (r + a)^3},

where the scale radius is related to the 3D half-mass radius by
:math:`a = r_e / (\sqrt{2} + 1)`, and the projected half-mass radius is
:math:`R_e \simeq 1.33\, r_e`. See :doc:`../profiles/index` and
:doc:`../mw/index` for the profiles provided.

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
to zero at infinity, as the potentials written by codes such as Gadget do.
Particles with :math:`\mathcal{E} \le 0` are unbound and have :math:`f_* = 0`.
The tracer density enters only through :math:`\nu(\Psi)`, which requires a
potential that decreases monotonically with radius.

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

1. To relate :math:`\nu` and :math:`\Psi`, and to compute :math:`g(E)`, use
   the **spherically averaged** potential
   :math:`\Phi(r) = \langle \Phi(\mathbf{r}) \rangle`.
2. To define each particle's energy
   :math:`E = \tfrac{1}{2}v^2 + \Phi(\mathbf{r})`, keep the **actual potential
   of the particle**.

Because the weight depends on the true energy, the stellar halo has
equidensity surfaces that follow the triaxial potential, not spheres.

.. _weights-limits:

Limits
------

- **Equilibrium.** Use only the virialised region of the halo (``rmax``
  :math:`\approx r_{vir}`). Haloes that are merging or strongly perturbed are
  not in equilibrium.
- **Isotropy.** :math:`f_* = f_*(E)` describes an isotropic tracer, and not
  every :math:`\nu` has a non-negative isotropic DF in a given potential. An &
  Evans (2006) showed that the inner slope of an isotropic tracer,
  :math:`-d\ln\nu/d\ln r`, must be at least half that of the potential's
  density cusp. A cored tracer (e.g. Plummer) in an NFW or Hernquist cusp
  therefore gives :math:`f_* < 0` at the most bound energies; those weights are
  set to zero, with a warning.
- **Equal-mass particles.** :math:`N(E)` is built from particle masses, and the
  weights assume all particles have the same mass.
- **Resolution.** The profile can only be followed where there are enough
  particles. Inside the innermost radial bin with ``min_per_bin`` particles
  (a few tenths of a kpc for 10\ :sup:`5` particles) the stellar halo is
  under-populated.

Computing the weights
=====================

A call to `~stracers.compute_weights` gives the weight of every particle. The
weights are normalised so that the stellar mass inside ``rmax`` is
``stellar_mass``::

    >>> res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8, rmax=250.0)
    >>> res.weights.shape
    (100000,)
    >>> print(f"{np.sum(res.weights * mass):.3g}")
    1e+08
    >>> bool(np.all(res.weights[r >= 250.0] == 0))   # particles outside rmax
    True

Stellar mass
------------

The stellar mass can be given in four ways:

========================================  ==================================================================
Arguments                                 Stellar mass carried by the particles inside ``rmax``
========================================  ==================================================================
``stellar_mass=M``                        :math:`M` (default ``within="rmax"``)
``stellar_fraction=f``                    :math:`f\, M_{DM}(<r_{max})`
``stellar_mass=M, within="total"``        :math:`M\, F(r_{max})`
``stellar_fraction=f, within="total"``    :math:`f\, M_{DM,\,all}\, F(r_{max})`
========================================  ==================================================================

:math:`F(r_{max})` is the fraction of the profile's mass inside ``rmax``
(`~stracers.mass_fraction_within`), and :math:`M_{DM,\,all}` the mass of all
the particles passed in. With ``within="total"`` the stars beyond ``rmax``
belong to the model but are carried by no particle; profiles whose total mass
diverges (outer slope −3 or shallower, e.g. NFW) are rejected. The stellar
mass assigned is in ``res.stellar_mass``::

    >>> res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_fraction=0.01)
    >>> m_dm = mass * np.sum((r > 0) & (r < 300))   # DM mass inside the default rmax
    >>> print(f"{res.stellar_mass / m_dm:.3f}")
    0.010
    >>> res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e9,
    ...                                within="total", rmax=100.0)
    >>> print(f"{res.stellar_mass:.4g}  {1e9 * (100 / 102)**2:.4g}")
    9.612e+08  9.612e+08

Truncating the stellar halo
---------------------------

A stellar halo can be given an outer edge :math:`r_t` in two ways.

**Energy truncation** (``r_trunc``): stars with :math:`E > \Phi(r_t)` get zero
weight. In a spherical potential such stars could reach :math:`r_t`, so no
star ever goes beyond :math:`r_t` and the truncated halo stays in
equilibrium. The density then falls to zero at :math:`r_t` and lies below
:math:`\nu` close to it, as in lowered (King-type) models::

    >>> res = stracers.compute_weights(pos, vel, pot, mass, nu, r_trunc=50.0)
    >>> float(np.sum(res.weights[r > 51.0]))
    0.0

**Tapered profile** (`~stracers.truncate`): the profile is multiplied by
:math:`e^{-(r/r_t)^p}` and the stars follow it, including the small mass it
has beyond :math:`r_t`. Tapering also gives a finite total mass to profiles
that lack one, so they can be used with ``within="total"``::

    >>> nu_t = stracers.truncate(nu, r_t=50.0, sharpness=4)
    >>> res = stracers.compute_weights(pos, vel, pot, mass, nu_t)

Choosing the parameters
-----------------------

``rmax``
   Radius of the region used. Set it to about the virial radius: the method
   assumes equilibrium, and the outskirts of real haloes are not in
   equilibrium. The default (300 kpc) comes from the original notebook.
``center``
   Halo centre, e.g. from a shrinking-sphere method. A wrong centre is the
   easiest way to get wrong weights without any warning.
``dlogr``, ``min_per_bin``
   Radial binning of the averaged potential. The defaults (0.1 dex, 20
   particles) work well from 10\ :sup:`5` to 10\ :sup:`7` particles. Use
   larger bins if :math:`\Psi(r)` is reported as non-monotonic.
``n_ebins``
   Number of energy bins (default 100). With fewer than about 10\ :sup:`5`
   particles, fewer bins reduce the noise in :math:`N(E)`.
``derivative``
   Keep ``"spline"`` unless you want to reproduce the original notebook
   (``"gradient"``).
``boundary_term``
   Keep `True`. Turning it off reproduces the notebook, and makes the stellar
   density fall off too early near ``rmax``.
``truncation_correction``
   Shifts the potential by :math:`-G M / r_{max}`, as in the original
   notebook. It is off by default (:math:`\Phi_0 = 0`, as in Laporte et al.).

Diagnostics
-----------

`stracers.plotting.plot_diagnostics` plots :math:`\Psi(r)`,
:math:`d^2\nu/d\Psi^2`, :math:`N(E)`, :math:`g(E)`, :math:`f_*(E)` and
:math:`w(E)` from a `~stracers.WeightsResult` (it needs matplotlib)::

    >>> from stracers.plotting import plot_diagnostics
    >>> fig = plot_diagnostics(res)  # doctest: +SKIP

The most direct check is to compare the stellar mass in radial shells with
that of the profile, as in :doc:`../getting_started`.

Many simulations
----------------

A call takes about 0.2 s per million particles, so for large sets of
simulations reading the snapshots usually takes longer than computing the
weights. Process the simulations in parallel and record the diagnostics::

    import warnings
    from concurrent.futures import ProcessPoolExecutor

    def process_one(path):
        pos, vel, pot, mass, ids, center = read_halo(path)      # your reader
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            res = stracers.compute_weights(pos, vel, pot, mass, nu,
                                           stellar_mass=1e8, rmax=250.0, center=center)
        stracers.save_weights(path + ".weights.h5", ids, res.weights,
                              stellar_mass=1e8, n_negative_f=res.n_negative_f)
        return path, res.n_negative_f, [str(w.message) for w in caught]

    with ProcessPoolExecutor(max_workers=16) as pool:
        for path, n_neg, msgs in pool.map(process_one, paths):
            print(path, n_neg, msgs)

How the code works
==================

`~stracers.compute_weights` calls one function per step of the method. Each
is listed in the API below.

1. **Selection.** The positions are shifted by ``center`` and only particles
   with :math:`0 < r < r_{max}` are kept. The energy of each particle uses its
   own potential, :math:`E_i = \tfrac{1}{2}|\mathbf{v}_i|^2 + \Phi_i`.
2. **Averaged potential** (`radial_potential`). Particles are binned in
   :math:`\log_{10} r` with width ``dlogr``; each bin gives the mean potential
   and its standard error. Bins with fewer than ``min_per_bin`` particles are
   dropped.
3. **Smoothing** (`smooth_potential`, ``derivative="spline"``). Eddington's
   formula needs :math:`d^2\nu/d\Psi^2`, and differentiating twice amplifies
   noise. A quartic smoothing spline is fitted to :math:`\Phi(\ln r)`,
   weighted by the inverse standard errors, and resampled on ``n_fine``
   logarithmic cells that cover exactly the radial range of the bins kept,
   from the inner edge of the first to the outer edge of the last (half a bin
   beyond the first and last bin centres). With ``derivative="gradient"`` the
   binned potential is used directly.
4. **Tracer derivatives** (`tracer_derivatives`). :math:`d\nu/d\Psi` and
   :math:`d^2\nu/d\Psi^2` by second-order finite differences on the radial
   grid.
5. **Energy grid.** One linear grid of ``n_ebins`` bins spans the energies of
   the selected particles. :math:`N(E)`, :math:`g(E)` and :math:`f_*(E)` are
   evaluated at the same bin centres; :math:`N(E)` is the particle mass in each
   bin divided by the bin width.
6. **Density of states** (`density_of_states`). Eq. :eq:`gE` as a sum over the
   radial grid, plus the sphere inside the first bin,

   .. math::

      g(E) \simeq (4\pi)^2 \left[\frac{r_{in}^3}{3}\sqrt{2(E-\Phi_1)}
      + \sum_{k:\,\Phi_k<E} r_k^2\, \Delta r_k \sqrt{2(E - \Phi_k)}\right].

   The sum stops at the outer edge of the outermost bin, not at :math:`r_E`:
   :math:`N(E)` only counts particles inside that radius (``res.r_range[1]``,
   the outermost selected particle unless sparse outer bins were dropped), so
   :math:`g` must be truncated in the same way for :math:`f_{DM} = N/g` to be
   unbiased. With ``derivative="gradient"`` the grid only has the bins kept,
   so particles in a dropped intermediate bin are left out of :math:`N(E)`
   too; the spline grid is continuous across such a bin.
7. **Eddington inversion** (`eddington`, `outer_eddington`). :math:`f_*` is
   evaluated at :math:`\mathcal{E} = -E` for each bin. Inside the data
   (:math:`\Psi \ge \Psi(r_{max})`), :math:`d^2\nu/d\Psi^2` is taken as constant
   on each segment and :math:`1/\sqrt{\mathcal{E}-\Psi}` is integrated
   exactly, which removes the noise from the integrable singularity. Beyond
   the data (``boundary_term=True``), :math:`\nu(\Psi)` is extrapolated as a
   power law, with the slope measured as a secant over the outermost
   ``dlogr`` of the grid, and the integral is done in closed form with the
   incomplete beta function, evaluated in log space so that steep profiles do
   not overflow. Where that slope differs from the derivative at the end of
   the grid, :math:`d\nu/d\Psi` jumps at :math:`\Psi_{min}`, and the term
   :math:`\Delta(d\nu/d\Psi)/\sqrt{\mathcal{E}-\Psi_{min}}` that the jump
   adds to Eq. :eq:`eddington` is included.
   Bins with :math:`f_* < 0` are set to zero and counted in
   ``n_negative_f``. With ``r_trunc``, :math:`f_*` and the weights of all
   particles with :math:`E > \Phi(r_t)` are set to zero.
8. **Weights.** :math:`w = f_* g / N` in each bin (zero where :math:`N = 0`);
   each particle gets the weight of its bin (:func:`numpy.digitize`), and the
   weights are scaled to the stellar mass.

Differences from the paper and the original notebook
----------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Item
     - Laporte et al. (2013)
     - Original notebook
     - ``stracers``
   * - :math:`\Phi_0`
     - 0
     - shifted by :math:`G M / r_{max}`
     - 0 (the shift is ``truncation_correction``)
   * - :math:`d^2\nu/d\Psi^2`
     - not specified
     - ``np.gradient`` on the binned potential
     - spline-smoothed potential (default), or ``"gradient"``
   * - Eddington integral
     - full, with the boundary term
     - Riemann sum over the data only
     - exact kernel per segment, plus a power-law extension beyond ``rmax``
   * - :math:`g(E)` lower limit
     - :math:`r = 0`
     - first radial bin
     - :math:`r = 0` (inner sphere added)
   * - :math:`f(\mathcal{E})` on the energy bins
     - :math:`\mathcal{E} = -E`
     - sorted :math:`f` (not equivalent)
     - evaluated at :math:`\mathcal{E} = -E`
   * - Normalisation
     - not specified
     - :math:`\sum w m = M_*`
     - :math:`\sum w m = M_*`

Validation and performance
--------------------------

The test suite samples an isotropic Hernquist halo (10\ :sup:`10` M☉,
:math:`a = 10` kpc, 2 × 10\ :sup:`5` particles) and checks that:

- :math:`f_*` for the halo's own profile matches the analytic Hernquist DF to
  ~1%, and :math:`g(E)` matches Eq. :eq:`gE` to within 3%;
- self-consistent weights (tracer = halo) are constant to within 5%;
- Hernquist, Plummer, Einasto, power-law and Milky Way tracers are recovered to
  within 10% between 0.5 and 100 kpc;
- truncation, the stellar mass options, recentring and masking behave as
  described.

The cost grows linearly with the number of particles: about 0.02 s for
10\ :sup:`5` particles, 0.17 s for 10\ :sup:`6` and 3 s for 10\ :sup:`7` on a
desktop CPU.

Troubleshooting
===============

"energy bins have f(E) < 0"
   A few negative bins at the most bound energies are normal: there are few
   particles there. Many negative bins mean the tracer has no isotropic DF in
   this potential, for example a cored profile in a cuspy halo (see
   :ref:`weights-limits`). Use a cuspier tracer, or accept that the centre will
   be under-populated.
"psi(r) is not strictly decreasing"
   The averaged potential is noisy or the centre is wrong. Check ``center``,
   increase ``dlogr``, or use ``derivative="spline"``.
"All weights are zero"
   Usually a potential with the wrong sign or zero point (it must be negative
   and go to zero at infinity), or a profile that is zero over the selected
   region.
"Non-finite weights"
   The profile is too steep or noisy for the radial grid. Use
   ``derivative="spline"`` or add a core.

References
==========

- An, J. H. & Evans, N. W. 2006, ApJ, 642, 752
- Binney, J. & Tremaine, S. 2008, *Galactic Dynamics*, 2nd ed., Princeton
  University Press
- Bullock, J. S. & Johnston, K. V. 2005, ApJ, 635, 931
- Hernquist, L. 1990, ApJ, 356, 359
- Laporte, C. F. P., White, S. D. M., Naab, T. & Gao, L. 2013, MNRAS, 435, 901,
  `arXiv:1301.5319 <https://arxiv.org/abs/1301.5319>`_
- Natarajan, P., Hjorth, J. & van Kampen, E. 1997, MNRAS, 286, 329

API
===

.. automodapi:: stracers.weights
   :no-heading:
   :no-main-docstr:
   :no-inheritance-diagram:
   :headings: =-

Algorithm steps
---------------

.. autofunction:: stracers.weights.radial_potential
.. autofunction:: stracers.weights.smooth_potential
.. autofunction:: stracers.weights.tracer_derivatives
.. autofunction:: stracers.weights.density_of_states
.. autofunction:: stracers.weights.eddington
.. autofunction:: stracers.weights.outer_eddington

Diagnostics
-----------

.. automodapi:: stracers.plotting
   :no-heading:
   :no-main-docstr:
   :headings: ^"
