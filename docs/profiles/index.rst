.. module:: stracers.profiles

****************************************
Tracer profiles (`stracers.profiles`)
****************************************

Introduction
============

`stracers.profiles` provides stellar density profiles :math:`\nu(r)` to use
as ``tracer_profile`` in `~stracers.compute_weights`. Any function of radius
[kpc] works; the normalisation of a profile does not matter for the weights.

========================  ===========================================  ==================================
Name                      Function                                     Parameters
========================  ===========================================  ==================================
``"plummer"``             `~stracers.plummer`                          ``a``
``"hernquist"``           `~stracers.hernquist`                        ``a``
``"nfw"``                 `~stracers.nfw`                              ``rs``, ``rvir``
``"einasto"``             `~stracers.einasto`                          ``n``, ``r_eff``
``"spl"``                 `~stracers.power_law`                        ``slope``, ``r_core``
``"bpl"``                 `~stracers.broken_power_law`                 ``r_break``, ``inner_slope``,
                                                                       ``outer_slope``, ``r_core``,
                                                                       ``smoothness``
``"dpl"``                 `~stracers.double_power_law`                 ``r_break``, ``inner_slope``,
                                                                       ``outer_slope``, ``alpha``,
                                                                       ``r_core``
========================  ===========================================  ==================================

`~stracers.get_profile` binds a profile to its parameters::

    >>> import numpy as np
    >>> import stracers
    >>> nu = stracers.get_profile("hernquist", a=2.0)
    >>> print(f"{float(nu(1.0)):.4g}  {float(nu(10.0)):.4g}")
    0.01179  1.842e-05

The Plummer, Hernquist, NFW and Einasto profiles are normalised to unit mass
(NFW inside ``rvir``); the power laws are relative densities. The observed
Milky Way stellar halo profiles are in :doc:`../mw/index`.

Analytic profiles
=================

- Plummer: :math:`\nu = \frac{3}{4\pi a^3}\left(1 + r^2/a^2\right)^{-5/2}`
- Hernquist (1990): :math:`\nu = \frac{a}{2\pi r (r + a)^3}`
- NFW: :math:`\nu \propto [x (1 + x)^2]^{-1}`, :math:`x = r/r_s`, normalised
  inside ``rvir``
- Einasto, with the :math:`d_n` approximation of Merritt et al. (2006),
  valid for :math:`n > 0.5`

Power laws
==========

Three power-law families take slopes and break radius set by the user. Slopes
are logarithmic slopes, :math:`d\ln\nu/d\ln r`, so they are negative, as in
the literature. All take an optional core radius,
:math:`r \to \sqrt{r^2 + r_{core}^2}`.

**Simple power law** (SPL): :math:`\nu = r^{\gamma}`::

    >>> spl = stracers.get_profile("spl", slope=-3.5, r_core=1.0)

**Broken power law** (BPL): slope :math:`\gamma_{in}` inside ``r_break`` and
:math:`\gamma_{out}` outside, with :math:`\nu(r_{break}) = 1`. The break is
sharp for ``smoothness=0``, and smoothed over about ``smoothness`` dex
otherwise::

    >>> bpl = stracers.get_profile("bpl", r_break=20.0, inner_slope=-2.5, outer_slope=-4.5)
    >>> float(bpl(20.0))
    1.0

**Double power law** (DPL, Zhao 1996), a smooth turnover that is sharper for
larger :math:`\alpha`:

.. math::

   \nu(r) = x^{\gamma_{in}} \left[\tfrac{1}{2}(1 + x^{\alpha})\right]^{(\gamma_{out} - \gamma_{in})/\alpha},
   \qquad x = r / r_{break}.

Hernquist (:math:`\gamma_{in}=-1, \gamma_{out}=-4, \alpha=1`), NFW
(−1, −3, 1) and Plummer (0, −5, 2) are special cases, and a smoothed BPL is a
DPL with :math:`\alpha = 1/\mathrm{smoothness}`::

    >>> r = np.logspace(-1, 2, 4)
    >>> dpl = stracers.double_power_law(r, 2.0, -1.0, -4.0, alpha=1.0)
    >>> ratio = dpl / stracers.hernquist(r, 2.0)
    >>> bool(np.allclose(ratio, ratio[0]))   # proportional to the Hernquist profile
    True

Two cautions:

- An inner slope of −3 or steeper needs a core, otherwise the mass diverges at
  the centre and most of the stellar mass ends up on the innermost particles.
- A sharp break (BPL with ``smoothness=0``, or a very large ``alpha``) leaves a
  numerical feature near the break, because the energy bins cannot resolve the
  narrow feature it puts into :math:`f_*(E)`.

Truncated profiles
==================

`~stracers.truncate` tapers any profile beyond a radius :math:`r_t`,
:math:`\nu_t(r) = \nu(r)\, e^{-(r/r_t)^p}` with :math:`p` = ``sharpness``.
The stars follow the tapered profile in `~stracers.compute_weights`. For an
equilibrium truncation in energy instead, see ``r_trunc`` in
:doc:`../weights/index`::

    >>> nu_t = stracers.truncate(nu, r_t=50.0, sharpness=4)

Enclosed mass
=============

`~stracers.mass_fraction_within` gives the fraction of a profile's total mass
inside a radius. It is used for ``within="total"`` in
`~stracers.compute_weights`, and needs a finite total mass: an outer slope
steeper than −3, and an inner slope shallower than −3 or a core. Tapering
makes the mass finite::

    >>> print(f"{stracers.mass_fraction_within(nu, 100.0):.6f}  {(100 / 102)**2:.6f}")
    0.961169  0.961169
    >>> nfw = stracers.get_profile("nfw", rs=20.0, rvir=200.0)
    >>> stracers.mass_fraction_within(nfw, 100.0)
    Traceback (most recent call last):
    ...
    ValueError: The total mass of the tracer profile does not converge (outer slope -3 or shallower?). Give the stellar mass inside rmax (within='rmax'), taper the profile with stracers.truncate(nu, r_t), or truncate the stellar halo with r_trunc.
    >>> print(f"{stracers.mass_fraction_within(stracers.truncate(nfw, 100.0), 100.0):.3f}")
    0.901

References
==========

- Hernquist, L. 1990, ApJ, 356, 359
- Merritt, D., Graham, A. W., Moore, B., Diemand, J. & Terzić, B. 2006, AJ, 132,
  2685
- Zhao, H. 1996, MNRAS, 278, 488

API
===

.. automodapi:: stracers.profiles
   :no-heading:
   :no-main-docstr:
   :headings: =-
