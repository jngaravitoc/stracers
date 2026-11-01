.. module:: stracers.nstars

**********************************
Star counts (`stracers.nstars`)
**********************************

Introduction
============

`stracers.nstars` gives the number of stars of a stellar halo in a range of
galactocentric radius,

.. math::

   N = \int_{r_{min}}^{r_{max}} r^2\, dr \int d\Omega\;
       n\!\left(r \sqrt{\sin^2\theta + \cos^2\theta / q^2}\right),

where :math:`n` is the number density at the ellipsoidal radius
:math:`r_q = \sqrt{R^2 + z^2/q^2}` and :math:`q` the flattening of the halo
(:math:`q < 1` is oblate). `~stracers.number_of_stars` takes any profile, or
the key of a profile in :doc:`../mw/index`, and `~stracers.number_of_stars_table`
runs it for all 31 Milky Way profiles. This is the calculation of Medina et
al. (2024, Sect. 5.2.1 and Fig. 7) for the number of RR Lyrae beyond a given
radius.

Normalisation
=============

The number density can be given:

- by the profile itself, in stars / kpc\ :sup:`3`, for a function of radius;
- by a number density ``n0`` at a radius ``r0`` in the plane, e.g. the local
  density of RR Lyrae at the solar radius;
- by a number of stars ``n_ref`` over the whole sky in a range ``ref_range``.

The Milky Way profiles have an arbitrary normalisation (:math:`\nu = 1` at the
break, or at 1 kpc), so ``n0`` or ``n_ref`` is required for them.

Counting RR Lyrae
=================

Medina et al. (2024) expect about 3000 RR Lyrae beyond 100 kpc over the whole
sky for their best fit, 0.073 per square degree, about 25 in their 350 deg²
survey, and about 350 beyond 230 kpc. Normalising their profile to 3000 stars
between 100 and 300 kpc reproduces these numbers::

    >>> import numpy as np
    >>> import stracers
    >>> from stracers.nstars import FULL_SKY_DEG2
    >>> kw = dict(n_ref=3000, ref_range=(100, 300))
    >>> n100, n230 = stracers.number_of_stars("medina2024", [100, 230], r_max=300, **kw)
    >>> print(f"{n100:.0f}  {n230:.0f}  {n100 / FULL_SKY_DEG2:.3f} per deg2")
    3000  356  0.073 per deg2
    >>> n_survey = stracers.number_of_stars("medina2024", 100, r_max=300, area_deg2=350, **kw)
    >>> print(f"{n_survey:.1f}")
    25.5

An array of inner radii gives cumulative counts N(>r), as in their Fig. 7.

The outer edge
--------------

The ~350 stars beyond 230 kpc require the halo to end at about 300 kpc
(``r_max``). An untruncated −4.47 power law gives 2.5 times more::

    >>> n100, n230 = stracers.number_of_stars("medina2024", [100, 230], r_max=np.inf,
    ...                                       n_ref=3000, ref_range=(100, np.inf))
    >>> print(f"{n230:.0f}")
    882

With ``r_max=np.inf`` the counts must converge: profiles with outer slopes of
−3 or shallower raise an error.

Flattening
----------

A flattened halo changes the counts, but for a power law only by a constant
factor, so ratios such as N(>230)/N(>100) do not change::

    >>> ratio = lambda q: stracers.number_of_stars("medina2024", 230, r_max=300, q=q, n_ref=1,
    ...                                            ref_range=(100, 300))
    >>> print(f"{ratio(1.0):.4f}  {ratio(0.6):.4f}")
    0.1186  0.1186

All the literature profiles
===========================

`~stracers.number_of_stars_table` applies the same normalisation to every
profile in :doc:`../mw/index`, to compare the predictions of different
slopes. Profiles whose counts diverge for the requested range get ``nan``,
with a warning::

    >>> table = stracers.number_of_stars_table(230, r_max=300, area_deg2=350, **kw)
    >>> print(f"{table['medina2024']:.2f}  {table['fukushima2019']:.2f}  {table['sesar2013']:.2f}")
    3.02  0.06  7.71

With ``rng``, the profile parameters are drawn within their uncertainties;
repeat the call to get the spread of the predicted counts.

Caveats
=======

- Survey areas are scaled isotropically (``area_deg2 / 41253``), as in Medina et
  al. (2024). This ignores the direction dependence of a flattened halo and
  the survey footprint.
- The literature fits are only valid over their radial ranges; counts at
  larger radii are extrapolations.

How the counts are computed
===========================

- The angular integral uses 64-point Gauss–Legendre quadrature in
  :math:`\mu = \cos\theta`, exact for a spherical halo.
- The radial integral uses :func:`scipy.integrate.quad` in :math:`\ln r`, with
  a short piece in linear :math:`r` near the centre when :math:`r_{min} = 0`.
- For ``r_max = np.inf`` the integral is taken to 10\ :sup:`8` × ``r_min``. It
  is rejected if radii beyond 10\ :sup:`6` × ``r_min`` contribute more than
  1%.

References
==========

- Medina, G. E., et al. 2024, MNRAS, 531, 4762,
  `doi:10.1093/mnras/stae1137 <https://doi.org/10.1093/mnras/stae1137>`_

API
===

.. automodapi:: stracers.nstars
   :no-heading:
   :no-main-docstr:
   :headings: =-
