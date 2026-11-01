.. module:: stracers.tracers

**********************************************
Mock tracer populations (`stracers.tracers`)
**********************************************

Introduction
============

`~stracers.paint_tracer` paints a population of a stellar tracer, such as RR
Lyrae, BHB stars or K giants, onto the dark matter particles of a halo, in a
given radial range:

1. it chooses a literature profile for that tracer from :doc:`../mw/index`;
2. it computes the weights of the **whole** stellar halo with that profile
   (`~stracers.compute_weights`);
3. it keeps the particles in the radial range, with the expected number of
   stars of that type on each one.

The result is a `~stracers.TracerPopulation`, with the velocity profile of
the population and mock star catalogues. Radii are in kpc, like the profiles
(see :doc:`../units`).

The radial range is applied **after** the weighting, as a selection, like an
observational cut. Truncating the profile to the range instead
(``r_trunc`` or `~stracers.truncate`) would change the distribution function
of the stars and therefore their kinematics.

Tracer names
============

Names are matched without regard to case, spaces, hyphens or a plural "s":

======================  ======================================  ==========================================
Name                    Accepted spellings (e.g.)               Fits with tracer
======================  ======================================  ==========================================
``"RR Lyrae"``          ``"RRLyr"``, ``"RRL"``, ``"rr lyrae"``  RRLs, RRab
``"RRab"``              ``"RRab"``                              RRab
``"BHB"``               ``"BHB"``, ``"bhb stars"``              BHB stars, BHB and BS stars
``"K giants"``          ``"K-giant"``, ``"k giants"``           K giants
``"M giants"``          ``"M-giant"``                           M giants
``"Giants"``            ``"giant"``                             Giant stars, K giants, M giants
``"MSTO"``              ``"MSTO"``                              MSTO stars, near MSTO stars
``"F stars"``           ``"F-star"``                            F stars
``"A-type stars"``      ``"A-type"``                            A-type stars
======================  ======================================  ==========================================

`~stracers.tracers.available_tracers` returns the keys of the fits for each
name, and `~stracers.tracers.tracer_fits` those of one tracer.

Choosing the fit
================

Most tracers have several published fits. By default
`~stracers.tracers.choose_fit` takes the **most recent fit whose fitted range
covers the requested range**. Ties go to broken over simple power laws, then
to the widest fitted range. If no fit covers the range, the one covering the
largest part of it is used, with a warning that the rest of the range uses the
extrapolated profile::

    >>> from stracers.tracers import choose_fit
    >>> choose_fit("RRLyr", (20, 100))      # Medina+2024 (7-145 kpc); Feng+2024 is also 2024, but an SPL
    'medina2024'
    >>> choose_fit("BHB", (20, 100))        # the only BHB fit covering 20-100 kpc
    'thomas2018'
    >>> choose_fit("K-giant", (20, 100))    # warns: no K-giant fit reaches 100 kpc
    'xue2015_bpl'

Pass ``key`` to `~stracers.paint_tracer` to use a specific fit instead.

Painting a population
=====================

On the toy halo of `stracers.halos` (in kpc, km/s, M☉), a population of 3000
RR Lyrae between 20 and 100 kpc::

    >>> import numpy as np
    >>> import stracers
    >>> from stracers.halos import sample_hernquist
    >>> pos, vel, pot = sample_hernquist(100_000, 1e10, 10.0, seed=1)
    >>> mass = 1e5
    >>> pop = stracers.paint_tracer(pos, vel, pot, mass, "RRLyr", r_range=(20, 100), n_stars=3000)
    >>> pop.key, pop.reference, pop.url
    ('medina2024', 'Medina et al. 2024', 'https://arxiv.org/abs/2402.14055')
    >>> print(f"{pop.total:.0f} stars on {pop.in_range.sum()} particles")
    3000 stars on 38085 particles

``pop.n_stars`` is the expected number of stars of the tracer on each
particle, and zero outside the range.

Normalisation
-------------

- ``n_stars``: the number of stars in ``r_range``.
- ``n_ref`` and ``ref_range``: the number of stars in another range, e.g. an
  observed count. The profile then sets the number in ``r_range``::

      >>> p = stracers.paint_tracer(pos, vel, pot, mass, "RRLyr", (20, 100),
      ...                           n_ref=3000, ref_range=(100, 300))
      >>> print(f"{p.total:.0f}")
      36120

- Neither: ``n_stars`` holds the fraction of the tracer population (inside
  ``rmax``) on each particle.

Kinematics
----------

`~stracers.TracerPopulation.velocity_profile` gives the velocity profile of
the population, weighted by its number of stars (see
:doc:`../kinematics/index`). An integer number of bins spans ``r_range``::

    >>> prof = pop.velocity_profile(bins=4)
    >>> print(np.round(prof.sigma_r, 1))
    [13.6 12.1 10.5  9.1]

A BHB population in the same range follows a different profile (Thomas et al.
2018), and so has different kinematics::

    >>> bhb = stracers.paint_tracer(pos, vel, pot, mass, "BHB", r_range=(20, 100), n_stars=3000)
    >>> print(bhb.key, np.round(bhb.velocity_profile(bins=4).sigma_r, 1))
    thomas2018 [14.  13.2 12.  10.3]

Mock catalogues
---------------

`~stracers.TracerPopulation.draw_stars` draws a Poisson number of stars on
each particle and returns the particle index of every star, so the stars'
positions and velocities are those of their particles::

    >>> idx = pop.draw_stars(rng=0)
    >>> stars_pos, stars_vel = pop.pos[idx], pop.vel[idx]
    >>> stars_pos.shape
    (2958, 3)

Other options
=============

- ``profile_kwargs`` is passed to `stracers.mw.mw_profile`: ``r_core``,
  ``smoothness``, and ``rng`` to draw the profile parameters within their
  uncertainties.
- Other keyword arguments go to `~stracers.compute_weights`. ``rmax`` defaults
  to ``max(300, r_range[1])`` kpc, and ``center`` is used both for the weights
  and for the radial selection.

API
===

.. automodapi:: stracers.tracers
   :no-heading:
   :no-main-docstr:
   :no-inheritance-diagram:
   :headings: =-
