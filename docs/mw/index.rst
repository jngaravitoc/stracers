.. module:: stracers.mw

**************************************************
Milky Way stellar halo profiles (`stracers.mw`)
**************************************************

Introduction
============

`stracers.mw` provides the 31 simple and broken power-law fits to the number
density of the Milky Way stellar halo compiled by Medina et al. (2024, MNRAS,
531, 4762, Table 6). The data are in
``src/stracers/data/mw_stellar_halo_profiles.csv``, with the uncertainties,
the radial range of each fit, the tracer and the reference.

The compilation
===============

Radii are galactocentric, in kpc; the range is the one used in each fit.
Slopes are logarithmic slopes of the number density, :math:`d\ln n/d\ln r`.

Simple power laws (``"spl"``), :math:`n \propto r^{\gamma}`:

.. csv-table::
   :file: spl.csv
   :header: "Key", "Slope", "Range", "Tracer", "Reference"
   :widths: auto
   :class: mw-table

Broken power laws (``"bpl"``), with slope :math:`\gamma_{in}` inside
:math:`r_{break}` and :math:`\gamma_{out}` outside:

.. csv-table::
   :file: bpl.csv
   :header: "Key", ":math:`r_{break}`", ":math:`\gamma_{in}`", ":math:`\gamma_{out}`", "Range", "Tracer", "Reference"
   :widths: auto
   :class: mw-table

Some entries need care: ``deason2014`` is an upper limit (slope < −6),
``vivas2006`` an approximate value, and ``thomas2018`` has an outer slope
shallower than its inner one.

Using the profiles
==================

`~stracers.mw.available` lists the keys and `~stracers.mw.get_entry` returns
one row of the table::

    >>> import numpy as np
    >>> from stracers import mw
    >>> len(mw.available()), len(mw.available("spl")), len(mw.available("bpl"))
    (31, 16, 15)
    >>> e = mw.get_entry("medina2024")
    >>> e["r_break"], e["inner_slope"], e["outer_slope"], e["tracer"]
    (18.1, -2.05, -4.47, 'RRLs')

`~stracers.mw.mw_profile` returns the profile as a function of radius, to use
in `~stracers.compute_weights` or `~stracers.number_of_stars`::

    >>> nu = mw.mw_profile("medina2024")
    >>> nu.params
    {'r_break': 18.1, 'inner_slope': -2.05, 'outer_slope': -4.47, 'r_core': 1.0, 'smoothness': 0.1}

Core and break
--------------

- None of the fits extends to the centre, and power laws steeper than −3 have
  a divergent mass there. `~stracers.mw.mw_profile` adds a core,
  :math:`r \to \sqrt{r^2 + r_{core}^2}`, with ``r_core=1`` kpc by default.
  For steep profiles the stellar mass inside a few ``r_core`` depends strongly
  on its value: without a core, the Feng et al. (2024) slope of −4.09 puts 91%
  of the stellar mass inside 1 kpc.
- Broken power laws are smoothed over ``smoothness=0.1`` dex around the break.
  The sharp published break puts a narrow feature into :math:`f_*(E)` that the
  energy bins cannot resolve: a ~10% dip at the break for Medina et al. (2024),
  and an excess of up to 2× for the −15 outer slope of Fukushima et al.
  (2019). ``smoothness=0`` gives the published form.

With these defaults, all 31 profiles are recovered to within 5% over their
fitted ranges by `~stracers.compute_weights` in the test halo.

Uncertainties
-------------

With ``rng``, each parameter is drawn within its uncertainties, read as the
median and the 16th and 84th percentiles (two half-normals joined at the
median). Parameters without uncertainties stay fixed. Correlations between
parameters, e.g. the break radius and the slopes, are not in the
compilation::

    >>> rng = np.random.default_rng(1)
    >>> r_break = [mw.mw_profile("medina2024", rng=rng).params["r_break"] for _ in range(4000)]
    >>> print(np.round(np.percentile(r_break, [16, 50, 84]), 1))   # 18.1 (-1.1, +2.1)
    [17.  18.1 20.1]

Total mass
----------

23 of the 31 profiles have a finite total mass and can be used with
``within="total"`` in `~stracers.compute_weights`. The other eight have outer
slopes of −3 or shallower (``wetterer1996``, ``vivas2006``, ``bell2008``,
``depropris2010``, ``deason2011_spl``, ``sesar2013``, ``iorio2018``,
``thomas2018``; for ``thomas2018``, −3.21, the mass converges only at absurdly
large radii). Taper them with `~stracers.truncate`, truncate them with
``r_trunc``, or give the stellar mass inside ``rmax``.

Caveats
=======

- The fits are **spherical**, while several of the papers fit flattened haloes;
  only the radial profile is kept.
- Each fit is only valid over its range (mostly ~5–100 kpc). Outside it the
  power laws are extrapolations.

References
==========

- Medina, G. E., et al. 2024, MNRAS, 531, 4762,
  `doi:10.1093/mnras/stae1137 <https://doi.org/10.1093/mnras/stae1137>`_. The
  references of the individual fits are listed in the table.

API
===

.. automodapi:: stracers.mw
   :no-heading:
   :no-main-docstr:
   :headings: =-
