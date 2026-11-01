.. module:: stracers.mw

**************************************************
Milky Way stellar halo profiles (`stracers.mw`)
**************************************************

Introduction
============

`stracers.mw` provides the 31 simple and broken power-law fits to the number
density of the Milky Way stellar halo compiled by Medina et al. (2024, MNRAS,
531, 4762, Table 6), and the power law with two breaks fitted by Li et al.
(2026) to DESI DR2 K giants, 32 profiles in all. The data are in
``src/stracers/data/mw_stellar_halo_profiles.csv``, with the uncertainties,
the radial range of each fit, the tracer and the reference.

The compilation
===============

Radii are galactocentric, in kpc; the range is the one used in each fit, and
each reference links to its paper on arXiv (or its DOI, for the two papers that
are not on arXiv).
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

A power law with two breaks (``"tpl"``), with slope :math:`\gamma_{in}`
inside :math:`r_{break}`, :math:`\gamma_{mid}` between the breaks and
:math:`\gamma_{out}` outside :math:`r_{break,2}`. It is from Li et al. (2026,
Table 1, whole halo), not from the compilation:

.. csv-table::
   :file: tpl.csv
   :header: "Key", ":math:`r_{break}`", ":math:`r_{break,2}`", ":math:`\gamma_{in}`", ":math:`\gamma_{mid}`", ":math:`\gamma_{out}`", "Range", "Tracer", "Reference"
   :widths: auto
   :class: mw-table

Some entries need care: ``deason2014`` is an upper limit (slope < −6),
``vivas2006`` an approximate value, and ``thomas2018`` has an outer slope
shallower than its inner one.

The profiles
============

The figures show the profiles in their published form (no core, sharp
breaks), each normalised to 1 at 30 kpc (dotted line). Each profile is a solid
line over the radial range of its fit, and a dashed line extrapolating it out
to 200 kpc. Fits quoted only with an outer radius (e.g. "<40 kpc") are drawn
from 1 kpc, and the profiles whose fits do not reach 30 kpc are normalised with
their extrapolation.

RR Lyrae
--------

.. plot::

   from stracers.plotting import plot_mw_profiles
   ax = plot_mw_profiles(tracers=["RRLs", "RRab"])
   ax.set_title(r"$\mathrm{RR\ Lyrae\ (RRLs\ and\ RRab)}$")

Blue horizontal branch stars
----------------------------

.. plot::

   from stracers.plotting import plot_mw_profiles
   ax = plot_mw_profiles(tracers=["BHB stars", "BHB and BS stars"])
   ax.set_title(r"$\mathrm{Blue\ horizontal\ branch\ stars}$")

Giants and other tracers
------------------------

K, M and other giants, main-sequence turn-off (MSTO), F and A-type stars:

.. plot::

   from stracers.plotting import plot_mw_profiles
   ax = plot_mw_profiles(tracers=["Giant stars", "K giants", "M giants", "MSTO stars",
                                  "near MSTO stars", "F stars", "A-type stars"])
   ax.set_title(r"$\mathrm{Giants\ and\ other\ tracers}$")

The figures are made with `stracers.plotting.plot_mw_profiles`, in LaTeX style
(Computer Modern, all text in math mode) with
`stracers.plotting.LATEX_STYLE`::

    import matplotlib.pyplot as plt
    from stracers.plotting import plot_mw_profiles, LATEX_STYLE

    with plt.rc_context(LATEX_STYLE):
        ax = plot_mw_profiles(tracers=["RRLs", "RRab"])
        ax.figure.savefig("mw_rrl.pdf", bbox_inches="tight")

Using the profiles
==================

`~stracers.mw.available` lists the keys and `~stracers.mw.get_entry` returns
one row of the table::

    >>> import numpy as np
    >>> from stracers import mw
    >>> len(mw.available()), len(mw.available("spl")), len(mw.available("bpl")), len(mw.available("tpl"))
    (32, 16, 15, 1)
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

With these defaults, all 32 profiles are recovered to within 5% over their
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

24 of the 32 profiles have a finite total mass and can be used with
``within="total"`` in `~stracers.compute_weights`. The other eight have outer
slopes of −3 or shallower (``wetterer1996``, ``vivas2006``, ``bell2008``,
``depropris2010``, ``deason2011_spl``, ``sesar2013``, ``iorio2018``,
``thomas2018``; for ``thomas2018``, −3.21, the mass converges only at absurdly
large radii). Taper them with `~stracers.truncate`, truncate them with
``r_trunc``, or give the stellar mass inside ``rmax``.

Caveats
=======

- The fits are **spherical**, while several of the papers fit flattened haloes;
  only the radial profile is kept. For Li et al. (2026), a triaxial halo
  (axis ratios 1 : 0.85 : 0.74, tilted by 44° from the disc), the radius of the
  profile is the ellipsoidal radius.
- Each fit is only valid over its range (mostly ~5–100 kpc). Outside it the
  power laws are extrapolations.

References
==========

The compilation (all but ``li2026``) is from

- Medina, G. E., et al. 2024, MNRAS, 531, 4762,
  `arXiv:2402.14055 <https://arxiv.org/abs/2402.14055>`__,
  `doi:10.1093/mnras/stae1137 <https://doi.org/10.1093/mnras/stae1137>`__

and the fits are from

.. include:: references.inc

The published reference of ``li2026`` is

- Li, S., Wang, W., Koposov, S. E., et al. 2026, ApJ, 999, 108,
  `arXiv:2512.01350 <https://arxiv.org/abs/2512.01350>`__,
  `doi:10.3847/1538-4357/ae41b9 <https://doi.org/10.3847/1538-4357/ae41b9>`__

API
===

.. automodapi:: stracers.mw
   :no-heading:
   :no-main-docstr:
   :headings: =-
