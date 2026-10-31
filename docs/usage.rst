Usage
=====

Installation
------------

.. code-block:: bash

   pip install -e ".[io,plot]"     # add "test" or "docs" for development

The only required dependencies are NumPy and SciPy. ``h5py`` is needed for
saving and loading weights, and ``matplotlib`` for the diagnostic plots.

Weights for one snapshot
------------------------

``stracers`` takes plain arrays. With a Gadget snapshot, for example:

.. code-block:: python

   import numpy as np
   import stracers
   from pygadgetreader import readsnap

   pos  = readsnap(snap, "pos", "dm")            # kpc
   vel  = readsnap(snap, "vel", "dm")            # km/s
   pot  = readsnap(snap, "pot", "dm")            # (km/s)^2
   ids  = readsnap(snap, "pid", "dm")
   mass = readsnap(snap, "mass", "dm") * 1e10    # Msun

   nu = stracers.get_profile("hernquist", a=2.0)
   res = stracers.compute_weights(
       pos, vel, pot, mass, nu,
       stellar_mass=1e8, rmax=250.0, center=halo_center,
   )

   m_star = res.weights * mass        # stellar mass per particle, sums to 1e8

The stellar mass can also be given as a fraction of the DM mass, or as the
total mass of the profile rather than the mass inside ``rmax``:

.. code-block:: python

   # 1% of the DM mass inside rmax
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_fraction=0.01, rmax=250.0)

   # total stellar mass 1e9 Msun; only the part inside rmax is carried by particles
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e9,
                                  within="total", rmax=250.0)
   res.stellar_mass        # mass assigned inside rmax = 1e9 * F(rmax)

   # 0.5% of the DM mass of all particles passed in, as a total stellar mass
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_fraction=0.005,
                                  within="total", rmax=250.0)

Milky Way stellar halo profiles
-------------------------------

:mod:`stracers.mw` provides the 31 simple and broken power-law fits to the
Milky Way stellar halo compiled by Medina et al. (2024, Table 6). The data are
in ``stracers/data/mw_stellar_halo_profiles.csv``.

.. code-block:: python

   from stracers import mw

   mw.available()              # all keys; mw.available("bpl") for broken power laws
   mw.get_entry("medina2024")  # parameters, uncertainties, fitted range, tracer, reference

   nu = mw.mw_profile("medina2024")        # r_core=1 kpc, smoothness=0.1 by default
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e9, rmax=250.0)

   # a different random draw of the parameters within their uncertainties
   rng = np.random.default_rng(42)
   nu = mw.mw_profile("medina2024", rng=rng)
   nu.params                    # the parameters actually used

Things to keep in mind:

- The fits are **spherical** and only valid over their fitted range
  (``r_min`` to ``r_max``, mostly ~5–100 kpc). Outside it the power laws are
  extrapolations.
- None of them is fitted at the centre, and slopes steeper than −3 have a
  divergent mass there. ``mw_profile`` therefore adds a core,
  :math:`r \to \sqrt{r^2 + r_{core}^2}`, with ``r_core=1`` kpc by default.
  For steep profiles the stellar mass inside a few ``r_core`` depends strongly
  on this choice. Without a core, the Feng et al. (2024) slope of −4.09 puts
  91% of the stellar mass inside 1 kpc.
- Broken power laws are smoothed over ``smoothness=0.1`` dex around the
  break. The sharp published break puts a narrow feature into f(E) that the
  energy bins cannot resolve. That shows up as a ~10% dip at the break for
  Medina et al. (2024), and as an excess of up to 2× for the −15 outer slope
  of Fukushima et al. (2019). Use ``smoothness=0`` to reproduce the published
  form exactly.
- With these defaults, all 31 profiles are recovered to within 5% over their
  fitted ranges in the test halo.
- 23 of the 31 profiles have a finite total mass and can be used with
  ``within="total"``. The other 8 have outer slopes of −3 or shallower
  (``wetterer1996``, ``vivas2006``, ``bell2008``, ``depropris2010``,
  ``deason2011_spl``, ``sesar2013``, ``iorio2018``, ``thomas2018``). For
  ``thomas2018`` (−3.21) the mass converges only at absurdly large radii.
  ``compute_weights`` raises an error for these 8. Taper them with
  :func:`stracers.truncate`, truncate them with ``r_trunc``, or give the mass
  inside ``rmax``.
- Uncertainties are read as median and 16th/84th percentiles, and each
  parameter is drawn independently. Correlations between parameters (e.g. the
  break radius and the slopes) are not available in the compilation.

Truncating the stellar halo
---------------------------

Two options, described in :doc:`method`:

.. code-block:: python

   # energy truncation: no star ever goes beyond 50 kpc; the density drops
   # below nu(r) near the edge
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8, r_trunc=50.0)

   # tapered profile: follows nu(r) exp(-(r/50)^4), with a little mass beyond 50 kpc
   nu_t = stracers.truncate(nu, r_t=50.0, sharpness=4)
   res = stracers.compute_weights(pos, vel, pot, mass, nu_t, stellar_mass=1e8)

The two can be combined. Tapering also gives profiles with divergent mass,
such as NFW, a finite total mass, so they can be used with
``within="total"``.

Power-law profiles with your own parameters
-------------------------------------------

Three power-law families take user-chosen slopes and break radius. All take an
optional core radius ``r_core`` (:math:`r \to \sqrt{r^2 + r_{core}^2}`).
Slopes are logarithmic slopes, :math:`d\ln\nu/d\ln r`, so they are negative,
as in the literature tables.

.. code-block:: python

   # simple power law (SPL): nu = r^slope
   nu = stracers.get_profile("spl", slope=-3.5, r_core=1.0)

   # broken power law (BPL): sharp break at r_break (or smoothed, smoothness in dex)
   nu = stracers.get_profile("bpl", r_break=20.0, inner_slope=-2.5, outer_slope=-4.5,
                             r_core=1.0, smoothness=0.1)

   # double power law (DPL, Zhao 1996): smooth turnover, sharper for larger alpha
   nu = stracers.get_profile("dpl", r_break=20.0, inner_slope=-2.5, outer_slope=-4.5,
                             alpha=2.0, r_core=1.0)

The same functions are available directly as :func:`stracers.power_law`,
:func:`stracers.broken_power_law` and :func:`stracers.double_power_law`.
The DPL is

.. math::

   \nu(r) = x^{\gamma_{in}} \left[\tfrac{1}{2}(1 + x^{\alpha})\right]^{(\gamma_{out} - \gamma_{in})/\alpha},
   \qquad x = r / r_{break},

normalised to :math:`\nu(r_{break}) = 1`. Hernquist
(:math:`\gamma_{in}=-1, \gamma_{out}=-4, \alpha=1`), NFW (−1, −3, 1) and
Plummer (0, −5, 2) are special cases. A smoothed BPL is a DPL with
:math:`\alpha = 1/\text{smoothness}`. Two cautions apply, as for the Milky Way
fits:

- An inner slope of −3 or steeper needs a core, otherwise the mass diverges at
  the centre.
- A sharp break (BPL with ``smoothness=0``, or a very large ``alpha``) leaves
  a numerical feature near the break.

The total stellar mass (``within="total"``) can be used with these profiles,
as with any other, when their total mass is finite. That requires an outer
slope steeper than −3, and an inner slope shallower than −3 or a core. For
example:

.. code-block:: python

   nu = stracers.get_profile("dpl", r_break=20.0, inner_slope=-2.5, outer_slope=-4.5, alpha=2.0)
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e9,
                                  within="total", rmax=250.0)
   res.stellar_mass            # part of the 1e9 Msun inside 250 kpc

   # outer slope -2.5: infinite total mass, so taper it (or use r_trunc) first
   nu = stracers.truncate(stracers.get_profile("spl", slope=-2.5, r_core=1.0), r_t=150.0)

Any function of radius works as a tracer profile:

.. code-block:: python

   res = stracers.compute_weights(pos, vel, pot, mass, lambda r: r**-1 * (1 + r / 5) ** -3)

Choosing the parameters
-----------------------

``stellar_mass``, ``stellar_fraction``, ``within``
   How much stellar mass to assign: an absolute mass or a fraction of the DM
   mass, either inside ``rmax`` (``within="rmax"``, default) or for the whole
   profile (``within="total"``). See :doc:`implementation` (step 8).
``r_trunc``
   Outer edge of the stellar halo, applied as an energy truncation (see
   above). Leave as ``None`` for an untruncated halo.
``rmax``
   Radius of the region used. Set it to about the virial radius: the method
   assumes equilibrium, and the outskirts of real haloes are not in
   equilibrium. The default (300 kpc) comes from the original notebook.
``center``
   Halo centre, e.g. from a shrinking-sphere method. A wrong centre is the
   easiest way to get wrong weights without any warning.
``dlogr``, ``min_per_bin``
   Radial binning of the averaged potential. The defaults (0.1 dex, 20
   particles) work well from 10⁵ to 10⁷ particles. Use larger bins if
   :math:`\Psi(r)` is reported as non-monotonic.
``n_ebins``
   Number of energy bins (default 100). With fewer than ~10⁵ particles, fewer
   bins reduce the noise in :math:`N(E)`.
``derivative``
   Keep ``"spline"`` unless you want to reproduce the original notebook
   (``"gradient"``).
``boundary_term``
   Keep ``True``. Turning it off reproduces the notebook, and makes the tracer
   density fall off too early near ``rmax``.

Following the stars through a simulation
----------------------------------------

Weights are computed once, normally on the initial conditions, and then
attached to the particle IDs:

.. code-block:: python

   stracers.save_weights("weights.h5", ids, res.weights,
                         stellar_mass=1e8, profile="hernquist", a=2.0, rmax=250.0)

   ids_w, w, attrs = stracers.load_weights("weights.h5")
   ids_later = readsnap(snap_later, "pid", "dm")
   w_later = stracers.match_weights(ids_w, w, ids_later)   # aligned with ids_later

Particles that are missing from the weights file get weight 0.

Many simulations
----------------

A single call takes about 0.2 s per million particles, so for large sets of
simulations reading the snapshots usually takes longer than computing the
weights. Process the simulations in parallel and record the diagnostics:

.. code-block:: python

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

Checking the result
-------------------

.. code-block:: python

   from stracers.plotting import plot_diagnostics
   fig = plot_diagnostics(res)

This plots :math:`\Psi(r)`, :math:`d^2\nu/d\Psi^2`, :math:`N(E)`,
:math:`g(E)`, :math:`f_*(E)` and :math:`w(E)`. The most direct check is to
compare the weighted density profile with the input :math:`\nu(r)`; see
``examples/hernquist_tracers.ipynb``.

Troubleshooting
---------------

.. _negative-f:

"energy bins have f(E) < 0"
   A few negative bins at the most bound energies are normal: there are few
   particles there. Many negative bins mean the tracer has no isotropic DF in
   this potential, for example a cored profile in a cuspy halo (see the
   An & Evans condition in :doc:`method`). Use a cuspier tracer, or accept that
   the centre will be under-populated.
"psi(r) is not strictly decreasing"
   The averaged potential is noisy or the centre is wrong. Check ``center``,
   increase ``dlogr``, or use ``derivative="spline"``.
"All weights are zero"
   Usually a potential with the wrong sign or zero point (it must be negative
   and go to zero at infinity), or a tracer profile that is zero over the
   selected region.
