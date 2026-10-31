How the code works
==================

Everything happens in :func:`stracers.compute_weights`. It calls a few small
functions in :mod:`stracers.weights`, one per step of the method. This page goes
through the steps in order, gives the function that implements each one, and
explains the numerical choices.

Inputs
------

The code works on plain arrays, so it does not depend on the snapshot format:

====================  ==============================================================
``pos``, ``vel``      (N, 3) positions [kpc] and velocities [km/s]
``pot``               (N,) potential of each particle [(km/s)²], zero at infinity
``mass``              particle mass [M☉] (scalar, or (N,) array of equal masses)
``tracer_profile``    a function ν(r); its normalisation does not matter
====================  ==============================================================

The gravitational constant is only needed for the optional
``truncation_correction``. It defaults to
:math:`G = 4.30071\times10^{-6}` kpc (km/s)² / M☉.

Pipeline
--------

1. Select particles
~~~~~~~~~~~~~~~~~~~

The positions are shifted by ``center`` (if given), and only particles with
:math:`0 < r < r_{max}` are kept. ``rmax`` should cover the virialised region.
The energy of each particle uses **its own potential**,
:math:`E_i = \tfrac{1}{2}|\mathbf{v}_i|^2 + \Phi_i`, which is the triaxial
prescription of Laporte et al.

2. Spherically averaged potential — :func:`~stracers.weights.radial_potential`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Particles are binned in :math:`\log_{10} r` with width ``dlogr`` (default 0.1
dex). In each bin the code computes the mean potential
:math:`\langle\Phi\rangle` and its standard error. Bins with fewer than
``min_per_bin`` particles (default 20, as in the original notebook) are
dropped. Bin centres are the midpoints in :math:`\log r`.

3. Smoothing — :func:`~stracers.weights.smooth_potential` (``derivative="spline"``, default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eddington's formula needs the *second* derivative :math:`d^2\nu/d\Psi^2`.
Differentiating twice amplifies any noise in the binned potential. The code
therefore fits a quartic smoothing spline to :math:`\Phi(\ln r)`, weighted by
the inverse standard errors with smoothing factor equal to the number of bins
(about one :math:`\chi^2` unit per bin). The spline is then resampled on
``n_fine`` (default 1000) logarithmically spaced radii. All later radial
quantities, including :math:`g(E)`, use this grid.

With ``derivative="gradient"`` the spline is skipped and the binned potential
is used directly, as in the original notebook.

4. Tracer derivatives — :func:`~stracers.weights.tracer_derivatives`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`\Psi = -\langle\Phi\rangle` and :math:`\nu` evaluated on the
radial grid, :math:`d\nu/d\Psi` and :math:`d^2\nu/d\Psi^2` are computed with
second-order finite differences (:func:`numpy.gradient`). If
:math:`\Psi(r)` is not strictly decreasing, the code warns.

5. Energy grid and :math:`N(E)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One linear grid of ``n_ebins`` bins (default 100) spans the energies of the
selected particles. :math:`N(E)`, :math:`g(E)` and :math:`f_*(E)` are all
evaluated at the centres of the same bins, so they line up exactly.
:math:`N(E)` is the particle mass in each bin divided by the bin width.

6. Density of states — :func:`~stracers.weights.density_of_states`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eq. :eq:`gE` is evaluated as a sum over the radial grid,

.. math::

   g(E) \simeq (4\pi)^2 \left[\frac{r_{in}^3}{3}\sqrt{2(E-\Phi_1)}
   + \sum_{k:\,\Phi_k<E} r_k^2\, \Delta r_k \sqrt{2(E - \Phi_k)}\right].

The first term is the sphere inside the first radial bin, taking
:math:`\Phi = \Phi_1` there. It matters for the most bound energies. The sum
stops at the outermost bin, not at :math:`r_E`. This is deliberate:
:math:`N(E)` only counts particles inside ``rmax``, so it equals
:math:`f_{DM}(E)` times the density of states *truncated at rmax*. Using the
same truncation for :math:`g` keeps :math:`f_{DM} = N/g` unbiased.

7. Eddington inversion — :func:`~stracers.weights.eddington` and :func:`~stracers.weights.outer_eddington`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`f_*` is evaluated at :math:`\mathcal{E} = -E` for each energy-bin
centre. The integral in Eq. :eq:`eddington` is split in two:

- **Inside the data**, :math:`\Psi_{min} \le \Psi < \mathcal{E}`, where
  :math:`\Psi_{min} = \Psi(r_{max})`. :math:`d^2\nu/d\Psi^2` is treated as
  constant on each segment between grid points, and
  :math:`1/\sqrt{\mathcal{E}-\Psi}` is integrated exactly on each segment.
  This removes the integrable singularity at :math:`\Psi = \mathcal{E}`, which
  makes a plain Riemann sum noisy.
- **Beyond the data**, :math:`0 \le \Psi < \Psi_{min}`, i.e.
  :math:`r > r_{max}`, plus the boundary term at :math:`\Psi = 0`. Here
  :math:`\nu(\Psi)` is extrapolated as a power law
  :math:`\nu \propto \Psi^k`, with :math:`k` measured at :math:`\Psi_{min}`,
  and the integral is done in closed form with the incomplete beta function.
  This is ``boundary_term=True`` (the default). It is needed to recover the
  tracer profile near ``rmax``, and it gives non-zero :math:`f_*` to weakly
  bound particles with :math:`\mathcal{E} < \Psi_{min}`.

Bins with a negative :math:`f_*` are set to zero. The count is stored in
``n_negative_f`` and a warning is issued (see :ref:`negative-f`).

Energy truncation (``r_trunc``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``r_trunc`` is given, :math:`\Phi_t = \langle\Phi\rangle(r_t)` is
interpolated from the radial grid. :math:`f_*` is set to zero in the energy
bins with :math:`E > \Phi_t`, and so is the weight of every *particle* with
:math:`E > \Phi_t`, so that the bin containing :math:`\Phi_t` is also cut
exactly. ``r_trunc`` must lie inside the radial range of the selected
particles (and so be :math:`\le` ``rmax``). With ``within="total"`` no
:math:`F(r_{max})` correction is applied, because all the stars are inside
``rmax``.

8. Weights
~~~~~~~~~~

:math:`w(E) = f_* g / N` is computed in each energy bin (zero where
:math:`N = 0`). Each particle gets the weight of its bin, found with
:func:`numpy.digitize`. The weights are then rescaled so that

.. math::

   \sum_i w_i\, m_i = M_*\ (\texttt{stellar\_mass}),

where :math:`M_*` is the stellar mass inside ``rmax``. Particles that were not
selected get :math:`w = 0`. :math:`M_*` can be given in four ways:

========================================  ==================================================================
``stellar_mass=M``                        :math:`M_* = M` inside ``rmax`` (default ``within="rmax"``)
``stellar_fraction=f``                    :math:`M_* = f\, M_{DM}(<r_{max})`
``stellar_mass=M, within="total"``        :math:`M_* = M\, F(r_{max})`
``stellar_fraction=f, within="total"``    :math:`M_* = f\, M_{DM,\,all}\, F(r_{max})`
========================================  ==================================================================

Here :math:`F(r_{max}) = M_\nu(<r_{max}) / M_\nu(<\infty)` is the fraction of
the tracer profile's mass inside ``rmax``
(:func:`~stracers.profiles.mass_fraction_within`), and :math:`M_{DM,\,all}` is
the mass of all particles passed in. With ``within="total"`` the stars beyond
``rmax`` exist in the model but are not carried by any particle. This option
fails for profiles whose mass diverges (e.g. NFW). The stellar mass actually
assigned is stored in ``WeightsResult.stellar_mass``. The result is a :class:`~stracers.WeightsResult`
holding the weights and every intermediate quantity, for diagnostics
(:func:`stracers.plotting.plot_diagnostics`).

Differences from the paper and the original notebook
----------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Item
     - Laporte+13
     - Original notebook
     - ``stracers``
   * - Φ₀
     - 0
     - shifted by :math:`G M / r_{max}`
     - 0. The shift is available as ``truncation_correction``
   * - :math:`d^2\nu/d\Psi^2`
     - not specified
     - ``np.gradient`` on binned Φ
     - spline-smoothed Φ (default), or ``"gradient"``
   * - Eddington integral
     - full, incl. boundary term
     - Riemann sum over the data only
     - exact kernel per segment, plus power-law extension beyond ``rmax``
   * - :math:`g(E)` lower limit
     - :math:`r = 0`
     - first radial bin
     - :math:`r = 0` (inner sphere added)
   * - Mapping :math:`f(\mathcal{E})` to E bins
     - :math:`\mathcal{E} = -E`
     - sorted :math:`f` (not equivalent)
     - evaluated at :math:`\mathcal{E} = -E_{bin}`
   * - Normalisation
     - not specified
     - :math:`\sum w m = M_*`
     - :math:`\sum w m = M_*`

The notebook's :math:`1/\sqrt{2(\mathcal{E}-\Psi)}` kernel differs from
Eddington's by a constant factor, which cancels in the normalisation.

Validation
----------

The test suite (``tests/``) samples an isotropic Hernquist halo (:math:`10^{10}`
M☉, :math:`a = 10` kpc, :math:`2\times10^5` particles) from its analytic DF,
using the analytic potential, and checks that:

- :math:`f_*` for the halo's own profile matches the analytic Hernquist DF to
  within 3% (in practice ~1%), and within 15% in the extrapolated outermost
  bins;
- :math:`g(E)` matches Eq. :eq:`gE` integrated numerically to within 3%;
- self-consistent weights (tracer = halo) are constant to within 5%;
- Hernquist (:math:`a=2` kpc), Plummer (:math:`a=3` kpc) and Einasto tracers
  are recovered to within 10% between 0.5 and 100 kpc;
- the energy uses each particle's own potential, recentring works, and the
  normalisation and masking behave as described.

Performance
-----------

The cost is dominated by a few passes over the particle arrays, so it grows
linearly with :math:`N`: about 0.02 s for :math:`10^5` particles, 0.17 s for
:math:`10^6` and 3 s for :math:`10^7` on a desktop CPU. The binned steps
(spline, Eddington, :math:`g(E)`) take about 1 ms in total.
