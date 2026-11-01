.. module:: stracers.kinematics

***************************************************
Kinematics of the tracers (`stracers.kinematics`)
***************************************************

Introduction
============

The stars painted onto dark matter particles move with them, so the
kinematics of the stellar halo are weighted averages over the particles, with
the weights :math:`w_i` of `~stracers.compute_weights`:

.. math::

   \bar v = \frac{\sum_i w_i v_i}{\sum_i w_i}, \qquad
   \sigma^2 = \frac{\sum_i w_i (v_i - \bar v)^2}{\sum_i w_i}.

The stellar masses ``weights * mass`` can be used instead: the normalisation
does not matter. `stracers.kinematics` provides:

=================================================  ===========================================================
Function                                           Returns
=================================================  ===========================================================
`~stracers.kinematics.spherical_velocities`        :math:`v_r`, :math:`v_\theta`, :math:`v_\phi` of each particle
`~stracers.kinematics.bulk_velocity`               the mean velocity vector of the stellar halo
`~stracers.kinematics.weighted_mean`,              weighted mean and dispersion of any quantity
`~stracers.kinematics.weighted_dispersion`
`~stracers.velocity_profile`                       radial profiles of mean velocities, dispersions and
                                                   :math:`\beta`
`~stracers.kinematics.anisotropy`                  :math:`\beta = 1 - (\sigma_\theta^2 + \sigma_\phi^2)/2\sigma_r^2`
=================================================  ===========================================================

Spherical components follow the physics convention: :math:`\theta` is the
polar angle from the :math:`z` axis and :math:`\phi` the azimuth. The
:doc:`../tutorials` include a worked example with plots.

The examples on this page use the stellar halo of :doc:`../getting_started`, a
Hernquist profile with a scale radius of 2 kpc painted onto an isotropic
Hernquist halo of 10\ :sup:`10` M☉::

    >>> import numpy as np
    >>> import stracers
    >>> from stracers.halos import sample_hernquist
    >>> pos, vel, pot = sample_hernquist(100_000, 1e10, 10.0, seed=1)
    >>> mass = 1e10 / 100_000
    >>> nu = stracers.get_profile("hernquist", a=2.0)
    >>> w = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8).weights

Velocities of the tracers
=========================

Each particle's stars share its velocity. `~stracers.kinematics.spherical_velocities`
gives its components, and `~stracers.kinematics.bulk_velocity` the mean
velocity of the stellar halo, here at rest to within the noise (a few tenths
of a km/s for dispersions of ~15 km/s)::

    >>> from stracers.kinematics import spherical_velocities, bulk_velocity, weighted_mean
    >>> v_r, v_theta, v_phi = spherical_velocities(pos, vel)
    >>> print(np.round(bulk_velocity(vel, w), 1))
    [-0.2  0.1 -0.1]
    >>> print(f"{weighted_mean(v_r, w):.1f}")
    0.1

Velocity dispersion profiles
============================

`~stracers.velocity_profile` bins the particles in radius and returns a
`~stracers.kinematics.VelocityProfile` with, in each bin, the weighted mean
velocities, the dispersions :math:`\sigma_r`, :math:`\sigma_\theta`,
:math:`\sigma_\phi` and :math:`\sigma_{3D}`, the anisotropy :math:`\beta`, the
number of particles and the effective number of particles,
:math:`N_{\rm eff} = (\sum w)^2/\sum w^2`. The statistical uncertainty of a
dispersion is about :math:`\sigma/\sqrt{2 N_{\rm eff}}`::

    >>> prof = stracers.velocity_profile(pos, vel, w, bins=[1, 3, 10, 30, 100])
    >>> print(np.round(prof.sigma_r, 1))
    [16.  17.4 15.7 11.6]
    >>> print(np.round(prof.sigma_r / np.sqrt(2 * prof.n_eff), 2))   # uncertainty
    [0.2  0.1  0.07 0.05]
    >>> print(np.round(prof.beta, 2) + 0.0)
    [-0.02  0.    0.    0.  ]

Without weights, the profile is that of the dark matter, which is hotter than
this compact stellar halo::

    >>> print(np.round(stracers.velocity_profile(pos, vel, bins=[1, 3, 10, 30, 100]).sigma_r, 1))
    [20.8 20.6 17.  11.9]

``bins`` can also be a number of logarithmic bins between ``rmin`` and
``rmax``. Pass the position and velocity of the halo centre as ``center`` and
``vcenter`` if the particles are not already centred.

Validation
==========

For an isotropic tracer in a spherical potential the Jeans equation gives
:math:`\sigma_r^2(r) = \nu(r)^{-1} \int_r^\infty \nu\, (d\Phi/dr')\, dr'`. In
the test halo the weighted dispersions of Hernquist, Plummer and Milky Way
(Medina et al. 2024) stellar haloes follow it to within a few percent from 1 to
100 kpc, with :math:`\beta \approx 0`: the weights give the stars the right
dynamics, not only the right density. The tests also check the velocity
components, the weighted statistics, the recovery of a rotation, and the
``center`` and ``vcenter`` options.

API
===

.. automodapi:: stracers.kinematics
   :no-heading:
   :no-main-docstr:
   :no-inheritance-diagram:
   :headings: =-
