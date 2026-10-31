stracers
========

``stracers`` assigns stellar-tracer weights to the dark matter (DM) particles of
an N-body halo, using the particle-tagging scheme of
`Laporte et al. (2013) <https://arxiv.org/abs/1301.5319>`_.

Every DM particle stands for both dark matter and stars, in different amounts.
The amount of stars is set by a weight :math:`w(E)` that depends only on the
particle's energy. The weights are chosen so that the weighted particles follow
a chosen stellar density profile :math:`\nu(r)`, for example a Hernquist
profile, while staying in equilibrium in the halo potential. Weights belong to
particle IDs, so the stars can be followed through any later snapshot of the
simulation.

.. code-block:: python

   import stracers

   nu = stracers.get_profile("hernquist", a=2.0)
   res = stracers.compute_weights(pos, vel, pot, mass, nu, stellar_mass=1e8, rmax=250)
   m_star = res.weights * mass          # stellar mass carried by each DM particle

.. toctree::
   :maxdepth: 2

   method
   implementation
   usage
   api
   references
