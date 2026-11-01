.. module:: stracers

********
stracers
********

``stracers`` paints stellar haloes onto the dark matter particles of N-body
simulations. It gives each dark matter particle a weight, the amount of stellar
mass it carries, following the particle-tagging scheme of `Laporte et al.
(2013) <https://arxiv.org/abs/1301.5319>`_. The weighted particles follow a
chosen stellar density profile while staying in equilibrium in the halo
potential, including in triaxial haloes. The weights belong to particle IDs, so
the stars can be followed through any later snapshot.

``stracers`` works on plain NumPy arrays (positions, velocities, potential and
mass), so it can be used with any simulation code. It also provides the Milky
Way stellar halo density profiles compiled by Medina et al. (2024), and star
counts for a halo or a survey.

``stracers`` is developed in a `public repository on GitHub
<https://github.com/jngaravitoc/tracers_triaxial>`_. Bug reports and feature
requests are welcome as `GitHub issues
<https://github.com/jngaravitoc/tracers_triaxial/issues>`_.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   install
   getting_started
   tutorials
   user_guide
   changes
