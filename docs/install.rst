************
Installation
************

``stracers`` needs Python 3.9 or later, NumPy and SciPy. Install it from the
GitHub repository::

    git clone https://github.com/jngaravitoc/stracers.git
    cd stracers
    python -m pip install .

Use ``python -m pip install -e .`` instead for an editable install, to work on
the code.

Optional dependencies are installed with extras:

=============  =====================================================================
Extra          Installs
=============  =====================================================================
``io``         h5py, to save and load weights
``plot``       matplotlib, for the diagnostic plots
``test``       pytest and h5py, to run the tests
``docs``       Sphinx and the extensions used to build this documentation
=============  =====================================================================

For example::

    python -m pip install -e ".[io,plot,test,docs]"

Running the tests
=================

From the top of the repository::

    pytest

Building the documentation
==========================

From the ``docs/`` folder::

    make html

The pages are written to ``docs/_build/html``; open ``index.html`` in a browser.
The examples in the pages are tested with::

    pytest --doctest-plus --doctest-rst docs/
