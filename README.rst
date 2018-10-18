Notebooks
=========

Here exists a repo containing Jupyter Notebooks that have been used in
the development of code for TCRM.

.. note:: We apply the `nbstripout <http://github.com/kynan/nbstripout>`_ package to clean output from
   all notebooks prior to committing to the repo. This makes it easier
   to manage the content of the notebooks, but does mean that we can't
   use the notebooks to display content on the GitHub repo.

Dependencies
------------

Some of these notebooks depend on modules in the TCRM
repository. Ensure that you have set the environment variables
correctly (see the `User
Guide <http://geoscienceaustralia.github.io/tcrm/docs/install.html#setting-the-environment>`_
for details).

External dependencies:
......................

* corner
* emcee
* ipython
* ipywidgets
* jupyter
* lmfit
* nbconvert
* notebook
* pandas
* pymc
* statsmodels

.. warning:: Some of these notebooks may be transitioned to Python 3.x
   in the near future. This may make them incompatible with
   Python 2. We will try to flag this at the top of each
   notebook as the changes are made. Once a critical mass is
   reached that have been transitioned, we will flag those
   that have not been converted.
