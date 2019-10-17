Notebooks
=========

Here exists a repo containing Jupyter Notebooks that have been used in
the development of code for TCRM.

.. note:: 
    We apply the `nbstripout <http://github.com/kynan/nbstripout>`_ package to clean output from
    all notebooks prior to committing to the repo. This makes it easier
    to manage the content of the notebooks, but does mean that we can't
    use the notebooks to display content on the GitHub repo.
   
Installing
----------

Clone this branch of the repository to a folder separate from your TCRM installation::

    git clone git@github.com/GeoscienceAustralia/tcrm.git --branch notebooks --single-branch

Alternatively, download the zip file containing the notebooks from `here <https://github.com/GeoscienceAustralia/tcrm/archive/notebooks.zip>`_ and extract the files.


Dependencies
------------

Some of these notebooks depend on modules in the TCRM
repository (`v2.1 <https://github.com/GeoscienceAustralia/tcrm/tree/v2.1>`_ ). Ensure that you have set the environment variables
correctly (see the `User
Guide <http://geoscienceaustralia.github.io/tcrm/docs/install.html#setting-the-environment>`_
for details).

Generally, the following commands will set the required environment variables::

    $ export PYTHONPATH=$PYTHONPATH:<path/to>/tcrm

where `<path/to>` is the path to the location where you installed TCRM. For example, I generally install TCRM directly into my home directory, so the following command works for me::

    $ export PYTHONPATH=$PYTHONPATH:$HOME/tcrm


External dependencies:
......................

External dependencies are listed in `requirements.txt`

We recommend using a virtual environment to run the notebooks, so the dependencies can be readily managed - there may be some conflicts with Python versions, expecially as we make the transition across to Python 3.6. 

To install the required modules::

    $ virtualenv notebooks
    $ source notebooks/bin/activate
    (notebooks)$ cd <path/to>/notebooks
    (notebooks)$ pip install -r requirements.txt



.. WARNING::
    Some of these notebooks may be transitioned to Python 3.x
    in the near future. This may make them incompatible with
    Python 2. We will try to flag this at the top of each
    notebook as the changes are made. Once a critical mass is
    reached that have been transitioned, we will flag those
    that have not been converted.
