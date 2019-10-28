Notebooks
=========

Here exists a repo containing Jupyter Notebooks that have been used in
the development of code for TCRM.

.. note:: 
    We apply the `nbstripout <http://github.com/kynan/nbstripout>`_ package to clean output from
    all notebooks prior to committing to the repo. This makes it easier
    to manage the content of the notebooks, but does mean that we can't
    use the notebooks to display content on the GitHub repo.



Dependencies
------------

External dependencies are listed in `requirements.txt`

We recommend using a virtual environment to run the notebooks, so the dependencies can be readily managed - there may be some conflicts with Python versions, expecially as we make the transition across to Python 3.6. 

Using Python3 virtualenv:
~~~~~~~~~~~~~~~~~~~~~~~~~

To install the required modules::

    $ virtualenv notebooks
    $ source notebooks/bin/activate
    (notebooks)$ cd <path/to>/notebooks/tcrm
    (notebooks)$ pip install -r requirements.txt
    

Using conda environments:
~~~~~~~~~~~~~~~~~~~~~~~~~

A conda environment file is also included (`notebooks.yml`). This can be used to build a conda environmnent (for Python3) using the following steps::

    $ conda env create -f notebooks.yml
    $ conda activate notebooks
    (notebooks)$ 
    
Note you may wish to modify the last line of `notebooks.yml` to install the environment in a location different to the default (currently set to the home directory).

See Anaconda's pages on `managing environments <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ for more details.

Installing
----------

Clone this branch of the repository to a folder separate from your TCRM installation::

    $ cd $HOME
    $ mkdir notebooks
    $ cd notebooks
    $ git clone git@github.com/GeoscienceAustralia/tcrm.git --branch notebooks --single-branch tcrm

Alternatively, download the zip file containing the notebooks from `here <https://github.com/GeoscienceAustralia/tcrm/archive/notebooks.zip>`_ and extract the files.

As some of the notebooks rely on modules in the TCRM code, you'll need to download the`[TCRM code] <https://github.com/GeoscienceAustralia/tcrm>`_ (check the branch, as some of the notebooks are Python3 code) and add the path to the code to the `PYTHONPATH`::

    $ cd $HOME
    $ git clone git@github.com/GeoscienceAustralia/tcrm.git --branch py3 --single-branch tcrm
    $ export PYTHONPATH=$PYTHONPATH:$HOME/tcrm
    
where `<path/to>` is the path to the location where you installed TCRM. For example, I generally install TCRM directly into my home directory, so the following command works for me::

    $ export PYTHONPATH=$PYTHONPATH:$HOME/tcrm 

Starting the notebook server
----------------------------

Starting the Jupyter notebook server:: 

    (notebooks)$ python <path/to>/anaconda3/Scripts/jupyter-notebook-script.py $HOME/notebooks
    
The `(notebooks)` before the prompt indicates that you are running in a virtual or conda environment.

See the `Jupyter documentation <https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/execute.html#>`_ for more details on starting the notebook server.

Warning
-------
Some of these notebooks have been transitioned to Python 3, and others will be transitioned in the near future. This may make them incompatible with Python 2. We will try to flag this at the top of each notebook as the changes are made. Once a critical mass is reached that have been transitioned, we will flag those that have not been converted.
