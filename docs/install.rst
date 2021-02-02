.. _installation:

Installation
============

Installing TCRM is intended to be a simple process, requiring only basic
understanding of command line operations. TCRM has been installed and 
tested on a range of unix-based systems, Windows and Mac OS/X systems.

.. _downloading:

Downloading
-----------

The TCRM code can be downloaded from the `Geoscience Australia GitHub
page <https://github.com/GeoscienceAustralia/tcrm>`_.

For users wanting to only run the code, a zip file or gzipped tar file
of the latest releases can be downloaded from the `Releases page
<https://github.com/GeoscienceAustralia/tcrm/releases>`_.

To have access to the latest updates, users should clone the repository, and
then regularly pull from the repository as updates are made. 

Those wanting to contribute to development can `fork
<https://github.com/GeoscienceAustralia/tcrm/fork>`_ the
repository. Submit a pull request to have your changes integrated into
TCRM. Read more about contribting_ to the TCRM code.

.. _dependencies:

Dependencies
------------

TCRM relies on a number of additional libraries that are not part of
the standard Pyhton library. There are several ways to obtain the required
libraries -- using Python's recommended tool `pip
<https://pip.readthedocs.org/en/latest/>`_, installing a distribution
such as `Python(x,y) package <http://code.google.com/p/pythonxy/>`_
(for Windows environments) or `Anaconda
<https://www.anaconda.com/distribution/#download-section>`_ (cross-platform), or
installing the libraries from source or binary installers
(pre-compiled binary Windows installer versions for all the libraries
(both 32-bit and 64-bit) can be obtained `here
<http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_).

For detailed instructions on installation of these dependencies,
please see the documentation for each individual library.

* `Python <https://www.python.org/>`_ - v3.5 or later
* `Numpy <http://www.numpy.org/>`_ - v1.6 or later
* `Scipy <http://www.scipy.org/>`_ - v0.12 or later
* `Matplotlib <http://matplotlib.org/>`_ v1.2 or later. 
* `Basemap <http://matplotlib.org/basemap/index.html>`_
* `netcdf4-python <https://code.google.com/p/netcdf4-python/>`_ -
  version 1.0.8 or later
* `Shapely <http://toblerity.org/shapely/index.html>`_ - v1.2.15 or later
* `statsmodels <http://statsmodels.sourceforge.net>`_ 
* `seaborn <http://seaborn.pydata.org>`_
* `pandas <http://pandas.pydata.org>`_
* `gitpython <http://gitpython.readthedocs.org>`_
* Parallel execution in multi-processor environments (with MPI
  installed) requires `mpi4py <https://mpi4py.readthedocs.io/>`_


Using Anaconda
~~~~~~~~~~~~~~

To install ``tcrm``, make a new environment:

.. code-block:: bash

    conda env create -f tcrmenv.yml

After creating the environment the user needs to move to that environment using the command

.. code-block:: bash

     conda activate tcrm

The bash promt will look like

.. code-block::

    (tcrm) user@server:~/tcrm$

Using pip
~~~~~~~~~

If you have `pip <https://pip.readthedocs.org/en/latest/>`_ installed,
the required modules can be installed using the following command,
executed in the main TCRM directory

.. code-block:: bash

    pip -v install -r requirements.txt

This will automatically build the required libraries (listed in the
``requirements.txt`` file) and any dependencies. ``pip`` must be on
the ``$PATH`` for this to work.



.. _environment:

Setting the environment
-----------------------

To enable TCRM to run correctly, you may need to change some
environment settings. The important variable to set is the
``PYTHONPATH`` variable. This should be set to the path where you have
extracted the contents of the zip file. In the examples below, change
``/path/to/tcrm`` to the location where you extracted the TCRM files.

A complete discussion on environment variables in Python is given in
the `Python documentation
<https://docs.python.org/2/using/cmdline.html#environment-variables>`_.


Windows
~~~~~~~

The Python documentation contains some simple instructions for setting
environment variables on Windows systems `here
<https://docs.python.org/2/using/windows.html>`_. See `this link
<http://www.computerhope.com/issues/ch000549.htm>`_ for setting the
variables on different Windows systems.

BASH shell
~~~~~~~~~~

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/path/to/tcrm:/path/to/tcrm/Utilities


CSH/TCSH shell
~~~~~~~~~~~~~~

.. code-block:: tcsh

    setenv PYTHONPATH $PYTHONPATH:/path/to/tcrm:/path/to/tcrm/Utilities


.. _compilation:

Windows
~~~~~~~

For Windows users, the code includes the ``compile.cmd`` script in the
main TCRM diretory that will build these extensions in place. By
default, TCRM uses the MinGW suite (http://www.mingw.org) for
compiling the extensions. Other Windows-based packages can also be
used (e.g. Cygwin). See the Python documentation on writing
configuration files for the :mod:`distutils` package for more details.

Notes
~~~~~

It is recommended to use a stand-alone Python installation for
compiling and running TCRM. Installations linked to other software
such as ArcGIS have resulted in compilation errors, as the required
:mod:`numpy` libraries are pre-compiled and packaged with such
installations.

.. _testing:

Testing the installation
------------------------

The model code includes a suite of unit tests that ensure elements of
the code base will work as expected, even if a user makes
modificaitons to the code.

The test suite can be run from the main directory. On Windows, run the
``run_test_all.cmd`` script from the main TCRM directory. On Unix, use
the command

.. code-block:: bash

    python ./tests/run.py

This should report no errors or failures. 

Special note for Windows systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On a Windows system, :func:`tests.test_files.testModulePath` may fail
due to the different path separators (``/`` versus ``\\``) used by the
Windows system. This test failure will appear as::

    ======================================================================
    FAIL: testModulePath (tests.test_files.TestModuleUtilities)
    Test flModulePath returns correct path, base & extension
    ----------------------------------------------------------------------
    Traceback (most recent call last):
      File "tcrm\tests\test_files.py", line 22, in testModulePath
        self.assertEqual(self.path, p)
    AssertionError: 'tcrm/tests' != 'tcrm\\tests'

    ---------------------------------------------------------------------- 
    Ran 111 tests in 92.513s

    FAILED (failures=1)

Such an error will not affect model execution.


Using Docker
------------

As an alternative way to install TCRM, you can use Docker.
Docker is a very convenient way to run containerized software which
avoids all the hassle with compilation or dependencies.

Prerequisites
~~~~~~~~~~~~~

Install `Docker Community Edition
<https://docs.docker.com/install/#supported-platforms>`_ for your
system.

Test the installation
~~~~~~~~~~~~~~~~~~~~~

Run this command

.. code-block:: bash

    docker run olivierdalang/tcrm nosetests --exe

The first time, this will take some time, as it needs to download the docker image.
If it works, you should see (after some time), something like ``OK (SKIP=1)``.
If not, you would see something like ``FAILED (SKIP=1, errors=1)``.

Normal usage
~~~~~~~~~~~~

To run TCRM though Docker, you need to mount a folders containing your
inputs and the output folder in the container.

This can be done like this (assuming you have a my_conf.ini file in
a folder)

.. code-block:: bash

    docker run -v /path_to/my_data_folder:/home/src/mount -v /path_to/my_output_folder:/home/src/output olivierdalang/tcrm python tcevent.py -v -c mount/my_conf.ini

Replace ``/path_to/my_data_folder`` and ``/path_to/my_output_folder``
by the folders you want to use on your system, and ``python tcevent.py 
v -c example/yasi.ini`` by the TCRM command you want to use.

The first time, the docker image will have to be downloaded which will
take some time.

Developement
~~~~~~~~~~~~

You can also use Docker when developping TCRM by mounting the source

.. code-block:: bash

    git checkout https://github.com/GeoscienceAustralia/tcrm.git
    cd tcrm
    docker run -v ${PWD}:/home/src olivierdalang/tcrm python tcevent.py -c example/yasi.ini

If you wish to make changes to the builds steps or dependencies, you need to rebuild the image locally

.. code-block:: bash

    docker build -t olivierdalang/tcrm .

Releases
~~~~~~~~

For users to be able to use the docker image out of the box without having to rebuild it locally,
the image must be pushed to the docker hub repository like this

.. code-block:: bash

    docker build -t olivierdalang/tcrm .
    docker login
    docker push olivierdalang/tcrm

This can be setup to be done automatically after pushes through docker hub.
