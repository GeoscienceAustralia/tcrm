.. _installation:

Installation
============

Installing TCRM is intended to be a simple process, requiring only a
small amount of compilation and basic understanding of command line
operations. TCRM has been installed and (lightly) tested on a range of
unix-based systems, Windows and Mac OS/X systems.

.. _downloading:

Downloading
-----------

The TCRM code can be downloaded from the `Geoscience Australia GitHub
page <https://github.com/GeoscienceAustralia/tcrm>`_.

For users wanting to only run the code, a zip file or gzipped tar file
of the latest releases can be downloaded from the `Releases page
<https://github.com/GeoscienceAustralia/tcrm/releases>`_.

Those wanting to contribute to development can `fork
<https://github.com/GeoscienceAustralia/tcrm/fork>`_ the
repository. Submit a pull request to have your changes integrated into
TCRM.

.. _environment:

Setting the environment
-----------------------

To enable TCRM to run flawlessly, you may need to change some environment settings. The important variable to set is the ``PYTHONPATH`` variable. This should be set to the path where you have extracted the contents of the zip file. In the examples below, change ``/path/to/tcrm`` to the location where you extracted the TCRM files.

A complete discussion on environment variables in Python is given in the `Python documentation <https://docs.python.org/2/using/cmdline.html#environment-variables>`_. 

Windows
~~~~~~~
The Python documentation contains some simple instructions for setting environment variables on Windows systems `here <https://docs.python.org/2/using/windows.html>`_. See `this link <http://www.computerhope.com/issues/ch000549.htm>`_ for setting the variables on different Windows systems.

BASH shell
~~~~~~~~~~

::

    export PYTHONPATH=$PYTHONPATH:/path/to/tcrm:/path/to/tcrm/Utilities


CSH/TCSH shell
~~~~~~~~~~~~~~

::

    setenv PYTHONPATH $PYTHONPATH:/path/to/tcrm:/path/to/tcrm/Utilities





.. _dependencies:

Dependencies
------------

TCRM relies on a number of additional libraries that are not part of
the standard library. There are several ways to obtain the required
libraries -- using Python's recommended tool `pip
<https://pip.readthedocs.org/en/latest/>`_, installing a distribution
such as `Python(x,y) package <http://code.google.com/p/pythonxy/>`_
(for Windows environments) or `Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ (cross-platform), or
installing the libraries from source or binary installers
(pre-compiled binary Windows installer versions for all the libraries
(both 32-bit and 64-bit) can be obtained `here
<http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_).

For detailed instructions on installation of these dependencies,
please see the documentation for each individual library.

* `Python <https://www.python.org/>`_ - v2.7 preferred
* `Numpy <http://www.numpy.org/>`_ - v1.6 or later
* `Scipy <http://www.scipy.org/>`_ - v0.12 or later
* `Matplotlib <http://matplotlib.org/>`_ v1.2 or later. 
* `Basemap <http://matplotlib.org/basemap/index.html>`_
* `netcdf4-python <https://code.google.com/p/netcdf4-python/>`_ -
  version 1.0.8 or later
* `Shapely <http://toblerity.org/shapely/index.html>`_ - v1.2.15 or later
* Parallel execution in multi-processor environments (with MPI
  installed) requires `Pypar <http://github.com/daleroberts/pypar>`_

Using pip
~~~~~~~~~

If you have `pip <https://pip.readthedocs.org/en/latest/>`_ installed,
the required modules can be installed using the following command,
executed in the main TCRM directory::

   pip -v install -r requirements.txt

This will automatically build the required libraries (listed in the
``requirements.txt`` file) and any dependencies. ``pip`` must be on
the ``$PATH`` for this to work.

.. _compilation:

Compiling the extensions
------------------------

The model requires a number of C extensions to be compiled before
execution. These can be built using Python's inbuilt :mod:`distutils`
module.


Unix
~~~~
From the base directory, execute the build process::

    python installer/setup.py build_ext -i

Ubuntu
~~~~~~
The github branch issue_25 (created from the v2.0 branch) had an environment created by `installing miniconda
<https://conda.io/docs/install/quick.html#linux-miniconda-install>`_ and executing the following commands.

        ~/miniconda2/bin/conda create --name tcrm
        ~/miniconda2/bin/source activate tcrm
        ~/miniconda2/bin/conda install numpy
        ~/miniconda2/bin/conda install scipy
        ~/miniconda2/bin/conda install matplotlib
        ~/miniconda2/bin/conda install basemap
        ~/miniconda2/bin/conda install netcdf4
        ~/miniconda2/bin/conda install shapely
        ~/miniconda2/bin/conda install Tornado
        ~/miniconda2/bin/conda install statsmodel
        ~/miniconda2/bin/conda install seaborn
        ~/miniconda2/bin/pip --proxy=http://localhost:3128 install simplejson


The following libraries were needed to compile the C extensions, and run the unit tests:
    sudo apt install libgl1-mesa-glx
    sudo apt-get install python-numpy-dev

The C extensions were compiled from the trcm directory with:
        (tcrm) user@server:~/tcrm$ python intaller/setup.py build_ext -i

An error occurred where the include file seems to have changed paths. It may be a one off,
or it may reoccur in another version of Linux. The error was in KPDF.c and the change was to
comment out one line and replace it with another.

        #include "numpy/arrayobject.h"
        /* #include "arrayobject.h" */

A requiremements file was created in the root directory called ``linux_v20.yml`` and should (it hasn't been tested)
replace the ``conda install`` commands above. The command to use this file is:

        conda env create -f linux_v20.yml

Activating the environment would be
        source activate linux_v20


Windows
~~~~~~~

For Windows users, the code includes the ``compile.cmd`` script in the
main TCRM diretory that will build these extensions in place. By default, TCRM uses the MinGW suite (http://www.mingw.org) for compiling the extensions. Other Windows-based packages can also be used (e.g. Cygwin). See the Python documentation on writing configuration files for the :mod:`distutils` package for more details.

Notes
~~~~~

It is recommended to use a stand-alone Python installation for compiling and running TCRM. Installations linked to other software such as ArcGIS have resulted in compilation errors, as the required :mod:`numpy` libraries are pre-compiled and packaged with such installations. 

.. _testing:

Testing the installation
------------------------

The model code includes a suite of unit tests that ensure elements of
the code base will work as expected, even if a user makes
modificaitons to the code.

The test suite can be run from the main directory. On Windows, run the
``run_test_all.cmd`` script from the main TCRM directory. On Unix, use
the command::

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

