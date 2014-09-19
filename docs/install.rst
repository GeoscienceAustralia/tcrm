.. _installation:
Installation
============

Installing TCRM is intended to be a simple process, requiring only a small amount of compilation and basic understanding of command line operations. TCRM has been installed and (lightly) tested on a range of unix-based systems, Windows and Mac OS/X systems.

.. _dependencies:
Dependencies
------------
For instructions on installation of these dependencies, please see the documentation for each package.

* `Python <https://www.python.org/>`_ - v2.7 preferred
* `Numpy <http://www.numpy.org/>`_ - v1.6 or later
* `Scipy <http://www.scipy.org/>`_ - v0.12 or later
* `Matplotlib <http://matplotlib.org/>`_ v1.2 or later. 
* `Basemap <http://matplotlib.org/basemap/index.html>`_
* `netcdf4-python <https://code.google.com/p/netcdf4-python/>`_ - version 1.0.8 or later
* Parallel execution in multi-processor environments (with MPI installed) requires `Pypar <http://github.com/daleroberts/pypar>`_ 

For Windows environments, we recommend using the `Python(x,y) package
<http://code.google.com/p/pythonxy/>`_, that includes all the required
modules (excluding Basemap and Pypar). Pre-compiled binary installer
versions can be obtained `here
<http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.

.. _compilation:
Compiling the extensions
------------------------

The model requires a number of C extensions to be compiled before
execution. These can be built using Python's inbuilt :mod:`distutils`
module. Copy the required files from the `installer` directory to the
base directory and then execute the build process::

    $ python intaller/setup.py build_ext -i

For Windows users, the model includes the ``compile.cmd`` script that
will build these extensions in place.

.. _testing:
Testing the installation
------------------------

The model code includes a suite of unit tests that ensure elements of
the code base will work as expected, even if a user makes
modificaitons to the code.

The test suite can be run from the main directory::

    $ python ./tests/run.py

This should report no errors. 

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

Such an error will not affect model execution.
