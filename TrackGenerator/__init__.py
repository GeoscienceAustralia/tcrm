"""
:mod:`TrackGenerator` -- Tropical cyclone track generation
==========================================================

This module contains the core objects for tropical cyclone track
generation.

Track generation can be run in parallel using MPI if the :term:`mpi4py`
library is found and TCRM is run using the :term:`mpirun` command. For
example, to run with 10 processors::

    $ mpirun -n 10 python main.py cairns.ini

:class:`TrackGenerator` can be correctly initialised and started by
calling :meth:`run` with the location of a configuration file::

    >>> import TrackGenerator
    >>> TrackGenerator.run('cairns.ini')

Alternatively, it can be run from the command line::

    $ python TrackGenerator.py cairns.ini

Note the absence of a '-c' flag (c.f. the normal TCRM command line arguments).
"""

import os
import sys

from .TrackGenerator import run
from Utilities.parallel import attemptParallel

if __name__ == "__main__":

    global MPI
    MPI = attemptParallel()
    import atexit
    atexit.register(MPI.Finalize)
    try:
        configFile = sys.argv[1]
    except IndexError:

        configFile = __file__.rstrip('.py') + '.ini'

        if not os.path.exists(configFile):
            errorMsg = 'No configuration file specified, please' + \
                       ' type: python main.py {config filename}.ini'
            raise IOError(errorMsg)

    if not os.path.exists(configFile):
        errorMsg = "Configuration file '" + configFile + "' not found"
        raise IOError(errorMsg)

    run(configFile)
    MPI.Finalize()
