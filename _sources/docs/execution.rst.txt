.. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

.. _execution:
===================
 Running the model
===================

The primary way to run TCRM is at the command line. Command line arguments 
are passed to the main ``tcrm.py`` script for defining the path 
to the configuration file, as well as enabling verbose output and/or debugging.

Command line arguments
======================

 -c file, --config file   Path to a configuration file
 -v, --verbose            If given, then logging messages will be printed to the console
 -d, --debug              In the case that execution results in an exception, allow the 
                          Python stack to call into the stack trace (through 
                          implementation of a custom hook script). 

Examples
========

Make sure ``python`` is in your system path, then from the base
directory, call the ``tcrm.py`` script, with the configuration file
option included. For example, to run the example simulation::

    python tcrm.py -c example/port_hedland.ini

The model will print a simple progress indicator to the console to
show that the model is working.

For Python3 users, some may have to use the ``python3`` command::

    python3 tcrm.py -c example/port_hedland.ini

If the ``-v`` option is included, then all logging messages will be
printed to the console. The level of logging detail is set in the
configuration file, using Python's inbuilt logging levels ('DEBUG', 'INFO',
'WARNING', etc.).

Running on a parallel system
============================

As a stochastic model, TCRM can generate massive numbers of synthetic
events, which implies run times can become very long. If TCRM is
installed on a multiprocessor system (either a shared memory or
distributed memory system), then the workload can be shared across the
workload. TCRM uses the :mod:`mpi4py` module to enable multi-threaded
processing. Instructions for installing :mod:`mpi4py` are given on the
`mpi4py website <https://mpi4py.readthedocs.io/>`_.

``mpi4py`` is built around the MPI library, and so uses the ``mpirun``
command to execute the model across multiple processors. As an example
the following command will execute across 16 processors::

    mpirun -np 16 python tcrm.py -c example/port_hedland.ini

Running across multiple processors means that logging messages from
each individual processor can get mixed up with others. To avoid this,
a separate log file is created for each thread, and output to the
console is suppressed (even if the ``-v`` or ``--verbose`` option is
given at the command line).
