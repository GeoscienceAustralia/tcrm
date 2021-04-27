"""
:mod:`parallel` -- base parallel processing functions
=====================================================

.. module:: parallel
    :synopsis: Provides a couple of base functions that set up a
               parallel processing environment. There are only a 
               small number of parallel processing functions 
               required for TCRM, so we only need to ensure 
               those functions are available, either as 
               the real thing or, if the required modules are 
               not available, dummy functions that pass straight 
               through. 
               We base our parallel processing on :term:`mpi4py`

.. moduleauthor: Craig Arthur, <craig.arthur@ga.gov.au>

"""

from functools import wraps

class DummyStatus(object):
    """
    A dummy `Status` class that provides a placeholder 
    for the methods that are used to control processing
    in parallel implementation

    """
    def __init__(self):
        self.source = 0
        self.tag = -1
        self.error = 0

    def __call__(self):
        pass

class DummyCommWorld(object):
    """
    A dummy COMM_WORLD class that provides the bare
    essential methods for running the code. This is used
    for basic parallelisation (task distribution).

    This is returned only if mpi4py raises an ImportError or 
    ModuleNotFoundError. 

    """
    
    def __init__(self):
        self._rank = 0
        self._size = 1
        self._name = 'DummyCommWorld'

    @property
    def name(self):
        return self._name
        
    @property
    def rank(self):
        return self._rank

    @property
    def size(self):
        return self._size

    def Get_size(self):
        return self._size

    def Get_rank(self):
        return self._rank

    def barrier(self):
        pass

    def finalize(self):
        pass

def attemptParallel():
    """
    Attempt to load `mpi4py.MPI` globally as `MPI`.  If mpi4py loads
    successfully, then a call to `MPI.Finalize` is registered to be
    called at exit of the Python interpreter. This is to ensure that
    MPI exits cleanly.

    If `mpi4py.MPI` cannot be loaded then a dummy `mpi4py.MPI` is created.

    :returns: An `mpi4py.MPI` object - either a dummy or the real thing

    """

    global MPI

    try:
        # load mpi4py for everyone

        from mpi4py import MPI
    except (ImportError, ModuleNotFoundError):

        # no mpi4py, create a dummy version of COMM_WORLD and
        # additional methods and attributes

        class DummyMPI(object):
            def __init__(self):
                self.COMM_WORLD = DummyCommWorld()
                self.Status = DummyStatus()
                self.ANY_SOURCE = -1
            
            def Init(self):
                pass

            def Finalize(self):
                pass

        MPI = DummyMPI()
    return MPI

 

def disableOnWorkers(f):
    """
    Decorator to disable function `f` calculation on workers.
    The function will only be evaluated on the master thread.

    :param f: Function to be wrapped
    :type f: function
    """

    @wraps(f)
    def wrap(*args, **kwargs):
        if MPI.COMM_WORLD.size > 1 and MPI.COMM_WORLD.rank > 0:
            return
        else:
            return f(*args, **kwargs)
    return wrap
