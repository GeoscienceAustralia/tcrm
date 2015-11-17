"""
:mod:`parallel` -- base parallel processing functions
=====================================================

.. module:: parallel
    :synopsis: Provides a couple of base functions that set up a
               parallel processing environment.

.. moduleauthor: Craig Arthur, <craig.arthur@ga.gov.au>

"""

from functools import wraps

def attemptParallel():
    """
    Attempt to load Pypar globally as `pp`.  If Pypar loads
    successfully, then a call to `pypar.finalize` is registered to be
    called at exit of the Python interpreter. This is to ensure that
    MPI exits cleanly.

    If pypar cannot be loaded then a dummy `pp` is created.

    :returns: A pypar object - either a dummy or the real thing

    """

    global pp

    try:
        # load pypar for everyone

        import pypar as pp

    except ImportError:

        # no pypar, create a dummy one

        class DummyPypar(object):

            def size(self):
                return 1

            def rank(self):
                return 0

            def barrier(self):
                pass

            def finalize(self):
                pass

        pp = DummyPypar()
    return pp

def disableOnWorkers(f):
    """
    Decorator to disable function `f` calculation on workers.
    The function will only be evaluated on the master thread.

    :param f: Function to be wrapped
    :type f: function
    """

    @wraps(f)
    def wrap(*args, **kwargs):
        if pp.size() > 1 and pp.rank() > 0:
            return
        else:
            return f(*args, **kwargs)
    return wrap
