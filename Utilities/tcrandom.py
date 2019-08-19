"""
:mod:`tcrandom` -- extended version of Python's :class:`random` library
=======================================================================

.. module: tcrandom
    :synopsis: Provides additional random variates beyond those in the `random` libray.
               - logisticvariate
               - cauchyvariate

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>
.. |mu| unicode:: U+003BC .. GREEK SMALL LETTER MU
.. |sigma|  unicode:: U+003C3 .. GREEK SMALL LETTER SIGMA
.. |gamma|  unicode:: U+003B3 .. GREEK SMALL LETTER GAMMA

"""
import math
from scipy.special import nctdtrit, ndtri

try:
    from numpy.random import Generator, Philox
except ImportError: # prior to numpy 1.18
    try:
        from randomgen import Generator, Philox
    except ImportError: # prior to randomgen 1.17
        from randomgen import RandomGenerator as Generator, Philox

#pylint: disable-msg=R0904

class Random:
    """
    Pseudorandom number generator.

    Expect each simulation to instantiate this class with the same
    seed integer, and with a unique stream integer (drawn from e.g. an
    enumeration of all the simulations)
    """
    # numpy 1.18 recommends Philox for independent pseudorandom streams
    def __init__(self, seed, stream):
        self.PRNG = Generator(Philox(key=seed + stream))
    def normalvariate(self, loc=0, scale=1, shape=None):
        return self.PRNG.normal(loc, scale, shape)
    def uniform(self, low=0, high=1, shape=None):
        return self.PRNG.uniform(low, high, shape)
    def random(self): # TODO: refactor elsewhere to call .uniform() directly
        return self.uniform()
#    TODO: migrate to use library implementations,
#          rather than custom implementations:
#    def logisticvariate(self, loc, sigma):
#        return self.PRNG.logistic(loc, sigma)
#    def lognormvariate:
#        return self.PRNG.lognormal(mean, sigma)

    def logisticvariate(self, mu, sigma):
        """
        Random variate from the logistic distribution.

        :param float mu: Location parameter (|mu| real).
        :param float sigma: Scale parameter (|sigma| > 0).

        :returns: A random variate from the logistic distribution.

        """
        u1 = self.random()
        if sigma <= 0.0:
            raise ValueError("Invalid input parameter: `sigma` must be positive")
        return mu + sigma * math.log(u1 / (1 - u1))

#    def cauchyvariate(self, mu, sigma):
#        """
#        Random variate from the Cauchy distribution.
#
#        :param float mu: Location parameter.
#        :param float sigma: Scale parameter (|sigma| > 0)
#
#        :returns: A random variate from the Cauchy distribution.
#
#        """
#        u1 = self.random()
#        if sigma <= 0.0:
#            raise ValueError("Invalid input parameter: `sigma` must be positive")
#        return mu + sigma * math.tan(math.pi * (u1 - 0.5))

#    def nctvariate(self, df, nc, mu=0.0, sigma=1.0):
#        """
#        Random variate from the non-central T distribution.
#
#        :param float df: degrees of freedom for the distribution.
#        :param float nc: non-centrality parameter.
#        :param float mu: Location parameter.
#        :param float sigma: Scale parameter.
#
#        :returns: A random variate from the non-central T distribution.
#        """
#        if df <= 0.0:
#            raise ValueError("Invalid input parameter: `df` must be positive")
#        if sigma <= 0.0:
#            raise ValueError("Invalid input parameter: `sigma` must be positive")
#
#        u1 = self.random()
#        return mu + sigma * nctdtrit(df, nc, u1)

    def lognormvariate(self, xi, mu=0.0, sigma=1.0):
        """
        Random variate from the lognormal distribution.

        :param float xi: Shape parameter
        :param float mu: Location parameter
        :param float sigma: Scale paramter (|sigma| > 0)

        :returns: A random variate from the lognormal distribution
        """
        if xi <= 0.0:
            raise ValueError("Invalid input parameter: `xi` must be positive")
        if sigma <= 0.0:
            raise ValueError("Invalid input parameter: `sigma` must be positive")

        u1 = self.random()
        return mu + sigma * math.exp(xi * ndtri(u1))
