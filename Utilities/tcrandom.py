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
import random
import math
from scipy.special import nctdtrit

#pylint: disable-msg=R0904

class Random(random.Random):
    """
    An extension of the standard :mod:`random` library to
    allow sampling from additional distributions.

    """

    def __init__(self, value=None):
        random.Random.__init__(self, value)

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

    def cauchyvariate(self, mu, sigma):
        """
        Random variate from the Cauchy distribution.

        :param float mu: Location parameter.
        :param float sigma: Scale parameter (|sigma| > 0)

        :returns: A random variate from the Cauchy distribution.

        """
        u1 = self.random()
        if sigma <= 0.0:
            raise ValueError("Invalid input parameter: `sigma` must be positive")
        return mu + sigma * math.tan(math.pi * (u1 - 0.5))

    def nctvariate(self, df, nc, mu=0.0, sigma=1.0):
        """
        Random variate from the non-central T distribution.

        :param float df: degrees of freedom for the distribution.
        :param float nc: non-centrality parameter.
        :param float mu: Location parameter.
        :param float sigma: Scale parameter.

        :returns: A random variate from the non-central T distribution.
        """
        if df <= 0.0:
            raise ValueError("Invalid input parameter: `df` must be positive")
        if sigma <= 0.0:
            raise ValueError("Invalid input parameter: `sigma` must be positive")

        u1 = self.random()
        return mu + sigma * nctdtrit(df, nc, u1)
