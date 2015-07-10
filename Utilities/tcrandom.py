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
        return mu + sigma * math.log(u1 / (1 - u1))

    def cauchyvariate(self, x0, gamma):
        """
        Random variate from the Cauchy distribution.

        :param float x0: Location parameter.
        :param float gamma: Scale parameter (|gamma| > 0)

        :returns: A random variate from the Cauchy distribution.

        """
        u1 = self.random()
        return x0 + gamma * math.tan(math.pi * (u1 - 0.5))
