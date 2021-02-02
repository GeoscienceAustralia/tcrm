"""
Test random number generation extensions

FIXME: Test additional variates available through the Random module
"""

import unittest
import numpy as np

from numpy.testing import assert_almost_equal
from Utilities.tcrandom import Random

class TestRandom(unittest.TestCase):

    def setUp(self):

        self.seed = 1
        self.prng = Random(self.seed, stream=0)

    @unittest.skip("Generator implementation specific")
    def testLogistic(self):
        """Testing logistic variates"""
        self.prng.seed(self.seed)
        result = self.prng.logisticvariate(0, 1)
        assert_almost_equal(result, -1.86290986)

    @unittest.skip("Generator implementation specific")
    def testNormal(self):
        """Testing normal variates"""
        self.prng.seed(self.seed)
        result = self.prng.normalvariate(0, 1)
        assert_almost_equal(result, 0.607455857)

    @unittest.skip("Generator implementation specific")
    def testCauchy(self):
        """Testing cauchy variates"""
        self.prng.seed(self.seed)
        result = self.prng.cauchyvariate(0, 1)
        assert_almost_equal(result, -2.22660116)

    @unittest.skip("Generator implementation specific")
    def testNCT(self):
        self.prng.seed(self.seed)
        result = self.prng.nctvariate(1, 0)
        assert_almost_equal(result, -2.22660116)
        self.prng.seed(1)
        result = self.prng.nctvariate(10, 0.5, 1, 0.5)
        assert_almost_equal(result, 0.68386525)

    def testLogisticInvalidParams(self):
        self.assertRaises(ValueError, self.prng.logisticvariate,
                          0, -1)

#    def testCauchyInvalidParams(self):
#        self.assertRaises(ValueError, self.prng.cauchyvariate,
#                          0, -1)
#
#    def testNCTInvalidParams(self):
#        self.assertRaises(ValueError, self.prng.nctvariate,
#                          -5, 0)
#        self.assertRaises(ValueError, self.prng.nctvariate,
#                          1, 0, 1, -1)

    @unittest.skip("Generator implementation specific")
    def testLognorm(self):
        """
        Testing lognorm variates
        """
        self.prng.seed(self.seed)
        result = self.prng.lognormvariate(1.)
        assert_almost_equal(result, 0.3308813420)

    def testLognormInvalidParams(self):
        self.assertRaises(ValueError, self.prng.lognormvariate,
                          -1., 0., 1.)
        self.assertRaises(ValueError, self.prng.lognormvariate,
                          1., 0., -1.)

if __name__ == '__main__':
    unittest.main()
