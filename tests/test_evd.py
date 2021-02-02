"""
Testing the extreme value distributions
"""

import unittest
import numpy as np

from numpy.testing import assert_almost_equal
from hazard.evd import gevfit, empfit


class TestEvd(unittest.TestCase):

    def setUp(self):
        self.v = np.array([0., 0., 38.55, 41.12, 59.29, 61.75, 74.79])
        self.years = np.array([25.0, 50.0, 100.0, 250.0, 500.0, 2000.0])
        self.w0 = np.array([ 47.65027995, 66.40860811, 79.83843278,
                             93.03113963, 100.67294472,111.81331626])

        self.loc0 = np.array([28.64508831])
        self.scale0 = np.array([31.21743669])
        self.shp0 = np.array([0.29781455])
        self.missingValue = np.array(-9999.0)

    def testEVD(self):
        """Testing extreme value distribution"""
        w, loc, scale, shp = gevfit(self.v,
                                    self.years,
                                    nodata=-9999,
                                    minrecords=3,
                                    yrspersim=10)

        assert_almost_equal(w, self.w0, decimal=5)
        assert_almost_equal(loc, self.loc0, decimal=5)
        assert_almost_equal(scale, self.scale0, decimal=5)
        assert_almost_equal(shp, self.shp0, decimal=5)

        w2, loc2, scale2, shp2 = gevfit(self.v,
                                        self.years,
                                        nodata=-9999,
                                        minrecords=50,
                                        yrspersim=10)

        assert_almost_equal(w2, np.ones(6) * self.missingValue, decimal=5)
        assert_almost_equal(loc2, self.missingValue, decimal=5)
        assert_almost_equal(scale2, self.missingValue, decimal=5)
        assert_almost_equal(shp2, self.missingValue, decimal=5)

class TestEmpiricalFit(unittest.TestCase):

    def setUp(self):
        self.v = np.array([ 24.01518592,  20.71014074,  21.72734664,  23.55090159,
                            22.55483303,  24.49634915,  24.21779962,  29.13960142,
                            24.01517175,  29.22021474,  22.40286699,  21.95521168,
                            25.37122473,  21.72289402,  22.99731093,  25.80831913,
                            20.68754605,  24.0068121 ,  40.67975816,  24.31328418])

        self.w = np.array([ 24.01515916,  25.80802   ,  29.22018164,  40.67662114])
        self.nodata = -9999.

    def testEmpFit(self):
        """Test calculating empirical return levels"""
        years = np.array([2, 5, 10, 20])
        w, loc, scale, shp = empfit(self.v, years, 20)
        assert_almost_equal(w, self.w)
        assert_almost_equal(loc, self.nodata)
        assert_almost_equal(shp, self.nodata)
        assert_almost_equal(scale, self.nodata)


if __name__ == "__main__":
    suite = unittest.makeSuite(TestEvd, 'test')
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.makeSuite(TestEmpiricalFit, 'test')
    unittest.TextTestRunner(verbosity=2).run(suite)
