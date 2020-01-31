import os
import sys
import unittest
from . import NumpyTestCase
import numpy as np

from wind.windmodels import *

try:
    from . import pathLocate
except:
    from tests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

class TestWindProfile(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.rMax = 40e3
        self.lat = -15.
        self.lon = 0.
        self.eP_hPa = 1010.
        self.cP_hPa = 970.
        self.eP_Pa = 101000.
        self.cP_Pa = 97000.
        self.beta = 1.9
        self.R = 1000.*np.linspace(1., 200.)

    def testHolland(self):
        """Test radial profile returns same value for different input units"""
        profile_hPa = HollandWindProfile(self.lat, self.lon, self.eP_hPa,
                                         self.cP_hPa, self.rMax, self.beta)
        V_hPa = profile_hPa.velocity(self.R)

        profile_Pa = HollandWindProfile(self.lat, self.lon, self.eP_Pa,
                                         self.cP_Pa, self.rMax, self.beta)
        V_Pa = profile_Pa.velocity(self.R)

        self.numpyAssertAlmostEqual(V_hPa, V_Pa, prec=1.0000000000000001e-004)

if __name__ == "__main__":
    testSuite = unittest.makeSuite(TestWindProfile, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)

        
