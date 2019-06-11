import os
import sys
import unittest
import numpy as np
from . import NumpyTestCase
import pickle
try:
    from . import pathLocate
except:
    from unittests import pathLocate

    
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

from TrackGenerator import trackSize

class TestRmaxModel(NumpyTestCase.NumpyTestCase):
    
    def setUp(self):
        np.random.seed(10)
        self.dparray = np.arange(10, 51, 5)
        self.latarray = np.arange(-23, -5, 2)
        self.rmaxout = np.array([66.21418972, 54.96078787, 45.7678418,
                                 39.33639152, 35.22023191, 32.67958816,
                                 31.07181726, 29.95463557, 29.06996118])

        self.rmaxoutcoeffs = np.array([57.74931727, 49.76638251, 42.67145839,
                                       37.24632781, 33.47199844, 30.98197273,
                                       29.3570741, 28.25288743, 27.43217413])


    def test_rmaxDefaults(self):
        """Test rmax returns correct value based on defaults"""
        dp = 25
        lat = -15
        eps = 0
        rmw = trackSize.rmax(dp, lat, eps)
        self.assertAlmostEqual(rmw, 39.21410711, places=1)

    def test_rmaxWrongLengths(self):
        """rmax raises exception when inputs are different lengths"""
        eps = 0
        latarray = np.arange(-23, 5, 2)
        self.assertRaises(Exception, trackSize.rmax, self.dparray, latarray, eps)

    def test_rmaxArrayInput(self):
        """Test rmax with array input"""
        eps = 0
        rmw = trackSize.rmax(self.dparray, self.latarray, eps)
        self.numpyAssertAlmostEqual(rmw, self.rmaxout, prec=0.1)

    def test_rmaxWithCoeffs(self):
        """Test rmax with user-defined coefficients"""
        eps = 0
        coeffs = [3.5, -0.004, 0.7, 0.002, .001]
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxoutcoeffs)

    def test_rmaxWithIncompleteCoeffs(self):
        """Test rmax falls back to default coefficients if not enough given"""
        coeffs = [4.45, -0.05, 0.0002]
        eps = 0
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs=coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxout, prec=0.1)

class TestFitRmax(NumpyTestCase.NumpyTestCase):
    """
    Note: this test uses a subset of the full dataset used to determine
    the default coefficients for the Rmw model used in `TrackGenerator`
    """
    def setUp(self):
        pklfile = open(os.path.join(unittest_dir, 'test_data', 'rmw.pkl'), 'rb')
        self.dp = pickle.load(pklfile)
        self.lat = pickle.load(pklfile)
        self.rmw = pickle.load(pklfile)
        self.params = [4.4650608902114888,
                       -0.042494641709203987,
                       0.00033723892839458182,
                       0.00021502458395316267,
                       0.35665997379737535]
        pklfile.close()

    def test_fitRmaxWrongLengths(self):
        """Test fitRmax raises exception when inputs are different lengths"""
        self.assertRaises(Exception, trackSize.fitRmax, self.rmw, self.dp[:-1], self.lat)


    def test_fitRmax(self):
        """Test fitRmax returns expected value"""
        params = trackSize.fitRmax(self.rmw, self.dp, self.lat)
        self.assertEqual(type(params), list)
        self.assertEqual(len(params), 5)
        self.numpyAssertAlmostEqual(np.array(params), np.array(self.params))
    
if __name__ == "__main__":
    unittest.main(verbosity=2)
