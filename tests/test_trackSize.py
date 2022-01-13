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
        self.rmaxout = np.array([58.84459583, 53.05345552, 47.8322453,
                                 43.124876, 38.8807784, 35.05436002,
                                 31.60451531, 28.49418411, 25.68995347])

        self.rmaxoutcoeffs = np.array([42.01387832, 33.98774437, 27.49488535,
                                       22.24239161, 17.9933096, 14.55595226,
                                       11.77525152,  9.52576278,  7.70600581])


    def test_rmaxDefaults(self):
        """Test rmax returns correct value based on defaults"""
        dp = 25
        lat = -15
        eps = 0
        rmw = trackSize.rmax(dp, lat, eps)
        self.assertAlmostEqual(rmw, 42.926957133432225, places=1)

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
        coeffs = [4.0, -0.04, 0.006]
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxoutcoeffs)

    def test_rmaxWithIncompleteCoeffs(self):
        """Test rmax falls back to default coefficients if not enough given"""
        coeffs = [4.45, -0.05]
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
        self.params = [-0.006778075190728221, 0.24194398102420178, 0.9410510407821646]
        pklfile.close()

    def test_fitRmaxWrongLengths(self):
        """Test fitRmax raises exception when inputs are different lengths"""
        self.assertRaises(Exception, trackSize.fitRmax, self.rmw, self.dp[:-1], self.lat)


    def test_fitRmax(self):
        """Test fitRmax returns expected value"""
        params = trackSize.fitRmax(self.rmw, self.dp, self.lat)
        self.assertEqual(type(params), list)
        self.assertEqual(len(params), 3)
        print(params)
        self.numpyAssertAlmostEqual(np.array(params), np.array(self.params))
        _ = trackSize.rmax(self.dp, self.lat, 0, params)


if __name__ == "__main__":
    unittest.main(verbosity=2)
