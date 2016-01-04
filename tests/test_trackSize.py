import os
import sys
import unittest
import numpy as np
import NumpyTestCase
import cPickle
try:
    import pathLocate
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
        self.rmaxout = np.array([104.08747056, 87.24373672, 74.36633603,
                                 64.46512563, 56.83025497, 50.94959054,
                                 46.45239542, 43.07069234, 40.61270477])
        self.rmaxoutcoeffs = np.array([109.4918411, 90.31016614,
                                       75.35799588, 63.61504735,
                                       54.32855894, 46.93905396,
                                       41.02780616, 36.27939815,
                                       32.45485465])


    def test_rmaxDefaults(self):
        """Test rmax returns correct value based on defaults"""
        dp = 25
        lat = -15
        eps = np.random.normal(0, scale=0.353)
        rmw = trackSize.rmax(dp, lat, eps)
        self.assertEqual(rmw, 63.867449833342235)

    def test_rmaxWrongLengths(self):
        """rmax raises exception when inputs are different lengths"""
        eps = np.random.normal(0, scale=0.353)
        latarray = np.arange(-23, 5, 2)
        self.assertRaises(Exception, trackSize.rmax, self.dparray, latarray, eps)

    def test_rmaxArrayInput(self):
        """Test rmax with array input"""
        eps = np.random.normal(0, scale=0.353)
        rmw = trackSize.rmax(self.dparray, self.latarray, eps)
        self.numpyAssertAlmostEqual(rmw, self.rmaxout)

    def test_rmaxWithCoeffs(self):
        """Test rmax with user-defined coefficients"""
        eps = np.random.normal(0, scale=0.353)
        coeffs = [4.5, -0.04, 0.0002, 0.0002]
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxoutcoeffs)

    def test_rmaxWithIncompleteCoeffs(self):
        """Test rmax falls back to default coefficients if not enough given"""
        coeffs = [4.45, -0.05, 0.0002]
        eps = np.random.normal(0, scale=0.353)
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs=coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxout)

class TestFitRmax(NumpyTestCase.NumpyTestCase):
    """
    Note: this test uses a subset of the full dataset used to determine
    the default coefficients for the Rmw model used in `TrackGenerator`
    """
    def setUp(self):
        pklfile = open(os.path.join(unittest_dir, 'test_data', 'rmw.pck'), 'r')
        self.dp = cPickle.load(pklfile)
        self.lat = cPickle.load(pklfile)
        self.rmw = cPickle.load(pklfile)
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
