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
        self.rmaxout = np.array([ 105.95720455, 87.69107765, 73.93504829,
                                  63.50608684, 55.57127783, 49.53993723,
                                  44.99150895, 41.62705669, 39.23655673])
        self.rmaxoutcoeffs = np.array([110.07658843, 90.79247266,
                                       75.76044949, 63.95478709,
                                       54.61870366, 47.18973462,
                                       41.24691747, 36.47315031,
                                       32.62818163])


    def test_rmaxDefaults(self):
        """Test rmax returns correct value based on defaults"""
        dp = 25
        lat = -15
        eps = np.random.normal(0, scale=0.357)
        rmw = trackSize.rmax(dp, lat, eps)
        self.assertEqual(rmw, 62.638227644545715)

    def test_rmaxWrongLengths(self):
        """rmax raises exception when inputs are different lengths"""
        eps = np.random.normal(0, scale=0.357)
        latarray = np.arange(-23, 5, 2)
        self.assertRaises(Exception, trackSize.rmax, self.dparray, latarray, eps)

    def test_rmaxArrayInput(self):
        """Test rmax with array input"""
        eps = np.random.normal(0, scale=0.357)
        rmw = trackSize.rmax(self.dparray, self.latarray, eps)
        self.numpyAssertAlmostEqual(rmw, self.rmaxout)

    def test_rmaxWithCoeffs(self):
        """Test rmax with user-defined coefficients"""
        eps = np.random.normal(0, scale=0.357)
        coeffs = [4.5, -0.04, 0.0002, 0.0002]
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxoutcoeffs)

    def test_rmaxWithIncompleteCoeffs(self):
        """Test rmax falls back to default coefficients if not enough given"""
        coeffs = [4.45, -0.05, 0.0002]
        eps = np.random.normal(0, scale=0.357)
        rmw = trackSize.rmax(self.dparray, self.latarray, eps, coeffs=coeffs)
        self.numpyAssertAlmostEqual(rmw, self.rmaxout)

class TestFitRmax(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pklfile = open(os.path.join(unittest_dir, 'test_data', 'rmw.pck'), 'r')
        self.dp = cPickle.load(pklfile)
        self.lat = cPickle.load(pklfile)
        self.rmw = cPickle.load(pklfile)
        self.params = [4.465060890211495, -0.042494641709204063,
                       0.0003372389283945545, 0.0002150245839531746,
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
