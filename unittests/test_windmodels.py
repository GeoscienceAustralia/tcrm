import os
import sys
import unittest
import cPickle
import NumpyTestCase

from WindfieldInterface.windmodels import *

try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

class TestWindVelocity(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pkl_file = open(os.path.join(
            unittest_dir, 'test_data', 'windProfile_testdata.pck'), 'rb')
        self.R = cPickle.load(pkl_file)
        self.pEnv = cPickle.load(pkl_file)
        self.pCentre = cPickle.load(pkl_file)
        self.rMax = cPickle.load(pkl_file)
        self.cLat = cPickle.load(pkl_file)
        self.cLon = cPickle.load(pkl_file)
        self.beta = cPickle.load(pkl_file)
        self.rMax2 = cPickle.load(pkl_file)
        self.beta1 = cPickle.load(pkl_file)
        self.beta2 = cPickle.load(pkl_file)
        self.test_wP_rankine = cPickle.load(pkl_file)
        self.test_wP_jelesnianski = cPickle.load(pkl_file)
        self.test_wP_holland = cPickle.load(pkl_file)
        self.test_wP_willoughby = cPickle.load(pkl_file)
        self.test_wP_doubleHolland = cPickle.load(pkl_file)
        self.test_wP_powell = cPickle.load(pkl_file)
        pkl_file.close()

    def testRankine(self):
        profile = RankineWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_rankine)

    def testJelesnianski(self):
        profile = JelesnianskiWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_jelesnianski)

    def testHolland(self):
        profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                     self.pCentre, self.rMax, self.beta)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_holland)

    def testWilloughby(self):
        profile = WilloughbyWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_willoughby)

    def testPowell(self):
        profile = PowellWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_powell)

    def testDoubleHolland(self):
        profile = DoubleHollandWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax,
            self.beta1, self.beta2, self.rMax2)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_doubleHolland)

class TestWindVorticity(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pkl_file = open(os.path.join(unittest_dir, 'test_data', 'vorticity_testdata.pck'), 'rb')
        self.R = cPickle.load(pkl_file)
        self.pEnv = cPickle.load(pkl_file)
        self.pCentre = cPickle.load(pkl_file)
        self.rMax = cPickle.load(pkl_file)
        self.cLat = cPickle.load(pkl_file)
        self.cLon = cPickle.load(pkl_file)
        self.beta = cPickle.load(pkl_file)
        self.rMax2 = cPickle.load(pkl_file)
        self.beta1 = cPickle.load(pkl_file)
        self.beta2 = cPickle.load(pkl_file)
        self.vMax = cPickle.load(pkl_file)
        self.test_vorticity_rankine = cPickle.load(pkl_file)
        self.test_vorticity_jelesnianski = cPickle.load(pkl_file)
        self.test_vorticity_holland = cPickle.load(pkl_file)
        self.test_vorticity_willoughby = cPickle.load(pkl_file)
        self.test_vorticity_doubleHolland = cPickle.load(pkl_file)
        self.test_vorticity_powell = cPickle.load(pkl_file)
        pkl_file.close()

    def testRankine(self):
        profile = RankineWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.vMax = self.vMax
        V = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_vorticity_rankine)

    def testJelesnianski(self):
        profile = JelesnianskiWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.vMax = self.vMax
        V = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_vorticity_jelesnianski)

    def testHolland(self):
        profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                     self.pCentre, self.rMax, self.beta)
        profile.vMax = self.vMax
        V = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_vorticity_holland)

    def testWilloughby(self):
        profile = WilloughbyWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        # Hack for testing as vMax needs to be set
        profile.vMax = self.vMax
        profile.beta = 1.0036 + 0.0173 * profile.vMax - 0.313 * np.log(self.rMax) + 0.0087 * np.abs(self.cLat)
        V = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_vorticity_willoughby)

    def testDoubleHolland(self):
        profile = DoubleHollandWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax,
            self.beta1, self.beta2, self.rMax2)
        profile.vMax = self.vMax
        V = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_vorticity_doubleHolland)

    def testPowell(self):
        profile = PowellWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.vMax = self.vMax
        V = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_vorticity_powell)

if __name__ == "__main__":
    testSuite = unittest.makeSuite(TestWindVelocity, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)

    testSuite = unittest.makeSuite(TestWindVorticity, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
