import os
import sys
import unittest
import pickle
from . import NumpyTestCase

from wind.windmodels import *

try:
    from . import pathLocate
except:
    from tests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

class TestWindVelocity(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pkl_file = open(os.path.join(
            unittest_dir, 'test_data', 'windProfileTestData.pkl'), 'r')
        self.R = pickle.load(pkl_file)
        self.pEnv = pickle.load(pkl_file)
        self.pCentre = pickle.load(pkl_file)
        self.rMax = pickle.load(pkl_file)
        self.cLat = pickle.load(pkl_file)
        self.cLon = pickle.load(pkl_file)
        self.beta = pickle.load(pkl_file)
        self.rMax2 = pickle.load(pkl_file)
        self.beta1 = pickle.load(pkl_file)
        self.beta2 = pickle.load(pkl_file)
        self.test_wP_rankine = pickle.load(pkl_file)
        self.test_wP_jelesnianski = pickle.load(pkl_file)
        self.test_wP_holland = pickle.load(pkl_file)
        self.test_wP_willoughby = pickle.load(pkl_file)
        self.test_wP_powell = pickle.load(pkl_file)
        self.test_wP_doubleHolland = pickle.load(pkl_file)
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
        self.numpyAssertAlmostEqual(V, self.test_wP_holland,
                                    prec=1.0000000000000001e-004)

    def testWilloughby(self):
        profile = WilloughbyWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_willoughby,
                                    prec=1.0000000000000001e-004)

    def testPowell(self):
        profile = PowellWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_powell,
                                    prec=1.0000000000000001e-004)

    def testDoubleHolland(self):
        profile = DoubleHollandWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax,
            self.beta1, self.beta2, self.rMax2)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_doubleHolland,
                                    prec=1.0000000000000001e-004)

class TestWindVorticity(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pkl_file = open(os.path.join(unittest_dir, 'test_data', 'vorticityTestData.pkl'), 'r')
        self.R = pickle.load(pkl_file)
        self.pEnv = pickle.load(pkl_file)
        self.pCentre = pickle.load(pkl_file)
        self.rMax = pickle.load(pkl_file)
        self.cLat = pickle.load(pkl_file)
        self.cLon = pickle.load(pkl_file)
        self.beta = pickle.load(pkl_file)
        self.rMax2 = pickle.load(pkl_file)
        self.beta1 = pickle.load(pkl_file)
        self.beta2 = pickle.load(pkl_file)
        self.vMax = pickle.load(pkl_file)
        self.test_vorticity_rankine = pickle.load(pkl_file)
        self.test_vorticity_jelesnianski = pickle.load(pkl_file)
        self.test_vorticity_holland = pickle.load(pkl_file)
        self.test_vorticity_willoughby = pickle.load(pkl_file)
        self.test_vorticity_doubleHolland = pickle.load(pkl_file)
        self.test_vorticity_powell = pickle.load(pkl_file)
        pkl_file.close()

    def testRankine(self):
        profile = RankineWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.vMax = self.vMax
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_rankine)

    def testJelesnianski(self):
        profile = JelesnianskiWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.vMax = self.vMax
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_jelesnianski)

    def testHolland(self):
        profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                     self.pCentre, self.rMax, self.beta)
        profile.vMax = self.vMax
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_holland)

    def testWilloughby(self):
        profile = WilloughbyWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        # Hack for testing as vMax needs to be set
        profile.vMax = self.vMax
        profile.beta = 1.0036 + 0.0173 * profile.vMax - 0.313 * np.log(self.rMax) + 0.0087 * np.abs(self.cLat)
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_willoughby)

    def testDoubleHolland(self):
        profile = DoubleHollandWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax,
            self.beta1, self.beta2, self.rMax2)
        profile.vMax = self.vMax
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_doubleHolland)

    def testPowell(self):
        profile = PowellWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.vMax = self.vMax
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_powell)


class TestWindField(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pkl_file = open(os.path.join(unittest_dir, 'test_data', 'windFieldTestData.pkl'), 'r')
        self.R = pickle.load(pkl_file)
        self.lam = pickle.load(pkl_file)
        self.rMax = pickle.load(pkl_file)
        self.f = pickle.load(pkl_file)
        self.V = pickle.load(pkl_file)
        self.Z = pickle.load(pkl_file)
        self.vFm = pickle.load(pkl_file)
        self.thetaFm = pickle.load(pkl_file)
        self.thetaMax = pickle.load(pkl_file)
        self.test_kepert_Ux = pickle.load(pkl_file)
        self.test_kepert_Vy = pickle.load(pkl_file)
        self.test_mcconochie_Ux = pickle.load(pkl_file)
        self.test_mcconochie_Vy = pickle.load(pkl_file)
        self.test_hubbert_Ux = pickle.load(pkl_file)
        self.test_hubbert_Vy = pickle.load(pkl_file)
        pkl_file.close()

    def test_Kepert(self):
        profile = WindProfileModel(-15., 0.0, 1000., 990., self.rMax, WindSpeedModel)
        profile.f = self.f
        windField = KepertWindField(profile)
        windField.V = self.V
        windField.Z = self.Z

        Ux, Vy = windField.field(self.R, self.lam, self.vFm, self.thetaFm, self.thetaMax)
        self.numpyAssertAlmostEqual(Ux, self.test_kepert_Ux)
        self.numpyAssertAlmostEqual(Vy, self.test_kepert_Vy)

    def test_McConochie(self):
        profile = WindProfileModel(-15., 0.0, 1000., 990., self.rMax, WindSpeedModel)
        profile.f = self.f
        windField = McConochieWindField(profile)
        windField.V = self.V
        windField.Z = self.Z
        Ux, Vy = windField.field(self.R, self.lam, self.vFm, self.thetaFm, self.thetaMax)
        self.numpyAssertAlmostEqual(Ux, self.test_mcconochie_Ux)
        self.numpyAssertAlmostEqual(Vy, self.test_mcconochie_Vy)

    def test_Hubbert(self):
        profile = WindProfileModel(-15., 0.0, 1000., 990., self.rMax, WindSpeedModel)
        profile.f = self.f
        windField = HubbertWindField(profile)
        windField.V = self.V
        windField.Z = self.Z
        Ux, Vy = windField.field(self.R, self.lam, self.vFm, self.thetaFm, self.thetaMax)
        self.numpyAssertAlmostEqual(Ux, self.test_hubbert_Ux)
        self.numpyAssertAlmostEqual(Vy, self.test_hubbert_Vy)

if __name__ == "__main__":
    testSuite = unittest.makeSuite(TestWindVelocity, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)

    testSuite = unittest.makeSuite(TestWindVorticity, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)

    testSuite = unittest.makeSuite(TestWindField, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
