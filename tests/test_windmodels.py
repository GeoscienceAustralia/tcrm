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
        self.cLon = 158.17
        self.cLat = -20.13
        self.pCentre = 95330.
        self.pEnv = 101445.0
        self.rMax = 50000.
        self.rMax2 = 90000.
        self.beta = 1.7
        self.beta1 = 1.7
        self.beta2 = 1.3
        self.vFm = 10.
        self.thetaFm = 70. * np.pi / 180.
        self.thetaMax = 70. * np.pi / 180.

        pkl_file = open(os.path.join(
            unittest_dir, 'test_data', 'windProfileTestData.pkl'), 'rb')
        self.R = pickle.load(pkl_file)
        self.lam = pickle.load(pkl_file)
        self.test_wP_rankine = pickle.load(pkl_file)
        self.test_wP_jelesnianski = pickle.load(pkl_file)
        self.test_wP_holland = pickle.load(pkl_file)
        self.test_wP_willoughby = pickle.load(pkl_file)
        self.test_wP_powell = pickle.load(pkl_file)
        self.test_wP_doubleHolland = pickle.load(pkl_file)
        pkl_file.close()

    def testRankine(self):
        """Test Rankine radial profile"""
        profile = RankineWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_rankine)

    def testJelesnianski(self):
        """Test Jelesnianski radial profile"""
        profile = JelesnianskiWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_jelesnianski)

    def testHolland(self):
        """Test Holland radial profile"""
        profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                     self.pCentre, self.rMax, self.beta)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_holland,
                                    prec=1.0000000000000001e-004)

    def testWilloughby(self):
        """Test Willoughby radial profile"""
        profile = WilloughbyWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_willoughby,
                                    prec=1.0000000000000001e-004)

    def testPowell(self):
        """Test Powell radial profile"""
        profile = PowellWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_powell,
                                    prec=1.0000000000000001e-004)

    def testDoubleHolland(self):
        """Test Double Holland radial profile"""
        profile = DoubleHollandWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax,
            self.beta1, self.beta2, self.rMax2)
        V = profile.velocity(self.R)
        self.numpyAssertAlmostEqual(V, self.test_wP_doubleHolland,
                                    prec=1.0000000000000001e-004)

class TestWindField(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.cLon = 158.17
        self.cLat = -20.13
        self.pCentre = 95330.
        self.pEnv = 101445.0
        self.rMax = 50000.
        self.rMax2 = 90000.
        self.beta = 1.7
        self.beta1 = 1.7
        self.beta2 = 1.3
        self.vFm = 10.
        self.thetaFm = 70. * np.pi / 180.
        self.thetaMax = 70. * np.pi / 180.
        self.profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                          self.pCentre, self.rMax, self.beta)
        pkl_file = open(os.path.join(unittest_dir, 'test_data',
                        'windFieldTestData.pkl'), 'rb')
        self.R = pickle.load(pkl_file)
        self.lam = pickle.load(pkl_file)
        self.test_kepert_Ux = pickle.load(pkl_file)
        self.test_kepert_Vy = pickle.load(pkl_file)
        self.test_mcconochie_Ux = pickle.load(pkl_file)
        self.test_mcconochie_Vy = pickle.load(pkl_file)
        self.test_hubbert_Ux = pickle.load(pkl_file)
        self.test_hubbert_Vy = pickle.load(pkl_file)
        pkl_file.close()

    def test_Kepert(self):
        windField = KepertWindField(self.profile)
        Ux, Vy = windField.field(self.R, self.lam, self.vFm, self.thetaFm,
                                 self.thetaMax)
        self.numpyAssertAlmostEqual(Ux, self.test_kepert_Ux)
        self.numpyAssertAlmostEqual(Vy, self.test_kepert_Vy)

    def test_McConochie(self):

        windField = McConochieWindField(self.profile)
        Ux, Vy = windField.field(self.R, self.lam, self.vFm, self.thetaFm,
                                 self.thetaMax)
        self.numpyAssertAlmostEqual(Ux, self.test_mcconochie_Ux)
        self.numpyAssertAlmostEqual(Vy, self.test_mcconochie_Vy)

    def test_Hubbert(self):

        windField = HubbertWindField(self.profile)
        Ux, Vy = windField.field(self.R, self.lam, self.vFm, self.thetaFm,
                                 self.thetaMax)
        self.numpyAssertAlmostEqual(Ux, self.test_hubbert_Ux)
        self.numpyAssertAlmostEqual(Vy, self.test_hubbert_Vy)

class TestWindVorticity(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.cLon = 158.17
        self.cLat = -20.13
        self.pCentre = 95330.
        self.pEnv = 101445.0
        self.rMax = 50000.
        self.rMax2 = 90000.
        self.beta = 1.7
        self.beta1 = 1.7
        self.beta2 = 1.3
        self.vFm = 10.
        self.thetaFm = 70. * np.pi / 180.
        self.thetaMax = 70. * np.pi / 180.
        self.profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                          self.pCentre, self.rMax, self.beta)

        pkl_file = open(os.path.join(unittest_dir, 'test_data',
                        'vorticityTestData.pkl'), 'rb')
        self.R = pickle.load(pkl_file)
        self.lam = pickle.load(pkl_file)
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
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_rankine)

    def testJelesnianski(self):
        profile = JelesnianskiWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_jelesnianski)

    def testHolland(self):
        profile = HollandWindProfile(self.cLat, self.cLon, self.pEnv,
                                     self.pCentre, self.rMax, self.beta)
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_holland)

    def testWilloughby(self):
        profile = WilloughbyWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        profile.beta = 1.0036 + 0.0173 * profile.vMax - \
            0.0313 * np.log(self.rMax) + 0.0087 * np.abs(self.cLat)
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_willoughby)

    def testDoubleHolland(self):
        profile = DoubleHollandWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax,
            self.beta1, self.beta2, self.rMax2)
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_doubleHolland)

    def testPowell(self):
        profile = PowellWindProfile(
            self.cLat, self.cLon, self.pEnv, self.pCentre, self.rMax)
        Z = profile.vorticity(self.R)
        self.numpyAssertAlmostEqual(Z, self.test_vorticity_powell)




if __name__ == "__main__":
    testSuite = unittest.makeSuite(TestWindVelocity, 'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)

    #testSuite = unittest.makeSuite(TestWindVorticity, 'test')
    #unittest.TextTestRunner(verbosity=2).run(testSuite)

    #testSuite = unittest.makeSuite(TestWindField, 'test')
    #unittest.TextTestRunner(verbosity=2).run(testSuite)
