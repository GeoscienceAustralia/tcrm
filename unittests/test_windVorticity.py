#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0
    Copyright (C) 2011  Geoscience Australia

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Title: testVorticity.py

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-12-04
Description: Unit test for the windVorticity class

ModifiedBy: Nicholas Summon, nicholas.summons@ga.gov.au
ModifiedDate: 2011-06-10
Modification: Replaced with more thorough test using picked input/output data

$Id: testVorticity.py 563 2007-10-24 02:52:40Z carthur $
"""
import os, sys
import unittest
import cPickle
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from WindfieldInterface import windVorticity
from Utilities.files import flStartLog


class TestVorticity(NumpyTestCase.NumpyTestCase):

    """Testing vorticity"""
    pkl_file = open(os.path.join(unittest_dir, 'test_data', 'vorticity_testdata.pck'), 'rb')
    R = cPickle.load(pkl_file)
    pEnv = cPickle.load(pkl_file)
    pCentre = cPickle.load(pkl_file)
    rMax = cPickle.load(pkl_file)
    cLat = cPickle.load(pkl_file)
    cLon = cPickle.load(pkl_file)
    beta = cPickle.load(pkl_file)
    rMax2 = cPickle.load(pkl_file)
    beta1 = cPickle.load(pkl_file)
    beta2 = cPickle.load(pkl_file)
    vMax = cPickle.load(pkl_file)
    test_vorticity_rankine = cPickle.load(pkl_file)
    test_vorticity_jelesnianski = cPickle.load(pkl_file)
    test_vorticity_holland = cPickle.load(pkl_file)
    test_vorticity_willoughby = cPickle.load(pkl_file)
    test_vorticity_doubleHolland = cPickle.load(pkl_file)
    test_vorticity_powell = cPickle.load(pkl_file)
    pkl_file.close()

    vorticity = windVorticity.WindVorticity(R, pEnv, pCentre, rMax, cLat, cLon, beta, rMax2, beta1, beta2)


    def test_rankine(self):
        """Testing rankine
        """
        vorticity_rankine = self.vorticity.rankine(self.vMax)
        self.numpyAssertAlmostEqual(vorticity_rankine, self.test_vorticity_rankine)

    def test_jelesnianski(self):
        """Testing jelesnianski
        """
        vorticity_jelesnianski = self.vorticity.jelesnianski(self.vMax)
        self.numpyAssertAlmostEqual(vorticity_jelesnianski, self.test_vorticity_jelesnianski)

    def test_holland(self):
        """Testing holland
        """
        vorticity_holland = self.vorticity.holland(self.vMax)
        self.numpyAssertAlmostEqual(vorticity_holland, self.test_vorticity_holland)
    
    def test_willoughby(self):
        """Testing willoughby
        """
        vorticity_willoughby = self.vorticity.willoughby(self.vMax)
        self.numpyAssertAlmostEqual(vorticity_willoughby, self.test_vorticity_willoughby)
        
    def test_doubleHolland(self):
        """Testing doubleHolland
        """
        vorticity_doubleHolland = self.vorticity.doubleHolland(self.vMax)
        self.numpyAssertAlmostEqual(vorticity_doubleHolland, self.test_vorticity_doubleHolland)

    def test_powell(self):
        """Testing powell
        """       
        vorticity_powell = self.vorticity.powell(self.vMax)
        self.numpyAssertAlmostEqual(vorticity_powell, self.test_vorticity_powell)      

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestVorticity,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)


# import sys
# import unittest
# from scipy import sign
# from Utilities import maputils
# from Utilities import metutils
# import WindfieldInterface.windProfile as windProfile
# import WindfieldInterface.windVorticity as windVorticity

# class TestVorticity(unittest.TestCase):
    # """Testing vorticity"""
    # cLat = -17
    # cLon = 150
    # R, lam = maputils.makeGrid(cLon,cLat)
    # f = metutils.coriolis(cLat)
    # xGrid,yGrid = maputils.meshLatLon(cLon,cLat)
    # rMax = 18
    # pEnv = 100700
    # pCentre = 95900
    # beta = 1.3
    # vorticity = windVorticity.WindVorticity(R,pEnv,pCentre,rMax,cLat,cLon,beta)
    # profile = windProfile.WindProfile(R,pEnv,pCentre,rMax,cLat,cLon,beta)

    # def test_All(self):
        # """testing all

        # Assumes there are corresponding methods in both windProfile and windVorticity
        # """
        # methodList = [method for method in dir(self.vorticity) if callable(getattr(self.vorticity, method))]
        # for method in methodList:
            # if method == '__doc__' or method == '__init__':
                # pass
            # else:
                # print method
                # vort = getattr(self.vorticity, method)
                # wind = getattr(self.profile, method)
                # V = wind()
                # Z = vort(abs(V).max())
                # self.assertEqual(sign(Z.max()), sign(Z.min()))

# mytoolTestSuite = unittest.TestLoader().loadTestsFromTestCase(TestVorticity)
# unittest.TextTestRunner(verbosity=2).run(mytoolTestSuite)
