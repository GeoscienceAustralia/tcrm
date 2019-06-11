"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

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


 Title: TestPressureProfile.py
 Author: Nicholas Summons, nicholas.summons@ga.gov.au
 CreationDate: 2011-06-09
 Description: Unit testing module for pressureProfile

 Version: $Rev$

 $Id$
"""
import os, sys
import pickle
import unittest
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from PressureInterface import pressureProfile
from Utilities.files import flStartLog

class TestPressureProfile(NumpyTestCase.NumpyTestCase):

    pkl_file = open(os.path.join(unittest_dir, 'test_data', 'pressureProfileTestData.pkl'), 'rb')
    R = pickle.load(pkl_file)
    pEnv = pickle.load(pkl_file)
    pCentre = pickle.load(pkl_file)
    rMax = pickle.load(pkl_file)
    cLat = pickle.load(pkl_file)
    cLon = pickle.load(pkl_file)
    beta = pickle.load(pkl_file)
    rMax2 = pickle.load(pkl_file)
    beta1 = pickle.load(pkl_file)
    beta2 = pickle.load(pkl_file)
    test_pHolland = pickle.load(pkl_file)
    test_pWilloughby = pickle.load(pkl_file)
    test_pdoubleHolland = pickle.load(pkl_file)
    test_pPowell = pickle.load(pkl_file)
    pkl_file.close()

    prP = pressureProfile.PrsProfile(R, pEnv, pCentre, rMax, cLat, cLon, beta, rMax2, beta1, beta2)

    def test_Holland(self):
        """Testing Holland profile """
        pHolland = self.prP.holland()
        self.numpyAssertAlmostEqual(pHolland, self.test_pHolland)

    def test_Willoughby(self):
        """Testing Willoughby profile """
        pWilloughby = self.prP.willoughby()
        self.numpyAssertAlmostEqual(pWilloughby, self.test_pWilloughby)

    def test_doubleHolland(self):
        """Testing Double Holland profile """
        pdoubleHolland = self.prP.doubleHolland()
        self.numpyAssertAlmostEqual(pdoubleHolland, self.test_pdoubleHolland)

    def test_Powell(self):
        """Testing Powell profile """
        pPowell = self.prP.powell()
        self.numpyAssertAlmostEqual(pPowell, self.test_pPowell)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestPressureProfile,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
