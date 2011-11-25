#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
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


 Title: TestWindField.py
 Author: Nicholas Summons, nicholas.summons@ga.gov.au
 CreationDate: 2011-06-09
 Description: Unit testing module for windField

 Version: $Rev$

 $Id$
"""
import os, sys
import cPickle
import unittest
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from WindfieldInterface import windField
from Utilities.files import flStartLog


class TestWindField(NumpyTestCase.NumpyTestCase):

    pkl_file = open(os.path.join(unittest_dir, 'test_data', 'windField_testdata.pck'), 'rb')
    R = cPickle.load(pkl_file)
    lam = cPickle.load(pkl_file)
    rMax = cPickle.load(pkl_file)
    f = cPickle.load(pkl_file)
    V = cPickle.load(pkl_file)
    Z = cPickle.load(pkl_file)
    vFm = cPickle.load(pkl_file)
    thetaFm = cPickle.load(pkl_file)
    thetaMax = cPickle.load(pkl_file)
    test_kepert_Ux = cPickle.load(pkl_file)
    test_kepert_Vy = cPickle.load(pkl_file)
    test_mcconochie_Ux = cPickle.load(pkl_file)
    test_mcconochie_Vy = cPickle.load(pkl_file)    
    test_hubbert_Ux = cPickle.load(pkl_file)
    test_hubbert_Vy = cPickle.load(pkl_file)    
    pkl_file.close()

    wf = windField.WindField(R, lam, rMax, f, V, Z, vFm, thetaFm, thetaMax)

    def test_Kepert(self):
        """Testing Kepert windfield"""
        kepert_Ux, kepert_Vy = self.wf.kepert()
        self.numpyAssertAlmostEqual(kepert_Ux, self.test_kepert_Ux)
        self.numpyAssertAlmostEqual(kepert_Vy, self.test_kepert_Vy)

    def test_McConochie(self):
        """Testing McConochie windfield"""
        mcconochie_Ux, mcconochie_Vy = self.wf.mcconochie()
        self.numpyAssertAlmostEqual(mcconochie_Ux, self.test_mcconochie_Ux)
        self.numpyAssertAlmostEqual(mcconochie_Vy, self.test_mcconochie_Vy)

    def test_Hubbert(self):
        """Testing Hubbert windfield"""
        hubbert_Ux, hubbert_Vy = self.wf.hubbert()
        self.numpyAssertAlmostEqual(hubbert_Ux, self.test_hubbert_Ux)
        self.numpyAssertAlmostEqual(hubbert_Vy, self.test_hubbert_Vy)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestWindField,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
