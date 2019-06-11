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


 Title: test_lmomentFit.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2012-11-21
 Description: Unit testing module for lmomentFit.py

 Version: $Rev: 737 $

 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: test_lmomentFit.py 737 2012-11-21 00:59:42Z carthur $
"""

import os
import sys
import unittest
import pickle
from . import NumpyTestCase
import numpy

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))
from Utilities import lmomentFit as lmom
from Utilities.files import flStartLog

class Testlmoments(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pfile = open(os.path.join(unittest_dir, 'test_data', 'testlmom.pkl'),'rb')
        self.values = pickle.load(pfile)
        self.moments = pickle.load(pfile)
        self.params = pickle.load(pfile)

        pfile.close()


    def test_samlmu_list(self):
        """Test samlmu works with list input"""
        moments = lmom.samlmu(list(self.values),3)
        self.numpyAssertAlmostEqual(moments,self.moments[0:3])

    def test_samlmu_mom(self):
        """Test samlmu returns correct values for moments"""
        moments = lmom.samlmu(self.values,5)
        self.numpyAssertAlmostEqual(moments,self.moments)

    def test_samlmu3(self):
        """Test samlmu3 returns same values as samlmu"""
        moments = lmom.samlmu(self.values,3)
        self.numpyAssertAlmostEqual(moments,self.moments[0:3])

    def test_pelgev(self):
        """Test pelgev returns correct GEV parameters"""
        l1 = self.moments[0]
        l2 = self.moments[1]
        t3 = self.moments[2]/self.moments[1]
        xmom = [l1,l2,t3]
        params = lmom.pelgev(xmom)
        self.numpyAssertAlmostEqual(params,self.params)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(Testlmoments,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
