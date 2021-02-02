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


 Title: test_GPD.py
 Author: Claire Krause, claire.krause@ga.gov.au
 CreationDate: 2017-08-30
 Description: Unit testing module for GPD.py

 Version:

 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id:
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
from hazard import GPD
from Utilities.files import flStartLog

class TestGPD(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        pfile = open(os.path.join(unittest_dir, 'test_data', 'testGPD.pkl'),'rb')
        self.values = pickle.load(pfile)
        self.mu = pickle.load(pfile)
        self.scale = pickle.load(pfile)
        self.shape = pickle.load(pfile)
        self.w = pickle.load(pfile)
        self.rate = pickle.load(pfile)

        self.years = numpy.array((2,5,10,20,25,50,100,200,250,500,1000,2000,5000,10000))

        pfile.close()


    def test_gpdfit(self):
        """Test gpdfit works with list input"""
        # set numSim = 100
        w, mu, scale, shape = GPD.gpdfit(self.values, self.years, 100)
        self.numpyAssertAlmostEqual(w, self.w)
        self.numpyAssertAlmostEqual(mu, self.mu)
        self.numpyAssertAlmostEqual(scale, self.scale)
        self.numpyAssertAlmostEqual(shape, self.shape)
        
    def test_gpdReturnLevel(self):
        """ Test the GDP return level calculations work """
        w = GPD.gpdReturnLevel(self.years, self.mu, self.shape, self.scale, self.rate)
        self.numpyAssertAlmostEqual(w, self.w)


if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestGPD,'test')
    unittest.TextTestRunner(verbosity=3).run(testSuite)
