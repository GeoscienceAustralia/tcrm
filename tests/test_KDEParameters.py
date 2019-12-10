"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
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


 Title: TestKDEParameters.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2007-05-07
 Description: Unit testing module for TestKDEParameters.py

 Version: $Rev: 480 $

 $Id: TestKDEParameters.py 480 2007-05-01 05:03:15Z nhabili $
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
from StatInterface import KDEParameters
from Utilities.files import flStartLog

class TestKDEParameters(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.pressure_rate = pickle.load(open(os.path.join(unittest_dir, 'test_data', 'kde_parameters_pressure_rate.pkl'), 'rb'))
        self.resultp = pickle.load(open(os.path.join(unittest_dir, 'test_data', 'kde_parameters_result.pkl'), 'rb'))
        kdeType = 'gau'
        self.k = KDEParameters.KDEParameters(kdeType)

    def test_GenerateKDE(self):
        """Testing GenerateKDE for 1-D data"""
        kdeStep = 0.1
        result = self.k.generateKDE(self.pressure_rate, kdeStep)
        self.numpyAssertAlmostEqual(self.resultp, result)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestKDEParameters,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
