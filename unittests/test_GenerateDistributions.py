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


 Title: TestGenerateDistributions.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2007-05-03
 Description: Unit testing module for GenerateDistributions

 Version: $Rev$

 $Id$
"""
import os, sys
import cPickle
import unittest
from scipy import random
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from StatInterface import GenerateDistributions
from Utilities.files import flStartLog

class TestGenerateDistributions(NumpyTestCase.NumpyTestCase):
    
    def setUp(self):
            gridLimit = {'xMin':70, 'xMax':180, 'yMin':-40, 'yMax':0}
            gridSpace = {'x':10, 'y':10}
            gridInc = {'x':2,'y':1}
            kdeType = "Epanechnikov"
            kdeStep = 0.1
            minSamplesCell = 10

            self.gDist = GenerateDistributions.GenerateDistributions(None, gridLimit, gridSpace, gridInc, kdeType, kdeStep, minSamplesCell)

    def test_GenerateSamples(self):
        """Testing GenerateSamples"""
        init_lon_lat = cPickle.load(open(os.path.join(unittest_dir, 'test_data', 'generate_distributions_init_lon_lat.pck')))
        init_pressure = cPickle.load(open(os.path.join(unittest_dir, 'test_data', 'generate_distributions_init_pressure.pck')))
        distp = cPickle.load(open(os.path.join(unittest_dir, 'test_data', 'generate_distributions_dist.pck')))
        
        dist = self.gDist.allDistributions(init_lon_lat, init_pressure)
        self.numpyAssertAlmostEqual(dist, distp)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestGenerateDistributions,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)