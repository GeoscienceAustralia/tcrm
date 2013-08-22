#!/usr/bin/env python
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


 Title: KDEOriginTestCase.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2006-12-13
 Description: Unit testing module for KDEOrigin.py

 Version: $Rev: 563 $

 $Id: TestKDEOrigin.py 563 2007-10-24 02:52:40Z carthur $
"""
import os, sys, pdb
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
from StatInterface import KDEOrigin
from Utilities.files import flStartLog

class TestKDEOrigin(NumpyTestCase.NumpyTestCase):
    
    def setUp(self):
        kdeType = 'Epanechnikov'
        gridLimit={'xMin':70, 'xMax':180, 'yMin':-36, 'yMax':0}
        kdeStep = 0.1
        lonLat = cPickle.load(open(os.path.join(unittest_dir, 'test_data', 'kde_origin_lonLat.pck')))
        self.kdeOrigin = KDEOrigin.KDEOrigin(None, kdeType, gridLimit, kdeStep, lonLat)

    def test_GenerateKDE(self):
        """Testing GenerateKDE for 2D data"""
        pkl_file = open(os.path.join(unittest_dir, 'test_data', 'kde_origin_xyz.pck'), 'rb')
        xp = cPickle.load(pkl_file)
        yp = cPickle.load(pkl_file)
        zp = cPickle.load(pkl_file)
        zp = zp.transpose()
        zp = zp / zp.sum()  # This line was added to match the normalisation added to KDEOrigin
        x, y, z = self.kdeOrigin.generateKDE()
        self.numpyAssertAlmostEqual(xp, x)
        self.numpyAssertAlmostEqual(yp, y)
        self.numpyAssertAlmostEqual(zp, z)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestKDEOrigin,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
