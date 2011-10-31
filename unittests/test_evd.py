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


 Title: TestEvd.py
 Author: Nicholas Summons, nicholas.summons@ga.gov.au
 CreationDate: 2011-06-08
 Description: Unit testing module for evd

 Version: $Rev$

 $Id$
"""
import os, sys
import cPickle
import unittest
import numpy
import NumpyTestCase

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))
from HazardInterface import evd
from Utilities.files import flStartLog

class TestEvd(NumpyTestCase.NumpyTestCase):

    v = numpy.array([0., 0., 38.55, 41.12, 59.29, 61.75, 74.79], dtype='float32')
    yrs = numpy.array([25.0, 50.0, 100.0, 250.0, 500.0, 2000.0])
    test_w = numpy.array([59.26156235, 69.34857941, 76.71388245, 84.10202789, 88.47135925, 95.00366974], dtype='float32')
    test_loc, test_scale, test_shp = numpy.array([49.2291362594, 16.3463688259, 0.272970209861], dtype='float32')
    missingValue = numpy.array(-9999, dtype='float32')

    def test_evd(self):
        """Testing evd"""
        w, loc, scale, shp = evd.estimate_EVD(self.v, self.yrs, missingValue=-9999, minRecords=3, yrspersim=10)
        self.numpyAssertAlmostEqual(w, self.test_w)
        self.numpyAssertAlmostEqual(loc, self.test_loc)
        self.numpyAssertAlmostEqual(scale, self.test_scale)
        self.numpyAssertAlmostEqual(shp, self.test_shp)
        w2, loc2, scale2, shp2 = evd.estimate_EVD(self.v, self.yrs, missingValue=-9999, minRecords=50, yrspersim=10)
        self.numpyAssertAlmostEqual(w2, numpy.ones(6, dtype='float32') * self.missingValue)
        self.numpyAssertAlmostEqual(loc2, self.missingValue)
        self.numpyAssertAlmostEqual(scale2, self.missingValue)
        self.numpyAssertAlmostEqual(shp2, self.missingValue)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestEvd,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)