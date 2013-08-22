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


 Title: TestMSLP.py
 Author: Nicholas Summons, nicholas.summons@ga.gov.au
 CreationDate: 2011-06-09
 Description: Unit testing module for mslp_seasonal_clim

 Version: $Rev$

 $Id$
"""
import os, sys
import cPickle
import unittest
import numpy
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from MSLP.mslp_seasonal_clim import MSLPGrid
from Utilities.files import flStartLog

class TestMSLP(NumpyTestCase.NumpyTestCase):

    mslpSeasonalAv = MSLPGrid([12,1,2])
    mslpSA_test = numpy.array([1009.8081665, 1010.480896, 1010.38543701, 1010.22454834,
                               1007.98876953, 1009.05102539, 1013.85418701], 'float')

    def test_mslp(self):
        """Testing mslp_seasonal_clim"""
        mslpSA = self.mslpSeasonalAv.sampleGrid([150, 160, 170, 180, 180, 180, 180], [-20, -20, -20, -20, 0, 10, 20])
        self.numpyAssertAlmostEqual(mslpSA, self.mslpSA_test)

if __name__ == "__main__":
    #flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestMSLP,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
