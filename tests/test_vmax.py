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


 Title: TestVMax.Py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2006-10-24
 Description: Unit test for vmax.py

 Version: $Rev: 563 $

 ModifiedBy: N. Habili
 ModifiedDate: 2007-04-27
 Modification: vMax and pDiff are now two separate tests

 ModifiedBy: N. Habili
 ModifiedDate: 2007-04-27
 Modification: Added tests for exceptions


 $Id: TestVMax.py 563 2007-10-24 02:52:40Z carthur $
"""
import os, sys
from scipy import arange
import unittest
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
import wind.vmax as vmax
from Utilities.files import flStartLog


class TestVMax(unittest.TestCase):

    pc = arange(88000, 101100, 1000)
    pe = 101000
    vmax = {"will":[71.283768,68.487229,65.571529,62.52,59.31168,55.919588,52.307985,48.427784,44.208316,39.54112,34.243614,27.959794,19.77056,0],
            "holl":[73.527058,70.642513,67.635057,64.487496,61.178211,57.67937,53.95411,49.9518,45.599546,40.785474,35.321256,28.839685,20.392737,0],
            "atkin":[69.863944556,66.3538774077,62.7379715373,59.0029306788,55.1322554425,51.1050167812,46.8939397787,42.4622486228,37.7581116559,
                    32.7039505937,27.1731303211,20.9284420943,13.3928678506,0.0]}

    def test_vMax(self):
        """Testing vMax"""

        #testing willoughby
        for i in range(self.pc.size):
            self.assertAlmostEqual(vmax.vmax(self.pc[i], self.pe, "willoughby"), self.vmax["will"][i], 3)

        #testing holland
        for i in range(self.pc.size):
            self.assertAlmostEqual(vmax.vmax(self.pc[i], self.pe, "holland"), self.vmax["holl"][i], 3)

        #testing atkinson
        for i in range(self.pc.size):
            self.assertAlmostEqual(vmax.vmax(self.pc[i], self.pe, "atkinson"), self.vmax["atkin"][i], 3)

        #testing exceptions
        self.assertRaises(ValueError, vmax.vmax, 1500.0, 1000.0)
        self.assertRaises(NotImplementedError, vmax.vmax, 1000.0, 1500.0, "will")


    def test_pDiff(self):
        """Testing pDiff"""

        #testing willoughby
        for i in range(self.pc.size):
            self.assertAlmostEqual(vmax.pDiff(self.vmax["will"][i], self.pe, "willoughby"), self.pc[i], 3)

        #testing holland
        for i in range(self.pc.size):
            self.assertAlmostEqual(vmax.pDiff(self.vmax["holl"][i], self.pe, "holland"), self.pc[i], 3)

        #testing atkinson
        for i in range(self.pc.size):
            self.assertAlmostEqual(vmax.pDiff(self.vmax["atkin"][i], self.pe, "atkinson"), self.pc[i], 3)

        #testing exceptions
        self.assertRaises(NotImplementedError, vmax.pDiff, 1000.0, 1500.0, "will")

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestVMax,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
