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


 Title: TestDerivatives.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2007-05-03
 Description: Unit testing module for Derivatives

 Version: $Rev: 717 $

 $Id: test_derivatives.py 717 2011-11-08 06:48:02Z nsummons $
"""
import os, sys
import unittest
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from WindfieldInterface import derivatives
from Utilities.files import flStartLog

class TestDerivatives(NumpyTestCase.NumpyTestCase):

    def test_Holland(self):
        """Testing Holland profile
        """
        d2Vm = derivatives.holland(1.0, 10, 1.3, 1000, 1.15)
        self.assertAlmostEqual(d2Vm, -0.156130084588)

    def test_doubleHolland(self):
        """Testing doubleHolland profile
        """
        d2Vm = derivatives.doubleHolland(1.0, 10, 10,1.3,1.3, 1000, 9000, 1.15)
        self.assertAlmostEqual(d2Vm, -0.242010340982)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestDerivatives,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)