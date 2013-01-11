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


 Title: test_files.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2009-11-27 9:46:AM
 Description: Unit test for files.py

 Throughout, I convert all directory separators to *NIX style
 separators (i.e. the Windows style '\\' become '/')

 PLEASE READ THE UNITTEST DOCUMENTATION AT
 http://docs.python.org/library/unittest.html
 before adding more tests

 Version :$Rev: 733 $

 $Id: test_GetType.py 733 2012-11-14 02:52:01Z carthur $
"""
import os
import sys
import pdb
import logging
import unittest
import numpy
import NumpyTestCase

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))
import GetType
from files import flStartLog

class Test_GetType(NumpyTestCase.NumpyTestCase):
    Types = GetType.GetType()
    typedict = {'age':'float',
                'bearing':'float',
                'date':'string',
                'day':'float',
                'direction':'float',
                'hour':'float',
                'index':'int',
                'intensity':'string',
                'lat':'float',
                'lon':'float',
                'maxspeed':'float',
                'minute':'float',
                'month':'float',
                'mpi':'float',
                'name':'string',
                'num':'int',
                'penv':'float',
                'pressure':'float',
                'rmax':'float',
                'season':'float',
                'second':'float',
                'shear':'float',
                'skip':'string',
                'speed':'float',
                'surfcode':'string',
                'tcserialno':'string',
                'unit':'string',
                'vmax':'float',
                'year':'float'}

    def test_getType(self):
        """Test the getType() function"""
        for t,v in self.typedict.iteritems():
            retval = GetType.GetType.getType(self.Types,t)
            self.assertEqual(retval,v)
        

if __name__ == "__main__":
    flStartLog('','CRITICAL',False)
    testSuite = unittest.makeSuite(Test_GetType,'test')
    
    unittest.TextTestRunner(verbosity=2).run(testSuite)
