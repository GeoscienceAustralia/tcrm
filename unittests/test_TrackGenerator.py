#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0
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


 Title: TestTrackGenerator.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2007-05-03
 Description: Unit testing module for TrackGenerator

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
import Utilities.my_tool as myutils
from TrackGenerator import TrackGenerator
from Utilities.files import flStartLog


class TestTrackGenerator(NumpyTestCase.NumpyTestCase):
    
    gridLimit = {'xMin':70, 'xMax':180, 'yMin':-40, 'yMax':0}
    gridSpace = {'x':10, 'y':10}
    numCyclones = 100
    dt = 1
    pEnv = 1010
    tsteps = 20
    timeOverflow = tsteps * dt
    distanceOverflow = 1500

    cfgFiles = myutils.loadConfig(open(os.path.join(unittest_dir, 'test_data', 'files.ini')))

    tracks = TrackGenerator.TrackGenerator(cfgFiles, numCyclones, dt, gridLimit, gridSpace, pEnv, timeOverflow, distanceOverflow, tsteps)

    def test_GeneratePath(self):
        """Testing GenerateSamples"""
        random.seed(10)
        stp = cPickle.load(open(os.path.join(unittest_dir, 'test_data', 'track_generation_sample_tracks.pck')))
        st = self.tracks.generatePath()
        self.numpyAssertAlmostEqual(stp, st)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestTrackGenerator,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
