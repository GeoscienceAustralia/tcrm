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


 Title: TestSamplingParameters.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2007-05-03
 Description: Unit testing module for SamplingParameters.py

 Version: $Rev$

 $Id$
"""
import os, sys
import pickle
import unittest
from scipy import random
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from StatInterface import SamplingParameters
from Utilities.files import flStartLog


class TestSamplingParameters(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.numberOfSamples = 1000
        pr = pickle.load(open(os.path.join(unittest_dir, 'test_data', 'sampling_parameters_xacy.pkl'), 'rb'))
        self.sampPar = SamplingParameters.SamplingParameters(pr)
        random.seed(10)

    def test_GenerateSamples(self):
        """Testing GenerateSamples"""
        samples = pickle.load(open(os.path.join(unittest_dir, 'test_data', 'sampling_parameters_samples.pkl'), 'rb'))
        samplesp = self.sampPar.generateSamples(self.numberOfSamples)
        self.numpyAssertAlmostEqual(samples, samplesp)

    def test_GenerateOneSample(self):
        """Testing GenerateOneSample"""
        sp = 0.2
        s = self.sampPar.generateOneSample()
        self.assertAlmostEqual(sp, s)

if __name__ == "__main__":
    #flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestSamplingParameters,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
