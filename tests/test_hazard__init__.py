"""
Testing the hazard module __init__ functions

 Title: test_hazard__init__.py
 Author: Claire Krause, claire.krause@ga.gov.au
 CreationDate: 2017-09-08
 Description: Unit testing for hazard module

"""
import os
import sys
import unittest
import pickle
from . import NumpyTestCase
import numpy

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))

from hazard import *
from numpy.testing import assert_raises

class TestloadFile(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.filepath = os.path.join(unittest_dir, 'test_data/folder_of_files/')
        self.filename = self.filepath + 'gust.000-00000.nc'
        self.expectedsubset = numpy.ma.core.MaskedArray(
            [[ 19.33163452,  19.36814117,  19.40473557,  19.44141579, 19.47817612], 
             [ 19.33328438,  19.37013245,  19.40707016,  19.4440937, 19.48120117],
             [19.33486748,  19.37206078,  19.40934563,  19.44671822, 19.48417473],
             [ 19.33638191,  19.37392426,  19.41156006,  19.44928551, 19.48709679],
             [ 19.33782578,  19.37572098,  19.41371346,  19.45179558, 19.48996544]],
            dtype='float32')
        pfile = open(os.path.join(unittest_dir, 'test_data', 'testDomain.pkl'),'rb')
        self.wf_lat = pickle.load(pfile)
        self.wf_lon = pickle.load(pfile)
        pfile.close()

    def testloadFile(self):
        """ test the loadFile function """
        limits = (500, 505, 500, 505)
        data_subset = loadFile(self.filename, limits)
        self.numpyAssertAlmostEqual(data_subset, self.expectedsubset)

        badfilename = self.filename + '1'
        assert_raises(IOError, loadFile, badfilename, limits)


    def testsetDomain(self):
        """ test setDomain """
        wf_lon, wf_lat = setDomain(self.filepath)
        self.numpyAssertAlmostEqual(wf_lon, self.wf_lon)
        self.numpyAssertAlmostEqual(wf_lat, self.wf_lat)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestloadFile,'test')
    unittest.TextTestRunner(verbosity=3).run(testSuite)
