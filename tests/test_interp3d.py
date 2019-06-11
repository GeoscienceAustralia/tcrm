"""
 Title: test_interp3d.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2013-01-17
 Description: Unit test suite for interp3d.py

 Version: $Rev$

 Id: $Id$
"""
import os
import sys
import unittest
import pickle
from . import NumpyTestCase
import numpy

unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))
from Utilities import nctools
from Utilities import interp3d
from Utilities.files import flStartLog

class TestInterp3d(NumpyTestCase.NumpyTestCase):
    """
    TestInterp3d:
    Description:
    Parameters:
    Members:
    Methods:
    Internal methods:
    """
    def setUp(self):
        self.filename = os.path.join(unittest_dir, 'test_data', 'mslp_ltm.nc')
        self.ncobj = nctools.ncLoadFile(self.filename)
        self.lon = nctools.ncGetDims(self.ncobj, 'lon')
        self.lat = nctools.ncGetDims(self.ncobj, 'lat')
        self.time = nctools.ncGetDims(self.ncobj, 'time')
        self.data = nctools.ncGetData(self.ncobj, 'mslp')
        self.ncobj.close()
        self.scale = [365.,-180.,360.]
        self.offset = [0., 90., 0.]

        # Load data
        pfile = open(os.path.join(unittest_dir, 'test_data',
                                  'testinterp3d.pkl'), 'rb')
        self.xlon = pickle.load(pfile)
        self.ylat = pickle.load(pfile)
        self.ztime = pickle.load(pfile)
        self.values = pickle.load(pfile)
        pfile.close()
        self.coords = numpy.array([self.ztime, self.ylat, self.xlon])


    def test_interp3d(self):
        """Test interp3d returns expected values"""
        output = interp3d.interp3d(self.data, self.coords,
                                   self.scale, self.offset)
        self.numpyAssertAlmostEqual(output, self.values)

    def test_interp3d_2dinput(self):
        """Test interp3d raises ValueError if input_array is 2D"""
        self.assertRaises(ValueError,
                          interp3d.interp3d,
                          self.data[1, :, :],
                          self.coords,
                          self.scale,
                          self.offset)



if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestInterp3d,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
