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
"""

import os, sys, pdb, logging
import unittest
import cPickle
import NumpyTestCase
import numpy
try:
    import pathLocate
except:
    from unittests import pathLocate

unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import nctools
from Utilities import grid
from Utilities.files import flStartLog

class TestSampleGrid(NumpyTestCase.NumpyTestCase):
    """test_SampleGrid:

    Description: Test that the range of methods to load gridded data
                produce expected results. Uses the 0.083 degree land-sea
                mask dataset as a test dataset.

    Parameters:
    Members:
    Methods:
    Internal methods:
    """
    def setUp(self):
        self.filename = os.path.join(unittest_dir,'test_data','landmask.nc')
        # Load the data using grid.grdRead:
        self.lslon,self.lslat,self.lsgrid = grid.grdRead(self.filename)
        # Load the data using nctools.ncLoadFile and nctools.ncGetData:
        ncobj = nctools.ncLoadFile(self.filename)
        self.nclon = nctools.ncGetDims(ncobj,'lon')
        self.nclat = nctools.ncGetDims(ncobj,'lat')
        self.ncgrid = nctools.ncGetData(ncobj,'landmask')
        ncobj.close()
        # Set up an instance of SampleGrid:
        self.sample = grid.SampleGrid(self.filename)
        # Sample the land-sea mask at these points around the globe:
        self.xlon = [100.,130.,180.,250.,300.]
        self.ylat = [-80.,-20.,40.]
        # Known point values of the land-sea mask data:
        self.ls = [3.,0.,3.,3.,3.,0.,3.,0.,0.,3.,0.,3.,3.,3.,0.]

    def test_gridNctools(self):
        """Test grid.grdRead produces same result as nctools.ncGetData"""
        self.numpyAssertEqual(self.lslon,self.nclon)
        self.numpyAssertEqual(self.lslat,self.nclat)
        self.numpyAssertEqual(self.lsgrid,self.ncgrid)

    def test_sampleGridArray(self):
        """Test grid.SampleGrid.grid produces flipped result of nctools.ncGetData"""
        self.numpyAssertEqual(self.sample.lon,self.nclon)
        self.numpyAssertEqual(self.sample.lat,self.nclat)
        self.numpyAssertEqual(self.sample.grid,numpy.flipud(self.ncgrid))

    def test_sampleGridPoints(self):
        """Test grid.SampleGrid returns correct values for known points"""
        i = 0
        for x in self.xlon:
            for y in self.ylat:
                self.assertEqual(self.sample.sampleGrid(x,y),self.ls[i])
                i +=1

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestSampleGrid,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
