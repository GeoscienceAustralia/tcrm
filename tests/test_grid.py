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


 Title: testGrid.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 07/15/08 1:46:PM
 Description: Unit test for grid.py

 Version :$Rev: 810 $

 $Id: test_grid.py 810 2012-02-21 07:52:50Z nsummons $
"""
import os, sys, pdb
import unittest
import pickle
from . import NumpyTestCase
import numpy
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import grid
from Utilities.files import flStartLog


class TestGrid(NumpyTestCase.NumpyTestCase):
    """TestReadGrid:

    Description:

    """
    def setUp(self):
        self.gridfile = os.path.join(unittest_dir, 'test_data', 'mslp.txt')
        self.ncfile = os.path.join(unittest_dir,'test_data','sfc_pres_temp.nc')
        self.gridobj = grid.SampleGrid(self.gridfile)
        self.ncgridobj = grid.SampleGrid(self.ncfile)
        self.ilon = 150.0
        self.ilat = -25.0

        self.nlons = 12; self.nlats = 6
        self.lon_check = 125.0 + 5.0*numpy.arange(self.nlons,dtype='float')
        self.lat_check = -25.0 + 5.0*numpy.arange(self.nlats,dtype='float')
        self.press_check = 900. + numpy.arange(self.nlats*self.nlons,dtype='float')
        self.press_check.shape = (self.nlats,self.nlons)

    def test_grdRead(self):
        """Test grid data is read correctly from ascii file"""
        pfile = open(os.path.join(unittest_dir, 'test_data', 'gridReadTestData.pkl'),'rb')
        pdata = pickle.load(pfile)
        plon = pickle.load(pfile)
        plat = pickle.load(pfile)
        pfile.close()
        lon, lat, data = grid.grdRead(self.gridfile)
        self.numpyAssertAlmostEqual(pdata, data)

    def test_grdReadNetCDF(self):
        """Test grid data is read correctly from netCDF file"""
        lons,lats,press = grid.grdRead(self.ncfile)
        self.numpyAssertAlmostEqual(self.lon_check,lons)
        self.numpyAssertAlmostEqual(self.lat_check,lats)
        self.numpyAssertAlmostEqual(self.press_check,press)

    def test_SampleGrid(self):
        """Test SampleGrid class using an ascii file as input"""
        pfile = open(os.path.join(unittest_dir, 'test_data', 'samplegrid.pkl'),'rb')
        pvalue = 1015.018
        pfile.close()
        value = self.gridobj.sampleGrid(self.ilon, self.ilat)
        self.numpyAssertAlmostEqual(numpy.array(pvalue), numpy.array(value))

    def test_SampleGridNetCDF(self):
        """Test SampleGrid class using a netCDF file as input"""
        pvalue = 965.0
        value = self.ncgridobj.sampleGrid(self.ilon,self.ilat)
        self.assertAlmostEqual(pvalue,value)

if __name__ == "__main__":
    unittest.main()
