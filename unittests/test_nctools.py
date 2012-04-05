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


 Title: testNetCDF.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 07/15/08 1:46:PM
 Description: Unit test for nctools.py

 Version :$Rev: 276 $

 $Id: testNetCDF.py 276 2010-04-16 02:24:00Z nsummons $
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

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import nctools
from Utilities.files import flStartLog

class TestNetCDF(NumpyTestCase.NumpyTestCase):
    """
    TestNetCDF: test the suite of netCDF tools
    """
    def setUp(self):
        self.ncfile=os.path.join(unittest_dir,'test_data','pres_temp_4D.nc')
        self.nrecs = 2; self.nlevs = 2; self.nlats = 6; self.nlons = 12
        self.lats_out = -25.0 + 5.0*numpy.arange(self.nlats,dtype='float')
        self.lons_out = 125.0 + 5.0*numpy.arange(self.nlons,dtype='float')
        self.latrange = [min(self.lats_out), max(self.lats_out)]
        self.lonrange = [min(self.lons_out), max(self.lons_out)]
        self.lat_atts = {'units':'degrees_north',
                         'long_name':'latitude',
                         'actual_range': self.latrange}
        self.lon_atts = {'units':'degrees_east',
                         'long_name':'longitude',
                         'actual_range': self.lonrange}

    def test_ncCreateFile(self):
        """Test nctools creates a file successfully"""
        ncobj = nctools.ncCreateFile(self.ncfile)
        # output data.
        press_out = 900. + numpy.arange(self.nlevs*self.nlats*self.nlons,dtype='f') # 1d array
        press_out.shape = (self.nlevs,self.nlats,self.nlons) # reshape to 2d array
        temp_out = 9. + numpy.arange(self.nlevs*self.nlats*self.nlons,dtype='f') # 1d array
        temp_out.shape = (self.nlevs,self.nlats,self.nlons) # reshape to 2d array
        # create the lat and lon dimensions.
        nctools.ncCreateDim(ncobj, 'lat', self.lats_out,'f',self.lat_atts)
        nctools.ncCreateDim(ncobj, 'lon', self.lons_out,'f',self.lon_atts)
        # create level dimension.
        nctools.ncCreateDim(ncobj, 'level', numpy.arange(self.nlevs),'f')
        # create time dimension (record, or unlimited dimension)
        nctools.ncCreateDim(ncobj, 'time', numpy.arange(self.nrecs),'f')
        dimensions = ('time','level','lat','lon')
        press = nctools.ncCreateVar(ncobj, 'pressure', dimensions, 'f')
        temp = nctools.ncCreateVar(ncobj, 'temperature', dimensions, 'f')
        press.units =  'hPa'
        temp.units = 'celsius'
        for nrec in range(self.nrecs):
            press[nrec,:,::] = press_out
            temp[nrec,:,::] = temp_out
        ncobj.close()

    def test_ncReadFile(self):
        """Test nctools functions for reading dimensions and variables"""
        ncobj = nctools.ncLoadFile(self.ncfile)
        lats_check = -25.0 + 5.0*numpy.arange(self.nlats,dtype='float')
        lons_check = 125.0 + 5.0*numpy.arange(self.nlons,dtype='float')
        press_check = 900. + numpy.arange(self.nlevs*self.nlats*self.nlons,dtype='float') # 1d array
        press_check.shape = (self.nlevs,self.nlats,self.nlons) # reshape to 2d array
        temp_check = 9. + numpy.arange(self.nlevs*self.nlats*self.nlons,dtype='float') # 1d array
        temp_check.shape = (self.nlevs,self.nlats,self.nlons) # reshape to 2d array
        lats = nctools.ncGetDims(ncobj,'lat')
        lons = nctools.ncGetDims(ncobj,'lon')
        self.numpyAssertEqual(lats_check,lats)
        self.numpyAssertEqual(lons_check,lons)
        press = nctools.ncGetData(ncobj, 'pressure')
        temp = nctools.ncGetData(ncobj, 'temperature')
        for nrec in range(self.nrecs):
            self.numpyAssertEqual(press_check,press[nrec])
            self.numpyAssertEqual(temp_check,temp[nrec])

    def test_ncLoadNonFile(self):
        """Test ncLoadFile raises IOError when attempting to load non-file"""
        self.assertRaises(IOError,nctools.ncLoadFile,'file')

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestNetCDF,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
