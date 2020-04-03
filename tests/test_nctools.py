"""
 Title: testNetCDF.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 07/15/08 1:46:PM
 Description: Unit test for nctools.py

 Version :$Rev: 276 $

 $Id: testNetCDF.py 276 2010-04-16 02:24:00Z nsummons $
"""
import os
import sys
from os.path import join as pjoin
import unittest
from . import NumpyTestCase
import numpy as np
import netCDF4
from datetime import datetime, timedelta

try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import nctools
from netCDF4 import Dataset

class TestNetCDF(NumpyTestCase.NumpyTestCase):
    """
    TestNetCDF: test the suite of netCDF tools
    """
    def setUp(self):
        self.ncfile = pjoin(unittest_dir, 'test_data', 'pres_temp_4D.nc')

        self.nrecs = 2; self.nlevs = 2; self.nlats = 6; self.nlons = 12
        self.lats_out = -25.0 + 5.0*np.arange(self.nlats, dtype='float')
        self.lons_out = 125.0 + 5.0*np.arange(self.nlons, dtype='float')
        self.latrange = [min(self.lats_out), max(self.lats_out)]
        self.lonrange = [min(self.lons_out), max(self.lons_out)]
        self.lat_atts = {'units': 'degrees_north',
                         'long_name': 'latitude',
                         'actual_range': self.latrange}
        self.lon_atts = {'units': 'degrees_east',
                         'long_name': 'longitude',
                         'actual_range': self.lonrange}
        self.time_atts = {'units': 'hours since 2000-01-01 00:00:00',
                          'calendar': 'standard'}

        press_out = 900. + np.arange(self.nrecs * self.nlevs * \
                                     self.nlats*self.nlons,
                                     dtype='f')
        press_out.shape = (self.nrecs, self.nlevs,
                           self.nlats, self.nlons)
        temp_out = 9. + np.arange(self.nrecs * self.nlevs * \
                                  self.nlats * self.nlons,
                                  dtype='f')
        temp_out.shape = (self.nrecs, self.nlevs,
                          self.nlats, self.nlons)


        self.dimensions = {
            0: {
                'name': 'time',
                'values': np.arange(self.nrecs),
                'dtype': int,
                'atts': self.time_atts
                },
            1: {
                'name': 'level',
                'values': np.arange(self.nlevs),
                'dtype': 'f',
                'atts': {}
                },
            2: {
                'name': 'lat',
                'values': self.lats_out,
                'dtype': 'f',
                'atts': self.lat_atts
                },
            3: {
                'name': 'lon',
                'values': self.lons_out,
                'dtype': 'f',
                'atts': self.lon_atts
                }
            }

        self.variables = {
            0: {
                'name': 'pressure',
                'values': press_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'hPa',
                    'standard_name': 'air_pressure'
                    }
                },
            1: {
                'name': 'temperature',
                'values': temp_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'degrees_celcius',
                    'standard_name': 'air_temperature'
                    }
                }
            }

        # remove a required key from one of the dimensions:
        self.missing_key_dims = {
            0: {
                'name': 'time',
                'values': self.nrecs,
                'atts': self.time_atts
                },
            1: {
                'name': 'level',
                'values': np.arange(self.nlevs),
                'dtype': 'f',
                'atts': {}
                },
            2: {
                'name': 'lat',
                'values': self.lats_out,
                'dtype': 'f',
                'atts': self.lat_atts
                },
            3: {
                'name': 'lon',
                'values': self.lons_out,
                'dtype': 'f',
                'atts': self.lon_atts
                }
            }

        self.missing_key_vars = {
            0: {
                'name': 'pressure',
                'values': press_out,
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'hPa',
                    'standard_name': 'air_pressure'
                    }
                },
            1: {
                'name': 'temperature',
                'values': temp_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'degrees_celcius',
                    'standard_name': 'air_temperature'
                    }
                }
            }

        self.misshapen_vars = {
            0: {
                'name': 'pressure',
                'values': press_out,
                'dtype': 'float64',
                'dims': ('time', 'lat', 'lon'),
                'atts': {
                    'units': 'hPa',
                    'standard_name': 'air_pressure'
                    }
                },
            1: {
                'name': 'temperature',
                'values': temp_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'degrees_celcius',
                    'standard_name': 'air_temperature'
                    }
                }
            }

        self.nullvalue_var = {
            0: {
                'name': 'pressure',
                'values': None,
                'dtype': 'float64',
                'dims': (),
                'atts': {
                    'units': 'hPa',
                    'standard_name': 'air_pressure'
                    }
                },
            1: {
                'name': 'temperature',
                'values': temp_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'degrees_celcius',
                    'standard_name': 'air_temperature'
                    }
                }
            }

    def tearDown(self):
        if os.path.exists(self.ncfile):
            os.unlink(self.ncfile)

    def test_ncCreateFile(self):
        """Test nctools creates a file successfully"""
        ncobj = netCDF4.Dataset(self.ncfile, 'w', format='NETCDF4',
                                clobber=True)
        # output data.
        press_out = 900. + np.arange(self.nlevs * self.nlats * \
                                     self.nlons, dtype='f')
        press_out.shape = (self.nlevs, self.nlats, self.nlons)
        temp_out = 9. + np.arange(self.nlevs * self.nlats * \
                                  self.nlons, dtype='f')
        temp_out.shape = (self.nlevs, self.nlats, self.nlons)
        # create the lat and lon dimensions.
        nctools.ncCreateDim(ncobj, 'lat', self.lats_out, 'f', self.lat_atts)
        nctools.ncCreateDim(ncobj, 'lon', self.lons_out, 'f', self.lon_atts)
        # create level dimension.
        nctools.ncCreateDim(ncobj, 'level', np.arange(self.nlevs), 'f')
        # create time dimension (record, or unlimited dimension)
        nctools.ncCreateDim(ncobj, 'time', np.arange(self.nrecs), 'f',
                            self.time_atts)
        dimensions = ('time', 'level', 'lat', 'lon')

        press = nctools.ncCreateVar(ncobj, 'pressure',
                                    dimensions, 'float64')
        temp = nctools.ncCreateVar(ncobj, 'temperature',
                                   dimensions, 'float64')

        press.units =  'hPa'
        temp.units = 'celsius'

        for nrec in range(self.nrecs):
            press[nrec,:,::] = press_out
            temp[nrec,:,::] = temp_out
        ncobj.close()

    def test_ncSaveGridExceptions(self):
        """Test ncSaveGrid raises correct exceptions"""
        self.assertRaises((KeyError,),
                          nctools.ncSaveGrid,
                          self.ncfile, self.missing_key_dims,
                          self.variables)
        self.assertRaises((KeyError,),
                          nctools.ncSaveGrid,
                          self.ncfile, self.dimensions,
                          self.missing_key_vars)
        self.assertRaises((ValueError,),
                          nctools.ncSaveGrid,
                          self.ncfile, self.dimensions,
                          self.misshapen_vars)


    def test_ncSaveGridOpenFile(self):
        """Test ncSaveGrid returns netCDF4.Dataset if keepfileopen=True"""
        ncobj = nctools.ncSaveGrid(self.ncfile,
                                   self.dimensions,
                                   self.variables,
                                   keepfileopen=True)

        self.assertEqual(type(ncobj), Dataset)
        ncobj.close()

    def test_ncSaveGridNullValue(self):
        """Test ncSaveGrid can save a variable with no values"""
        nctools.ncSaveGrid(self.ncfile,
                           self.dimensions,
                           self.nullvalue_var)


class TestNCReading(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.ncfile = pjoin(unittest_dir, 'test_data', 'pres_temp_4D.nc')

        self.nrecs = 2; self.nlevs = 2; self.nlats = 6; self.nlons = 12
        self.lats_out = -25.0 + 5.0*np.arange(self.nlats, dtype='float')
        self.lons_out = 125.0 + 5.0*np.arange(self.nlons, dtype='float')
        self.latrange = [min(self.lats_out), max(self.lats_out)]
        self.lonrange = [min(self.lons_out), max(self.lons_out)]
        self.lat_atts = {'units': 'degrees_north',
                         'long_name': 'latitude',
                         'actual_range': self.latrange}
        self.lon_atts = {'units': 'degrees_east',
                         'long_name': 'longitude',
                         'actual_range': self.lonrange}
        self.time_atts = {'units': 'hours since 2000-01-01 00:00:00',
                          'calendar': 'standard'}

        press_out = 900. + np.arange(self.nrecs * self.nlevs * \
                                     self.nlats * self.nlons,
                                     dtype='f')
        press_out.shape = (self.nrecs, self.nlevs,
                           self.nlats, self.nlons)
        temp_out = 9. + np.arange(self.nrecs * self.nlevs * \
                                  self.nlats * self.nlons,
                                  dtype='f')
        temp_out.shape = (self.nrecs, self.nlevs,
                          self.nlats, self.nlons)


        self.dimensions = {
            0: {
                'name': 'time',
                'values': np.arange(self.nrecs),
                'dtype': int,
                'atts': self.time_atts
                },
            1: {
                'name': 'level',
                'values': np.arange(self.nlevs),
                'dtype': 'f',
                'atts': {}
                },
            2: {
                'name': 'lat',
                'values': self.lats_out,
                'dtype': 'f',
                'atts': self.lat_atts
                },
            3: {
                'name': 'lon',
                'values': self.lons_out,
                'dtype': 'f',
                'atts': self.lon_atts
                }
            }

        self.variables = {
            0: {
                'name': 'pressure',
                'values': press_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'hPa',
                    'standard_name': 'air_pressure'
                    }
                },
            1: {
                'name': 'temperature',
                'values': temp_out,
                'dtype': 'float64',
                'dims': ('time', 'level', 'lat', 'lon'),
                'atts': {
                    'units': 'degrees_celcius',
                    'standard_name': 'air_temperature'
                    }
                }
            }

        ncobj = netCDF4.Dataset(self.ncfile, 'w', format='NETCDF4',
                                clobber=True)
        # output data.
        press_out = 900. + np.arange(self.nlevs * self.nlats * \
                                     self.nlons, dtype='f')
        press_out.shape = (self.nlevs, self.nlats, self.nlons)
        temp_out = 9. + np.arange(self.nlevs * self.nlats * \
                                  self.nlons, dtype='f')
        temp_out.shape = (self.nlevs, self.nlats, self.nlons)
        # create the lat and lon dimensions.
        nctools.ncCreateDim(ncobj, 'lat', self.lats_out, 'f',
                            self.lat_atts)
        nctools.ncCreateDim(ncobj, 'lon', self.lons_out, 'f',
                            self.lon_atts)
        # create level dimension.
        nctools.ncCreateDim(ncobj, 'level', np.arange(self.nlevs), 'f')
        # create time dimension (record, or unlimited dimension)
        nctools.ncCreateDim(ncobj, 'time', np.arange(self.nrecs), 'f',
                            self.time_atts)
        dimensions = ('time', 'level', 'lat', 'lon')

        press = nctools.ncCreateVar(ncobj, 'pressure', dimensions, 'float64')
        temp = nctools.ncCreateVar(ncobj, 'temperature', dimensions, 'float64')

        press.units =  'hPa'
        temp.units = 'celsius'

        for nrec in range(self.nrecs):
            press[nrec,:,::] = press_out
            temp[nrec,:,::] = temp_out
        ncobj.close()

    def test_ncReadFile(self):
        """Test nctools functions for reading dimensions and variables"""
        ncobj = nctools.ncLoadFile(self.ncfile)
        lats_check = -25.0 + 5.0*np.arange(self.nlats, dtype='float')
        lons_check = 125.0 + 5.0*np.arange(self.nlons, dtype='float')
        press_check = 900. + np.arange(self.nlevs * self.nlats * \
                                       self.nlons, dtype='float64')
        press_check.shape = (self.nlevs, self.nlats, self.nlons)
        temp_check = 9. + np.arange(self.nlevs * self.nlats * \
                                    self.nlons, dtype='float64')
        temp_check.shape = (self.nlevs, self.nlats, self.nlons)
        lats = nctools.ncGetDims(ncobj, 'lat')
        lons = nctools.ncGetDims(ncobj, 'lon')
        self.numpyAssertAlmostEqual(lats_check, lats)
        self.numpyAssertAlmostEqual(lons_check, lons)
        press = nctools.ncGetData(ncobj, 'pressure')
        temp = nctools.ncGetData(ncobj, 'temperature')
        for nrec in range(self.nrecs):
            self.numpyAssertEqual(press_check, press[nrec])
            self.numpyAssertEqual(temp_check, temp[nrec])
        ncobj.close()

    def test_ncGetTimes(self):
        """Test ncGetTimes returns datetime objects"""
        ncobj = netCDF4.Dataset(self.ncfile)
        times = nctools.ncGetTimes(ncobj)
        ncobj.close()
        #self.assertEqual(type(times[0]), datetime)
        #  Note: cftype.real_datetime inherits from datetime.datetime
        print(type(times[0]))
        self.assertTrue(issubclass(type(times[0]), datetime))

    def test_ncGetTimeValues(self):
        """Test ncGetTimes returns correct time values"""
        ncobj = netCDF4.Dataset(self.ncfile)
        times = nctools.ncGetTimes(ncobj)
        start = datetime.strptime(ncobj.variables['time'].units,
                                  'hours since %Y-%m-%d %H:%M:%S')
        t = np.array([start + timedelta(hours=t) for t in range(self.nrecs)])
        self.numpyAssertEqual(t, times)
        ncobj.close()

    def tearDown(self):
        os.unlink(self.ncfile)

if __name__ == "__main__":
    unittest.main()
