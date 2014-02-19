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

 Title: grid.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 07/01/08 3:28:PM
 Description: Provide basic functions to read, write and sample ascii
 grid files. These files are assumed to have an ArcGIS GRID format:

 ncols        nx
 nrows        ny
 xllcorner    xll
 yllcorner    yll
 cellsize     dx
 NODATA_value -9999
 data data data ...
  :    :    :
  :    :    :
 Constraints:
 SeeAlso: files.py, config.py
 Version: 74

Version: $Rev: 685 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-04-30 1:09:PM
Modification: Standardised function names

$Id: grid.py 685 2012-03-29 04:22:32Z carthur $
"""
import os, sys
import logging as log

import numpy
import threading
from lat_long_UTM_conversion import LLtoUTM, UTMtoLL
import metutils
import nctools

__version__ = '$Id: grid.py 685 2012-03-29 04:22:32Z carthur $'


def grdSave(filename, data, lon, lat, delta, delimiter=' ', nodata=-9999,
            fmt='%.10e', coords='latlon'):
    """
    Save formatted data to an ascii grid format file.
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    grdSave(filename, data, lon, lat, delta, delimiter=' ',
            nodata=-9999, fmt='%.10e, coords='latlon')
    """

    if filename and os.path.isfile(filename):
        if filename.endswith('.gz'):
            import gzip
            fh = gzip.open(filename, 'wb')
        else:
            fh = file(filename, 'w')
    elif hasattr(filename, 'seek'):
        fh = filename
    else:
        try:
            fh = open(filename,'w')
        except:
            raise ValueError('Filename must be a string or file handle')

    if coords == 'UTM':
        zone, xllcorner, yllcorner = LLtoUTM(lat.min(),lon.min())
        delta = metutils.convert(delta, "deg", "m")
    else:
        # Assume geographic coordinates
        xllcorner = lon.min()
        yllcorner = lat.min()

    fh.write('ncols         '+str(len(lon))+'\n')
    fh.write('nrows         '+str(len(lat))+'\n')
    fh.write('xllcorner     '+str(xllcorner)+'\n')
    fh.write('yllcorner     '+str(yllcorner)+'\n')
    fh.write('cellsize      '+str(delta)+'\n')
    fh.write('NODATA_value  '+str(nodata)+'\n')
    X = numpy.array(data)
    origShape = None
    if len(X.shape) == 1:
        origShape = X.shape
        X.shape = len(X), 1
    for row in X:
        fh.write(delimiter.join([fmt%val for val in row]) + '\n')
    fh.close()
    if origShape is not None:
        X.shape = origShape

def grdReadFromNetcdf(filename):
    """
    Read formatted data from an ascii grid format file.
    Returns the longitude and latitude of the grid and the data values
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    longitude, latitude, data = grdRead(filename, [delimiter])
    """
    from netCDF4 import Dataset

    lon = None
    lat = None
    data = None

    if filename.endswith('nc'):
        ncdf = Dataset(filename, 'r')
        lon = ncdf.variables['lon'][:]
        lat = ncdf.variables['lat'][:]
        dataname = (set(ncdf.variables.keys()) -
                set(ncdf.dimensions.keys())).pop()
        data = ncdf.variables[dataname][:]
        ncdf.close()

    log.debug('Loaded filename %s, memory usage (bytes): lon %i lat %i data %i' % (filename, lon.nbytes, lat.nbytes, data.nbytes))

    return lon, lat, data


def grdRead(filename, delimiter=None):
    """
    Read formatted data from an ascii grid format file.
    Returns the longitude and latitude of the grid and the data values
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    longitude, latitude, data = grdRead(filename, [delimiter])
    """

    fileext = filename.rsplit('.')[-1]

    # If file extention is '.nc' then load as netcdf file
    # Otherwise load with grdRead
    if fileext == 'nc':
        nc_obj = nctools.ncLoadFile(filename)
        lon = numpy.array(nctools.ncGetDims(nc_obj, 'lon'),dtype=float)
        lat = numpy.array(nctools.ncGetDims(nc_obj, 'lat'),dtype=float)
        #lat = numpy.flipud(lat)
        data_varname = set.difference(set(nc_obj.variables.keys()),
                                      set(nc_obj.dimensions.keys()))
        if len(data_varname) != 1:
            raise IOError, 'Cannot resolve data variable in netcdf file: ' + filename
        data = numpy.array(nctools.ncGetData(nc_obj, data_varname.pop()),dtype=float)
        nc_obj.close()
    else:
        try:
            fh = open(filename, 'r')
        except:
            #g_logger.flLog("Cannot open %s"%filename)
            raise IOError, "Cannot open %s"%filename
            return

        metadata = {}
        metadata["ncols"] = []
        metadata["nrows"] = []
        metadata["xllcorner"] = []
        metadata["yllcorner"] = []
        metadata["cellsize"] = []
        metadata["NODATA_value"] = []

        for i in xrange(0,6):
            line = fh.readline()
            contents = line.split()
            label = contents[0]
            metadata[label] = float(contents[1])

        lon0 = metadata["xllcorner"]
        lon = numpy.array(range(int(metadata["ncols"])), dtype=float)
        lon = lon*metadata["cellsize"]+lon0
        lat0 = metadata["yllcorner"]
        lat = numpy.array(range(int(metadata["nrows"])), dtype=float)
        lat = lat*metadata["cellsize"]+lat0
        lat = numpy.flipud(lat)

        data = numpy.zeros([metadata["nrows"], metadata["ncols"]], dtype=float)

        for i in xrange(int(metadata["nrows"])):
            row = numpy.zeros([metadata["ncols"]], dtype=float)
            line = fh.readline()
            for j, val in enumerate(line.split(delimiter)):
                value = float(val)
                if value == metadata["NODATA_value"]:
                    value = Nan
                row[j] = value
            data[i,:] = row
        fh.close()

    log.debug('filename %s mem:: lon %i lat %i data %i' % (filename, lon.nbytes, lat.nbytes, data.nbytes))

    return lon, lat, data

class SampleGrid:
    """
    Description: Sample data from an ascii grid file
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    grid = SampleGrid(filename)
    value = grid.sampleGrid(lon, lat)

    Parameters: file (filename) containing gridded data in ascii grd
                format.

    Members:
    sampleGrid(cLon, cLat): sample a value from the grid at the given cLon, cLat
        At this time, does not interplolate from teh input grid to the given location

    Methods:
    sampleGrid(cLon, cLat): sample a value from the grid at the given cLon, cLat
        At this time, does not interplolate from the input grid to the given location.

    Internal Methods:

    """

    def __init__(self, filename):
        """
        Read in the data and ensure it's the right way around.
        """
        if filename.endswith('nc'):
            self.lon, self.lat, self.grid = grdReadFromNetcdf(filename)
        else:
            self.lon, self.lat, self.grid = grdRead(filename)
        self.grid = numpy.flipud(self.grid)

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Sample data from an ascii grid file \
    The files have 6 header lines describing the data, followed by the \
    data in a gridded format.\
    \
    Headers:\
    ncols\
    nrows\
    xllcorner\
    yllcorner\
    cellsize\
    NODATA_value\
    \
    Usage:\
    grid = SampleGrid(filename)\
    value = grid.sampleGrid(lon, lat)'

    def sampleGrid(self, lon, lat):
        """sampleGrid(self, lon, lat):
        Sample a value from the grid at the given cLon, cLat
        At this time, does not interplolate from the input grid to the given location
        
        Input: lon - longitude of the point to sample
               lat - latitude of the point to sample
        Output: value of the nearest grid point to the given lon/lat point
        
        Example:
        
        grid = SampleGrid( '/foo/bar/grid.nc' )
        value = grid.sampleGrid( 100., -25. )

        """
        indi = self.lon.searchsorted(lon)
        indj = self.lat.searchsorted(lat)

        return self.grid[indj, indi]
