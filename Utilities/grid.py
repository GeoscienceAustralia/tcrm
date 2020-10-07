"""
:mod:`grid` -- read, write and sample ascii grid files
======================================================

Provide functions to read, write and sample ascii grid files. The
ascii files are assumend to have and ArcGIS GIRD format::

    ncols        nx
    nrows        ny
    xllcorner    xll
    yllcorner    yll
    cellsize     dx
    NODATA_value -9999
    data data data ...
     :    :    :
     :    :    :


.. module:: grid
    :synopsis: Read, write and sample ascii grid files.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import logging as log

import numpy

from .lat_long_UTM_conversion import LLtoUTM, UTMtoLL
from . import metutils
from . import nctools

# pylint: disable=R0913, R0914
def grdSave(filename, data, lon, lat, delta, delimiter=' ', nodata=-9999,
            fmt='%.10e', coords='latlon'):
    """
    Save formatted data to an ascii grid format file.
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    :param str filename: Path to the file to be written.
    :param data: 2-d array of data values to store.
    :param lon: Array of longitudes corresponding to data points.
    :param lat: Array of latitudes corresponding to data points.
    :param float delta: Spacing between grid points.
    :param str delimiter: Delimiter to put between data points (default ' ').
    :param float nodata: Value to indicate missing values (default -9999).
    :param str fmt: String format statement.
    :param str coords: Optionally store the data in UTM
                       coordinates. Default is to store the data in
                       geographic coordinates (``coords='latlon'``).
                       If ``coords='UTM'``, then the latitude &
                       longitudes are converted to the local UTM
                       coordinate system.

    :raises ValueError: If the ``filename`` is not a string of file handle.

    Usage::

     >>> grdSave(filename, data, lon, lat, delta, delimiter=' ',
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
            fh = open(filename, 'w')
        except:
            raise ValueError('Filename must be a string or file handle')

    if coords == 'UTM':
        zone, xllcorner, yllcorner = LLtoUTM(lat.min(), lon.min())
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
    Read formatted data from a netcdf file.
    Returns the longitude and latitude of the grid and the data
    values. Assumes that there is only one (non-coordinate) variable
    in the file, which is only 2 dimensional.

    :param str filename: Path to a netcdf file to read.

    :returns: longitude, latitude and grid data.

    Usage:
    longitude, latitude, data = grdReadFromNetcdf(filename)
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

    log.debug('Loaded filename %s, memory use (bytes): lon %i lat %i data %i' %
              (filename, lon.nbytes, lat.nbytes, data.nbytes))

    return lon, lat, data


def grdRead(filename, delimiter=None):
    """
    Read formatted data from an ascii grid format file.
    Returns the longitude and latitude of the grid and the data values

    :param str filename: Path to an ascii grid format or netcdf file.
    :param delimiter: Delimiter for the ascii format file (optional).

    :returns: longitude, latitude, data
    :rtype: :class:`numpy.ndarray`

    Usage:
    longitude, latitude, data = grdRead(filename, [delimiter])
    """

    fileext = filename.rsplit('.')[-1]

    # If file extention is '.nc' then load as netcdf file
    # Otherwise load with grdRead
    if fileext == 'nc':
        nc_obj = nctools.ncLoadFile(filename)
        lon = numpy.array(nctools.ncGetDims(nc_obj, 'lon'), dtype=float)
        lat = numpy.array(nctools.ncGetDims(nc_obj, 'lat'), dtype=float)
        data_varname = set.difference(set(nc_obj.variables.keys()),
                                      set(nc_obj.dimensions.keys()))
        if len(data_varname) != 1:
            raise IOError('Cannot resolve data variable in netcdf file: '
                          + filename)
        data = numpy.array(nctools.ncGetData(nc_obj, data_varname.pop()),
                           dtype=float)
        nc_obj.close()
    else:
        try:
            fh = open(filename, 'r')
        except:
            raise IOError("Cannot open %s"%filename)
            return

        metadata = {}
        metadata["ncols"] = []
        metadata["nrows"] = []
        metadata["xllcorner"] = []
        metadata["yllcorner"] = []
        metadata["cellsize"] = []
        metadata["NODATA_value"] = []

        for i in range(0, 6):
            line = fh.readline()
            contents = line.split()
            label = contents[0]
            metadata[label] = float(contents[1])

        lon0 = metadata["xllcorner"]
        lon = numpy.array(list(range(int(metadata["ncols"]))), dtype=float)
        lon = lon * metadata["cellsize"] + lon0
        lat0 = metadata["yllcorner"]
        lat = numpy.array(list(range(int(metadata["nrows"]))), dtype=float)
        lat = lat * metadata["cellsize"] + lat0
        lat = numpy.flipud(lat)

        data = numpy.zeros([int(metadata["nrows"]),
                            int(metadata["ncols"])], 
                           dtype=float)

        for i in range(int(metadata["nrows"])):
            row = numpy.zeros([int(metadata["ncols"])], dtype=float)
            line = fh.readline()
            for j, val in enumerate(line.split(delimiter)):
                value = float(val)
                if value == metadata["NODATA_value"]:
                    value = numpy.nan
                row[j] = value
            data[i, :] = row
        fh.close()

    log.debug('filename %s mem:: lon %i lat %i data %i' %
              (filename, lon.nbytes, lat.nbytes, data.nbytes))

    return lon, lat, data

class SampleGrid(object):
    """
    Sample data from a gridded data file. The class is instantiated
    with a gridded data file (either an ascii file or a netcdf file
    with single 2-d variable), and the :meth:`SampleGrid.sampleGrid`
    method returns the value of the grid point closest to the given
    longitude and latitude.

    :param str filename: Path to a file containing gridded data.

    Example::

          >>> grid = SampleGrid( '/foo/bar/grid.nc' )
          >>> value = grid.sampleGrid( 100., -25. )

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

    def sampleGrid(self, lon, lat):
        """
        Sample a value from the grid at the given cLon, cLat
        At this time, does not interplolate from the input grid to the
        given location.

        :param float lon: Longitude of the point to sample.
        :param float lat: Latitude of the point to sample.

        :returns: Value of the nearest grid point to the given lon/lat point.

        """
        indi = self.lon.searchsorted(numpy.mod(lon, 360.))
        indj = self.lat.searchsorted(lat)

        try:
            value = self.grid[indj, indi]
        except IndexError:
            log.exception(f"Index is out of bounds for point ({lon}, {lat})")
            log.exception(f"Bounds of grid are ({self.lon.min()} - {self.lon.max()},"
                          f"{self.lat.min()} - {self.lat.max()}")
            raise
        return value
