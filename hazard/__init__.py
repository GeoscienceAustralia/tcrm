"""
:mod:`hazard` -- Hazard calculation
===================================

.. module:: hazard
.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

This module contains the core objects for the return period hazard
calculation.

Hazard calculations can be run in parallel using MPI if the
:term:`mpi4py` library is found and TCRM is run using the
:term:`mpirun` command. For example, to run with 10 processors::

    mpirun -n 10 python tcrm.py cairns.ini

:class:`hazard` can be correctly initialised and started by
calling the :meth: `run` with the location of a *configFile*::

    import hazard
    hazard.run('cairns.ini')

"""

import os
import sys
import numpy as np
import logging
import random

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile

from Utilities.files import flProgramVersion
from Utilities.config import ConfigParser
from Utilities.parallel import attemptParallel, disableOnWorkers
import Utilities.nctools as nctools
from .evd import EVFUNCS
from . import GPD

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

def setDomain(inputPath):
    """
    Establish the full extent of input wind field files

    :param str inputPath: path of folder containing wind field files

    :return:  Longitudes and latitudes of the wind field grid.
    :rtype: `numpy.ndarray`

    """

    fileList = os.listdir(inputPath)
    inputFile = pjoin(inputPath, fileList[0])
    ncobj = nctools.ncLoadFile(inputFile)
    wf_lon = nctools.ncGetDims(ncobj, 'lon')
    wf_lat = nctools.ncGetDims(ncobj, 'lat')
    ncobj.close()
    return wf_lon, wf_lat

class Tile(object):
    """
    Tile object.

    Because there is a buffer region around the outer edge of
    the :class:`dict` gridLimit, the indices where data is pulled from
    (the wind field files) are different from those where the data is
    stored (the output hazard array).

    This object holds the index ranges for the input array and output array,
    to indicate this relationship.

    """

    def __init__(self, number, input_limits, output_limits):
        self.input_limits = input_limits
        self.output_limits = output_limits
        self.number = number

class TileGrid(object):
    """
    Tiling to minimise MemoryErrors and enable parallelisation.

    """

    def __init__(self, gridLimit, wf_lon, wf_lat, xstep=100, ystep=100):
        """
        Initialise the tile grid for dividing up the domain

        :param gridLimit: :class:`dict` describing the domain where the
                          tracks will be generated.
                          The :class:`dict` should contain the keys :attr:`xMin`,
                          :attr:`xMax`, :attr:`yMin` and :attr:`yMax`. The *y*
                          variable bounds the latitude and the *x* variable
                          bounds the longitude.

        :param wf_lon: `numpy.ndarray` of longitudes in thw wind field.
        :param wf_lat: `numpy.ndarray` of latitudes of the wind field.
        :param xstep: `int` size of the tile in the x-direction.
        :param ystep: `int` size of the tile in the y-direction.

        """


        ii, = ((np.round(wf_lon*1000).astype(int) >=
                int(gridLimit['xMin']*1000)) &
               (np.round(wf_lon*1000).astype(int) <=
                int(gridLimit['xMax']*1000))).nonzero()


        jj, = ((np.round(wf_lat*1000).astype(int) >=
                int(gridLimit['yMin']*1000)) &
               (np.round(wf_lat*1000).astype(int) <=
                int(gridLimit['yMax']*1000))).nonzero()

        self.imin = ii[0]
        self.imax = ii[-1]
        self.jmin = jj[0]
        self.jmax = jj[-1]

        self.xdim = self.imax - self.imin + 1
        self.ydim = self.jmax - self.jmin + 1

        self.xstep = xstep
        self.ystep = ystep
        self.wf_lon = wf_lon
        self.wf_lat = wf_lat

        self.tileGrid()

    def tileGrid(self):
        """
        Defines the indices required to subset a 2D array into smaller
        rectangular 2D arrays (of dimension x_step * y_step).
        """

        subset_maxcols = int(np.ceil(self.xdim / float(self.xstep)))
        subset_maxrows = int(np.ceil(self.ydim / float(self.ystep)))
        self.num_tiles = subset_maxcols * subset_maxrows
        self.x_start = np.zeros(self.num_tiles, 'i')
        self.x_end = np.zeros(self.num_tiles, 'i')
        self.y_start = np.zeros(self.num_tiles, 'i')
        self.y_end = np.zeros(self.num_tiles, 'i')
        k = 0

        for i in range(subset_maxcols):
            for j in range(subset_maxrows):
                self.x_start[k] = i * self.xstep + self.imin
                self.x_end[k] = min((i + 1) * self.xstep + self.imin,
                                    self.xdim + self.imin) - 1
                self.y_start[k] = j * self.ystep + self.jmin
                self.y_end[k] = min((j + 1) * self.ystep + self.jmin,
                                    self.ydim + self.jmin) - 1
                k += 1


    def getGridLimit(self, k):
        """
        Return the limits for tile `k`. x-indices correspond to the
        east-west coordinate, y-indices correspond to the north-south
        coordinate.

        :param int k: tile number

        :return x1: minimum x-index for tile `k`
        :return x2: maximum x-index for tile `k`
        :return y1: minimum y-index for tile `k`
        :return y2: maximum y-index for tile `k`

        """

        x1 = int(self.x_start[k])
        x2 = int(self.x_end[k] + 1)
        y1 = int(self.y_start[k])
        y2 = int(self.y_end[k] + 1)

        return x1, x2, y1, y2

    def getDomainExtent(self):
        """
        Return the longitude and latitude values that lie within
        the modelled domain

        :return lon: :class:`numpy.ndarray` containing longitude values
        :return lat: :class:`numpy.ndarray` containing latitude values

        """

        lon = self.wf_lon[self.imin:self.imax + 1]
        lat = self.wf_lat[self.jmin:self.jmax + 1]

        return lon, lat



class HazardCalculator(object):
    """
    Calculate return period wind speeds using GEV fitting

    """

    def __init__(self, configFile, tilegrid, numSim, minRecords, yrsPerSim,
                 calcCI=False, evd='GEV'):
        """
        Initialise HazardCalculator object.

        :param str configFile: path to TCRM configuration file.
        :param tilegrid: :class:`TileGrid` instance
        :param int numSim: number of simulations created.
        :param int minRecords: minimum number of valid wind speed values required
                               to do fitting.
        :param int yrsPerSim:
        :param boolean calcCI:
        :param str extreme_value_distribution: evd to use. Options so far are GEV
                                               and GPD.
        """
        config = ConfigParser()
        config.read(configFile)

        self.nodata = -9999.
        self.years = np.array(config.get('Hazard',
                                         'Years').split(',')).astype('f')
        self.outputPath = pjoin(config.get('Output', 'Path'), 'hazard')
        self.inputPath = pjoin(config.get('Output', 'Path'), 'windfield')
        gridLimit = config.geteval('Region', 'gridLimit')

        self.numSim = numSim
        self.minRecords = minRecords
        self.yrsPerSim = yrsPerSim
        self.calcCI = calcCI
        if self.calcCI:
            log.debug("Bootstrap confidence intervals will be calculated")
            self.sample_size = config.getint('Hazard', 'SampleSize')
            self.prange = config.getint('Hazard', 'PercentileRange')
        self.evd = evd

        self.tilegrid = tilegrid
        lon, lat = self.tilegrid.getDomainExtent()

        # Create arrays for storing output data:
        self.loc = np.zeros((len(lat), len(lon)), dtype='f')
        self.shp = np.zeros((len(lat), len(lon)), dtype='f')
        self.scale = np.zeros((len(lat), len(lon)), dtype='f')
        self.Rp = np.zeros((len(self.years), len(lat), len(lon)), dtype='f')

        self.RPupper = np.zeros((len(self.years), len(lat), len(lon)), dtype='f')
        self.RPlower = np.zeros((len(self.years), len(lat), len(lon)), dtype='f')

        self.global_atts = {'title': ('TCRM hazard simulation - '
                            'return period wind speeds'),
                            'tcrm_version': flProgramVersion(),
                            'python_version': sys.version}


        # Add configuration settings to global attributes:
        for section in config.sections():
            for option in config.options(section):
                key = "{0}_{1}".format(section, option)
                value = config.get(section, option)
                self.global_atts[key] = value

    def calculateHazard(self, tilelimits):
        """
        Load input hazard data and then calculate the return period and
        distribution parameters for a given tile. The extreme value distribution
        used in the calculation can be set in the config file. The default
        distribution is set to GEV.

        :param tilelimits: `tuple` of tile limits       

        Returns:
        --------

        :param Rp: `numpy.ndarray` of return period wind speed values for each lat/lon
        :param loc: `numpy.ndarray` of location parameters for each lat/lon
        :param scale: `numpy.ndarray` of scale parameters for each lat/lon
        :param shp: `numpy.ndarray` of shape parameters for each lat/lon
        :param RpUpper: Upper CI return period wind speed values for each lat/lon
        :param RpLower: Lower CI return period wind speed values for each lat/lon

        """

        #try:
        #    evfunc = EVFUNCS['{0}fit'.format(self.evd.lower())]
        #except KeyError:
        #    log.exception("{0} distribution not implemented for hazard calculation".format(self.evd))
        #    raise

        log.info("Using {0} distribution for the hazard curves".format(self.evd))
        #Vr = loadFilesFromPath(self.inputPath, tilelimits)
        #Rp, loc, scale, shp = evfunc(Vr, self.years, self.numSim, self.nodata,
        #self.minRecords)
        if self.evd not in ["GPD", "GEV", "power", "emp"]:
            msg = (f"Invalid extreme value distribution function: {self.evd} \n"
                   "Set 'Hazard--ExtremeValueDistribution' to one of the following: \n"
                   "GPD, GEV, power or emp")
            raise ValueError (msg)

        if self.evd == 'GPD':
            log.info("Using the GPD distribution for the hazard curves")
            Vr = loadFilesFromPath(self.inputPath, tilelimits) 
            Rp, loc, scale, shp = calculateGPD(Vr, self.years, self.numSim, self.nodata,
                                               self.minRecords, self.yrsPerSim)
        elif self.evd == 'GEV':
            log.info("Using the GEV distribution for the hazard curves")
            Vr = aggregateWindFields(self.inputPath, self.numSim, tilelimits)
            Rp, loc, scale, shp = calculateGEV(Vr, self.years, self.nodata,
                                               self.minRecords, self.yrsPerSim)
        elif self.evd == 'power':
            log.info("Using the power law function for the hazard curves")
            Vr = loadFilesFromPath(self.inputPath, tilelimits)
            Rp, loc, scale, shp = calculatePower(Vr, self.years, self.numSim, self.nodata,
                                                 self.minRecords, self.yrsPerSim)
        elif self.evd == 'emp':
            log.info("Using empirical hazard curve")
            Vr = loadFilesFromPath(self.inputPath, tilelimits)
            Rp, loc, scale, shp = calculateEMP(Vr, self.years, self.numSim, self.nodata,
                                               self.minRecords, self.yrsPerSim)
        if self.calcCI: # set in config
            RpUpper, RpLower = calculateCI(Vr, self.years, self.nodata,
                                           self.minRecords, self.yrsPerSim,
                                           self.sample_size, self.prange)
 
            return (tilelimits, Rp, loc, scale, shp, RpUpper, RpLower)
        else:
            return (tilelimits, Rp, loc, scale, shp)

    def dumpHazardFromTiles(self, tiles, progressCallback=None):
        """
        Iterate over tiles to calculate return period hazard levels

        :param tileiter: generator that yields tuples of tile dimensions.

        """
        status = MPI.Status()
        work_tag = 0
        result_tag = 1
        if (comm.rank == 0) and (comm.size > 1):
            w = 0
            p = comm.size - 1
            for d in range(1, comm.size):
                if w < len(tiles):
                    comm.send(tiles[w], dest=d, tag=work_tag)
                    log.debug("Processing tile %d of %d" % (w, len(tiles)))
                    w += 1
                else:
                    comm.send(None, dest=d, tag=work_tag)
                    p = w

            terminated = 0
            while(terminated < p):

                result = comm.recv(source=MPI.ANY_SOURCE, status=status, tag=MPI.ANY_TAG)

                if self.calcCI:
                    limits, Rp, loc, scale, shp, RPupper, RPlower = result
                else:
                    limits, Rp, loc, scale, shp = result

                # Reset the min/max bounds for the output array:
                (xmin, xmax, ymin, ymax) = limits
                xmin -= self.tilegrid.imin
                xmax -= self.tilegrid.imin
                ymin -= self.tilegrid.jmin
                ymax -= self.tilegrid.jmin

                self.loc[ymin:ymax, xmin:xmax] = loc
                self.scale[ymin:ymax, xmin:xmax] = scale
                self.shp[ymin:ymax, xmin:xmax] = shp
                self.Rp[:, ymin:ymax, xmin:xmax] = Rp[:, :, :]

                if self.calcCI:
                    self.RPupper[:, ymin:ymax, xmin:xmax] = RPupper[:, :, :]
                    self.RPlower[:, ymin:ymax, xmin:xmax] = RPlower[:, :, :]

                d = status.source

                if w < len(tiles):
                    comm.send(tiles[w], dest=d, tag=status.tag)
                    log.debug("Processing tile %d of %d" % (w, len(tiles)))
                    w += 1
                else:
                    comm.send(None, dest=d, tag=status.tag)
                    terminated += 1

                log.debug("Number of terminated threads is %d"%terminated)
                if progressCallback:
                    progressCallback(w)

        elif (comm.size > 1) and (comm.rank != 0):
            status = MPI.Status()
            W = None
            while(True):
                W = comm.recv(source=0, tag=work_tag, status=status)
                if W is None:
                    log.debug("No work to be done on this processor: {0}".format(comm.rank))
                    break
                results = self.calculateHazard(W)
                comm.send(results, dest=0, tag=status.tag)

        elif comm.size == 1 and comm.rank == 0:
            # Assumed no Pypar - helps avoid the need to extend DummyPypar()
            for i, tile in enumerate(tiles):
                log.debug("Processing tile %d of %d" % (i, len(tiles)))
                result = self.calculateHazard(tile)
                if self.calcCI:
                    limits, Rp, loc, scale, shp, RPupper, RPlower = result
                else:
                    limits, Rp, loc, scale, shp = result

                # Reset the min/max bounds for the output array:
                (xmin, xmax, ymin, ymax) = limits
                xmin -= self.tilegrid.imin
                xmax -= self.tilegrid.imin
                ymin -= self.tilegrid.jmin
                ymax -= self.tilegrid.jmin

                self.loc[ymin:ymax, xmin:xmax] = loc
                self.scale[ymin:ymax, xmin:xmax] = scale
                self.shp[ymin:ymax, xmin:xmax] = shp
                self.Rp[:, ymin:ymax, xmin:xmax] = Rp[:, :, :]
                if self.calcCI:
                    self.RPupper[:, ymin:ymax, xmin:xmax] = RPupper[:, :, :]
                    self.RPlower[:, ymin:ymax, xmin:xmax] = RPlower[:, :, :]

                if progressCallback:
                    progressCallback(i)


    @disableOnWorkers
    def saveHazard(self):
        """
        Save hazard data to a netCDF file.

        """

        log.info("Saving hazard data file")
        # FIXME: need to ensure CF-1.6 and OGC compliance in output files.
        lon, lat = self.tilegrid.getDomainExtent()

        dimensions = {
            0: {
                'name': 'ari',
                'values': self.years,
                'dtype': 'f',
                'atts': {
                    'long_name' : 'Average recurrence interval',
                    'units' : 'years',
                    'axis' : 'Z'
                }
            },
            1: {
                'name': 'lat',
                'values': lat,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Latitude',
                    'standard_name': 'latitude',
                    'units': 'degrees_north',
                    'axis': 'Y'
                }
            },
            2: {
                'name': 'lon',
                'values': lon,
                'dtype': 'd',
                'atts': {
                    'long_name': 'Longitude',
                    'standard_name': 'longitude',
                    'units': 'degrees_east',
                    'axis': 'X'
                }
            }
        }

        # Create variables:
        variables = {
            0: {
                'name': 'loc',
                'dims': ('lat', 'lon'),
                'values': self.loc,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Location parameter for GEV distribution',
                    'units': 'm/s',
                    'actual_range': (np.min(self.loc), np.max(self.loc)),
                    'valid_range': (0.0, 200.),
                    'grid_mapping': 'crs'
                }
            },
            1: {
                'name': 'scale',
                'dims': ('lat', 'lon'),
                'values': self.scale,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Scale parameter for GEV distribution',
                    'units': '',
                    'grid_mapping': 'crs'
                }
            },
            2: {
                'name': 'shp',
                'dims': ('lat', 'lon'),
                'values': self.shp,
                'dtype': 'f',
                'least_significant_digit': 5,
                'atts': {
                    'long_name': 'Shape parameter for GEV distribution',
                    'units': '',
                    'grid_mapping': 'crs'
                }
            },
            3: {
                'name': 'wspd',
                'dims': ('ari', 'lat', 'lon'),
                'values': self.Rp,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Average recurrence interval wind speed',
                    'units': 'm/s',
                    'actual_range': (np.min(self.Rp), np.max(self.Rp)),
                    'valid_range': (0.0, 200.),
                    'grid_mapping': 'crs'
                }
            },
            4: {
                'name': 'wspdupper',
                'dims': ('ari', 'lat', 'lon'),
                'values': self.RPupper,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Upper percentile ARI wind speed',
                    'units': 'm/s',
                    'percentile': 95,
                    'valid_range': (0.0, 200.),
                    'grid_mapping': 'crs'
                }
            },
            5: {
                'name': 'wspdlower',
                'dims': ('ari', 'lat', 'lon'),
                'values': self.RPlower,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lower percentile ARI wind speed',
                    'units': 'm/s',
                    'percentile': 5,
                    'valid_range': (0.0, 200.),
                    'grid_mapping': 'crs'
                }
            },
            6: {
                'name': 'crs',
                'dims': (),
                'values': None,
                'dtype': 'i',
                'atts': {
                    'grid_mapping_name': 'latitude_longitude',
                    'semi_major_axis': 6378137.0,
                    'inverse_flattening': 298.257222101,
                    'longitude_of_prime_meridian': 0.0
                }
            }
        }

        # Create output file for return-period gust wind speeds and
        # GEV parameters
        nctools.ncSaveGrid(pjoin(self.outputPath, 'hazard.nc'),
                           dimensions, variables,
                           nodata=self.nodata,
                           datatitle='TCRM hazard simulation',
                           gatts=self.global_atts, writedata=True,
                           keepfileopen=False)


def calculateGEV(Vr, years, nodata, minRecords, yrsPerSim):
    """
    Fit a GEV to the wind speed records for a 2-D extent of
    wind speed values

    :param Vr: `numpy.ndarray` of wind speeds (3-D - event, lat, lon)
               block maxima processed with aggregateWindRecords
    :param years: `numpy.ndarray` of years for which to evaluate
                  return period values
    :param float nodata: missing data value.
    :param int minRecords: minimum number of valid wind speed values required
                           to fit distribution.
    :param int yrsPerSim: Taken from the config file

    Returns:
    --------
    
    GEV fit parameters and return period wind speeds for each grid cell in
    simulation domain

    :param Rp: `numpy.ndarray` of return period wind speed values
    :param loc: `numpy.ndarray` of location parameters in the domain of `Vr`
    :param scale: `numpy.ndarray` of scale parameters in the domain of `Vr`
    :param shp: `numpy.ndarray` of shape parameters in the domain of `Vr`

    """

    Vr.sort(axis=0) #axis 0 = year
    Rp = np.zeros((len(years),) + Vr.shape[1:], dtype='f') # Rp = years x lat x lon
    loc = np.zeros(Vr.shape[1:], dtype='f') # loc = lat x lon
    scale = np.zeros(Vr.shape[1:], dtype='f') # scale = lat x lon
    shp = np.zeros(Vr.shape[1:], dtype='f') # shp = lat x lon

    for i in range(Vr.shape[1]): # lat
        for j in range(Vr.shape[2]): # lon
            if Vr[:,i,j].max() > 0.0: # all years at one lat/lon
                w, l, sc, sh = evd.gevfit(Vr[:,i,j], years, nodata,
                                          minRecords, yrsPerSim)
                # w = array of return period wind speed values
                # l = location parameter of fit
                # sc = scale parameter of fit
                # sh = shape parameter of fit

                # Put the returned values back into the lat/lon grid
                Rp[:, i, j] = w
                loc[i, j] = l
                scale[i, j] = sc
                shp[i, j] = sh

    return Rp, loc, scale, shp

def calculateEMP(Vr, years, numsim, nodata, minRecords, yrsPerSim):
    """
    Calculate empirical return levels the wind speed records for a 2-D extent of
    wind speed values

    :param Vr: `numpy.ndarray` of wind speeds (3-D - event, lat, lon)
               block maxima processed with aggregateWindRecords
    :param years: `numpy.ndarray` of years for which to evaluate
                  return period values
    :param float nodata: missing data value.
    :param int minRecords: minimum number of valid wind speed values required
                           to fit distribution.
    :param int yrsPerSim: Taken from the config file

    Returns:
    --------
    
    GEV fit parameters and return period wind speeds for each grid cell in
    simulation domain

    :param Rp: `numpy.ndarray` of return period wind speed values
    :param loc: `numpy.ndarray` of location parameters in the domain of `Vr`
    :param scale: `numpy.ndarray` of scale parameters in the domain of `Vr`
    :param shp: `numpy.ndarray` of shape parameters in the domain of `Vr`

    """

    Vr.sort(axis=0) #axis 0 = year
    Rp = np.zeros((len(years),) + Vr.shape[1:], dtype='f') # Rp = years x lat x lon
    loc = np.zeros(Vr.shape[1:], dtype='f') # loc = lat x lon
    scale = np.zeros(Vr.shape[1:], dtype='f') # scale = lat x lon
    shp = np.zeros(Vr.shape[1:], dtype='f') # shp = lat x lon

    for i in range(Vr.shape[1]): # lat
        for j in range(Vr.shape[2]): # lon
            if Vr[:,i,j].max() > 0.0: # all years at one lat/lon
                w, l, sc, sh = evd.empfit(Vr[:,i,j], years, numsim, nodata,
                                          minRecords)
                # w = array of return period wind speed values
                # l = location parameter of fit
                # sc = scale parameter of fit
                # sh = shape parameter of fit

                # Put the returned values back into the lat/lon grid
                Rp[:, i, j] = w
                loc[i, j] = l
                scale[i, j] = sc
                shp[i, j] = sh

    return Rp, loc, scale, shp



def calculateGPD(Vr, years, numsim, nodata, minRecords, yrsPerSim):
    """
    Fit a GPD to the wind speed records for a 2-D extent of
    wind speed values

    :param inputPath: path to individual wind field files.
    :param tuple tilelimits: tuple of index limits of a tile.
    :param years: `numpy.ndarray` of years for which to evaluate
                  return period values
    :param int numsim: number of simulations created.
    :param float nodata: missing data value.
    :param int minRecords: minimum number of valid wind speed values required
                           to fit distribution.
    :param int yrsPerSim: Taken from the config file

    Returns:
    --------
    
    GPD fit parameters and return period wind speeds for each grid cell in
    simulation domain

    :param Rp: `numpy.ndarray` of return period wind speed values
    :param loc: `numpy.ndarray` of location parameters in the domain of `Vr`
    :param scale: `numpy.ndarray` of scale parameters in the domain of `Vr`
    :param shp: `numpy.ndarray` of shape parameters in the domain of `Vr`

    """

    Vr.sort(axis=0) #axis 0 = year
    Rp = np.zeros((len(years),) + Vr.shape[1:], dtype='f') # Rp = years x lat x lon
    loc = np.zeros(Vr.shape[1:], dtype='f') # loc = lat x lon
    scale = np.zeros(Vr.shape[1:], dtype='f') # scale = lat x lon
    shp = np.zeros(Vr.shape[1:], dtype='f') # shp = lat x lon

    for i in range(Vr.shape[1]): # lat
        for j in range(Vr.shape[2]): # lon
            if Vr[:,i,j].max() > 0.0: # all years at one lat/lon
                log.debug("lat: {0}, lon: {1}".format(i, j))
                w, l, sc, sh = GPD.gpdfit(Vr[:,i,j],
                                          years,
                                          numsim,
                                          nodata,
                                          minRecords,
                                          threshold=99.5)
                # w = array of return period wind speed values
                # l = location parameter of fit
                # sc = scale parameter of fit
                # sh = shape parameter of fit

                # Put the returned values back into the lat/lon grid
                Rp[:, i, j] = w
                loc[i, j] = l
                scale[i, j] = sc
                shp[i, j] = sh

    return Rp, loc, scale, shp

def calculatePower(Vr, years, numsim, nodata, minRecords, yrsPerSim):
    """
    Fit a GPD to the wind speed records for a 2-D extent of
    wind speed values

    :param inputPath: path to individual wind field files.
    :param tuple tilelimits: tuple of index limits of a tile.
    :param years: `numpy.ndarray` of years for which to evaluate
                  return period values
    :param int numSim: number of simulations created.
    :param float nodata: missing data value.
    :param int minRecords: minimum number of valid wind speed values required
                           to fit distribution.
    :param int yrsPerSim: Taken from the config file

    Returns:
    --------
    
    GPD fit parameters and return period wind speeds for each grid cell in
    simulation domain

    :param Rp: `numpy.ndarray` of return period wind speed values
    :param loc: `numpy.ndarray` of location parameters in the domain of `Vr`
    :param scale: `numpy.ndarray` of scale parameters in the domain of `Vr`
    :param shp: `numpy.ndarray` of shape parameters in the domain of `Vr`

    """

    Vr.sort(axis=0) #axis 0 = year
    Rp = np.zeros((len(years),) + Vr.shape[1:], dtype='f') # Rp = years x lat x lon
    loc = np.zeros(Vr.shape[1:], dtype='f') # loc = lat x lon
    scale = np.zeros(Vr.shape[1:], dtype='f') # scale = lat x lon
    shp = np.zeros(Vr.shape[1:], dtype='f') # shp = lat x lon

    for i in range(Vr.shape[1]): # lat
        for j in range(Vr.shape[2]): # lon
            if Vr[:,i,j].max() > 0.0: # all years at one lat/lon
                log.debug("lat: {0}, lon: {1}".format(i, j))
                w, l, sc, sh = evd.powerfit(Vr[:,i,j],
                                            years,
                                            numsim,
                                            nodata,
                                            minRecords)
                # w = array of return period wind speed values
                # l = location parameter of fit
                # sc = scale parameter of fit
                # sh = shape parameter of fit

                # Put the returned values back into the lat/lon grid
                Rp[:, i, j] = w
                loc[i, j] = l
                scale[i, j] = sc
                shp[i, j] = sh

    return Rp, loc, scale, shp

def calculateCI(Vr, years, nodata, minRecords, yrsPerSim=1,
                sample_size=50, prange=90):
    """
    Fit a GEV to the wind speed records for a 2-D extent of
    wind speed values, providing a confidence range by resampling at
    random from the input values.

    :param Vr: `numpy.ndarray` of wind speeds (3-D - event, lat, lon)
    :param years: `numpy.ndarray` of years for which to evaluate
                  return period values.
    :param float nodata: missing data value.
    :param int minRecords: minimum number of valid wind speed values required
                           to fit distribution.
    :param int yrsPerSim: Values represent block maxima - this value indicates
                          the time span of the block (default 1).
    :param int sample_size: number of records to randomly sample for calculating
                            confidence interval of the fit.
    :param float prange: percentile range.


    Return:
    -------

    :param RpUpper: Upper CI return period wind speed values for each lat/lon
    :param RpLower: Lower CI return period wind speed values for each lat/lon

    """

    lower = (100 - prange) / 2. # 5th percentile default
    upper = 100. - lower # 95th percentile default

    nrecords = Vr.shape[0] # number of years (since we have aggregated into 1/yr)
    nsamples = nrecords / sample_size # number of iterations to perform
    
    # RpUpper/RpLower = years x lat x lon
    RpUpper = nodata*np.ones((len(years), Vr.shape[1], Vr.shape[2]), dtype='f')
    RpLower = nodata*np.ones((len(years), Vr.shape[1], Vr.shape[2]), dtype='f')

    # w: years x number of iterations
    w = np.zeros((len(years), nsamples), dtype='f')
    wUpper = np.zeros((len(years)), dtype='f')
    wLower = np.zeros((len(years)), dtype='f')

    for i in range(Vr.shape[1]): # lat
        for j in range(Vr.shape[2]): # lon
            if Vr[:, i, j].max() > 0.0: # check for valid data
                random.shuffle(Vr[:, i, j]) # shuffle the years
                for n in range(nsamples): # iterate through fitting of random samples
                    nstart = n*sample_size
                    nend  = (n + 1)*sample_size - 1
                    vsub = Vr[nstart:nend, i, j] # select random 50(default) events

                    vsub.sort()
                    if vsub.max( ) > 0.:
                        # Perform the fitting on a random subset of samples
                        w[:, n], loc, scale, shp = evd.gevfit(vsub, years, nodata,
                                                              minRecords/10, yrsPerSim)

                # Pull out the upper and lower percentiles from the random sample fits
                for n in range(len(years)):
                    wUpper[n] = percentile(w[n,:], upper)
                    wLower[n] = percentile(w[n,:], lower)

                # Store upper and lower percentiles for each return period, for each grid cell
                RpUpper[:, i, j] = wUpper
                RpLower[:, i, j] = wLower

    return RpUpper, RpLower

def aggregateWindFields(inputPath, numSimulations, tilelimits):
    """
    Aggregate wind field data into annual maxima for use in fitting
    extreme value distributions.

    :param str inputPath: path to individual wind field files.
    :param int numSimulations: Number of simulated years of activity.

    """
    from glob import glob
    log.info("Aggregating individual events to annual maxima")
    ysize = tilelimits[3] - tilelimits[2]
    xsize = tilelimits[1] - tilelimits[0]
    Vm = np.zeros((numSimulations, ysize, xsize), dtype='f')

    for year in range(numSimulations):
        filespec = pjoin(inputPath, "gust.*-%05d.nc"%year)
        fileList = glob(filespec)
        if len(fileList) == 0:
            log.debug("No files for year: {0}".format(year))
            Vm[year, :, :] = np.zeros((ysize, xsize), dtype='f')
            continue

        Va = np.zeros((len(fileList), ysize, xsize), dtype='f')
        for n, f in enumerate(fileList):
            Va[n, :, :] = loadFile(f, tilelimits)

        Vm[year, :, :] = np.max(Va, axis=0)

    return Vm

def loadFilesFromPath(inputPath, tilelimits):
    """
    Load wind field data for each subset into a 3-D array.

    :param str inputPath: str path to wind field files.

    :param tuple tilelimits: tuple of index limits of a tile.

    :returns: 3-D `numpy.narray` of wind field records.

    """

    fileList = os.listdir(inputPath)
    files = [pjoin(inputPath, f) for f in fileList]
    files = [f for f in files if os.path.isfile(f)]
    log.debug("Loading data from %d files" % (len(files)))

    ysize = tilelimits[3] - tilelimits[2]
    xsize = tilelimits[1] - tilelimits[0]
    Vr = np.empty((len(files), ysize, xsize), dtype='f')

    for n, f in enumerate(sorted(files)):
        Vr[n,:,:] = loadFile(f, tilelimits)

    return Vr

def loadFile(filename, limits):
    """
    Load a subset of the data from the given file, with the extent
    of the subset specified in the `limits` tuple

    :param str filename: str full path to file to load.

    :param tuple limits: tuple of index limits of a tile.

    :returns: 2-D `numpy.ndarray` of wind speed values.

    """

    (xmin, xmax, ymin, ymax) = limits

    try:
        ncobj = nctools.ncLoadFile(filename)
        ncobj_vmax = nctools.ncGetVar(ncobj, 'vmax')
        data_subset = ncobj_vmax[ymin:ymax, xmin:xmax]
        ncobj.close()

        if xmax < xmin or ymax < ymin:
            log.debug("max tile limits are not smaller than min")
            
        return data_subset

    except IOError:
        log.debug('{0} file does not exist'.format(filename))
        raise

def getTiles(tilegrid):
    """
    Helper to obtain a generator that yields tile numbers

    :param tilegrid: :class:`TileGrid` instance
    """

    tilenums = list(range(tilegrid.num_tiles))
    return getTileLimits(tilegrid, tilenums)

def getTileLimits(tilegrid, tilenums):
    """
    Generate a list of tuples of the x- and y- limits of a tile

    :param tilegrid: :class:`TileGrid` instance
    :param tilenums: list of tile numbers (must be sequential)

    :return: list of tuples of tile limits

    """

    tilelimits = [tilegrid.getGridLimit(t) for t in tilenums]
    return tilelimits


def run(configFile, callback=None):
    """
    Run the hazard calculations.

    This will attempt to run the calculation in parallel by tiling the
    domain, but also provides a sane fallback mechanism to execute
    in serial.

    :param str configFile: path to configuration file

    """

    log.info("Loading hazard calculation settings")

    config = ConfigParser()
    config.read(configFile)

    outputPath = config.get('Output', 'Path')
    inputPath = pjoin(outputPath, 'windfield')
    gridLimit = config.geteval('Region', 'gridLimit')
    numsimulations = config.getint('TrackGenerator', 'NumSimulations')
    yrsPerSim = config.getint('TrackGenerator', 'YearsPerSimulation')
    minRecords = config.getint('Hazard', 'MinimumRecords')
    calculate_confidence = config.getboolean('Hazard', 'CalculateCI')
    extreme_value_distribution = config.get('Hazard', 'ExtremeValueDistribution')

    wf_lon, wf_lat = setDomain(inputPath)

    global MPI, comm
    MPI = attemptParallel()
    comm = MPI.COMM_WORLD
    log.info("Running hazard calculations")
    TG = TileGrid(gridLimit, wf_lon, wf_lat)
    tiles = getTiles(TG)

    #def progress(i):
    #    callback(i, len(tiles))

    comm.barrier()
    hc = HazardCalculator(configFile, TG,
                          numsimulations,
                          minRecords,
                          yrsPerSim,
                          calculate_confidence,
                          extreme_value_distribution
                          )




    hc.dumpHazardFromTiles(tiles)
    log.debug("Finished hazard calculations")
    comm.barrier()

    hc.saveHazard()

    log.info("Completed hazard calculation")


