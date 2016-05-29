"""
:mod:`hazard` -- Hazard calculation
===================================

.. module:: hazard
.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

This module contains the core objects for the return period hazard
calculation.

Hazard calculations can be run in parallel using MPI if the
:term:`pypar` library is found and TCRM is run using the
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
import evd

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

        for i in xrange(subset_maxcols):
            for j in xrange(subset_maxrows):
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
                 calcCI=False):
        """
        Initialise HazardCalculator object.

        :param str configFile: path to TCRM configuration file.
        :param tilegrid: :class:`TileGrid` instance
        :param int numSim: number of simulations created.
        :param int minRecords: minimum number of valid wind speed values required
                               to do fitting.
        :param int yrsPerSim:
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
        distribution parameters for a given tile.

        :param tilelimits: `tuple` of tile limits
        """
        #Vr = loadFilesFromPath(self.inputPath, tilelimits)
        Vr = aggregateWindFields(self.inputPath, self.numSim, tilelimits)
        Rp, loc, scale, shp = calculate(Vr, self.years, self.nodata,
                                        self.minRecords, self.yrsPerSim)

        if self.calcCI:
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

        work_tag = 0
        result_tag = 1
        if (pp.rank() == 0) and (pp.size() > 1):
            w = 0
            p = pp.size() - 1
            for d in range(1, pp.size()):
                if w < len(tiles):
                    pp.send(tiles[w], destination=d, tag=work_tag)
                    log.debug("Processing tile %d of %d" % (w, len(tiles)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    p = w


            terminated = 0

            while(terminated < p):

                result, status = pp.receive(pp.any_source, tag=result_tag,
                                             return_status=True)

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
                    pp.send(tiles[w], destination=d, tag=work_tag)
                    log.debug("Processing tile %d of %d" % (w, len(tiles)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    terminated += 1

                log.debug("Number of terminated threads is %d"%terminated)
                if progressCallback:
                    progressCallback(w)

        elif (pp.size() > 1) and (pp.rank() != 0):
            while(True):
                W = pp.receive(source=0, tag=work_tag)
                if W is None:
                    break
                results = self.calculateHazard(W)
                pp.send(results, destination=0, tag=result_tag)

        elif pp.size() == 1 and pp.rank() == 0:
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
                'name': 'years',
                'values': self.years,
                'dtype': 'f',
                'atts': {
                    'long_name' : 'Return period',
                    'units' : 'years'
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
                'dims': ('years', 'lat', 'lon'),
                'values': self.Rp,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Return period wind speed',
                    'units': 'm/s',
                    'actual_range': (np.min(self.Rp), np.max(self.Rp)),
                    'valid_range': (0.0, 200.),
                    'grid_mapping': 'crs'
                }
            },
            4: {
                'name': 'wspdupper',
                'dims': ('years', 'lat', 'lon'),
                'values': self.RPupper,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Upper percentile return period wind speed',
                    'units': 'm/s',
                    'percentile': 95,
                    'valid_range': (0.0, 200.),
                    'grid_mapping': 'crs'
                }
            },
            5: {
                'name': 'wspdlower',
                'dims': ('years', 'lat', 'lon'),
                'values': self.RPlower,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lower percentile return period wind speed',
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


def calculate(Vr, years, nodata, minRecords, yrsPerSim):
    """
    Fit a GEV to the wind speed records for a 2-D extent of
    wind speed values

    :param Vr: `numpy.ndarray` of wind speeds (3-D - event, lat, lon)
    :param years: `numpy.ndarray` of years for which to evaluate
                  return period values

    Returns:
    --------

    :param Rp: `numpy.ndarray` of return period wind speed values
    :param loc: `numpy.ndarray` of location parameters in the domain of `Vr`
    :param scale: `numpy.ndarray` of scale parameters in the domain of `Vr`
    :param shp: `numpy.ndarray` of shape parameters in the domain of `Vr`

    """

    Vr.sort(axis=0)
    Rp = np.zeros((len(years),) + Vr.shape[1:], dtype='f')
    loc = np.zeros(Vr.shape[1:], dtype='f')
    scale = np.zeros(Vr.shape[1:], dtype='f')
    shp = np.zeros(Vr.shape[1:], dtype='f')

    for i in xrange(Vr.shape[1]):
        for j in xrange(Vr.shape[2]):
            if Vr[:,i,j].max() > 0.0:
                w, l, sc, sh = evd.estimateEVD(Vr[:,i,j],
                                               years,
                                               nodata,
                                               minRecords,
                                               yrsPerSim)

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


    :return: `numpy.ndarray` of return period wind speed values

    """

    lower = (100 - prange) / 2.
    upper = 100. - lower

    nrecords = Vr.shape[0]
    nsamples = nrecords / sample_size
    RpUpper = nodata*np.ones((len(years), Vr.shape[1], Vr.shape[2]), dtype='f')
    RpLower = nodata*np.ones((len(years), Vr.shape[1], Vr.shape[2]), dtype='f')

    w = np.zeros((len(years), nsamples), dtype='f')
    wUpper = np.zeros((len(years)), dtype='f')
    wLower = np.zeros((len(years)), dtype='f')

    for i in xrange(Vr.shape[1]):
        for j in xrange(Vr.shape[2]):
            if Vr[:, i, j].max() > 0.0:
                random.shuffle(Vr[:, i, j])
                for n in xrange(nsamples):
                    nstart = n*sample_size
                    nend  = (n + 1)*sample_size - 1
                    vsub = Vr[nstart:nend, i, j]

                    vsub.sort()
                    if vsub.max( ) > 0.:
                        w[:, n], loc, scale, shp = evd.estimateEVD(vsub, years, nodata,
                                                                   minRecords/10, yrsPerSim)

                for n in range(len(years)):
                    wUpper[n] = percentile(w[n,:], upper)
                    wLower[n] = percentile(w[n,:], lower)

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

    for year in xrange(numSimulations):
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

    ncobj = nctools.ncLoadFile(filename)
    ncobj_vmax = nctools.ncGetVar(ncobj, 'vmax')
    data_subset = ncobj_vmax[ymin:ymax, xmin:xmax]
    ncobj.close()
    return data_subset

def getTiles(tilegrid):
    """
    Helper to obtain a generator that yields tile numbers

    :param tilegrid: :class:`TileGrid` instance
    """

    tilenums = range(tilegrid.num_tiles)
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

    wf_lon, wf_lat = setDomain(inputPath)

    global pp
    pp = attemptParallel()

    log.info("Running hazard calculations")
    TG = TileGrid(gridLimit, wf_lon, wf_lat)
    tiles = getTiles(TG)

    #def progress(i):
    #    callback(i, len(tiles))

    pp.barrier()
    hc = HazardCalculator(configFile, TG,
                          numsimulations,
                          minRecords,
                          yrsPerSim,
                          calculate_confidence)




    hc.dumpHazardFromTiles(tiles)

    pp.barrier()

    hc.saveHazard()

    log.info("Completed hazard calculation")


