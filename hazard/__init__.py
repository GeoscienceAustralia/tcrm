"""
:mod:`hazard` -- Hazard calculation
===================================

This module contains the core objects for the return period 
hazard calculation.

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
import itertools
import numpy as np
import logging as log

from os.path import join as pjoin
from functools import wraps

from Utilities.files import flProgramVersion
from Utilities.config import ConfigParser
import Utilities.nctools as nctools
import evd


def setDomain(inputPath):
    """
    Establish the full extent of input wind field files
    
    Parameters:
    -----------

    :param inputPath: `str` path of folder containing wind field files

    Returns:
    --------

    :param wf_lon: `numpy.ndarray` of longitudes of the wind field

    :param wf_lat: `numpy.ndarray` of latitudes of the wind field

    """

    fileList = os.listdir(inputPath)
    inputFile = pjoin(inputPath, fileList[0])
    ncobj = nctools.ncLoadFile(inputFile)
    wf_lon = nctools.ncGetDims(ncobj, 'lon')
    wf_lat = nctools.ncGetDims(ncobj, 'lat')
    ncobj.close()
    return wf_lon, wf_lat

def disableOnWorkers(f):
    """
    Disable function calculation on workers. Function will
    only be evaluated on the master.
    """
    @wraps(f)
    def wrap(*args, **kwargs):
        if pp.size() > 1 and pp.rank() > 0:
            return
        else:
            return f(*args, **kwargs)
    return wrap

class TileGrid(object):
    """
    Tiling to minimise MemoryErrors and enable parallelisation.

    Parameters:
    -----------
    
    :param gridLimit: :class:`dict` the domain where the hazard will 
                      be calculated. The :class:`dict` should contain 
                      the keys :attr:`xMin`, :attr:`xMax`, 
                      :attr:`yMin` and :attr:`yMax`. The *x* variable 
                      bounds the longitude and the *y* variable bounds
                      the latitude.

    :param wf_lon: `numpy.ndarray` of longitudes of the wind field

    :param wf_lat: `numpy.ndarray` of latitudes of the wind field

    :param xstep: `int` size of the tile in the x-direction.

    :param ystep: `int` size of the tile in the y-direction.

    """

    def __init__(self, gridLimit, wf_lon, wf_lat, xstep=100, ystep=100):
        """
        Initialise the tile grid for dividing up the domain
        
        Parameters:
        -----------

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
        
        ii, = ((wf_lon >= gridLimit['xMin']) & 
               (wf_lon <= gridLimit['xMax'])).nonzero()
        

        jj, = ((wf_lat >= gridLimit['yMin']) & 
               (wf_lat <= gridLimit['yMax'])).nonzero()

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
                self.x_start[k] = i * self.xstep
                self.x_end[k] = min((i + 1) * self.xstep, self.xdim) - 1
                self.y_start[k] = j * self.ystep
                self.y_end[k] = min((j + 1) * self.ystep, self.ydim) - 1
                k += 1
    

    def getGridLimit(self, k):
        """
        Return the limits for tile `k`. x-indices correspond to the 
        east-west coordinate, y-indices correspond to the north-south
        coordinate.

        Parameters:
        -----------

        :param k: `int` tile number

        Returns:
        --------

        :param x1: minimum x-index for tile `k`
        :param x2: maximum x-index for tile `k`
        :param y1: minimum y-index for tile `k`
        :param y2: maximum y-index for tile `k`

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
        
        Returns:
        --------

        :param lon: :class:`numpy.ndarray` containing longitude values
        
        :param lat: :class:`numpy.ndarray` containing latitude values

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
        Initialise HazardCalculator object
        """
        config = ConfigParser()
        config.read(configFile)
        
        self.nodata = -9999.
        self.years = np.array(config.get('HazardInterface', 
                                         'Years').split(',')).astype('f')
        self.outputPath = pjoin(config.get('Output', 'Path'), 'hazard')
        self.inputPath = pjoin(config.get('Output', 'Path'), 'windfield')
        gridLimit = config.geteval('Region', 'gridLimit')
        
        self.numSim = numSim
        self.minRecords = minRecords
        self.yrsPerSim = yrsPerSim
        self.calcCI = calcCI
        
        self.tilegrid = tilegrid
        lon, lat = self.tilegrid.getDomainExtent()
        
        # Create arrays for storing output data:
        self.loc = np.zeros((len(lat), len(lon)), dtype='f')
        self.shp = np.zeros((len(lat), len(lon)), dtype='f')
        self.scale = np.zeros((len(lat), len(lon)), dtype='f')
        self.Rp = np.zeros((len(self.years), len(lat), len(lon)), dtype='f')

        self.global_atts = {'history': ('TCRM hazard simulation - '
                            'return period wind speeds'),
                 'version': flProgramVersion(),
                 'Python_ver': sys.version}
                 

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
        
        Parameters:
        -----------
        
        :param tilelimits: `tuple` of tile limits 
        
        """

        Vr = loadFilesFromPath(self.inputPath, tilelimits)
        Rp, loc, scale, shp = calculate(Vr, self.years, self.nodata, 
                                        self.minRecords, self.yrsPerSim)

        return (tilelimits, Rp, loc, scale, shp)

    def dumpHazardFromTiles(self, tileiter):
        """
        Iterate over tiles to calculate return period hazard levels

        Parameters:
        -----------
        
        :param tileiter: `generator` that yields tuples of tile dimensions.

        """
        
        work_tag = 0 
        result_tag = 1
        if (pp.rank() == 0) and (pp.size() > 1):
            w = 0
            for d in range(1, pp.size()):
                pp.send(tileiter[w], destination=d, tag = work_tag)
                w += 1

            terminated = 0
            while(terminated < pp.size() - 1):

                result, status = pp.receive(pp.any_source, tag=result_tag, 
                                             return_status=True)
                limits, Rp, loc, scale, shp = result

                (xmin, xmax, ymin, ymax) = limits
                self.loc[ymin:ymax, xmin:xmax] = loc
                self.scale[ymin:ymax, xmin:xmax] = scale
                self.shp[ymin:ymax, xmin:xmax] = shp
                self.Rp[:, ymin:ymax, xmin:xmax] = Rp[:, :, :]
                #if self.calcCI:
                #    self.RPupper[:, ymin:ymax, xmin:xmax] = RpUpper[:, :, :]
                #    self.RPlower[:, ymin:ymax, xmin:xmax] = RpLower[:, :, :]

                d = status.source

                if w < len(tileiter):
                    pp.send(tileiter[w], destination=d, tag=work_tag)
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    terminated += 1

        elif (pp.size() > 1) and (pp.rank() != 0):
            while(True):
                W = pp.receive(source=0, tag=work_tag)
                if W is None:
                    break
                results = self.calculateHazard(W)
                pp.send(results, destination=0, tag=result_tag)

        elif pp.size() == 1 and pp.rank() == 0:
            # Assumed no Pypar - helps avoid the need to extend DummyPypar()
            for tile in tileiter:
                result = self.calculateHazard(tile)
                limits, Rp, loc, scale, shp = result

                (xmin, xmax, ymin, ymax) = limits
                self.loc[ymin:ymax, xmin:xmax] = loc
                self.scale[ymin:ymax, xmin:xmax] = scale
                self.shp[ymin:ymax, xmin:xmax] = shp
                self.Rp[:, ymin:ymax, xmin:xmax] = Rp[:, :, :]


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
                    'valid_range': (0.0, 200.)
                }
            },
            1: {
                'name': 'scale',
                'dims': ('lat', 'lon'),
                'values': self.scale, 
                'dtype': 'f',
                'atts': {
                    'long_name': 'Scale parameter for GEV distribution',
                    'units': ''
                }
            },
            2: {
                'name': 'shp',
                'dims': ('lat', 'lon'),
                'values': self.shp,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Shape parameter for GEV distribution',
                    'units': ''
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
                    'valid_range': (0.0, 200.)
                } 
            },
            4: {
                'name': 'wspdupper',
                'dims': ('years', 'lat', 'lon'),
                'values': self.nodata * np.ones((len(self.years), 
                                                 len(lat), 
                                                 len(lon))), 
                'dtype': 'f',
                'atts': {
                    'long_name': 'Upper percentile return period wind speed',
                    'units': 'm/s',
                    'percentile': 95
                } 
            },
            5: {
                'name': 'wspdlower', 
                'dims': ('years', 'lat', 'lon'),
                'values': self.nodata * np.ones((len(self.years), 
                                                 len(lat), 
                                                 len(lon))), 
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lower percentile return period wind speed',
                    'units': 'm/s',
                    'percentile': 5
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

    Parameters:
    -----------

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

def loadFilesFromPath(inputPath, tilelimits):
    """
    Load wind field data for each subset into a 3-D array.

    Parameters:
    -----------

    :param inputPath: str path to wind field files.
    
    :param tilelimits: tuple of index limits of a tile.

    Returns:
    --------

    :param Vr: 3-D `numpy.narray` of wind field records.

    """

    fileList = os.listdir(inputPath)
    files = [pjoin(inputPath, f) for f in fileList]
    files = [f for f in files if os.path.isfile(f)]

    ysize = tilelimits[3] - tilelimits[2]
    xsize = tilelimits[1] - tilelimits[0]
    Vr = np.empty((len(files), ysize, xsize), dtype='f')

    for n, f in enumerate(files):
        Vr[n,:,:] = loadFile(f, tilelimits)

    return Vr

def loadFile(filename, limits):
    """
    Load a subset of the data from the given file, with the extent 
    of the subset specified in the `limits` tuple
    
    Parameters:
    -----------
    
    :param filename: str full path to file to load.

    :param limits: tuple of index limits of a tile.

    Returns:
    --------

    :param data_subset: 2-D `numpy.ndarray` of wind speed values.

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

    Parameters:
    -----------

    :param tilegrid: :class:`TileGrid` instance
    
    :param tilenums: list of tile numbers (must be sequential)

    Returns:
    --------
    
    :param tilelimits: list of tuples of tile imits

    """

    tilelimits = [tilegrid.getGridLimit(t) for t in tilenums]
    return tilelimits

def _getTileLimits(tilegrid, tilenums):
    """
    Generator that yields a tuple of the x- and y-limits of a tile
    
    Parameters:
    -----------
    :param tilegrid: :class:`TileGrid` instance


    """

    for tilenum in balanced(tilenums):
        limits = tilegrid.getGridLimit(tilenum)
        yield limits
    

def balanced(iterable):
    """
    Balance an iterator across processors.

    This partitions the work evenly across processors. However, it requires
    the iterator to have been generated on all processors before hand. This is
    only some magical slicing of the iterator, i.e., a poor man version of
    scattering.

    """

    P, p = pp.size(), pp.rank()
    return itertools.islice(iterable, p, None, P)

def balance(N):
    """
    Compute p'th interval when N is distributed over P bins
    """

    P, p = pp.size(), pp.rank()

    L = int(np.floor(float(N) / P))
    K = N - P * L
    if p < K:
        Nlo = p * L + p
        Nhi = Nlo + L + 1
    else:
        Nlo = p * L + K
        Nhi = Nlo + L

    return Nlo, Nhi


def attemptParallel():
    """
    Attempt to load Pypar globally as `pp`.  If pypar cannot be loaded then a
    dummy `pp` is created.

    """

    global pp

    try:
        # load pypar for everyone

        import pypar as pp

    except ImportError:

        # no pypar, create a dummy one
        
        class DummyPypar(object):

            def size(self):
                return 1

            def rank(self):
                return 0

            def barrier(self):
                pass

        pp = DummyPypar()


def run(configFile, callback=None):
    """
    Run the hazard calculations.

    This will attempt to run the calculation in parallel by tiling the
    domain, but also provides a sane fallback mechanism to execute 
    in serial.

    :param configFile: str

    """
    
    log.info("Loading hazard calculation settings")
    

    config = ConfigParser()
    config.read(configFile)
    
    outputPath = config.get('Output', 'Path')
    gridLimit = config.geteval('Region', 'gridLimit')
    numsimulations = config.getint('HazardInterface', 'NumSim')
    calculate_confidence = config.getboolean('HazardInterface', 'CalculateCI')
    
    inputPath = pjoin(outputPath, 'windfield')

    wf_lon, wf_lat = setDomain(inputPath)

    outputPath = pjoin(config.get('Output', 'Path'), 'hazard')

    years = np.array(config.get('HazardInterface', 
                                'Years').split(',')).astype('f')

    minRecords = config.get('HazardInterface', 'MinimumRecords')
    yrsPerSim = config.get('HazardInterface', 'YearsPerSimulation')
    minRecords=10 

    attemptParallel()

    log.info("Running hazard calculations")
    TG = TileGrid(gridLimit, wf_lon, wf_lat)
    tiles = getTiles(TG)
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
    pp.finalize()
    
    
