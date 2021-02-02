"""
:mod:`TrackGenerator` -- Tropical cyclone track generation
==========================================================

This module contains the core objects for tropical cyclone track
generation.

Track generation can be run in parallel using MPI if the :term:`mpi4py`
library is found and TCRM is run using the :term:`mpirun` command. For
example, to run with 10 processors::

    $ mpirun -n 10 python main.py cairns.ini

:class:`TrackGenerator` can be correctly initialised and started by
calling :meth:`run` with the location of a configuration file::

    >>> import TrackGenerator
    >>> TrackGenerator.run('cairns.ini')

Alternatively, it can be run from the command line::

    $ python TrackGenerator.py cairns.ini


Track format
------------

Tracks are stored in netCDF4 format files, making use of the compound
variable types and heirarchical group structure available in this
version of netCDF. Each event is contained as a separate group in the
file, and each file contains a simulation 'year' (though in the case
of multi-year simulations, multiple years are contained in the file).

Here's an example that shows the group structure in a track file. The
track file contains five tracks, each with only two variables - 'time'
and 'track'::

    >>> from netCDF4 import Dataset
    >>> rootgrp = Dataset("tracks.00001.nc", "r")
    >>> print rootgrp.groups
    OrderedDict([(u'tracks', <netCDF4.Group object at 0x04B92E88>)])
    >>> def walktree(top):
    >>>  values = top.groups.values()
    >>>  yield values
    >>>  for value in top.groups.values():
    >>>   for children in walktree(value):
    >>>    yield children

    >>> print rootgrp
    <type 'netCDF4.Dataset'>
    root group (NETCDF4 file format):
        created_on: 2015-09-23 14:12:20
        created_by: u12161
        dimensions:
        variables:
        groups: tracks

    >>> for children in walktree(rootgrp):
    >>>  for child in children:
    >>>   print child

    <type 'netCDF4.Group'>
    group /tracks:
        dimensions:
        variables:
        groups: tracks-0000, tracks-0001, tracks-0002, tracks-0003, tracks-0004

    <type 'netCDF4.Group'>
    group /tracks\tracks-0000:
        dimensions: time
        variables: time, track
        groups:

    <type 'netCDF4.Group'>
    group /tracks\tracks-0001:
        dimensions: time
        variables: time, track
        groups:

    <type 'netCDF4.Group'>
    group /tracks\tracks-0002:
        dimensions: time
        variables: time, track
        groups:

    <type 'netCDF4.Group'>
    group /tracks\tracks-0003:
        dimensions: time
        variables: time, track
        groups:

    <type 'netCDF4.Group'>
    group /tracks\tracks-0004:
        dimensions: time
        variables: time, track
        groups:

The 'time' variable is an unlimited dimension representing the length
of the track in hours since genesis.

The 'track' variable is a compound variable akin to a
:class:`numpy.recarray`. Here's an example of the summary of a 'track'
variable (taking the last track from the previous example)::

    >>> print child.variables['track']
    <type 'netCDF4.Variable'>
    compound track(time)
        long_name: Tropical cyclone track data
        time_units: hours since 1900-01-01 00:00
        calendar: julian
        lon_units: degrees east
        lat_units: degrees north
        pressure_units: hPa
        speed_units: m/s
        length_units: km
        trackId: (5, 1)
    compound data type: {'names':[u'CycloneNumber',u'Datetime',u'TimeElapsed',
                                  u'Longitude',u'Latitude',u'Speed',
                                  u'Bearing',u'CentralPressure',
                                  u'EnvPressure',u'rMax'],
                         'formats':['<i4','<f8','<f4','<f8','<f8','<f8',
                                    '<f8','<f8','<f8','<f8'],
                         'offsets':[0,8,16,24,32,40,48,56,64,72],
                         'itemsize':80}
    path = /tracks\tracks-0004
    unlimited dimensions: time
    current shape = (220,)

A thorough explanation of the netcdf4-python API is given on the
`Unidata Github page <http://unidata.github.io/netcdf4-python/>`_.

"""

import os
import sys
import logging as log
import math
import time
import itertools
from datetime import datetime, timedelta
from os.path import join as pjoin

import numpy as np


import Utilities.stats as stats
from . import trackLandfall
from . import trackSize
import Utilities.nctools as nctools
import Utilities.maputils as maputils
import Utilities.metutils as metutils
import Utilities.tcrandom as random
from netCDF4 import Dataset as netcdf_file
from scipy.ndimage.interpolation import spline_filter

from StatInterface.generateStats import GenerateStats
from StatInterface.SamplingOrigin import SamplingOrigin
from Utilities.files import flLoadFile, flSaveFile

from DataProcess.CalcFrequency import CalcFrequency
from DataProcess.CalcTrackDomain import CalcTrackDomain
from Utilities.config import ConfigParser
from Utilities.interp3d import interp3d
from Utilities.parallel import attemptParallel
from Utilities.track import ncSaveTracks, Track, TCRM_COLS, TCRM_FMTS
from Utilities.loadData import getPoci

class SamplePressure(object):
    """
    Provide a method to get a 3-d interpolated mean sea level
    pressure at a given location

    :param str mslp_file: path to a 3-d (time, lat, lon) MSLP
                          netcdf file.
    :param str var: Variable name (assumed 'slp')

    """
    def __init__(self, mslp_file, var='slp'):
        ncobj = nctools.ncLoadFile(mslp_file)
        data = nctools.ncGetData(ncobj, var)
        slpunits = getattr(ncobj.variables[var], 'units')
        data = metutils.convert(data, slpunits, 'hPa')
        self.data = spline_filter(data)

    def get_pressure(self, coords):
        """
        Interpolate daily long term mean sea level pressure at
        the given coordinates.

        :param list coords: [day of year, latitude, longitude]
        :rtype: float
        :return: long term MSLP

        """

        scale = [365., 180., 360.]
        offset = [0., -90., 0.]
        mslp = interp3d(self.data, coords, scale, offset, prefilter=False)
        return mslp

class TrackGenerator(object):

    """
    Generates tropical cyclone tracks based on the empirical
    probability distributions of speed, bearing, pressure and maximum
    radius.


    :type  processPath: string
    :param processPath: the location of the empirical probability
                        distribution files.

    :type  gridLimit: :class:`dict`
    :param gridLimit: the domain where the tracks will be generated.
                      The :class:`dict` should contain the keys
                      :attr:`xMin`, :attr:`xMax`, :attr:`yMin` and
                      :attr:`yMax`. The *x* variable bounds the
                      longitude and the *y* variable bounds the
                      latitude.

    :type  gridSpace: :class:`dict`
    :param gridSpace: the grid spacing of the domain.
                      This is used to chop the domain into cells. The
                      :class:`dict` should contain the keys :attr:`x`
                      and :attr:`y`.

    :type  gridInc: :class:`dict`
    :param gridInc: the increments for the grid to be used by
                    :class:`StatInterface.GenerateStats` when
                    insufficient observations exist in the cell being
                    analysed.

    :type  landfall: :class:`LandfallDecay`
    :param landfall: the object that calculates the decay rate of a
                     tropical cyclone after it has made landfall. The
                     object should have the methods :meth:`onLand`
                     and :meth:`pChange`.

    :type  innerGridLimit: :class:`dict`
    :param innerGridLimit: an optional domain limit that can be used to
                           cull tropical cyclone paths. The
                           :class:`dict` should contain the keys
                           :attr:`xMin`, :attr:`xMax`, :attr:`yMin` and
                           :attr:`yMax`. The *x* variable bounds the
                           longitude and the *y* variable bounds the
                           latitude.

    :type  dt: float
    :param dt: the time step used the the track simulation.

    :type  maxTimeSteps: int
    :param maxTimeSteps: the maximum number of tropical cyclone time
                         steps that will be simulated.

    :type  sizeMean: float (default: 46.5.0)
    :param sizeMean: the fallback average tropical cyclone size to use
                     when the empirical distribution data cannot be
                     loaded from file.

    :type  sizeStdDev: float (default: 0.5)
    :param sizeStdDev: the fallback standard deviation of the tropical
                       cyclone size to use when the empirical
                       distribution data cannot be loaded from file.

    """

    def __init__(self, processPath, gridLimit, gridSpace, gridInc, mslp,
                 landfall, innerGridLimit=None, dt=1.0, maxTimeSteps=360,
                 sizeMean=46.5, sizeStdDev=0.5):
        self.processPath = processPath
        self.gridLimit = gridLimit
        self.gridSpace = gridSpace
        self.gridInc = gridInc
        self.mslp = mslp
        self.landfall = landfall
        self.innerGridLimit = innerGridLimit
        self.dt = dt
        self.maxTimeSteps = maxTimeSteps
        self.sizeMean = sizeMean
        self.sizeStdDev = sizeStdDev
        self.timeOverflow = dt * maxTimeSteps
        self.missingValue = sys.maxsize  # FIXME: remove
        self.progressbar = None  # FIXME: remove
        self.allCDFInitBearing = None
        self.allCDFInitSpeed = None
        self.allCDFInitPressure = None
        self.allCDFInitSize = None
        self.cdfSize = None
        self.vStats = None
        self.pStats = None
        self.bStats = None
        self.dpStats = None

        self.dpChi = None
        self.dsChi = None
        self.pChi = None
        self.bChi = None
        self.vChi = None
        self.offshorePressure = None
        self.theta = None
        self.ds = None
        self.dp = None
        self.sChi = None
        self.v = None

        originDistFile = pjoin(processPath, 'originPDF.nc')
        self.originSampler = SamplingOrigin(originDistFile, None, None)
        log.debug("Track domain: {0}".format(repr(self.gridLimit)))
        log.debug("Inner gridLimit: {0}".format(repr(self.innerGridLimit)))

    def loadInitialConditionDistributions(self):
        """
        Load the tropical cyclone empirical distribution data for
        initial bearing, initial speed, and initial pressure. The
        method will try to load the files

            all_cell_cdf_init_bearing.nc
            all_cell_cdf_init_speed.nc
            all_cell_cdf_init_pressure.nc

        from the :attr:`processPath` directory. If it can't find those
        files then it will fallback to the csv versions of those files.
        Note: the csv file names do *not* end with `.csv`.

        The empirical distribution data will be used to initialise
        tropical cyclone tracks with random values when
        :attr:`initBearing`,  :attr:`initSpeed`, and
        :attr:`initPressure` are not provided to
        :meth:`generateTracks`.
        """

        def load(filename):
            """
            Helper function that loads the data from a file.
            """

            # try to load the netcdf version of the file

            if os.path.isfile(filename + '.nc'):
                log.debug('Loading data from %s.nc', filename)
                ncdf = netcdf_file(filename + '.nc', 'r')
                i = ncdf.variables['cell'][:]
                x = ncdf.variables['x'][:]
                y = ncdf.variables['CDF'][:]
                ncdf.close()
                return np.vstack((i, x, y)).T

            # otherwise, revert to old csv format

            log.info('Could not load %s.nc, reverting to old format.',
                     filename)

            try:
                return flLoadFile(filename, '%', ',')
            except IOError:
                log.critical('CDF file %s does not exist!',
                             filename)
                log.critical('Run AllDistribution option in main to' +
                             ' generate those files.')
                raise

        # Load the files

        log.debug('Loading the cyclone parameter data')

        path = self.processPath

        self.allCDFInitBearing = \
            load(pjoin(path, 'all_cell_cdf_init_bearing'))

        self.allCDFInitSpeed = \
            load(pjoin(path, 'all_cell_cdf_init_speed'))

        self.allCDFInitPressure = \
            load(pjoin(path, 'all_cell_cdf_init_pressure'))

        self.allCDFInitDay = \
            load(pjoin(path, 'all_cell_cdf_init_day'))

        try:
            self.allCDFInitSize = load(pjoin(path,
                                             'all_cell_cdf_init_rmax'))

        except IOError:
            log.warning(('RMW distribution file does not exist. '
                         'Using RMW model instead.'))
            self.cdfSize = np.array(stats.rMaxDist(self.sizeMean,
                                                   self.sizeStdDev,
                                                   maxrad=120.0)).T

    def loadCellStatistics(self):
        """
        Load the cell statistics for speed, bearing, pressure, and pressure
        rate of change from netcdf files.

        This method loads the statistics from the files

            speed_stats.nc
            pressure_stats.nc
            bearing_stats.nc
            pressure_rate_stats.nc

        in the :attr:`processPath` directory.
        """

        def init(filename, angular=False):
            """
            Helper function to initialise :class:`GenerateStats`.
            """
            return GenerateStats(
                pjoin(self.processPath, filename),
                pjoin(self.processPath, 'all_lon_lat'),
                self.gridLimit,
                self.gridSpace,
                self.gridInc,
                angular=angular,
                calculateLater=True)

        log.debug('Loading cell statistics for speed from netcdf file')
        self.vStats = init('all_speed')
        self.vStats.load(pjoin(self.processPath, 'speed_stats.nc'))

        log.debug('Loading cell statistics for pressure from netcdf file')
        self.pStats = init('all_pressure')
        self.pStats.load(pjoin(self.processPath, 'pressure_stats.nc'))

        log.debug('Loading cell statistics for bearing from netcdf file')
        self.bStats = init('all_bearing', angular=True)
        self.bStats.load(pjoin(self.processPath, 'bearing_stats.nc'))

        log.debug('Loading cell statistics for pressure_rate from netcdf file')
        self.dpStats = init('pressure_rate')
        self.dpStats.load(pjoin(self.processPath, 'pressure_rate_stats.nc'))

    def generateTracks(self, nTracks, simId, initLon=None, initLat=None,
                       initSpeed=None, initBearing=None,
                       initPressure=None, initEnvPressure=None,
                       initRmax=None, initDay=None):
        """
        Generate tropical cyclone tracks from a single genesis point.

        If the initial conditions for speed, pressure, and bearing are
        not provided then they will be drawn randomly from the
        empirical distributions that were calculated from historical
        data (see :meth:`loadInitialConditionDistributions`).

        If the genesis point (initLon, initLat) is not provided, then
        an origin will be randomly chosen using the empirical genesis
        point distribution calculated by
        :class:`StatInterface.SamplingOrigin`. However, if this random
        point (or the initial point given) falls outside the domain
        defined by :attr:`gridLimit` then no tracks will be generated.

        :type  nTracks: int
        :param nTracks: the number of tracks to generate from the
                        genesis point.

        :type  initLon: float
        :param initLon: the longitude of the genesis point.

        :type  initLat: float
        :param initLat: the latitude of the genesis point.

        :type  initSpeed: float
        :param initSpeed: the initial speed of the tropical cyclone.

        :type  initBearing: float
        :param initBearing: the initial bearing of the tropical
                            cyclone.

        :type  initPressure: float
        :param initPressure: the initial pressure of the tropical
                             cyclone.

        :type  initEnvPressure: float
        :param initEnvPressure: the initial environment pressure.

        :type  initRmax: float
        :param initRmax: the initial maximum radius of the tropical
                         cyclone.

        :rtype :class:`numpy.array`
        :return: the tracks generated.
        """
        # Define some filter functions

        def empty(track):
            """
            :return: True if the track is empty. False, otherwise.
            """
            return len(track.Longitude) == 0

        def diedEarly(track, minAge=12):
            """
            :return: True if the track dies before `minAge`. False,
            otherwise.
            """
            return track.TimeElapsed[-1] < minAge

        def insideDomain(track):
            """
            :return: True if the track stays inside the domain. False,
            otherwise.
            """
            inside = [track.Longitude[k] > self.innerGridLimit['xMin'] and
                      track.Longitude[k] < self.innerGridLimit['xMax'] and
                      track.Latitude[k] > self.innerGridLimit['yMin'] and
                      track.Latitude[k] < self.innerGridLimit['yMax']
                      for k in range(len(track.Longitude))]
            if not any(inside):
                log.debug("Track is outside inner grid domain")
                log.debug("Track id: {0}-{1}".format(*track.trackId))

            return all(inside)

        def validPressures(track):
            """
            :return: True if a valid pressure. False, otherwise.
            """
            if all(np.round(track.CentralPressure, 2) <\
                       np.round(track.EnvPressure, 2)):
                return True
            else:
                log.debug(("Invalid pressure values in track "
                           "{0}-{1}".format(*track.trackId)))
                return False

        def validSize(track):
            """
            :return: True if all rmax values are > 5.0 km, False otherwise.
            """

            if any(track.rMax > 500.):
                log.debug("Track {0}-{1} has rMax > 500 km"\
                              .format(*track.trackId))
            if any(track.rMax < 5.):
                log.debug("Track {0}-{1} has rMax < 5 km"\
                              .format(*track.trackId))
            return (all(track.rMax > 5.) and all(track.rMax < 500.))

        def validInitSize(track):
            return track.rMax[0] < 100.

        def validInitPressure(track):
            if track.EnvPressure[0] - track.CentralPressure[0] < 1.0:
                log.critical("Initial pressure difference too small")
            return (track.EnvPressure[0] - track.CentralPressure[0] > 1.0)

        log.debug('Generating %d tropical cyclone tracks', nTracks)
        genesisYear = int(uniform(1900, 9998))
        results = []
        j = 0
        while j < nTracks:

            if not (initLon and initLat):
                log.debug('Cyclone origin not given, sampling a' +
                          ' random one instead.')
                genesisLon, genesisLat = \
                    self.originSampler.ppf(uniform(), uniform())
            else:
                log.debug('Using prescribed initial position' +
                          ' ({0:.2f}, {1:.2f})'.format(initLon, initLat))
                genesisLon = initLon
                genesisLat = initLat

            # Get the initial grid cell

            initCellNum = stats.getCellNum(genesisLon, genesisLat,
                                           self.gridLimit,
                                           self.gridSpace)

            log.debug('Cyclones origin: (%6.2f, %6.2f) Cell: %i' +
                      ' Grid: %s', genesisLon, genesisLat, initCellNum,
                      self.gridLimit)

            # Sample an initial bearing if none is provided

            if not initBearing:
                ind = self.allCDFInitBearing[:, 0] == initCellNum
                cdfInitBearing = self.allCDFInitBearing[ind, 1:3]
                genesisBearing = ppf(uniform(), cdfInitBearing)
            else:
                genesisBearing = initBearing

            # Sample an initial speed if none is provided

            if not initSpeed:
                ind = self.allCDFInitSpeed[:, 0] == initCellNum
                cdfInitSpeed = self.allCDFInitSpeed[ind, 1:3]
                genesisSpeed = ppf(uniform(), cdfInitSpeed)
            else:
                genesisSpeed = initSpeed


            # Sample an initial day if none is provided

            if not initDay:
                ind = self.allCDFInitDay[:, 0] == initCellNum
                cdfInitDay = self.allCDFInitDay[ind, 1:3]
                genesisDay = ppf(uniform(), cdfInitDay)
            else:
                genesisDay = initDay

            genesisHour = int(uniform(0, 24))

            initTimeStr = "%04d-%03d %d:00" % (genesisYear,
                                               genesisDay,
                                               genesisHour)
            genesisTime = datetime.strptime(initTimeStr, "%Y-%j %H:%M")

            # Sample an initial environment pressure if none is
            # provided - dependent on initial day of year:

            if not initEnvPressure:
                initEnvPressure = \
                    self.mslp.get_pressure(np.array([[genesisDay],
                                                     [genesisLat],
                                                     [genesisLon]]))

            # Sample an initial pressure if none is provided

            if not initPressure:
                # Sample subject to the constraint initPressure <
                # initEnvPressure
                ind = self.allCDFInitPressure[:, 0] == initCellNum
                cdfInitPressure = self.allCDFInitPressure[ind, 1:3]
                ix = cdfInitPressure[:, 0].searchsorted(initEnvPressure)
                upperProb = cdfInitPressure[ix - 1, 1]
                genesisPressure = ppf(uniform(0.0, upperProb),
                                      cdfInitPressure)
            else:
                genesisPressure = initPressure

            # Sample an initial maximum radius if none is provided

            if not initRmax:
                if not self.allCDFInitSize:
                    cdfSize = self.cdfSize[:, [0, 2]]
                else:
                    ind = self.allCDFInitSize[:, 0] == initCellNum
                    cdfSize = self.allCDFInitSize[ind, 1:3]

                dp = initEnvPressure - genesisPressure
                self.rmwEps = np.random.normal(0, scale=0.335)
                genesisRmax = trackSize.rmax(dp, genesisLat, self.rmwEps)

                # Censor the initial Rmax to be < 100 km.
                if genesisRmax > 100.:
                    while genesisRmax > 100.:
                        self.rmwEps = np.random.normal(0, scale=0.335)
                        genesisRmax = trackSize.rmax(dp,
                                                     genesisLat,
                                                     self.rmwEps)

            else:
                genesisRmax = initRmax
            # Do not generate tracks from this genesis point if we are
            # going to exit the domain on the first step

            nextLon, nextLat = \
                maputils.bear2LatLon(genesisBearing,
                                     self.dt * genesisSpeed,
                                     genesisLon, genesisLat)

            log.debug('initBearing: %.2f initSpeed: %.2f' +
                      ' initEnvPressure: %.2f initPressure: %.2f' +
                      ' initRmax: %.2f',
                      genesisBearing, genesisSpeed, initEnvPressure,
                      genesisPressure, genesisRmax)

            xMin = self.gridLimit['xMin']
            xMax = self.gridLimit['xMax']
            yMin = self.gridLimit['yMin']
            yMax = self.gridLimit['yMax']

            if not ((xMin <= nextLon <= xMax) and
                    (yMin <= nextLat <= yMax)):
                log.debug('Tracks will exit domain immediately' +
                          ' for this genesis point.')
                continue

            if (initEnvPressure - genesisPressure) < 2:
                log.debug(("Track does not start with sufficient "
                           "pressure deficit"))
                continue

            log.debug('** Generating track %i from point (%.2f,%.2f)',
                      j, genesisLon, genesisLat)

            data = self._singleTrack(j, genesisLon, genesisLat,
                                     genesisSpeed, genesisBearing,
                                     genesisPressure, initEnvPressure,
                                     genesisRmax, genesisTime)

            track_dtype = np.dtype({'names':TCRM_COLS, 'formats':TCRM_FMTS})
            data = np.array(data)
            data = np.core.records.fromarrays(data, dtype=track_dtype)
            track = Track(data)
            track.trackId = (j, simId)

            if not (empty(track) or diedEarly(track))\
                    and validInitPressure(track) \
               and validPressures(track) and validSize(track) \
               and validInitSize(track):
                if self.innerGridLimit and not insideDomain(track):
                    log.debug("Track exits inner grid limit - rejecting")
                    continue
                else:
                    results.append(track)
                    log.debug("Completed track {0:03d}-{1:04d}".\
                              format(*track.trackId))
                    j += 1
            else:
                log.debug("Eliminated invalid track")

        return results

    def generateTracksToFile(self, outputFile, nTracks, simId, initLon=None,
                             initLat=None, initSpeed=None,
                             initBearing=None, initPressure=None,
                             initEnvPressure=None, initRmax=None,
                             initDay=None):
        """
        Generate tropical cyclone tracks from a single genesis point
        and save the tracks to a file.

        This is a helper function that calls :meth:`generateTracks`.

        :type  outputFile: str
        :param outputFile: the filename of the file where the tracks
                           will be saved. If `outputFile` has the `shp`
                           extension then it will be saved to a shp
                           file. Otherwise, the tracks will be saved in
                           csv format.
        """

        results = self.generateTracks(
            nTracks, simId, initLon=initLon, initLat=initLat,
            initSpeed=initSpeed, initBearing=initBearing,
            initPressure=initPressure,
            initEnvPressure=initEnvPressure,
            initRmax=initRmax, initDay=initDay)

        if outputFile.endswith("shp"):
            from Utilities.shptools import shpSaveTrackFile
            from Utilities.AsyncRun import AsyncRun

            log.debug('Outputting data into %s', outputFile)

            fields = {}

            fields['Index'] = {
                'Type': 1,
                'Length': 5,
                'Precision': 0,
                'Data': results[:, 0]
            }

            fields['Time'] = {
                'Type': 2,
                'Length': 7,
                'Precision': 1,
                'Data': results[:, 1]
            }

            fields['Longitude'] = {
                'Type': 2,
                'Length': 7,
                'Precision': 2,
                'Data': results[:, 2]
            }

            fields['Latitude'] = {
                'Type': 2,
                'Length': 7,
                'Precision': 2,
                'Data': results[:, 3]
            }

            fields['Speed'] = {
                'Type': 2,
                'Length': 6,
                'Precision': 1,
                'Data': results[:, 4]
            }

            fields['Bearing'] = {
                'Type': 2,
                'Length': 6,
                'Precision': 1,
                'Data': results[:, 5]
            }

            fields['Pressure'] = {
                'Type': 2,
                'Length': 6,
                'Precision': 1,
                'Data': results[:, 6]
            }

            fields['pEnv'] = {
                'Type': 2,
                'Length': 6,
                'Precision': 1,
                'Data': results[:, 7]
            }

            fields['rMax'] = {
                'Type': 2,
                'Length': 5,
                'Precision': 1,
                'Data': results[:, 8]
            }

            args = {
                'filename': outputFile,
                'lon': results[:, 2],
                'lat': results[:, 3],
                'fields': fields
            }

            thr = AsyncRun(shpSaveTrackFile, args)
            try:
                thr.start()
            except:
                raise
        else:
            log.debug('Outputting data into %s', outputFile)

            header = 'CycloneNumber,TimeElapsed(hr),' + \
                     'Longitude(degree),Latitude(degree),' + \
                     'Speed(km/hr),Bearing(degrees),' + \
                     'CentralPressure(hPa),EnvPressure(hPa),rMax(km)'

            args = {
                'filename': outputFile,
                'data': results,
                'header': header,
                'delimiter': ',',
                'fmt': '%7.2f'
            }

            fl = AsyncRun(flSaveFile, args)
            fl.start()

    def _singleTrack(self, cycloneNumber, initLon, initLat, initSpeed,
                     initBearing, initPressure, initEnvPressure,
                     initRmax, initTime):
        """
        Generate a single tropical cyclone track from a genesis point.

        :type  cycloneNumber: int
        :param cycloneNumer: the tropical cyclone index.

        :type  initLon: float
        :param initLon: the longitude of the genesis point.

        :type  initLat: float
        :param initLat: the latitude of the genesis point.

        :type  initSpeed: float
        :param initSpeed: the initial speed of the tropical cyclone.

        :type  initBearing: float
        :param initBearing: the initial bearing of the tropical
                            cyclone.

        :type  initPressure: float
        :param initPressure: the initial pressure of the tropical
                             cyclone.

        :type  initEnvPressure: float
        :param initEnvPressure: the initial environment pressure.

        :type  initRmax: float
        :param initRmax: the initial maximum radius of the tropical
                         cyclone.

        :type  initDay: float
        :param initDay: the initial day of year of the tropical cyclone.

        :return: a tuple of :class:`numpy.ndarray`'s
                 The tuple consists of::

                      index - the tropical cyclone index
                      age - age of the tropical cyclone
                      lon - longitude
                      lat - latitude
                      speed
                      bearing
                      pressure
                      penv - environment pressure
                      rmax - maximum radius
        """

        index = np.ones(self.maxTimeSteps, 'f') * cycloneNumber
        dates = np.empty(self.maxTimeSteps, dtype=datetime)
        age = np.empty(self.maxTimeSteps, 'f')
        jday = np.empty(self.maxTimeSteps, 'f')
        lon = np.empty(self.maxTimeSteps, 'f')
        lat = np.empty(self.maxTimeSteps, 'f')
        speed = np.empty(self.maxTimeSteps, 'f')
        bearing = np.empty(self.maxTimeSteps, 'f')
        pressure = np.empty(self.maxTimeSteps, 'f')
        poci = np.empty(self.maxTimeSteps, 'f')
        rmax = np.empty(self.maxTimeSteps, 'f')
        land = np.empty(self.maxTimeSteps, 'i')
        dist = np.empty(self.maxTimeSteps, 'f')

        # Initialise the track
        poci_eps = normal(0., 2.5717)
        lfeps = lognorm(0.69527, -0.06146, 0.0471)

        age[0] = 0
        dates[0] = initTime
        jday[0] = int(initTime.strftime("%j")) + initTime.hour/24.
        lon[0] = initLon
        lat[0] = initLat
        speed[0] = initSpeed
        bearing[0] = initBearing
        pressure[0] = initPressure
        poci[0] = getPoci(initEnvPressure, initPressure,
                          initLat, jday[0], poci_eps)
        rmax[0] = initRmax
        land[0] = 0
        dist[0] = self.dt * speed[0]

        timestep = timedelta(self.dt/24.)

        # Initialise variables that will be used when performing a step
        self.offshorePressure = initPressure
        self.offshorePoci = poci[0]
        self.landfallSpeed = initSpeed
        self.theta = initBearing
        self.v = initSpeed
        self.vChi = 0.0
        self.bChi = 0.0
        self.pChi = 0.0
        self.sChi = 0.0
        self.dpChi = 0.0
        self.dsChi = 0.0
        self.dp = 0.0
        self.ds = 0.0

        # Initialise the landfall time over land (`tol`)
        tol = 0.0

        # Generate the track
        for i in range(1, self.maxTimeSteps):

            # Get the new latitude and longitude from bearing and
            # distance
            lon[i], lat[i] = maputils.bear2LatLon(bearing[i - 1],
                                                  dist[i - 1],
                                                  lon[i - 1],
                                                  lat[i - 1])

            age[i] = age[i - 1] + self.dt
            dates[i] = dates[i - 1] + timestep
            jday[i] = jday[i - 1] + self.dt/24.
            jday[i] = np.mod(jday[i], 365)

            # Sample the environment pressure
            penv = self.mslp.get_pressure(np.array([[jday[i]],
                                                    [lat[i]],
                                                    [lon[i]]]))

            # Terminate and return the track if it steps out of the
            # domain
            if (lon[i] < self.gridLimit['xMin'] or
                    lon[i] >= self.gridLimit['xMax'] or
                    lat[i] <= self.gridLimit['yMin'] or
                    lat[i] > self.gridLimit['yMax']):

                log.debug('TC exited domain at point ' +
                          '(%.2f %.2f) and time %i', lon[i], lat[i], i)

                return (index[:i], dates[:i], age[:i], lon[:i], lat[:i],
                        speed[:i], bearing[:i], pressure[:i],
                        poci[:i], rmax[:i])

            cellNum = stats.getCellNum(lon[i], lat[i],
                                       self.gridLimit, self.gridSpace)
            onLand = self.landfall.onLand(lon[i], lat[i])

            land[i] = onLand

            # Do the real work: generate a step of the model
            self._stepPressureChange(cellNum, i, onLand)
            self._stepBearing(cellNum, i, onLand)
            self._stepSpeed(cellNum, i, onLand)

            # Update bearing, speed and pressure

            bearing[i] = self.theta
            speed[i] = abs(self.v)  # reflect negative speeds

            # Calculate the central pressure

            if onLand:
                tol += float(self.dt)
                deltaP = self.offshorePoci - self.offshorePressure
                alpha = 0.03515 + 0.000435 * deltaP +\
                        0.002865 * self.landfallSpeed + lfeps
                pressure[i] = poci[i - 1] - deltaP * np.exp(-alpha * tol)
                poci[i] = getPoci(penv, pressure[i], lat[i], jday[i], poci_eps)
                log.debug('alpha value for landfall decay: {0}'.format(alpha))
                log.debug(('Central pressure {0} hours after landfall:'
                           '{1:.2f}'.format(tol, pressure[i])))
            else:
                pstat = self.pStats.coeffs
                pressure[i] = pressure[i - 1] + self.dp * self.dt

                # If the central pressure of the synthetic storm is
                # more than 4 std deviations lower than the minimum
                # observed central pressure, automatically start
                # raising the central pressure.

                if (pressure[i] < (pstat.min[cellNum] -
                                   4. * pstat.sig[cellNum])):
                    log.debug('Recalculting pressure as extremely low')
                    pressure[i] = (pressure[i - 1] +
                                   abs(self.dp) * self.dt)

                self.offshorePressure = pressure[i]
                self.landfallSpeed = speed[i]

                poci[i] = getPoci(penv, pressure[i], lat[i], jday[i], poci_eps)
                self.offshorePoci = poci[i]

            # If the empirical distribution of tropical cyclone size is
            # loaded then sample and update the maximum radius.
            # Otherwise, keep the maximum radius constant.

            if self.allCDFInitSize:
                self._stepSizeChange(cellNum, i, onLand)
                rmax[i] = rmax[i - 1] + self.ds * self.dt
                # if the radius goes below 1.0, then do an
                # antithetic increment instead
                if rmax[i] <= 1.0:
                    rmax[i] = rmax[i - 1] - self.ds * self.dt
            else:
                dp = poci[i] - pressure[i]
                rmax[i] = trackSize.rmax(dp, lat[i], self.rmwEps)

            # Update the distance and the age of the cyclone

            dist[i] = self.dt * speed[i]

            # Terminate the track if it doesn't satisfy certain criteria

            if self._notValidTrackStep(pressure[i], poci[i], age[i],
                                       lon[0], lat[0], lon[i], lat[i]):
                log.debug('Track no longer satisfies criteria, ' +
                          'terminating at time %i.', i)

                return (index[:i], dates[:i], age[:i], lon[:i], lat[:i],
                        speed[:i], bearing[:i], pressure[:i], poci[:i],
                        rmax[:i])

        return (index, dates, age, lon, lat, speed, bearing, pressure,
                poci, rmax)

    def _stepPressureChange(self, c, i, onLand):
        """
        Take one step of the pressure change model.

        This updates :attr:`self.dpChi` and :attr:`self.dp` based on an
        (inhomogeneous) AR(1) model.

        :type  c: int
        :param c: a valid cell index in the domain

        :type  i: int
        :param i: the step number (i.e., time)

        :type  onLand: bool
        :param onLand: True if the tropical cyclone is currently over
                       land.
        """

        # Change the parameter set accordingly

        if onLand:
            alpha = self.dpStats.coeffs.lalpha
            phi = self.dpStats.coeffs.lphi
            mu = self.dpStats.coeffs.lmu
            sigma = self.dpStats.coeffs.lsig
        else:
            alpha = self.dpStats.coeffs.alpha
            phi = self.dpStats.coeffs.phi
            mu = self.dpStats.coeffs.mu
            sigma = self.dpStats.coeffs.sig

        # Do the step

        self.dpChi = alpha[c] * self.dpChi + phi[c] * logistic()

        if i == 1:
            self.dp += sigma[c] * self.dpChi
        else:
            self.dp = mu[c] + sigma[c] * self.dpChi

    def _stepBearing(self, c, i, onLand):
        """
        Take one step of the bearing model.

        This updates :attr:`self.bChi` and :attr:`self.theta` based on
        an (inhomogeneous) AR(1) model.

        :type  c: int
        :param c: a valid cell index in the domain

        :type  t: int
        :param t: the step number (i.e., time)

        :type  onLand: bool
        :param onLand: True if the tropical cyclone is currently over
                       land.
        """

        # Change the parameter set accordingly

        if onLand:
            alpha = self.bStats.coeffs.lalpha
            phi = self.bStats.coeffs.lphi
            mu = self.bStats.coeffs.lmu
            sigma = self.bStats.coeffs.lsig
        else:
            alpha = self.bStats.coeffs.alpha
            phi = self.bStats.coeffs.phi
            mu = self.bStats.coeffs.mu
            sigma = self.bStats.coeffs.sig

        # Do the step

        self.bChi = alpha[c] * self.bChi + phi[c] * logistic()

        # Update the bearing

        if i == 1:
            self.theta += math.degrees(sigma[c] * self.bChi)
        else:
            self.theta = math.degrees(mu[c] + sigma[c] * self.bChi)

        self.theta = np.mod(self.theta, 360.)

    def _stepSpeed(self, c, i, onLand):
        """
        Take one step of the speed model.

        This updates :attr:`self.vChi` and :attr:`self.v` based on an
        (inhomogeneous) AR(1) model.

        :type  c: int
        :param c: a valid cell index in the domain

        :type  t: int
        :param t: the step number (i.e., time)

        :type  onLand: bool
        :param onLand: True if the tropical cyclone is currently over
                       land.
        """

        # Change the parameter set accordingly

        if onLand:
            alpha = self.vStats.coeffs.lalpha
            phi = self.vStats.coeffs.lphi
            mu = self.vStats.coeffs.lmu
            sigma = self.vStats.coeffs.lsig
        else:
            alpha = self.vStats.coeffs.alpha
            phi = self.vStats.coeffs.phi
            mu = self.vStats.coeffs.mu
            sigma = self.vStats.coeffs.sig

        # Do the step

        self.vChi = alpha[c] * self.vChi + phi[c] * logistic()

        # Update the speed

        if i == 1:
            self.v += abs(sigma[c] * self.vChi)
        else:
            self.v = abs(mu[c] + sigma[c] * self.vChi)

    def _stepSizeChange(self, c, i, onLand):
        """
        Take one step of the size change model.

        This updates :attr:`self.vChi` and :attr:`self.v` based on an
        (inhomogeneous) AR(1) model.

        :type  c: int
        :param c: a valid cell index in the domain

        :type  t: int
        :param t: the step number (i.e., time)

        :type  onLand: bool
        :param onLand: True if the tropical cyclone is currently over
                       land.
        """

        # Change the parameter set accordingly

        if onLand:
            alpha = self.dsStats.coeffs.lalpha
            phi = self.dsStats.coeffs.lphi
            mu = self.dsStats.coeffs.lmu
            sigma = self.dsStats.coeffs.lsig
        else:
            alpha = self.dsStats.coeffs.alpha
            phi = self.dsStats.coeffs.phi
            mu = self.dsStats.coeffs.mu
            sigma = self.dsStats.coeffs.sig

        # Do the step

        self.dsChi = alpha[c] * self.dsChi + phi[c] * logistic()

        # Update the size change

        if i == 1:
            self.ds += sigma[c] * self.dsChi
        else:
            self.ds = mu[c] + sigma[c] * self.dsChi

    def _notValidTrackStep(self, pressure, poci, age, lon0, lat0,
                           nextlon, nextlat):
        """
        This is called to check if a tropical cyclone track meets
        certain conditions.
        """
        if np.isnan(poci):
            return True
        if age > 12 and ((poci - pressure) < 5.0):
            log.debug('Pressure difference < 5.0' +
                      ' (penv: %f pressure: %f)', poci, pressure)
            return True
        elif age <= 12 and ((poci - pressure) < 2.0):
            log.debug('Pressure difference < 2.0' +
                      ' (penv: %f pressure: %f)', poci, pressure)
            return True

        return False

    def dumpAllCellCoefficients(self):
        """
        Dump all cell coefficients to a netcdf file to permit further
        analysis.
        """
        lon = np.arange(self.gridLimit['xMin'],
                        self.gridLimit['xMax'],
                        self.gridSpace['x'])

        lat = np.arange(self.gridLimit['yMax'],
                        self.gridLimit['yMin'],
                        -1 * self.gridSpace['y'])

        nx = len(lon)
        ny = len(lat)

        dimensions = {
            0: {
                'name': 'lat',
                'values': lat,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Latitude',
                    'units': 'degrees_north'
                }
            },
            1: {
                'name': 'lon',
                'values': lon,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Longitude',
                    'units': 'degrees_east'
                }
            }
        }

        variables = {
            0: {
                'name': 'vmu',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.mu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean forward speed',
                    'units': 'km/h'
                }
            },
            1: {
                'name': 'valpha',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.alpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of forward' +
                                 ' speed',
                    'units': ''
                }
            },
            2: {
                'name': 'vsig',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.sig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation forward speed',
                    'units': 'km/h'
                }
            },
            3: {
                'name': 'vmin',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.min.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum forward speed',
                    'units': 'km/h'
                }
            },
            4: {
                'name': 'vlmu',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.lmu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean forward speed (over land)',
                    'units': 'km/h'
                }
            },
            5: {
                'name': 'vlalpha',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.lalpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of forward' +
                                 ' speed (over land)',
                    'units': ''
                }
            },
            6: {
                'name': 'vlsig',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.lsig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of forward' +
                                 ' speed (over land)',
                    'units': 'km/h'
                }
            },
            7: {
                'name': 'vlmin',
                'dims': ('lat', 'lon'),
                'values': self.vStats.coeffs.lmin.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum forward speed (over land)',
                    'units': 'km/h'
                }
            },
            8: {
                'name': 'bmu',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.mu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean bearing',
                    'units': 'degrees'
                }
            },
            9: {
                'name': 'balpha',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.alpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of bearing',
                    'units': ''
                }
            },
            10: {
                'name': 'bsig',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.sig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of bearing',
                    'units': 'degrees'
                }
            },
            11: {
                'name': 'bmin',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.min.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum bearing',
                    'units': 'degrees'
                }
            },
            12: {
                'name': 'blmu',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.lmu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean bearing(over land)',
                    'units': 'degrees'
                }
            },
            13: {
                'name': 'blalpha',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.lalpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of bearing' +
                                 ' (over land)',
                    'units': ''
                }
            },
            14: {
                'name': 'blsig',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.lsig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of bearing' +
                                 ' (over land)',
                    'units': 'degrees'
                }
            },
            15: {
                'name': 'blmin',
                'dims': ('lat', 'lon'),
                'values': self.bStats.coeffs.lmin.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum bearing (over land)',
                    'units': 'degrees'
                }
            },
            16: {
                'name': 'pmu',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.mu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean central pressure',
                    'units': 'hPa'
                }
            },
            17: {
                'name': 'palpha',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.alpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of central' +
                                 ' pressure',
                    'units': ''
                }
            },
            18: {
                'name': 'psig',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.sig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of central' +
                                 ' pressure',
                    'units': 'hPa'
                }
            },
            19: {
                'name': 'pmin',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.min.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum central pressure',
                    'units': 'hPa'
                }
            },
            20: {
                'name': 'plmu',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.lmu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean central pressure (over land)',
                    'units': 'hPa'
                }
            },
            21: {
                'name': 'plalpha',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.lalpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of central' +
                                 ' pressure (over land)',
                    'units': ''
                }
            },
            22: {
                'name': 'plsig',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.lsig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of central' +
                                 ' pressure (over land)',
                    'units': 'hPa'
                }
            },
            23: {
                'name': 'plmin',
                'dims': ('lat', 'lon'),
                'values': self.pStats.coeffs.lmin.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum central pressure (over land)',
                    'units': 'hPa'
                }
            },
            24: {
                'name': 'dpmu',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.mu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean rate of pressure change',
                    'units': 'hPa/h'
                }
            },
            25: {
                'name': 'dpalpha',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.alpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of rate of' +
                                 ' pressure change',
                    'units': ''
                }
            },
            26: {
                'name': 'dpsig',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.sig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of rate of' +
                                 ' pressure change',
                    'units': 'hPa/h'
                }
            },
            27: {
                'name': 'dpmin',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.min.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum rate of pressure change',
                    'units': 'hPa/h'
                }
            },
            28: {
                'name': 'dplmu',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.lmu.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean rate of pressure change' +
                                 ' (over land)',
                    'units': 'hPa/h'
                }
            },
            29: {
                'name': 'dplalpha',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.lalpha.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Lag-1 autocorrelation of rate of' +
                                 ' pressure change (over land)',
                    'units': ''
                }
            },
            30: {
                'name': 'dplsig',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.lsig.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Standard deviation of rate of' +
                                 ' pressure change (over land)',
                    'units': 'hPa/h'
                }
            },
            31: {
                'name': 'dplmin',
                'dims': ('lat', 'lon'),
                'values': self.dpStats.coeffs.lmin.reshape((ny, nx)),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Minimum rate of pressure' +
                                 ' change (over land)',
                    'units': 'hPa/h'
                }
            }
        }

        outputFile = pjoin(self.processPath, 'coefficients.nc')
        nctools.ncSaveGrid(outputFile, dimensions, variables,
                           nodata=self.missingValue, datatitle=None,
                           writedata=True, keepfileopen=False)

# Define a global pseudo-random number generator. This is done to
# ensure we are sampling correctly across processors when performing
# the simulation in parallel. We use the inbuilt Python `random`
# library as it provides the ability to `jumpahead` in the stream (as
# opposed to `numpy.random`).

PRNG = random.Random(seed=1234, stream=0)


def normal(mean=0.0, stddev=1.0):
    """
    Sample from a Normal distribution.
    """
    return PRNG.normalvariate(mean, stddev)


def uniform(a=0.0, b=1.0):
    """
    Sample from a uniform distribution.
    """
    return PRNG.uniform(a, b)

def logistic(loc=0., scale=1.0):
    """
    Sample from a logistic distribution.
    """
    return PRNG.logisticvariate(loc, scale)

def nct(df, nc, loc=0.0, scale=1.0):
    """
    Sample from a non-central T distribution.
    """
    return PRNG.nctvariate(df, nc, loc, scale)

def lognorm(xi, loc=0., scale=1.0):
    """
    Sample from a lognormal distribution.
    """
    return PRNG.lognormvariate(xi, loc, scale)

def ppf(q, cdf):
    """
    Percentage point function (aka. inverse CDF, quantile) of
    an empirical CDF.

    This is used to sample from an empirical distribution.
    """
    i = cdf[:, 1].searchsorted(q)
    return cdf[i, 0]


def balanced(iterable):
    """
    Balance an iterator across processors.

    This partitions the work evenly across processors. However, it
    requires the iterator to have been generated on all processors
    before hand. This is only some magical slicing of the iterator,
    i.e., a poor man version of scattering.
    """
    P, p = MPI.COMM_WORLD.size, MPI.COMM_WORLD.rank
    return itertools.islice(iterable, p, None, P)


class Simulation(object):

    """
    Simulation parameters.

    This is used to set the PRNG state before `ntracks` are simulated.

    :type  index: int
    :param index: the simulation index number.

    :type  seed: int
    :param seed: the initial seed used for the PRNG.

    :type  jumpahead: int
    :param jumpahead: the amount to jump ahead from the initial seed in
                      the PRNG stream.

    :type  ntracks: int
    :param ntracks: the number of tracks to be generated during the
                    simulation.

    :type  outfile: str
    :param outfile: the filename where the tracks will be saved to.
    """

    def __init__(self, index, seed, jumpahead, ntracks, outfile):
        self.index = index
        self.seed = seed
        self.jumpahead = jumpahead
        self.ntracks = ntracks
        self.outfile = outfile


def run(configFile, callback=None):
    """
    Run the tropical cyclone track generation.

    This will attempt to perform the simulation in parallel but also
    provides a sane fallback mechanism.

    :type  configFile: str
    :param configFile: the filename of the configuration file to load
                       the track generation configuration from.
    """

    log.info('Loading track generation settings')

    # Get configuration

    config = ConfigParser()
    config.read(configFile)

    outputPath = config.get('Output', 'Path')
    nSimulations = config.getint('TrackGenerator', 'NumSimulations')
    maxTimeSteps = config.getint('TrackGenerator', 'NumTimeSteps')
    dt = config.getfloat('TrackGenerator', 'TimeStep')
    fmt = config.get('TrackGenerator', 'Format')
    gridSpace = config.geteval('Region', 'GridSpace')
    gridInc = config.geteval('Region', 'GridInc')
    gridLimit = config.geteval('Region', 'gridLimit')
    mslpFile = config.get('Input', 'MSLPFile')
    mslpVar = config.get('Input', 'MSLPVariableName')
    seasonSeed = None
    trackSeed = None
    trackPath = pjoin(outputPath, 'tracks')
    processPath = pjoin(outputPath, 'process')
    trackFilename = 'tracks.%05i-%%04i.' + fmt

    if config.has_option('TrackGenerator', 'gridLimit'):
        gridLimit = config.geteval('TrackGenerator', 'gridLimit')
    else:
        CalcTD = CalcTrackDomain(configFile)
        gridLimit = CalcTD.calcDomainFromFile()

    if config.has_option('TrackGenerator', 'Frequency'):
        meanFreq = config.getfloat('TrackGenerator', 'Frequency')
    else:
        log.info('No genesis frequency specified: auto-calculating')
        CalcF = CalcFrequency(configFile, gridLimit)
        meanFreq = CalcF.calc()
        log.info('Estimated annual genesis frequency for domain: %s',
                 meanFreq)

    if config.has_option('TrackGenerator', 'SeasonSeed'):
        seasonSeed = config.getint('TrackGenerator', 'SeasonSeed')
    else:
        log.info("Setting seasonSeed")
        seasonSeed = int(time.time()/1314000) # Days since epoch
        config.set('TrackGenerator', 'SeasonSeed', seasonSeed)

    if config.has_option('TrackGenerator', 'TrackSeed'):
        trackSeed = config.getint('TrackGenerator', 'TrackSeed')
    else:
        log.info("Setting trackSeed")
        trackSeed = int(time.time())  # Seconds since epoch
        config.set('TrackGenerator', 'TrackSeed', trackSeed)

    if config.has_option('TrackGenerator', 'YearsPerSimulation'):
        yrsPerSim = config.getint('TrackGenerator', 'YearsPerSimulation')
    else:
        yrsPerSim = 1

    # Attempt to start the track generator in parallel
    global MPI
    MPI = attemptParallel()
    comm = MPI.COMM_WORLD

    if comm.size > 1 and (not seasonSeed or not trackSeed):
        log.critical('TrackSeed and GenesisSeed are needed' +
                     ' for parallel runs!')
        sys.exit(1)

    mslp = SamplePressure(mslpFile, var=mslpVar)

    # Initialise the landfall tracking

    landfall = trackLandfall.LandfallDecay(configFile, dt)

    # Wait for configuration to be loaded by all processors

    comm.barrier()

    # Seed the numpy PRNG. We use this PRNG to sample the number of
    # tropical cyclone tracks to simulate for each season. The inbuilt
    # Python `random` library does not provide a function to sample
    # from the Poisson distribution.

    if seasonSeed:
        np.random.seed(seasonSeed)

    # Do the first stage of the simulation (i.e., sample the number of
    # tracks to simulate at each genesis point) on all processors
    # simultaneously. Since the same seed is set on all processors,
    # they will all get exactly the same simulation outcome. This also
    # behaves correctly when not done in parallel.

    nCyclones = np.random.poisson(yrsPerSim * meanFreq, nSimulations)

    # Estimate the maximum number of random values to be drawn from the
    # PRNG for each track and calculate how much each track simulation
    # should jump ahead in the PRNG stream to ensure that it is
    # independent of all other simulations.

    maxRvsPerTrack = 8 * (maxTimeSteps + 2)
    jumpAhead = np.hstack([[0],
                           np.cumsum(nCyclones * maxRvsPerTrack)[:-1]])

    log.info('Generating %i total events for %i simulations',
             sum(nCyclones), nSimulations)

    # Setup the simulation parameters

    sims = []
    for i, n in enumerate(nCyclones):
        sims.append(Simulation(i, trackSeed, jumpAhead[i], n,
                               trackFilename % i))

    # Load the track generator

    tg = TrackGenerator(processPath, gridLimit, gridSpace, gridInc,
                        mslp, landfall, dt=dt,
                        maxTimeSteps=maxTimeSteps)

    tg.loadInitialConditionDistributions()
    tg.loadCellStatistics()

    # Hold until all processors are ready

    comm.barrier()

    N = sims[-1].index

    # Balance the simulations over the number of processors and do it

    for sim in balanced(sims):
        log.debug('Simulating tropical cyclone tracks:' +
                  ' %3.0f percent complete' % (sim.index / float(N)
                                               * 100.))
        if callback is not None:
            callback(sim.index, N)

        if sim.seed:
            global PRNG #TODO: explicitly-pass rather than mutate global state
            PRNG = random.Random(seed=sim.seed, stream=sim.index)
            log.debug('seed %i stream %i', sim.seed, sim.index)

        trackFile = pjoin(trackPath, sim.outfile)
        tracks = tg.generateTracks(sim.ntracks, sim.index)
        ncTrackFile = pjoin(trackPath, "tracks.{0:05d}.nc".format(sim.index))
        ncSaveTracks(ncTrackFile, tracks, calendar='julian')

    log.info('Simulating tropical cyclone tracks:' +
             ' 100 percent complete')

    comm.barrier()
