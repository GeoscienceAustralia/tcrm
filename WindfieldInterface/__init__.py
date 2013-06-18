"""
:mod:`WindfieldInterface` -- Wind field calculation
===================================================

This module contains the core object for the wind field calculations.

Wind field calculations can be run in parallel using MPI if the :term:`pypar`
library is found and TCRM is run using the :term:`mpirun` command. For example,
to run with 10 processors::

    mpirun -n 10 python main.py cairns.ini

:class:`WindfieldInterface` can be correctly initialised and started by
calling the :meth: `run` with the location of a *configFile*::

    import WindfieldInterface
    WindfieldInterface.run('cairns.ini')

"""

import numpy as np
import logging as log
import itertools
import math
import time
import os
import sys
import Utilities.metutils as metutils
import Utilities.maputils as maputils
import Utilities.nctools as nctools

from os.path import join as pjoin, dirname
from collections import defaultdict
from config import ConfigParser
from PressureInterface.pressureProfile import PrsProfile as PressureProfile
from windProfile import WindProfile
from windVorticity import WindVorticity
from windField import WindField
from Utilities.config import cnfGetIniValue
from Utilities.progressbar import ProgressBar
from Utilities.files import flModDate, flProgramVersion, flGetStat


"""
Trackfile .csv format.
"""

TRACKFILE_COLS = ('CycloneNumber', 'TimeElapsed', 'Longitude', 'Latitude',
                  'Speed', 'Bearing', 'CentralPressure', 'EnvPressure', 'rMax')
TRACKFILE_UNIT = ('', 'hr', 'degree',
                  'degree', 'kph', 'degrees', 'hPa', 'hPa', 'km')
TRACKFILE_FMTS = ('i', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f')
TRACKFILE_CNVT = {
    0: lambda s: int(float(s.strip() or 0)),
    4: lambda s: metutils.convert(float(s.strip() or 0), TRACKFILE_UNIT[4], 'mps'),
    5: lambda s: maputils.bearing2theta(float(s.strip() or 0)*np.pi/180.),
    6: lambda s: metutils.convert(float(s.strip() or 0), TRACKFILE_UNIT[6], 'Pa'),
    7: lambda s: metutils.convert(float(s.strip() or 0), TRACKFILE_UNIT[7], 'Pa'),
}


class Track(object):
    """
    A single tropical cyclone track.

    The object exposes the track data through the object attributes. For
    example, If `data` contains the tropical cyclone track data (`numpy.array`)
    loaded with the :meth:`readTrackData` function, then the central pressure
    column can be printed out with the code::

        t = Track(data)
        print(t.CentralPressure)


    :type  data: numpy.ndarray
    :param data: the tropical cyclone track data.
    """

    def __init__(self, data):
        self.data = data

    def __getattr__(self, key):
        """
        Get the `key` from the `data` object.

        :type  key: str
        :param key: the key to lookup in the `data` object.
        """
        if key.startswith('__') and key.endswith('__'):
            return super(Track, self).__getattr__(key)
        return self.data[key]

    def inRegion(self, gridLimit):
        """
        Check if the tropical cyclone track falls within a region.
    
        :type  gridLimit: :class:`dict`
        :param gridLimit: the region to check.
                          The :class:`dict` should contain the keys :attr:`xMin`,
                          :attr:`xMax`, :attr:`yMin` and :attr:`yMax`. The *x*
                          variable bounds the latitude and the *y* variable bounds
                          the longitude.
        """
        xMin = gridLimit['xMin']
        xMax = gridLimit['xMax']
        yMin = gridLimit['yMin']
        yMax = gridLimit['yMax']

        return (xMin <= np.min(self.Longitude)) and \
               (np.max(self.Latitude) <= xMax) and \
               (yMin <= np.min(self.Latitude)) and \
               (np.max(self.Latitude) <= yMax)



class WindfieldAroundTrack(object):
    """
    The windfield around the tropical cyclone track.


    :type  track: :class:`Track`
    :param track: the tropical cyclone track.

    :type  profileType: str
    :param profileType: the wind profile type.

    :type  windFieldType: str
    :param windFieldType: the wind field type.

    :type  beta: float
    :param beta: wind field parameter.

    :type  beta1: float
    :param beta1: wind field parameter.

    :type  beta2: float
    :param beta2: wind field parameter.

    :type  thetaMax: float
    :param thetaMax: 

    :type  margin: float
    :param margin: 

    :type  resolution: float
    :param resolution: 

    :type  gustFactor: float
    :param gustFactor: 

    :type  gridLimit: :class:`dict`
    :param gridLimit: the domain where the tracks will be generated.
                      The :class:`dict` should contain the keys :attr:`xMin`,
                      :attr:`xMax`, :attr:`yMin` and :attr:`yMax`. The *x*
                      variable bounds the latitude and the *y* variable bounds
                      the longitude.

    """

    def __init__(self, track, profileType='powell', windFieldType='kepert',
                 beta=1.5, beta1=1.5, beta2=1.4, thetaMax=70.0, margin=2.0,
                 resolution=0.05, gustFactor=1.38, gridLimit=None):
        self.track = track
        self.profileType = profileType
        self.windFieldType = windFieldType
        self.beta = beta
        self.beta1 = beta1
        self.beta2 = beta2
        self.thetaMax = math.radians(thetaMax)
        self.margin = margin
        self.resolution = resolution
        self.gustFactor = gustFactor
        self.gridLimit = gridLimit

    def polarGridAroundEye(self, i):
        """
        Generate a polar coordinate grid around the eye of the
        tropical cyclone at time i.

        :type  i: int
        :param i: the time.
        """
        R, theta = maputils.makeGrid(self.track.Longitude[i],
                                     self.track.Latitude[i],
                                     self.margin, self.resolution)
        return R, theta

    def windProfile(self, i, R):
        """
        Calculate the wind profile at time i at the radiuses R
        around the tropical cyclone.
        """
        p = WindProfile(R, self.track.EnvPressure[i],
                        self.track.CentralPressure[i],
                        self.track.rMax[i],
                        self.track.Latitude[i],
                        self.track.Longitude[i],
                        self.beta, beta1=self.beta1,
                        beta2=self.beta2)
        try:
            wind = getattr(p, self.profileType)
        except AttributeError:
            log.exception('%s not implemented in windProfile'
                              % self.profileType)
        V = wind()

        return V

    def windVorticity(self, i, R, maxAbsWindSpeed):
        """
        Calculate the wind vorticity at time `i` at the radiuses `R`
        around the tropical cyclone that has maximum absolute wind
        speed given by `maxAbsWindSpeed`.


        :type  i: int
        :param i: the time.

        :type  R: float
        :param R: the radius around the tropical cyclone.

        :type  maxAbsWindSpeed: float
        :param maxAbsWindSpeed: the maximum absolute wind speed.
        """
        v = WindVorticity(R, self.track.EnvPressure[i],
                          self.track.CentralPressure[i],
                          self.track.rMax[i],
                          self.track.Latitude[i],
                          self.track.Longitude[i],
                          self.beta, beta1=self.beta1,
                          beta2=self.beta2)
        try:
            vort = getattr(v, self.profileType)
        except AttributeError:
            log.exception('%s not implemented in windVorticity'
                              % self.profileType)
        Z = vort(maxAbsWindSpeed)
        return Z

    def pressureProfile(self, i, R):
        """
        Calculate the pressure profile at time `i` at the radiuses `R`
        around the tropical cyclone.


        :type  i: int
        :param i: the time.

        :type  R: :class:`numpy.ndarray`
        :param R: the radiuses around the tropical cyclone.
        """
        p = PressureProfile(R, self.track.EnvPressure[i],
                            self.track.CentralPressure[i],
                            self.track.rMax[i],
                            self.track.Latitude[i],
                            self.track.Longitude[i],
                            self.beta, beta1=self.beta1,
                            beta2=self.beta2)
        try:
            pressure = getattr(p, self.profileType)
        except AttributeError:
            log.exception('%s not implemented in pressureProfile'
                              % self.profileType)
        return pressure()

    def localWindField(self, i):
        """
        Calculate the local wind field at time `i` around the
        tropical cyclone.

        :type  i: int
        :param i: the time.
        """
        R, theta = self.polarGridAroundEye(i)
        V = self.windProfile(i, R)
        Z = self.windVorticity(i, R, np.absolute(V).max())
        P = self.pressureProfile(i, R)
        f = metutils.coriolis(self.track.Latitude[i])

        wf = WindField(R, theta, self.track.rMax[i], f, V, Z,
                       self.track.Speed[i],
                       self.track.Bearing[i],
                       self.thetaMax)
        try:
            field = getattr(wf, self.windFieldType)
        except AttributeError:
            log.exception('%s not implemented in windField'
                              % self.windFieldType)

        Ux, Vy = field()
        return (Ux, Vy, P)

    def regionalExtremes(self, gridLimit, timeStepCallback=None):
        """
        Calculate the maximum potential wind gust and minimum
        pressure over the region throughout the life of the
        tropical cyclone.


        :type  gridLimit: :class:`dict`
        :param gridLimit: the domain where the tracks will be considered.
                          The :class:`dict` should contain the keys 
                          :attr:`xMin`, :attr:`xMax`, :attr:`yMin` and 
                          :attr:`yMax`. The *x* variable bounds the 
                          latitude and the *y* variable bounds the longitude.

        :type  timeStepCallback: function
        :param timeStepCallback: the function to be called on each time step.
        """

        if len(self.track.data) > 0:
            envPressure = self.track.EnvPressure[0]
        else:
            envPressure = np.NaN

        # Get the limits of the region

        xMin = gridLimit['xMin']
        xMax = gridLimit['xMax']
        yMin = gridLimit['yMin']
        yMax = gridLimit['yMax']

        # Setup a 'centidegree' integer grid for the region

        gridMargin = int(100.*self.margin)
        gridStep = int(100.*self.resolution)

        minLat = int(100.*yMin)-gridMargin
        maxLat = int(100.*yMax)+gridMargin
        minLon = int(100.*xMin)-gridMargin
        maxLon = int(100.*xMax)+gridMargin

        latGrid = np.arange(minLat, maxLat+gridStep, gridStep, dtype=int)
        lonGrid = np.arange(minLon, maxLon+gridStep, gridStep, dtype=int)

        [cGridX, cGridY] = np.meshgrid(lonGrid, latGrid)

        # Initialise the region

        UU = np.zeros_like(cGridX, dtype='f')
        VV = np.zeros_like(cGridY, dtype='f')
        bearing = np.zeros_like(cGridX, dtype='f')
        gust = np.zeros_like(cGridX, dtype='f')
        pressure = np.ones_like(cGridX, dtype='f')*envPressure

        lonCDegree = np.array(100.*self.track.Longitude, dtype=int)
        latCDegree = np.array(100.*self.track.Latitude, dtype=int)

        # We only consider the times when the TC track falls in the region

        timesInRegion = np.where((xMin <= self.track.Longitude) &
                                (self.track.Longitude <= xMax) &
                                (yMin <= self.track.Latitude) &
                                (self.track.Latitude <= yMax))[0]

        for i in timesInRegion:

            # Map the local grid to the regional grid

            jmin = int((latCDegree[i]-minLat-gridMargin)/gridStep)
            jmax = int((latCDegree[i]-minLat+gridMargin)/gridStep)+1
            imin = int((lonCDegree[i]-minLon-gridMargin)/gridStep)
            imax = int((lonCDegree[i]-minLon+gridMargin)/gridStep)+1

            # Calculate the local wind speeds and pressure at time i

            Ux, Vy, P = self.localWindField(i)

            # Calculate the local wind gust and bearing

            Ux *= self.gustFactor
            Vy *= self.gustFactor

            localGust = np.sqrt(Ux**2 + Vy**2)
            localBearing = ((np.arctan2(-Ux, -Vy))*180./np.pi)

            # Handover this time step to a callback if required

            if timeStepCallback:
                timeStepCallback(i, localGust, Ux, Vy, P)

            # Retain when there is a new maximum gust

            mask = localGust > gust[jmin:jmax, imin:imax]

            gust[jmin:jmax, imin:imax] = np.where(
                mask, localGust, gust[jmin:jmax, imin:imax])
            bearing[jmin:jmax, imin:imax] = np.where(
                mask, localBearing, bearing[jmin:jmax, imin:imax])
            UU[jmin:jmax, imin:imax] = np.where(
                mask, Ux, UU[jmin:jmax, imin:imax])
            VV[jmin:jmax, imin:imax] = np.where(
                mask, Vy, VV[jmin:jmax, imin:imax])

            # Retain the lowest pressure

            pressure[jmin:jmax, imin:imax] = np.where(
                P < pressure[jmin:jmax, imin:imax],
                P, pressure[jmin:jmax, imin:imax])

        return gust, bearing, UU, VV, pressure, lonGrid/100., latGrid/100.


class WindfieldGenerator(object):
    """
    The wind field generator.


    :type  margin: float
    :param margin: 

    :type  resolution: float
    :param resolution: 

    :type  profileType: str
    :param profileType: the wind profile type.

    :type  windFieldType: str
    :param windFieldType: the wind field type.

    :type  beta: float
    :param beta: wind field parameter.

    :type  beta1: float
    :param beta1: wind field parameter.

    :type  beta2: float
    :param beta2: wind field parameter.

    :type  thetaMax: float
    :param thetaMax: 

    :type  gridLimit: :class:`dict`
    :param gridLimit: the domain where the tracks will be generated.
                      The :class:`dict` should contain the keys :attr:`xMin`,
                      :attr:`xMax`, :attr:`yMin` and :attr:`yMax`. The *x*
                      variable bounds the latitude and the *y* variable bounds
                      the longitude.

    """

    def __init__(self, margin=2.0, resolution=0.05, profileType='powell', 
            windFieldType='kepert', beta=1.5, beta1=1.5, beta2=1.4, thetaMax=70.0, 
            gridLimit=None):
        self.margin = margin
        self.resolution = resolution
        self.profileType = profileType
        self.windFieldType = windFieldType
        self.beta = beta
        self.beta1 = beta1
        self.beta2 = beta2
        self.thetaMax = thetaMax
        self.gridLimit = gridLimit

    def calculateExtremesFromTrack(self, track):
        """
        Calculate the wind extremes given a single tropical cyclone track.


        :type  track: :class:`Track`
        :param track: the tropical cyclone track.
        """
        wt = WindfieldAroundTrack(track, 
                profileType   = self.profileType,
                windFieldType = self.windFieldType,
                beta          = self.beta,
                beta1         = self.beta1,
                beta2         = self.beta2,
                thetaMax      = self.thetaMax,
                margin        = self.margin,
                resolution    = self.resolution)
        return track, wt.regionalExtremes(self.gridLimit)

    def calculateExtremesFromTrackfile(self, trackfile):
        """
        Calculate the wind extremes from a `trackfile` that might contain a
        number of tropical cyclone tracks. The wind extremes are calculated
        over the tracks, i.e., the maximum gusts and minimum pressures over all
        tracks are retained.


        :type  trackfile: str
        :param trackfile: the file name of the trackfile.
        """
        trackiter = loadTracks(trackfile)
        f = lambda track: self.calculateExtremesFromTrack(track)[1]

        results = map(f, trackiter)

        gust, bearing, Vx, Vy, P, lon, lat = results.next()

        for result in results:
            gust1, bearing1, Vx1, Vy1, P1, lon1, lat1 = result
            gust = np.where(gust1 > gust, gust1, gust)
            P = np.where(P1 < P, P1, P)

        return (gust, bearing, Vx, Vy, P, lon, lat)

    def dumpExtremesFromTrackfile(self, trackfile, dumpfile):
        """
        Helper method to calculate the wind extremes from a `trackfile` and
        save them to a file called `dumpfile`.


        :type  trackfile: str
        :param trackfile: the file name of the trackfile.

        :type  dumpfile: str
        :param dumpfile: the file name where to save the wind extremes.
        """
        result = self.calculateExtremesFromTrackfile(trackfile)
        self._saveExtremesToFile(result, dumpfile)

    def plotExtremesFromTrackfile(self, trackfile, windfieldfile, pressurefile):
        """
        Helper method to calculate the wind extremes from a `trackfile` and
        generate image files for the wind field and the pressure.


        :type  trackfile: str
        :param trackfile: the file name of the trackfile.

        :type  windfieldfile: str
        :param windfieldfile: the file name of the windfield image file to
                              write.

        :type  pressurefile: str
        :param pressurefile: the file name of the pressure image file to write.
        """
        result = self.calculateExtremesFromTrackfile(trackfile)
        self._plotExtremesToFile(result, winfieldfile, pressurefile)

    def _plotExtremesToFile(self, result, windfieldfile, pressurefile):
        """
        
        """
        from PlotInterface.plotWindfield import plotWindfield, plotPressurefield
        gust, bearing, Vx, Vy, P, lon, lat = result
        [gridX, gridY] = meshgrid(lon, lat)
        plotWindfield(gridX, gridY, speed, title="Windfield", 
                fileName=windfieldfile)
        plotPressurefield(gridX, gridY, pressure, title="Pressure field", 
                fileName=pressurefile)

    def dumpGustsFromTracks(self, trackiter, windfieldPath, filenameFormat, 
                            progressCallback=None):
        """
        Dump the maximum wind speeds (gusts) observed over a region to 
        netcdf files. One file is created for every track file.
        """
        results = itertools.imap(self.calculateExtremesFromTrack, trackiter)

        gusts = {}
        done = defaultdict(list)

        i = 0
        for track, result in results:
            gust, bearing, Vx, Vy, P, lon, lat = result

            if track.trackfile in gusts:
                gust1, bearing1, Vx1, Vy1, P1, lon1, lat1 = gusts[track.trackfile]
                gust = np.where(gust > gust1, gust, gust1)

            gusts[track.trackfile] = (gust, bearing, Vx, Vy, P, lon, lat)
            done[track.trackfile] += [track.id]

            if len(done[track.trackfile]) >= done[track.trackfile][0][1]:
                dumpfile = pjoin(windfieldPath, filenameFormat % i)
                self._saveGustToFile(track.trackfile, (lat, lon, gust), dumpfile)
                del done[track.trackfile]
                del gusts[track.trackfile]
                i += 1
                if progressCallback:
                    progressCallback((track.trackfile, dumpfile))

    def _saveGustToFile(self, trackfile, result, filename):
        """

        """
        lat, lon, speed = result

        inputFileDate = flModDate(trackfile)

        gatts = {'history'        : 'TCRM hazard simulation - synthetic event wind field',
                 'version'        : flProgramVersion(),
                 'Python_ver'     : sys.version,
                 'track_file'     : '%s (modified %s)' % (trackfile, inputFileDate),
                 'radial_profile' : self.profileType,
                 'boundary_layer' : self.windFieldType,
                 'beta'           : self.beta}

        dimensions = {0: {'name'  : 'lat', 
                          'values': lat, 
                          'dtype' : 'f',
                          'atts'  : {'long_name': 'Latitude',
                                     'units': 'degrees_north'}},
                      1: {'name'  : 'lon', 
                          'values': lon, 
                          'dtype' : 'f',
                          'atts'  : {'long_name': 'Longitude',
                                   'units': 'degrees_east'}}}

        variables =  {0: {'name'  : 'vmax', 
                          'dims'  : ('lat', 'lon'),
                          'values': speed, 
                          'dtype' : 'f',
                          'atts'  : {'long_name': 'Maximum 3-second gust wind speed', 
                                     'units': 'm/s'}}}

        nctools.ncSaveGrid(filename, dimensions, variables, gatts=gatts)

    def dumpGustsFromTrackfiles(self, trackfiles, windfieldPath, 
                                filenameFormat='gust-%04i.nc', 
                                progressCallback=None):
        """
        Helper method to dump the maximum wind speeds (gusts) observed over a
        region to netcdf files. One file is created for every track file.

        :type  trackfiles: list of str
        :param trackfiles: a list of track file filenames.

        :type  windfieldPath: str
        :param windfieldPath: the path where to store the gust output files.

        :type  filenameFormat: str
        :param filenameFormat: the format string for the output file names. The
                               default is set to 'gust-%04i.nc'.

        :type  progressCallback: function
        :param progressCallback: optional function to be called after a file is 
                                 saved. This can be used to track progress.
        """
        tracks = loadTracksFromFiles(sorted(trackfiles))
        self.dumpGustsFromTracks(tracks, windfieldPath, filenameFormat, 
                                 progressCallback=progressCallback)

    def _saveExtremesToFile(self, result, filename):
        gust, bearing, Vx, Vy, P, lon, lat = result

        dimensions = {0:{'name':'lat', 'values':lat, 'dtype':'f',
                             'atts':{'long_name':'Latitude',
                                         'units':'degrees_north'} },
                      1:{'name':'lon', 'values':lon, 'dtype':'f',
                             'atts':{'long_name':'Longitude',
                                         'units':'degrees_east'} } }

        variables = {0:{'name':'vmax', 'dims':('lat', 'lon'),
                        'values':gust, 'dtype':'f',
                        'atts':{'long_name':'Maximum 3-second gust wind speed',
                                'units':'m/s'} },
                     1:{'name':'ua', 'dims':('lat', 'lon'),
                        'values':Vx, 'dtype':'f',
                        'atts':{'long_name':'Maximum eastward wind',
                                'units':'m/s'} },
                     2:{'name':'va', 'dims':('lat', 'lon'),
                        'values':Vy, 'dtype':'f',
                        'atts':{'long_name':'Maximum northward wind',
                                'units':'m/s'} },
                     3:{'name':'psl', 'dims':('lat', 'lon'),
                        'values':P, 'dtype':'f',
                        'atts':{'long_name':'Minimum air pressure at sea level',
                                'units':'hPa'} } }

        nctools.ncSaveGrid(filename, dimensions, variables)


def readTrackData(trackfile):
    """
    Read a track .csv file into a numpy.ndarray.

    The track format and converters are specified with the global variables

        TRACKFILE_COLS -- The column names
        TRACKFILE_FMTS -- The entry formats
        TRACKFILE_CNVT -- The column converters

    :type  trackfile: str
    :param trackfile: the track data filename.
    """
    try:
        return np.loadtxt(trackfile,
                          comments='%',
                          delimiter=',',
                          dtype={
                          'names': TRACKFILE_COLS,
                          'formats': TRACKFILE_FMTS},
                          converters=TRACKFILE_CNVT)
    except:
        # return an empty array with the appropriate `dtype` field names
        return np.empty(0, dtype={
                        'names': TRACKFILE_COLS,
                        'formats': TRACKFILE_FMTS})


def readMultipleTrackData(trackfile):
    """
    Reads all the track datas from a .csv file into a list of numpy.ndarrays.
    The tracks are seperated based in their cyclone id. This function calls
    `readTrackData` to read the data from the file.
    
    :type  trackfile: str
    :param trackfile: the track data filename.
    """
    datas = []
    data = readTrackData(trackfile)
    if len(data) > 0:
        cycloneId = data['CycloneNumber']
        for i in range(1, np.max(cycloneId)+1):
            datas.append(data[cycloneId == i])
    else:
        datas.append(data)
    return datas


def loadTracksFromFiles(trackfiles):
    """
    Generator that yields :class:`Track` objects from a list of track
    filenames. 
    
    When run in parallel, the list `trackfiles` is distributed across the MPI
    processors using the `balanced` function. Track files are loaded in a lazy
    fashion to reduce memory consumption. The generator returns individual
    tracks (recall: a trackfile can contain multiple tracks) and only moves on
    to the next file once all the tracks from the current file have been
    returned.

    :type  trackfiles: list of strings
    :param trackfiles: list of track filenames. The filenames must include the
                       path to the file.
    """
    for f in balanced(trackfiles):
        log.info('Calculating wind fields for tracks in %s' % f)
        tracks = loadTracks(f)
        for track in tracks:
            yield track


def loadTracks(trackfile):
    """
    Read tracks from a track .csv file and return a list of :class:`Track`
    objects.

    This calls the function `readMultipleTrackData` to parse the track .csv
    file.

    :type  trackfile: str
    :param trackfile: the track data filename.
    """
    tracks = []
    datas = readMultipleTrackData(trackfile)
    n = len(datas)
    for i, data in enumerate(datas):
        track = Track(data)
        track.trackfile = trackfile
        track.id = (i, n)
        tracks.append(track)
    return tracks


def loadTracksFromPath(path):
    """
    Helper function to obtain a generator that yields :class:`Track` objects
    from a directory containing track .csv files.

    This function calls `loadTracksFromFiles` to obtain the generator and track
    filenames are processed in alphabetical order.

    :type  path: str
    :param path: the directory path.
    """
    files = os.listdir(path)
    trackfiles = [pjoin(path, f) for f in files if f.startswith('tracks')]
    log.info('Processing %d track files in %s' % (len(trackfiles), path))
    return loadTracksFromFiles(sorted(trackfiles))


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
            def size(self): return 1
            def rank(self): return 0
            def barrier(self): pass

        pp = DummyPypar()


def run(configFile):
    """
    Run the wind field calculations.
    """

    log.info('Loading wind field calculation settings')

    # Get configuration

    config = ConfigParser()
    config.read(configFile)

    outputPath      = config.get('Output', 'Path')
    profileType     = config.get('WindfieldInterface', 'profileType')
    windFieldType   = config.get('WindfieldInterface', 'windFieldType')
    beta            = config.getfloat('WindfieldInterface', 'beta')
    beta1           = config.getfloat('WindfieldInterface', 'beta1')
    beta2           = config.getfloat('WindfieldInterface', 'beta2')
    thetaMax        = config.getfloat('WindfieldInterface', 'thetaMax')
    margin          = config.getfloat('WindfieldInterface', 'Margin')
    resolution      = config.getfloat('WindfieldInterface', 'Resolution')
    gridLimit       = config.geteval('Region', 'gridLimit')
    plotOutput      = config.getboolean('WindfieldInterface', 'PlotOutput')
    showProgessBar  = config.getboolean('Logging', 'ProgressBar')
    windfieldPath   = pjoin(outputPath, 'windfield')
    trackPath       = pjoin(outputPath, 'tracks')
    windfieldFormat = 'gust-%04d.nc'

    if config.has_option('WindfieldInterface', 'gridLimit'):
        gridLimit = config.geteval('WindfieldInterface', 'gridLimit')

    if config.has_option('WindfieldInterface', 'TrackPath'):
        trackPath = config.get('WindfieldInterface', 'TrackPath') 

    thetaMax = math.radians(thetaMax)

    # Attempt to start the track generator in parallel
    
    attemptParallel()

    log.info('Running windfield generator')

    wfg = WindfieldGenerator(margin = margin,
                      resolution    = resolution,
                      profileType   = profileType,
                      windFieldType = windFieldType,
                      beta          = beta,
                      beta1         = beta1,
                      beta2         = beta2,
                      thetaMax      = thetaMax,
                      gridLimit     = gridLimit)

    log.info('Dumping gusts to %s' % windfieldPath)

    # Get the trackfile names and count

    files = os.listdir(trackPath)
    trackfiles = [pjoin(trackPath, f) for f in files if f.startswith('tracks')]
    nfiles = len(trackfiles)

    log.info('Processing %d track files in %s' % (nfiles, trackPath))

    # Do the work
    
    pp.barrier()

    wfg.dumpGustsFromTrackfiles(trackfiles, windfieldPath, windfieldFormat)

    pp.barrier()

    log.info('Completed windfield generator')

if __name__ == '__main__':
    try:
        configFile = sys.argv[1]
    except IndexError:

        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'

        # If no filename is specified and default filename doesn't exist =>
        # raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg

    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile + "' not found"
        raise IOError, error_msg

    run(configFile)

    # Finalise MPI
    
    pp.finalize()
