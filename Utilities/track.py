"""
:mod:`track` - track-related attributes and functions
=====================================================

.. module:: tracks
    :synopsis: This module contains funcitons for reading/writing
               track data from/to csv and netCDF formats.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import time
import logging
import getpass
import numpy as np
from datetime import datetime

from os.path import join as pjoin
import re
from shapely.geometry import Point, LineString

from Utilities.metutils import convert
from Utilities.maputils import bearing2theta

from netCDF4 import Dataset, date2num, num2date
from cftime import num2pydate

try:
    from exceptions import WindowsError
except:
    class WindowsError(IOError): pass

#if not getattr(__builtins__, "WindowsError", None):
#    class WindowsError(IOError):
#        pass

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

trackFields = ('Indicator', 'CycloneNumber', 'Year', 'Month',
               'Day', 'Hour', 'Minute', 'TimeElapsed', 'Datetime', 'Longitude',
               'Latitude', 'Speed', 'Bearing', 'CentralPressure',
               'WindSpeed', 'rMax', 'EnvPressure')

trackTypes = ('i', 'i', 'i', 'i',
              'i', 'i', 'i', 'f', datetime,
              'f', 'f', 'f', 'f', 'f',
              'f', 'f', 'f')

trackFormats = ('%i, %i, %i, %i,'
                '%i, %i, %i, %5.1f,' '%s',
                '%8.3f, %8.3f, %6.2f, %6.2f, %7.2f,'
                '%6.2f, %6.2f, %7.2f')

PATTERN = re.compile(r'\d+')

class Track(object):
    """
    A single tropical cyclone track.

    The object exposes the track data through the object attributes.
    For example, If `data` contains the tropical cyclone track data
    (`numpy.array`) loaded with the :meth:`readTrackData` function,
    then the central pressure column can be printed out with the
    code::

        t = Track(data)
        print(t.CentralPressure)

    :type  data: numpy.ndarray
    :param data: the tropical cyclone track data.
    """

    def __init__(self, data):
        """
        :type  data: numpy.ndarray
        :param data: the tropical cyclone track data.
        """
        self.data = data
        self.trackId = None
        self.trackfile = None
        if (len(data) > 0) and ('CentralPressure' in data.dtype.names):
            self.trackMinPressure = np.min(data['CentralPressure'])
        else:
            self.trackMinPressure = None
        if (len(data) > 0) and ('WindSpeed' in data.dtype.names):
            self.trackMaxWind = np.max(data['WindSpeed'])
        else:
            self.trackMaxWind = None

    def __getattr__(self, key):
        """
        Get the `key` from the `data` object.

        :type  key: str
        :param key: the key to lookup in the `data` object.
        """
        if (key.startswith('__') and key.endswith('__')) or (key == 'data'):
            return super(Track, self).__getattr__(key)

        return self.data[key]

    def __repr__(self):
        return "<Track of dtype [{}]>".format(", ".join(self.data.dtype.names))

    def inRegion(self, gridLimit):
        """
        Check if the tropical cyclone track starts within a region.

        :type  gridLimit: :class:`dict`
        :param gridLimit: the region to check.
                          The :class:`dict` should contain the keys
                          :attr:`xMin`, :attr:`xMax`, :attr:`yMin` and
                          :attr:`yMax`. The *x* variable bounds the
                          latitude and the *y* variable bounds the
                          longitude.

        """
        xMin = gridLimit['xMin']
        xMax = gridLimit['xMax']
        yMin = gridLimit['yMin']
        yMax = gridLimit['yMax']

        return ((xMin <= self.Longitude[0]) and
                (self.Longitude[0] <= xMax) and
                (yMin <= self.Latitude[0]) and
                (self.Latitude[0] <= yMax))

    def minimumDistance(self, points):
        """
        Calculate the minimum distance between a track and a
        collection of :class:`shapely.geometry.Point` points. Assumes
        the points and the :attr:`Longitude` and :attr:`Latitude`
        attributes share the same coordinate system (presumed to be
        geographic coordinates).

        :param points: sequence of :class:`shapely.geometry.Point` objects.

        :returns: :class:`numpy.ndarray` of minimum distances between
                  the set of points and the line features (in km).
        """
        coords = [(x, y) for x, y in zip(self.Longitude, self.Latitude)]

        if len(coords) == 1:
            point_feature = Point(self.Longitude, self.Latitude)
            distances = [point_feature.distance(point) for point in points]
        else:
            line_feature = LineString(coords)
            distances = [line_feature.distance(point) for point in points]

        return convert(distances, 'deg', 'km')


# Define format for TCRM output track files:
ISO_FORMAT = "%Y-%m-%d %H:%M:%S"
TCRM_COLS = ('CycloneNumber', 'Datetime', 'TimeElapsed', 'Longitude',
             'Latitude', 'Speed', 'Bearing', 'CentralPressure',
             'EnvPressure', 'rMax')

TCRM_UNIT = ('', '', 'hr', 'degree', 'degree', 'kph', 'degrees',
                  'hPa', 'hPa', 'km')

TCRM_FMTS = ('i', 'object', 'f', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')

TCRM_CNVT = {
    0: lambda s: int(float(s.strip() or 0)),
    1: lambda s: datetime.strptime(s.strip(), ISO_FORMAT),
    5: lambda s: convert(float(s.strip() or 0), TCRM_UNIT[5], 'mps'),
    6: lambda s: bearing2theta(float(s.strip() or 0) * np.pi / 180.),
    7: lambda s: convert(float(s.strip() or 0), TCRM_UNIT[7], 'Pa'),
    8: lambda s: convert(float(s.strip() or 0), TCRM_UNIT[8], 'Pa'),
}

TRACK_DT_ERR = "Track data does not have required \
attributes to convert to datetime object"

TRACK_EMPTY_GROUP = """No track groups in this netcdf file: {0}"""

def ncReadTrackData(trackfile):
    """
    Read a netcdf-format track file into a collection of
    :class:`Track` objects. The returned :class:`Track` objects *must*
    have all attributes accessed by the `__getattr__` method.

    :param str trackfile: track data filename (netCDF4 format).

    :return: track data
    :rtype: list of :class:`Track` objects

    """

    track_dtype = np.dtype({'names':TCRM_COLS,
                            'formats':TCRM_FMTS})
    try:
        ncobj = Dataset(trackfile, mode='r')
    except (IOError, RuntimeError):
        log.exception("Cannot open {0}".format(trackfile))
        raise IOError("Cannot open {0}".format(trackfile))

    g = ncobj.groups
    if not bool(g):
        # We have a track file that stores data in separate variables
        log.debug(f"Reading data from a single track file")
        dt = ncobj.variables['Datetime']
        units = ncobj.getncattr('time_units')
        calendar = ncobj.getncattr('calendar')
        dtt = num2date(dt[:], units, calendar)
        # Convert to true python datetimes
        dtconversion = [datetime.strptime(d.strftime(), "%Y-%m-%d %H:%M:%S") for d in dtt]
        newtd = np.zeros(len(dtt), dtype=track_dtype)
        for f in ncobj.variables.keys():
            if f != 'Datetime' and f in track_dtype.names:
                newtd[f] = ncobj.variables[f][:]
        newtd['Datetime'] = dtconversion
        track = Track(newtd)
        track.trackfile = trackfile
        track.trackId = eval(ncobj.trackId)

        return [track]

    tracks = []
    if 'tracks' in g:
        tgroup = g['tracks'].groups
        ntracks = len(tgroup)
        for i, (t, data) in enumerate(tgroup.items()):
            log.debug("Loading data for {0}".format(t))
            track_data = data.variables['track'][:]

            try: 
                dt = num2date(track_data['Datetime'],
                              data.variables['time'].units,
                              data.variables['time'].calendar)
            except AttributeError:
                log.exception(TRACK_DT_ERR)
                raise AttributeError

            newtd = np.zeros(len(track_data), dtype=track_dtype)
            for f in track_data.dtype.names:
                if f != 'Datetime' and f in track_dtype.names:
                    newtd[f] = track_data[f]
            dtconversion = [datetime.strptime(d.strftime(), "%Y-%m-%d %H:%M:%S") for d in dt]
            newtd['Datetime'] = dtconversion

            track = Track(newtd)
            track.trackfile = trackfile
            if hasattr(data.variables['track'], "trackId"):
                track.trackId = eval(data.variables['track'].trackId)
            else:
                track.trackId = (i+1, ntracks)
            tracks.append(track)

    else:
        log.warn(TRACK_EMPTY_GROUP.format(trackfile))

    ncobj.close()
    return tracks

def ncSaveTracks(trackfile, tracks,
                 timeunits='hours since 1900-01-01 00:00',
                 calendar='standard', attributes={}):
    """
    Save a collection of :class:`Track` objects to a netCDF file. This
    makes use of netCDF4 compound data types to store the data as a
    structure akin to a :class:`numpy.recarray`. Each track in the
    collection is stored as a separate :class:`netCDF4.Group` instance
    in the output file.

    The :class:`Track` objects hold datetime information as an array
    of :class:`datetime.datetime` objects - there is no equivalent
    data type in netCDF4, so the datetime information is converted to
    floats using the :class:`netCDF4.date2num` function.

    :param str trackfile: Path to the file to save data to.
    :param list tracks: Collection of :class:`Track` objects.
    :param str timeunits: A string of the form '*time units* since
                          *reference time*' describing the time units.
                          Default is 'hours since 1900-01-01 00:00'.
    :param str calendar: Calendar used for time calculations. Valid calendars
                         are 'standard', 'gregorian', 'proleptic_gregorian',
                         'noleap', '365_day', '360_day', 'julian', 'all_leap',
                         '366_day'. Default is 'standard', which is a mixed
                         Julian/Gregorian calendar.
    :param dict attributes: Global attributes to add to the file.

    """

    if len(tracks) == 0:
        log.info("No tracks to be stored in track file: {0}".format(trackfile))
        return

    try:
        ncobj = Dataset(trackfile, "w", format="NETCDF4", clobber=True)
    except IOError:
        log.exception("Cannot open {0} for writing".format(trackfile))
        raise IOError("Cannot open {0} for writing".format(trackfile))

    tgroup = ncobj.createGroup('tracks')

    # Fidget with the dtype to convert :class:`datetime` objects to floats:
    track_dtype = np.dtype(tracks[0].data.dtype)
    dtidx = track_dtype.names.index('Datetime')
    track_dtype = track_dtype.descr
    track_dtype[dtidx] = ('Datetime', 'f8')
    track_dtype = np.dtype(track_dtype)

    for n, t in enumerate(tracks):
        if len(t.data) == 0: # Empty track
            continue
        tname = "tracks-{:04d}".format(n)
        tdata = tgroup.createGroup(tname)
        tdtype = tdata.createCompoundType(track_dtype, 'track_dtype')

        dims = tdata.createDimension('time', None)
        times = tdata.createVariable('time', 'f8', ('time',),
                                     zlib=True, complevel=8, shuffle=True)
        tvar = tdata.createVariable('track', tdtype, ('time',),
                                    zlib=True, complevel=8, shuffle=True)
        t.data['Datetime'] = date2num(t.data['Datetime'], timeunits, calendar)
        times[:] = t.data['Datetime']
        times.units = 'hours since 1900-01-01 00:00'
        times.calendar = calendar
        tvar[:] = t.data.astype(track_dtype)
        tvar.long_name = "Tropical cyclone track data"
        tvar.time_units = 'hours since 1900-01-01 00:00'
        tvar.calendar = calendar
        tvar.lon_units = 'degrees east'
        tvar.lat_units = 'degrees north'
        tvar.pressure_units = 'hPa'
        tvar.speed_units = 'km/h'
        tvar.length_units = 'km'
        tvar.trackId = repr(t.trackId)

    attributes['created_on'] = time.strftime(ISO_FORMAT, time.localtime())
    attributes['created_by'] = getpass.getuser()
    ncobj.setncatts(attributes)
    ncobj.close()

    return

def readTrackData(trackfile):
    """
    Read a track .csv file into a numpy.ndarray.

    The track format and converters are specified with the global variables

        TRACKFILE_COLS -- The column names
        TRACKFILE_FMTS -- The entry formats
        TRACKFILE_CNVT -- The column converters

    :param str trackfile: the track data filename.

    :return: track data
    :rtype: :class:`numpy.ndarray`

    """

    try:
        return np.loadtxt(trackfile,
                          comments='%',
                          delimiter=',',
                          dtype={
                          'names': TCRM_COLS,
                          'formats': TCRM_FMTS},
                          converters=TCRM_CNVT)
    except ValueError:
        # return an empty array with the appropriate `dtype` field names
        return np.empty(0, dtype={
                        'names': TCRM_COLS,
                        'formats': TCRM_FMTS})

def readMultipleTrackData(trackfile):
    """
    Reads all the track datas from a .csv file into a list of numpy.ndarrays.
    The tracks are seperated based in their cyclone id. This function calls
    `readTrackData` to read the data from the file.

    :param str trackfile: the track data filename.

    :return: a collection of :class:`Track` objects

    """

    datas = []
    data = readTrackData(trackfile)
    if len(data) > 0:
        cycloneId = data['CycloneNumber']
        for i in range(1, np.max(cycloneId) + 1):
            datas.append(data[cycloneId == i])
    else:
        datas.append(data)
    return datas

def loadTracks(trackfile):
    """
    Read tracks from a track .nc file and return a list of :class:`Track`
    objects.

    This calls the function `ncReadTrackData` to parse the track .nc file.

    :param str trackfile: the track data filename.

    :return: list of :class:`Track` objects.

    """
    if not isinstance(trackfile, str):
        raise TypeError("Track file name is not a string: {0}".\
                        format(trackfile))
    if os.path.exists(trackfile):
        tracks = ncReadTrackData(trackfile)
        return tracks
    else:
        raise IOError("Track file doesn't exist: {0}".format(trackfile))

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

    :raises: TypeError if input argument is not a list.
    """
    if not isinstance(trackfiles, list):
        raise TypeError("Input argument is not a list")

    for f in trackfiles:
        msg = "Loading tracks in {0}".format(f)
        log.debug(msg)
        tracks = loadTracks(f)
        for track in tracks:
            yield track
    
def loadTracksFromPath(path):
    """
    Helper function to obtain a generator that yields :class:`Track` objects
    from a directory containing track .csv files.

    This function calls `loadTracksFromFiles` to obtain the generator and track
    filenames are processed in alphabetical order.

    :type  path: str
    :param path: the directory path.

    :raises: IOError if the path does not exist.
    """
    try: 
        files = os.listdir(path)
        trackfiles = [pjoin(path, f) for f in files if f.startswith('tracks')]
        msg = "Loading {0} track files in {1}".format(len(trackfiles), path)
        log.info(msg)
        return loadTracksFromFiles(sorted(trackfiles))
    except (IOError, OSError, WindowsError):
        raise IOError("Path {0} does not exist".format(path))
