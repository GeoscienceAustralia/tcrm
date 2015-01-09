"""
Track-related attributes
"""

import os
import logging
import numpy as np
from datetime import datetime
import re
from shapely.geometry import Point, LineString

from Utilities.metutils import convert
from Utilities.maputils import bearing2theta

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
                
pattern = re.compile('\d+')

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
        self.trackMinPressure = None
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

    def inRegion(self, gridLimit):
        """
        Check if the tropical cyclone track falls within a region.

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

        return ((xMin <= np.min(self.Longitude)) and
                (np.max(self.Latitude) <= xMax) and
                (yMin <= np.min(self.Latitude)) and
                (np.max(self.Latitude) <= yMax))
                
    def minimumDistance(self, points):
        """
        Calculate the minimum distance between a track and a
        collection of :class:`shapely.geometry.Point` points. Assumes
        the points and the :attr:`Longitude` and :attr:`Latitude`
        attributes share the same coordinate system.

        :param points: sequence of :class:`shapely.geometry.Point` objects.

        :returns: :class:`numpy.ndarray` of minimum distances between
                  the set of points and the line features (in km).
        """
        coords = [(x,y) for x, y in zip(self.Longitude, self.Latitude)]
        
        if len(coords) == 1:
            point_feature = Point(self.Longitude, self.Latitude)
            distances = [point_feature.distance(point) for point in points]
        else:
            line_feature = LineString(coords)
            distances = [line_feature.distance(point) for point in points]

        return convert(distances, 'deg', 'km')

    
# Define format for TCRM output track files:
DATEFORMAT = "%Y-%m-%d %H:%M:%S"
TCRM_COLS = ('CycloneNumber', 'Datetime', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'EnvPressure', 'rMax')

TCRM_UNIT = ('', '', 'hr', 'degree', 'degree', 'kph', 'degrees',
                  'hPa', 'hPa', 'km')

TCRM_FMTS = ('i', 'object', 'f', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')

TCRM_CNVT = {
    0: lambda s: int(float(s.strip() or 0)),
    1: lambda s: datetime.strptime(s.strip(), DATEFORMAT),
    5: lambda s: convert(float(s.strip() or 0), TCRM_UNIT[5], 'mps'),
    6: lambda s: bearing2theta(float(s.strip() or 0) * np.pi / 180.),
    7: lambda s: convert(float(s.strip() or 0), TCRM_UNIT[7], 'Pa'),
    8: lambda s: convert(float(s.strip() or 0), TCRM_UNIT[8], 'Pa'),
}

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
    Read tracks from a track .csv file and return a list of :class:`Track`
    objects.

    This calls the function `readMultipleTrackData` to parse the track .csv
    file.

    :param str trackfile: the track data filename.

    :return: list of :class:`Track` objects.

    """

    tracks = []
    datas = readMultipleTrackData(trackfile)
    n = len(datas)

    sim, num = pattern.findall(os.path.basename(trackfile))
    
    data = readTrackData(trackfile)
    track = Track(data)
    track.trackfile = trackfile
    track.trackId = (sim, num)
    
    return [track]

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
    for f in trackfiles:
        msg = 'Loading tracks in %s' % f
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
    """
    files = os.listdir(path)
    trackfiles = [pjoin(path, f) for f in files if f.startswith('tracks')]
    msg = 'Loading %d track files in %s' % (len(trackfiles), path)
    log.info(msg)
    return loadTracksFromFiles(sorted(trackfiles))
