"""
:mod:`tracks2shp` -- save a track instance as a shape file
===========================================================

.. module:: tracks2shp
    :synopsis: Writes track instances to a shape file as a
               collection of line, ployline or point features.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import Utilities.shapefile as shapefile

import numpy as np
import logging

LOG = logging.getLogger(__name__)

# For all observation points/line segments:
OBSFIELD_NAMES = ('Indicator', 'TCID', 'Year', 'Month',
                  'Day', 'Hour', 'Minute', 'TElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'Pcentre',
                  'MaxWind', 'rMax', 'Penv', 'Category')
OBSFIELD_TYPES = ('N',)*17
OBSFIELD_WIDTH = (1, 6, 4, 2, 2, 2, 2, 6, 7, 7, 6, 6, 7, 6, 6, 7, 1)
OBSFIELD_PREC =  (0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0)

OBSFIELDS = [[n, t, w, p] for n, t, w, p in zip(OBSFIELD_NAMES,
                                                OBSFIELD_TYPES,
                                                OBSFIELD_WIDTH,
                                                OBSFIELD_PREC)]

TCRM_FIELD_NAMES = ('CycloneNumber', 'TimeElapsed', 'Longitude', 'Latitude',
                    'Speed', 'Bearing', 'CentralPressure', 'EnvPressure',
                    'rMax','Category')
TCRM_FIELD_TYPES = ('N',) * 10
TCRM_FIELD_WIDTH = (6, 7, 9, 9, 8, 8, 8, 8, 8, 1)
TCRM_FIELD_PREC =  (0, 2, 4, 4, 4, 4, 3, 3, 4, 0)

TCRM_FIELDS = [[n, t, w, p] for n, t, w, p in zip(TCRM_FIELD_NAMES,
                                                  TCRM_FIELD_TYPES,
                                                  TCRM_FIELD_WIDTH,
                                                  TCRM_FIELD_PREC)]

# For storing events as a single polyline:
EVENTFIELD_NAMES = ('TCID', 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Age',
                    'MinCP', 'MaxWind' )
EVENTFIELD_TYPES = ('N',)*9
EVENTFIELD_WIDTH = (6, 4, 2, 2, 2, 2, 6, 7, 7)
EVENTFIELD_PREC =  (0, 0, 0, 0, 0, 0, 2, 2, 1)

EVENTFIELDS = [[n, t, w, p] for n, t, w, p in zip(EVENTFIELD_NAMES,
                                                  EVENTFIELD_TYPES,
                                                  EVENTFIELD_WIDTH,
                                                  EVENTFIELD_PREC)]

def recdropfields(rec, names):
    """
    Return a new numpy record array with fields in *names* dropped.
    From http://matplotlib.org/api/mlab_api.html#matplotlib.mlab.rec_drop_fields

    :param rec: :class:`numpy.recarray` containing a number of fields.
    :param list names: List of names to drop from the record array.

    :returns: A new :class:`numpy.recarray` with the fields in *names*
              dropped.

    """

    names = set(names)

    newdtype = np.dtype([(name, rec.dtype[name]) for name in rec.dtype.names
                       if name not in names])

    newrec = np.recarray(rec.shape, dtype=newdtype)
    for field in newdtype.names:
        newrec[field] = rec[field]

    return newrec


def add_field(a, descr):
    """
    Add a field to the description of a track.
    """
    if a.dtype.fields is None:
        raise ValueError("`A' must be a structured numpy array")
    b = np.empty(a.shape, dtype=a.dtype.descr + descr)
    for name in a.dtype.names:
        b[name] = a[name]
    return b


def add_category(tracks):
    """
    Add a category field (for central pressure) to the tracks.
    """
    for track in tracks:
        track.data = add_field(track.data, [('Category', int)])

        for rec in track.data:
            if rec["CentralPressure"] < 930:
                rec["Category"] = 5
            elif rec["CentralPressure"] < 955:
                rec["Category"] = 4
            elif rec["CentralPressure"] < 970:
                rec["Category"] = 3
            elif rec["CentralPressure"] < 985:
                rec["Category"] = 2
            elif rec["CentralPressure"] < 999:
                rec["Category"] = 1
            else:
                rec["Category"] = 0


def tracks2point(tracks, outputFile, netcdf_format=False):
    """
    Writes tracks to a shapefile as a collection of point features.

    :type  tracks: list of :class:`Track` objects
    :param tracks: :class:`Track` features to store in a shape file

    :param str outputFile: Path to output file destination
    :param bool netcdf_format: Whether tracks are in TCRM format

    :raises: :mod:`shapefile.ShapefileException` if there is an error
             when attempting to save the file.

    """
    LOG.info("Writing point shape file: {0}".format(outputFile))
    sf = shapefile.Writer(shapefile.POINT)
    if netcdf_format:
        sf.fields = TCRM_FIELDS
    else:
        sf.fields = OBSFIELDS

    LOG.debug("Processing {0} tracks".format(len(tracks)))

    for track in tracks:
        track.data = recdropfields(track.data, ['Datetime'])
        for lon, lat, rec in zip(track.Longitude, track.Latitude, track.data):
            sf.point(lon, lat)
            sf.record(*rec)

    try:
        sf.save(outputFile)
    except shapefile.ShapefileException:
        LOG.exception("Cannot save shape file: {0}".format(outputFile))
        raise

    return

def tracks2line(tracks, outputFile, dissolve=False, netcdf_format=False):
    """
    Writes tracks to a shapefile as a collection of line features

    If dissolve==True, then each track feature is written as a
    single polyline feature, otherwise each track segment is
    stored as a separate feature.

    :type  tracks: list of :class:`Track` objects
    :param tracks: :class:`Track` features to store in a shape file

    :type  outputFile: str
    :param outputFile: Path to output file destination

    :type  dissolve: boolean
    :param dissolve: Store track features or track segments.

    :param bool netcdf_format: Whether tracks are in TCRM format

    :raises: :mod:`shapefile.ShapefileException` if there is an error
             when attempting to save the file.
    """
    LOG.info("Writing line shape file: {0}".format(outputFile))
    sf = shapefile.Writer(shapefile.POLYLINE)
    if netcdf_format:
        sf.fields = TCRM_FIELDS
    elif dissolve:
        sf.fields = EVENTFIELDS
    else:
        sf.fields = OBSFIELDS

    LOG.debug("Processing {0} tracks".format(len(tracks)))

    for track in tracks:
        track.data = recdropfields(track.data, ['Datetime'])

        if dissolve:
            if len(track.data) > 1:
                dlon = np.diff(track.Longitude)
                if dlon.min() < -180:
                    # Track crosses 0E longitude - split track
                    # into multiple parts:
                    idx = np.argmin(dlon)
                    parts = []
                    lines = zip(track.Longitude[:idx],
                                 track.Latitude[:idx])

                    parts.append(lines)
                    lines = zip(track.Longitude[idx+1:],
                                 track.Latitude[idx+1:])

                    parts.append(lines)
                    sf.line(parts)
                else:
                    lines = zip(track.Longitude, track.Latitude)
                    sf.line([lines])
            else:
                lines = zip(track.Longitude, track.Latitude)
                sf.line([lines])


            if netcdf_format:
                sf.record(*track.data[0])
            else:
                minPressure = track.trackMinPressure
                maxWind = track.trackMaxWind

                age = track.TimeElapsed.max()

                startYear = track.Year[0]
                startMonth = track.Month[0]
                startDay = track.Day[0]
                startHour = track.Hour[0]
                startMin = track.Minute[0]
                record = [track.CycloneNumber[0], startYear, startMonth, startDay,
                          startHour, startMin, age, minPressure, maxWind]
                sf.record(*record)

        else:
            if len(track.data) == 1:
                line = [[[track.Longitude, track.Latitude],
                        [track.Longitude, track.Latitude]]]
                sf.line(line)
                sf.record(*track.data[0])
            else:
                for n in range(len(track.data) - 1):
                    dlon = track.Longitude[n + 1] - track.Longitude[n]
                    if dlon < -180.:
                        # case where the track crosses 0E:
                        segment = [[[track.Longitude[n], track.Latitude[n]],
                                    [track.Longitude[n], track.Latitude[n]]]]
                    else:
                        segment = [[[track.Longitude[n],
                                     track.Latitude[n]],
                                    [track.Longitude[n + 1],
                                     track.Latitude[n + 1]]]]
                    sf.line(segment)
                    sf.record(*track.data[n])

                # Last point in the track:
                sf.line([[[track.Longitude[n + 1],
                           track.Latitude[n + 1]],
                              [track.Longitude[n + 1],
                               track.Latitude[n + 1]]]])
                sf.record(*track.data[n+1])

    try:
        sf.save(outputFile)
    except shapefile.ShapefileException:
        LOG.exception("Cannot save shape file: {0}".format(outputFile))
        raise

if __name__ == '__main__':

    from Utilities.loadData import loadTrackFile
    from Utilities.config import ConfigParser
    from Utilities.files import flStartLog

    import argparse
    import os
    from os.path import join as pjoin, dirname, realpath, splitext, isdir

    # pylint: disable=C0103

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', help='Input configuration file')
    parser.add_argument('-f', '--file', help='Input TC track file')
    parser.add_argument('-s', '--source',
                        help='Input TC track file source format')
    parser.add_argument('-v', '--verbose',
                        help='Print log messages to STDOUT',
                        action='store_true')
    args = parser.parse_args()

    config_file = args.config_file
    config = ConfigParser()
    config.read(config_file)

    logfile = config.get('Logging', 'LogFile')
    logdir = dirname(realpath(logfile))

    # If log file directory does not exist, create it
    if not isdir(logdir):
        try:
            os.makedirs(logdir)
        except OSError:
            logfile = pjoin(os.getcwd(), 'tracks2shp.log')

    logLevel = config.get('Logging', 'LogLevel')
    verbose = config.getboolean('Logging', 'Verbose')
    datestamp = config.getboolean('Logging', 'Datestamp')

    if args.verbose:
        verbose = True

    flStartLog(logfile, logLevel, verbose, datestamp)

    if args.file:
        track_file = args.file
    else:
        track_file = config.get('DataProcess', 'InputFile')

    if args.source:
        source = args.source
    else:
        source = config.get('DataProcess', 'Source')

    output_path = dirname(realpath(track_file))
    filename, ext = splitext(track_file)
    pt_output_file = filename + '_pt.shp'
    line_output_file = filename + '_line.shp'
    dissolve_output_file = filename + '_dissolve.shp'

    if track_file.endswith(".nc"):

        from Utilities.track import ncReadTrackData
        tracks = ncReadTrackData(track_file)
        netcdf_format = True

    elif track_file.endswith(".csv"):
        tracks = loadTrackFile(config_file, track_file, source,
                               calculateWindSpeed=True)
        netcdf_format = False

    else:
        raise ValueError("format of {} is not recognizable".format(track_file))

    add_category(tracks)
    tracks2point(tracks, pt_output_file, netcdf_format=netcdf_format)
    tracks2line(tracks, line_output_file, netcdf_format=netcdf_format)
    tracks2line(tracks, dissolve_output_file, dissolve=True, netcdf_format=netcdf_format)
    LOG.info("Completed tracks2shp")

