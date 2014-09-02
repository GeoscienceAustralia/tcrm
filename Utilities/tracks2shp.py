"""
:mod:`tracks2shp` -- save a track instance as a shape file
===========================================================

.. module:: tracks2shp
    :synopsis: Writes track instances to a shape file as a
               collection of line, ployline or point features.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import shapefile
from itertools import izip
import numpy as np


# For all observation points/line segments:
OBSFIELD_NAMES = ('Indicator', 'TCID', 'Year', 'Month',
                  'Day', 'Hour', 'Minute', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'WindSpeed', 'rMax', 'EnvPressure')
OBSFIELD_TYPES = ('N',)*16
OBSFIELD_WIDTH = (1, 6, 4, 2, 2, 2, 2, 6, 7, 7, 6, 6, 7, 6, 6, 7)
OBSFIELD_PREC =  (0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 1, 1, 1, 1, 1)

OBSFIELDS = [[n, t, w, p] for n, t, w, p in zip(OBSFIELD_NAMES,
                                                OBSFIELD_TYPES,
                                                OBSFIELD_WIDTH,
                                                OBSFIELD_PREC)]

# For storing events as a single polyline:
EVENTFIELD_NAMES = ('TCID', 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Age',
                    'MinPressure', 'MaxWindSpeed' )
EVENTFIELD_TYPES = ('N',)*9
EVENTFIELD_WIDTH = (6, 4, 2, 2, 2, 2, 6, 7, 7)
EVENTFIELD_PREC =  (0, 0, 0, 0, 0, 0, 2, 2, 1)

EVENTFIELDS = [[n, t, w, p] for n, t, w, p in zip(EVENTFIELD_NAMES,
                                                  EVENTFIELD_TYPES,
                                                  EVENTFIELD_WIDTH,
                                                  EVENTFIELD_PREC)]


def tracks2point(tracks, outputFile):
    """
    Writes tracks to a shapefile as a collection of point features

    :type  tracks: list of :class:`Track` objects
    :param tracks: :class:`Track` features to store in a shape file

    :param str outputFile: Path to output file destination

    :raises: :mod:`shapefile.ShapefileException` if there is an error
             when attempting to save the file.

    """
    sf = shapefile.Writer(shapefile.POINT)
    sf.fields = OBSFIELDS
    
    for track in tracks:
        for x, y, rec in zip(track.Longitude, track.Latitude, track.data):
            sf.point(x, y)
            sf.record(*rec)

    try:
        sf.save(outputFile)
    except shapefile.ShapefileException:
        raise

    return

def tracks2line(tracks, outputFile, dissolve=False):
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

    :raises: :mod:`shapefile.ShapefileException` if there is an error
             when attempting to save the file.
    """

    sf = shapefile.Writer(shapefile.POLYLINE)
    if dissolve:
        sf.fields = EVENTFIELDS
    else:
        sf.fields = OBSFIELDS

    for track in tracks:
        if dissolve:
            if len(track.data) > 1:
                dlon = np.diff(track.Longitude)
                if dlon.min() < -180:
                    # Track crosses 0E longitude - split track
                    # into multiple parts:
                    idx = np.argmin(dlon)
                    parts = []
                    lines = izip(track.Longitude[:idx],
                                 track.Latitude[:idx])
                    
                    parts.append(lines)
                    lines = izip(track.Longitude[idx+1:],
                                 track.Latitude[idx+1:])
                    
                    parts.append(lines)
                    sf.line(parts)
                else:
                    lines = izip(track.Longitude, track.Latitude)
                    sf.line([lines])
            else:
                lines = izip(track.Longitude, track.Latitude)
                sf.line([lines])


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
                        segment =[[[track.Longitude[n], track.Latitude[n]],
                                   [track.Longitude[n], track.Latitude[n]]]]
                    else:
                        segment = [[[track.Longitude[n], track.Latitude[n]],
                                    [track.Longitude[n + 1], track.Latitude[n + 1]]]]
                    sf.line(segment) 
                    sf.record(*track.data[n])
                    
                # Last point in the track:
                sf.line([[[track.Longitude[n + 1], track.Latitude[n + 1]],
                              [track.Longitude[n + 1], track.Latitude[n + 1]]]])
                sf.record(*track.data[n+1])

    try:
        sf.save(outputFile)
    except shapefile.ShapefileException:
        raise

