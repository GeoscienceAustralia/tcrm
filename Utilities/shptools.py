"""
:mod:`shptools` - helper functions for manipulating shape files
===============================================================

.. module: shptools
    :synopsis: A collection of useful functions to manipulate shapefiles. Uses the `shapefile
               library <http://code.google.com/p/pyshp/>`_



"""

import logging

from . import shapefile
import numpy as np


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


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




def parseData(data):
    """
    Parse a dict of dicts to generate a list of lists describing
    the fields, and an array of the corresponding data records

    :param dict data: a dict of dicts with field names as keys, and each
                      sub-dict containing keys of 'Type', 'Length',
                      'Precision' and 'Data'.

    :returns: fields, records
    :rtype: list

    """

    fields = []
    records = []
    for k in list(data.keys()):
        t = data[k]['Type']
        l = data[k]['Length']
        if t == 0:
            nt = 'C'
            p = 0
        elif t == 1:
            nt = 'N'
            p = 0
        elif t == 2:
            nt = 'N'
            p = data[k]['Precision']
        else:
            nt = t
            p = 0

        fields.append([k, nt, l, p])
        records.append([data[k]['Data']])

    return fields, records


def shpCreateFile(outputFile, shapes, data, shpType):
    """
    Create a shapefile of the give type, containing the given fields.

    :param str outputFile: full path (excluding extension!) to the shapefile
                           to create.
    :param shapes: Collection of shape objects representing the geometry of
                   features.
    :type shapes: :class:`shapefile._Shape`

    :param int shptype: :class:`shapefile` object type (these are integer
                        values, but you can also use the shapelib.SHPT value).
    :param dict fields: a dictionary of dictionaries with field names as
                        keys, and each sub-dictionary containing keys of 'Type',
                        'Length','Precision' and 'Data':
                        'Type' must be one of the following integer values:
                        0 - strings
                        1 - integers
                        2 - doubles
                        4 - Invalid

    :raises: :mod:`shapefile.ShapefileException` if unable to write the file.

    """
    fields, records = parseData(data)

    log.info("Writing data to {0}.shp".format(outputFile))
    w = shapefile.Writer(shpType)
    # Set the fields:
    for f in fields:
        w.field(*f)

    # Add shapes and records:
    for r, s in zip(records, shapes):
        if len(s.parts) > 1:
            start = 0
            new_parts = []
            for p in s.parts[1:]:
                new_parts.append(list(s.points[start:p - 1]))
                start = p
            new_parts.append(list(s.points[p:-1]))
            w.poly(parts=new_parts)
        else:
            w.poly(parts=[s.points])
        w.record(*r)

    try:
        w.save(outputFile)
    except shapefile.ShapefileException:
        log.exception("Unable to save data to {0}".format(outputFile) )

    return

def shpSaveTrackFile(outputFile, tracks, fmt="point"):
    """
    Save track data to shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.

    :param str outputFile: name for the output shapefile, excluding extension.
    :param tracks: collection of track features.
    :type tracks: :class:`Track` object
    :param format: Type of features to save. "point" will save each record as a
                   single point. "lines" will save each individual TC track as
                   a single (POLYLINE) feature. "segments" will save each
                   segment of a track to a line feature between consecutive
                   observations.

    """

    if fmt == "points":
        tracks2point(tracks, outputFile)
    elif fmt == "lines":
        tracks2line(tracks, outputFile, dissolve=True)
    elif fmt == "segments":
        tracks2line(tracks, outputFile, dissolve=False)
    else:
        log.critical(("Unknown output format - must be one of 'points', "
                        "'lines' or 'segments'"))
        return



def shpWriteShapeFile(outputFile, shpType, fields, shapes, records):
    """
    Save data to a shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.

    :param outputFile: A dbf file object created by shpCreateFile
    :param shpType: The type of features to be created.
    :param dict data: A dictionary of dictionaries with field names as
                      keys, and each sub-dictionary containing keys of 'Type',
                      'Length','Precision' and 'Data'
                      'Type' must be one of the following integer values:
                      0 - strings
                      1 - integers
                      2 - doubles
                      4 - Invalid

    :raises: :mod:`shapefile.ShapefileException` if unable to write the file.

    """
    log.info("Writing data to {0}.shp".format(outputFile))
    w = shapefile.Writer(shpType)

    # Set the fields:
    for f in fields:
        w.field(*f)

    # Add shapes and records:
    for rec, shp in zip(records, shapes):
        if len(shp.parts) > 1:
            start = 0
            new_parts = []
            for p in shp.parts[1:]:
                new_parts.append(list(shp.points[start:p-1]))
                start = p
            new_parts.append(list(shp.points[p:-1]))
            w.poly(parts=new_parts)
        else:
            w.poly(parts=[shp.points])
        w.record(*rec)

    try:
        w.save(outputFile)
    except shapefile.ShapefileException:
        log.exception("Unable to save data to {0}".format(outputFile) )

    return

def shpGetVertices(shape_file, key_name=None):
    """
    Returns a dictionary of arrays containing coordinate pairs
    representing vertices contained in the shapefile.

    Dictionary keys can either be a field contained within the
    corresponding dbf file (through the optional kwarg keyName), or if
    no key name is provided, the object id number (i.e. the record
    number) is used.

    WARNING:
    If any records share the same value for the chosen key, then only
    one record will be retained in the returned values.

    Input:
    :type  shape_file: str
    :param shape_file: path to a shape file, excluding the extension

    :type  key_name: optional str
    :param key_name: name of a field in the shape file that
                     acts as a key for the dictionary of returned
                     vertices

    :return: dict keyed by object id number or optionally a field
            name, with values being arrays of vertices of the
            corresponding shape object.
    :rtype: dict

    Example: vertices = shpGetVertices('/foo/bar/baz/shp', 'FIELD1')

    This function is retained for backwards compatibility.
    We recommend using the shapefile interface directly for extracting
    shapes and records.

    """

    try:
        sf = shapefile.Reader(shape_file, "rb")

    except shapefile.ShapefileException:
        log.exception("Cannot open {0} for reading".format(shape_file))
        raise

    else:
        fields = sf.fields[1:] # Drop first field (DeletionFlag)

        field_names = [fields[i][0] for i in range(len(fields))]
        if key_name and (key_name in field_names):
            keyIndex = field_names.index(key_name)

        vertices = {}

        for i, shprec in enumerate(sf.shapeRecords()):
            if key_name and (key_name in field_names):
                recValue = shprec.record[keyIndex]
                vertices[recValue] = shprec.shape.points
            else:
                vertices[i] = shprec.shape.points

    return vertices

def shpGetField(shape_file, field_name, dtype=float):
    """
    Extract from the records the value of the field corresponding
    to fieldname.

    :type  shpFile: str
    :param shpFile: path to a valid shape file

    :type  fieldname: str
    :param fieldname: name of a field in the attribute table of the
                      shape file (.dbf)

    :type  dtype: `dtype`
    :param dtype: type of values in the requested field. Default is
                  float

    :return: the value of the given field for each feature
    :rtype: array or list (if `dtype` is a string)

    """

    log.debug("Extracting {0} from records".format(field_name))

    try:
        sf = shapefile.Reader(shape_file, "rb")

    except shapefile.ShapefileException:
        log.exception("Cannot read {0}for ".format(shape_file))
        raise

    else:
        fields = sf.fields[1:] # Drop first field (DeletionFlag)

        field_names = [fields[i][0] for i in range(len(fields))]
        if field_name not in field_names:
            log.warn("No field '{0}' in the list of fieldnames" .
                    format(field_name))
            log.warn("Unable to proceed with processing")
            raise ValueError

    records = sf.records()
    nrecords = len(records)

    # Get the index of the required field name:
    idx = field_names.index(field_name)

    if dtype != str:
        # For non-string data, return a numpy array:
        output = np.array([records[rec][idx] for rec in range(nrecords)],
                              dtype=dtype)

    else:
        # Otherwise, return a list:
        output = [records[rec][idx] for rec in range(nrecords)]


    return output

def shpReadShapeFile(shape_file):
    """
    Return the vertices and records for the given shape file

    :param str shape_file: path of input shape file.

    :return: vertices
    :rtype: dict
    :return: records
    :rtype: dict

    """

    try:
        sf = shapefile.Reader(shape_file,"rb")
    except shapefile.ShapefileException:
        log.exception("Cannot read {0}".format(shape_file))
        raise

    else:
        records = sf.records()
        vertices = shpGetVertices(shape_file)

    return vertices, records

def tracks2point(tracks, outputFile):
    """
    Writes tracks to a shapefile as a collection of point features

    :type  tracks: list of :class:`Track` objects
    :param tracks: :class:`Track` features to store in a shape file

    :param str outputFile: Path to output file destination

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
        raise
