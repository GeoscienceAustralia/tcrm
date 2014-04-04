"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience
    Australia)

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see
    <http://www.gnu.org/licenses/>.

 Title: shptools.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 05/30/08 8:54:AM
 Description: A collection of useful functions to handle shapefiles.
 Uses the shapelib library, normally provided with matplotlib v0.99
 and higher.

 Version :$Rev: 804 $

 $Id: shptools.py 804 2013-08-13 22:12:35Z carthur $
"""
import os, sys

import logging as log

filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
import shapelib
import dbflib
import shapefile
import numpy as np

__version__ = '$Id: shptools.py 804 2013-08-13 22:12:35Z carthur $'


def _shpGetField(shpFile, fieldName=None):
    """
Return an array containing the specified field (if it is contained
within the corresponding dbf file). If no field name is specified,
return all fields.
    """

    try:
        dbffh = dbflib.open(shpFile)
    except IOError:
        log.warn("Cannot open %s.dbf"%shpFile)
        raise
    nrecs = dbffh.record_count()
    nfields = dbffh.field_count()
    data = {}
    for f in xrange(nfields):
        fname = dbffh.field_info(f)[1]
        data[fname]=[]
    for rec in xrange(nrecs):
        recdata = dbffh.read_record(rec)
        for key in recdata.keys():
            data[key].append(recdata[key])
    if fieldName:
        if data.has_key(fieldName):
            returnData = data[fieldName]
        else:
            log.warn("%s.dbf does not contain a field called %s - returning all data" \
                         % (shpFile,fieldName))
            returnData = data
    else:
        returnData = data
    return returnData

def _shpType(shpFile):
    """
Return the type of objects contained in the shapefile as a string
    """
    try:
        shpfh = shapelib.open(shpFile)
    except IOError:
        log.warn("Cannot open %s"%shpFile)
        raise

    shp_type = shapelib.type_name(shpfh.info()[1])
    shpfh.close()
    return shp_type

def _shpExtents(shpFile):
    """
    Return the extents of the objects in the shapefile
    """
    pass

def _shpCreatePoint(lon, lat):
    """
    Create shape object of point type.
    Input: lon - longitude of point
           lat - latitude of point
    Output: shpObj shapelib point object
    """
    coord = lon, lat
    shp_coord = [coord]
    shp_point = [shp_coord]
    shpObj = shapelib.SHPObject(shapelib.SHPT_POINT, 1, shp_point)

    return shpObj

def _shpCreateLine(lons, lats):
    """
    Create a shape object of line type.
    Input: lons - array of longitudes of line segment
           lats - array of latitudes of line segment
    Output: shpObj shapelib line object.
    """
    coords = []
    for lon, lat in zip(lons, lats):
        coords.append( (lon, lat) )
    shp_coords = [coords]
    shpObj = shapelib.SHPObject( shapelib.SHPT_ARC, 1, shp_coords )
    return shpObj

def _shpCreateFile(fileName, shptype, fields):
    """
    Create a shapefile (and a corresponding dbf file) of the give type,
    containing the given fields.
    Input: fileName - full path (excluding extension!) to the shapefile
           to create.
           shptype - shapelib object type (these are integer
           values, but you can also use the shapelib.SHPT_ value).
           fields - a dictionary of dictionaries with field names as
           keys, and each sub-dictionary containing keys of 'Type',
           'Length','Precision' and 'Data':
           'Type' must be one of the following integer values:
                    0 - strings
                    1 - integers
                    2 - doubles
                    4 - Invalid
    Output: shapefile and dbffile objects
    """
    try:
        fshp = shapelib.create(fileName,shptype)
    except IOError:
        log.critical("Failed to create shapefile: %s.shp"%fileName)
        raise IOError
    try:
        fdbf = dbflib.create(fileName)
    except IOError:
        log.critical("Failed to create dbffile: %s.dbf"%fileName)
        raise IOError
    fieldNames = fields.keys()
    for f in sorted(fieldNames):
        fieldType = fields[f]['Type']
        fieldLength = fields[f]['Length']
        # Force the precision to be zero unless the field is a double
        if fieldType==2:
            fieldPrec = fields[f]['Precision']
        else:
            fieldPrec = 0
        fdbf.add_field(f, fieldType, fieldLength, fieldPrec)

    return fshp, fdbf

def _shpSaveTrackFile(filename, lon, lat, fields):
    """
    Save track data to shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.
    Input: filename - name for the output shapefile, excluding extension
           lon - array of longitudes of TC track data
           lat - array of latitudes of TC track data
           fields - a dictionary of dictionaries with field names as
                    keys, and each sub-dictionary containing keys of
                    'Type','Length','Precision' and 'Data'
                    'Type' must be one of the following integer values:
                    0 - strings
                    1 - integers
                    2 - doubles
                    4 - Invalid
    Output: None
    """
    fshp, fdbf = shpCreateFile(filename, 3, fields)
    for i in xrange(len(lon)-1):
        if fields['Index']['Data'][i] == fields['Index']['Data'][i+1]:
            obj = shpCreateLine([lon[i], lon[i+1]], [lat[i], lat[i+1]])
        else:
            obj = shpCreateLine([lon[i], lon[i]], [lat[i], lat[i]])
        rec = fshp.write_object(-1, obj)
        recdata = ()
        fieldNames = fields.keys()
        for f in sorted(fieldNames):
            recdata += (fields[f]['Data'][i],)
        fdbf.write_record(rec, recdata)
    fshp.close()
    fdbf.close()
    return

def _shpReadShapeFile(shapeFile):
    """
    Read in shape file, return vertices and information
    Based on the mpl_toolkit.basemap.Basemap.readshapefile() function
    Changes include removal of the optional kwargs to draw bounds on an
    axes object, and no conversion from lat/lon coords to map (x/y)
    coords.

    Input: shapeFile - path to shapefile components
    Output: coords - array of shapefile vertices or points (lon,lat)
            attributes - list of dictionaries, one for each shape,
            containing attributes of each shape from the corresponding
            dbf file.
    """
    try:
        shp = shapelib.ShapeFile(shapeFile)
    except:
        raise IOError, "Error reading shapefile %s.shp" % shapeFile
    try:
        dbf = dbflib.open(shapeFile)
    except:
        raise IOError, "Error reading dbffile %s.dbf" % shapeFile

    info = shp.info()
    if info[1] not in [1,3,5,8]:
        raise ValueError, "shpReadShapeFile can only handle 2D shape types"

    msg = """Shapefile must have lat/lon vertices - it appears this one
    has vertices in map projection coordinates. Convert the shapefile to
    geographic coordinates using the shpproj utility from the shapelib
    tools"""
    if info[1] in [1,8]:
        coords = []
        nelem = shp.info()[0]
        for elem in xrange(nelem):
            shpObj = shp.read_object(elem)
            verts = shpObj.vertices()
            lons, lats = zip(*verts)
            if max(lons) > 721. or min(lons) < -721. or max(lats) > 91. or \
               min(lats) < -91.:
                raise ValueError, msg
            if len(verts) > 1:
                coords.append(zip(lons, lats))
            else:
                coords.append((lons[0], lats[0]))
        attributes = [dbf.read_record(i) for i in xrange(nelem)]

    else:
        coords = []
        attributes = []
        for npoly in xrange(shp.info()[0]):
            shpObj = shp.read_object(npoly)
            verts = shpObj.vertices()
            rings = len(verts)
            for ring in xrange(rings):
                lons, lats = zip(*verts[ring])
                if max(lons) > 721. or min(lons) < -721. or max(lats) > 91. or \
                   min(lats) < -91.:
                    raise ValueError, msg
                coords.append(zip(lons, lats))
                if ring == 0:
                    shapedict = dbf.read_record(npoly)
                shapedict['RINGNUM'] = ring+1
                shapedict['SHAPENUM'] = npoly+1
                attributes.append(shapedict)
    shp.close()
    dbf.close()
    return coords, attributes

def _shpCreateDBFFile(fileName, fields):
    """
    Create a dbf file, containing the given fields
    Input: fileName - full path (excluding extension!) to the shapefile
                      to create
           fields - a dictionary of dictionaries with field names as
                    keys, and each sub-dictionary containing keys of
                    'Type','Length','Precision' and 'Data'
                    'Type' must be one of the following integer values:
                    0 - strings
                    1 - integers
                    2 - doubles
                    4 - Invalid
    Output: dbffile object
    """
    try:
        fdbf = dbflib.create(fileName)
    except IOError:
        print "Failed to create dbffile: %s.dbf"%fileName
        return None
    fieldNames = fields.keys()
    for f in sorted(fieldNames):
        fieldType = fields[f]['Type']
        fieldLength = fields[f]['Length']
        # Force the precision to be zero unless the field is a double
        if fieldType == 2:
            fieldPrec = fields[f]['Precision']
        else:
            fieldPrec = 0
        fdbf.add_field(f, fieldType, fieldLength, fieldPrec)

    return fdbf

def _shpSaveDBFFile(fileName, fields, nrecords):
    """
    Save track data to shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.
    Input: fileName - dbf file object created by shpCreateFile
           fields - a dictionary of dictionaries with field names as
           keys, and each sub-dictionary containing keys of 'Type',
           'Length','Precision' and 'Data'
                    'Type' must be one of the following integer values:
                    0 - strings
                    1 - integers
                    2 - doubles
                    4 - Invalid
           nrecords - integer number of records to be inserted into the
                      dbf file
    Output: None
    """
    fdbf = shpCreateDBFFile(fileName, fields)
    for rec in xrange(nrecords):
        recdata = ()
        fieldNames = fields.keys()
        for f in sorted(fieldNames):
            recdata += (fields[f]['Data'][rec],)
        fdbf.write_record(rec, recdata)
    fdbf.close()
    return 1

def parseData( data ):
    """
    Parse a dict of dicts to generate a list of lists describing
    the fields, and an array of the corresponding data records
    """
    fields = []
    records = []
    for k in data.keys():
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


def shpCreateFile(outputFile, shapes, data):
    """
    Create a shapefile (and a corresponding dbf file) of the give type,
    containing the given fields.
    Input: fileName - full path (excluding extension!) to the shapefile
           to create.
           shptype - shapelib object type (these are integer
           values, but you can also use the shapelib.SHPT_ value).
           fields - a dictionary of dictionaries with field names as
           keys, and each sub-dictionary containing keys of 'Type',
           'Length','Precision' and 'Data':
           'Type' must be one of the following integer values:
                    0 - strings
                    1 - integers
                    2 - doubles
                    4 - Invalid
    Output: None.
    """
    fields, records = parseData(data)

    log.info("Writing data to {0}.shp".format(outputFile))
    w = shapefile.Writer(shpType)
    # Set the fields:
    for f in fields:
        w.field(*f)

    # Add shapes and records:
    for r,s in zip(records, shapes):
        if len(s.parts) > 1:
            start = 0
            new_parts = []
            for p in s.parts[1:]:
                new_parts.append(list(s.points[start:p-1]))
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

def shpSaveTrackFile(outputFile, lon, lat, index, data):
    """
    Save track data to shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.
    Input: filename - name for the output shapefile, excluding extension
           lon - array of longitudes of TC track data
           lat - array of latitudes of TC track data
           index - array of 1's and 0's, where 1 indicates the start of a
                   new TC event.
           data - a dictionary of dictionaries with field names as
                    keys, and each sub-dictionary containing keys of
                    'Type','Length','Precision' and 'Data'
                    'Type' must be one of the following integer values:
                    0 - strings
                    1 - integers
                    2 - doubles
                    4 - Invalid
    Output: None
    """

    fields, records = parseData(data)
    logger.info("Writing data to {0}.shp".format(outputFile))
    w = shapefile.Writer(shapefile.POLYLINE)

    # Set the fields:
    for f in fields:
        w.field(*f)

    for i in xrange(len(lon)-1):
        if index[i] == index[i + 1]:
            obj = [[lon[i], lat[i]], [lon[i + 1], lat[i + 1]]]
        else:
            obj = [[lon[i], lat[i]], [lon[i], lat[i]]]

        w.line([obj])
        w.record(*records[i])

    # Append the last point
    obj = [[lon[-1], lat[-1]], [lon[-1], lat[-1]]]
    w.line([obj])
    w.record(*records[-1])

    try:
        w.save(outputFile)
    except shapefile.ShapefileException:
        logger.exception("Failed to write to {0}".
                         format(os.path.abspath(outputFile)))
        raise
    except AssertionError:
        logger.exception("Problem in length and precision of fields")
        raise
    return


def shpWriteShapeFile(outputFile, shpType, fields, shapes, records):
    """
    Save data to a shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.
    Input: outputFile - dbf file object created by shpCreateFile
           shpType -
           data - a dictionary of dictionaries with field names as
           keys, and each sub-dictionary containing keys of 'Type',
           'Length','Precision' and 'Data'
           'Type' must be one of the following integer values:
           0 - strings
           1 - integers
           2 - doubles
           4 - Invalid


    Output: None
    """
    log.info("Writing data to {0}.shp".format(outputFile))
    w = shapefile.Writer(shpType)

    fields = []
    for key in data.keys():
        fieldName = key
        fieldType = data[key]['Type']
        fieldLength = data[key]['Length']
        fieldPrec = data[key]['Precision']

        field = [fieldName, fieldType, fieldLength, fieldPrecision]
        fields.append(field)

    # Set the fields:
    for f in fields:
        w.field(*f)

    # Add shapes and records:
    for rec,shp in zip(records, shapes):
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
        logger.exception("Unable to save data to {0}".format(outputFile) )

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

    Output: dictionary keyed by object id number or optionally a field
            name, with values being arrays of vertices of the
            corresponding shape object.

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
        shapes = sf.shapes()

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

    """

    log.debug("Extracting {0} from records".format(field_name))

    try:
        sf = shapefile.Reader(shape_file,"rb")

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
        output = np.array([records[rec][idx] for rec in xrange(nrecords)],
                              dtype=dtype)

    else:
        # Otherwise, return a list:
        output = [records[rec][idx] for rec in xrange(nrecords)]


    return output

def shpReadShapeFile(shape_file):
    """
    Return the vertices and records for the given shape file

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
