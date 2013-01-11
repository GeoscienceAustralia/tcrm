#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Title: shptools.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 05/30/08 8:54:AM
 Description: A collection of useful functions to handle shapefiles.
 Uses the shapelib library,  normally provided with matplotlib v0.99 and
 higher.

 Version :$Rev: 686 $

 $Id: shptools.py 686 2012-03-29 04:24:59Z carthur $
"""
import os, sys, pdb, logging

filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
import shapelib
import dbflib

__version__ = '$Id: shptools.py 686 2012-03-29 04:24:59Z carthur $'

logger = logging.getLogger()

def shpGetVertices(shpFile, keyName=None):
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

    As yet untested on MULTIPOINT shapefiles (polylines, polygons).

    Input: shpFile - path to a shapefile, excluding the extension
           (shapelib/dbflib automatically append the correct extension)
           keyName - (optional) name of a field in the dbf file that
                     acts as a key for the dictionary of returned
                     vertices
    Output: dictionary keyed by object id number or optionally a field
            name, with values being arrays of vertices of the
            corresponding shape object.
    Example: vertices = shpGetVertices(shpFile,keyName)
    """
    try:
        shpfh = shapelib.open(shpFile)
    except IOError:
        logger.warn("Cannot open %s"%shpFile)
        raise
    try:
        dbffh = dbflib.open(shpFile)
    except IOError:
        logger.warn("Cannot open %s.dbf"%shpFile)
        raise
    vertices = {}
    nshapes = shpfh.info()[0]
    nfields = dbffh.field_count()

    for oid in xrange(nshapes):
        shpdata = shpfh.read_object(oid)
        dbfdata = dbffh.read_record(oid)
        nparts = len(shpdata.vertices())
        v = []
        for p in xrange(nparts):
            v.append(shpdata.vertices()[p])
        if keyName and dbfdata.has_key(keyName):
            vertices[dbfdata[keyName]] = v
        else:
            vertices[oid] = v

    shpfh.close()
    dbffh.close()
    return vertices

def shpGetField(shpFile, fieldName=None):
    """
Return an array containing the specified field (if it is contained
within the corresponding dbf file). If no field name is specified,
return all fields.
    """

    try:
        dbffh = dbflib.open(shpFile)
    except IOError:
        logger.warn("Cannot open %s.dbf"%shpFile)
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
            logger.warn("%s.dbf does not contain a field called %s - returning all data" \
                         % (shpFile,fieldName))
            returnData = data
    else:
        returnData = data
    return returnData

def shpType(shpFile):
    """
Return the type of objects contained in the shapefile as a string
    """
    try:
        shpfh = shapelib.open(shpFile)
    except IOError:
        logger.warn("Cannot open %s"%shpFile)
        raise

    shp_type = shapelib.type_name(shpfh.info()[1])
    shpfh.close()
    return shp_type

def shpExtents(shpFile):
    """
    Return the extents of the objects in the shapefile
    """
    pass

def shpCreatePoint(lon, lat):
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

def shpCreateLine(lons, lats):
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

def shpCreateFile(fileName, shptype, fields):
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
        logger.critical("Failed to create shapefile: %s.shp"%fileName)
        raise IOError
    try:
        fdbf = dbflib.create(fileName)
    except IOError:
        logger.critical("Failed to create dbffile: %s.dbf"%fileName)
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

def shpSaveTrackFile(filename, lon, lat, fields):
    """
    Save track data to shapefile. The fields are sorted using the same
    function as in shpCreateFile, so the fields should be in the correct
    order.
    Input: fshp - shapefile object created by shpCreateFile
           fdbf - dbf file object created by shpCreateFile
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

def shpReadShapeFile(shapeFile):
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

def shpCreateDBFFile(fileName, fields):
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

def shpSaveDBFFile(fileName, fields, nrecords):
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
