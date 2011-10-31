#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

 Title: my_tool.py
 Author: Geoff Xu
 CreationDate: 2005-12-12
 Description:
 Defines a class for miscellaneous tools. Includes methods for data
 loading, mathematical and plotting functions. This is a parent of
 GridTool, StatTool and WindfieldTool.

 Version: $Rev: 512 $

 ModifiedBy: C. Arthur
 ModifiedDate: 2006-10-24
 Modification: Added descriptive headers and metadata

 ModifiedBy: C. Arthur
 ModifiedDate: 2006-11-03
 Modification: Changed array<atan2, cos, sin> to use list comprehension rather than for loops. Modified removeNan.

 ModifiedBy: C. Arthur
 ModifiedDate: 2006-11-06
 Modification: Changed from a class instance to a module that contains only functions.

 ModifiedBy: N. Habili, nariman.habili@ga.gov.au
 ModifiedDate: 2007-01-17
 Modification: Removed arrayatan2, arraycos, and arraysin.

 Version: 351
 ModifiedBy: C. Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-01-31
 Modification: Updated textread() (removed need for separate file textread.py).

 Version: 403
 ModifiedBy: C.Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-03-01
 Modification: Added readgrid() to read in ASCII grid file

 Version: 493
 ModifiedBy: C.Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-03-28
 Modification: textread takes optional second argument nullValue to replace missing values in the data being read.

 Version: 536
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-07-06
 Modification: Added SampleGrid class to read an ascii grd file and return values for a given lon,lat

 Version: 541
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-10-03 11:31:AM
 Modification: Added helper functions to create output paths.

 Version: 544
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-10-19 9:50:AM
 Modification: Added savegrid() function

 Version: $Rev: 512 $
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-10-29 10:20:AM
 Modification: Upgraded savegrid to convert to UTM coords if requested

 $Id: my_tool.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
from scipy import array, asarray, randn, shape
from lat_long_UTM_conversion import *
import numpy
import math
import pylab
import sys
import csv
import datetime
import ConfigParser
__version__ = '$Id: my_tool.py 512 2011-10-31 07:20:38Z nsummons $'
""" Functions
    ---------
    load(fname,comments='%',delimiter=None) : N column array of float
                                             (N depend on number of
                                              columns in the file)
        Load ASCII data from fname into an array and return the array.
        The data must be regular, same number of values in every row

    save(fname, X, header='', delimiter=' ', fmt='%.18e')
        Save the data in X to file fname using fmt string to convert the
        data to strings

    textread(data) : tuple of 1D arrays
        read data from csv file and store in arrays representing each
        column

    removeNan(1D array or list): 1D list
        Removes all elements in an array or list for which value is 0

    get3DAusBndy(aus_blat,aut_blon): 3D array of float
        Sub-function that converts two 1D arrays of Australian coastline
        coordinates aus_blat and aus_blon into 3D coordinates

"""

def load(fname, comments='%', delimiter=None):
    """
    Load ASCII data from fname into an array and return the array.
    The data must be regular, same number of values in every row
    code originated from pylab.py in matplotlib
    by John D. Hunter <jdhunter@ace.bsd.uhicago.edu>
    modified by Geoff Xu, 2005
    """
    if pylab.is_string_like(fname):
        if fname.endswith('.gz'):
            import gzip
            fh = gzip.open(fname)
        else:
            fh = file(fname)
    elif hasattr(fname, 'seek'):
        fh = fname
    else:
        raise ValueError('fname must be a string or file handle')
    X = []
    numCols = None
    for line in fh:
        line = line[:line.find(comments)].strip()
        if not len(line): continue
        row = []
        for val in line.split(delimiter):
            if val == 'NaN': val = Nan
            row.append(float(val))
        thisLen = len(row)
        if numCols is not None and thisLen != numCols:
            raise ValueError('All rows must have the same number of columns')
        X.append(row)
    X = array(X)
    r,c = X.shape
    if r == 1 or c == 1:
        X.shape = max([r,c]),
    return X

def save(fname, X, header='', delimiter=' ', fmt='%.18e'):
    """
    Save the data in X to file fname using fmt string to convert the
    data to strings, originated from pylab.py in matplotlib
    by John D. Hunter <jdhunter@ace.bsd.uhicago.edu>
    modified by Geoff Xu, 2006
    """

    if pylab.is_string_like(fname):
        if fname.endswith('.gz'):
            import gzip
            fh = gzip.open(fname, 'wb')
        else:
            fh = file(fname, 'w')
    elif hasattr(fname, 'seek'):
        fh = fname
    else:
        raise ValueError('fname must be a string or file handle')

    fh.write('%'+header+'\n')
    X = asarray(X)
    origShape = None
    if len(X.shape) == 1:
        origShape = X.shape
        X.shape = len(X), 1
    for row in X:
        fh.write(delimiter.join([fmt%val for val in row]) + '\n')

    if origShape is not None:
        X.shape = origShape


def removeNum(a, Num=sys.maxint):
    """removeNum(a [, Num]):
    a function that removes all elements for which value is Num in an
    array (Num is by default sys.maxint)
    """
    if shape(a) == ():
        raise ValueError, "condition must be 1-d array"
    return a.compress(a < Num)

def get3DAusBndy(lat, lon):
    """get3DAusBndy(lat,lon)
    sub function that converts the two 1D arrays of australia coastline
    coordinates aus_blat and aus_blon into 3D coordinates
    """
    if len(lon) != len(lat):
        raise ValueError, "Array sizes do not match"
    return array([[lat[i], lon[i], 0] for i in range(len(lat))])

def rnorm(n=1, m=0, s=1):
    """rnorm(n=1, m=0, s=1):
    A method that generates normal random numbers
    p = rnorm(Number,Mean,StandardDeviation)
    extracted fron rnorm.m by  Anders Holtsberg, 18-11-93
    Modified by Geoff Xu 2006
    2006-11-06 - C. Arthur - changed to optional kwargs
    """

    return (randn(n)*s + m)[0]

def convert_animation(conversion_cmd, animation_path, model_name):
    """convert_animation(conversion_cmd,animation_path,model_name):
    Converts the current simulation frames into an mpg animation
    Uses ImageMagick's 'convert' command
    """
    #os.mkdir
    print "Initiating Image to Animation Conversion"
    #time.sleep(.2)
    cmd_convert = conversion_cmd + ' -normalize -quality 80 -delay 10 ' + \
                  model_name + '*.png "' + animation_path + model_name + \
                  '.mpg"'
    if sys.platform == 'win32':
        cmd_remove = 'del '+ model_name+ '*.png'
    else:
        cmd_remove = 'rm '+ model_name+ '*.png'
    os.system(cmd_convert)
    os.system(cmd_remove)
    print "Completed Image to Animation Conversion"

def loadConfig(fileName, config={}):
    """loadConfig(filename, config={}):
    Returns a dictionary with keys of the form
    <section>.<option> and the values
    """
    config = config.copy()
    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(fileName)
    for sec in cp.sections():
        name = sec
        for opt in cp.options(sec):
            try:
                config[name + "." + opt] = cp.getint(sec, opt)
            except ValueError:
                try:
                    config[name + "." + opt] = cp.getfloat(sec, opt)
                except ValueError:
                    try:
                        config[name + "." + opt] = cp.getboolean(sec, opt)
                    except ValueError:
                        config[name + "." + opt] = cp.get(sec, opt)
    return config

def textread(config_file, nullValue=sys.maxint):
    """textread(config_file, nullValue=sys.maxint):
    Loads a csv file containing 'column' data into a dictionary of
    (numpy) arrays with keys labelled by 'fields'. Values of the
    dictionary are arrays containing the data in the columns in the
    file. See the comments following the code for an example layout
    of the configuration file.

    Usage:
    data = textread('configfile.ini')
    """
    configFile = loadConfig(open(config_file))
    dataFile = configFile['File.Filename']

    source = configFile['File.Source']
    delimiter = configFile[source+'.FieldDelimiter']
    columns = array(configFile[source+'.Columns'].split(delimiter))
    fields = array(configFile[source+'.Fields'].split(delimiter))
    headingLine = configFile[source+'.HeadingLine']

    # Check the data file is readable:
    try:
        f = open(dataFile, 'r')
    except IOError:
        raise IOError, ("Cannot open %s"%dataFile)
    else:
        numRecords = len(f.readlines())
        f.close()

    # A dictionary of column names to field names:
    mappings = {}
    for colname in columns:
        if colname not in fields:
            mappings[colname]=configFile['Mappings.' + colname]
        else:
            mappings[colname] = colname

    # Create a temporary dictionary of *lists*:
    temp = {}
    for colname in columns:
        temp[colname] = []

    data = csv.DictReader(file(dataFile), columns, restval=nullValue,
                          delimiter=str(delimiter))
    for record in data:
        for key in record.keys():
            if key in temp.keys():
                if configFile['Types.' + key] == 'int':
                    try:
                        temp[key].append(int(record[key]))
                    except ValueError:
                        temp[key].append(nullValue)

                elif configFile['Types.' + key] == 'float':
                    try:
                        temp[key].append(float(record[key]))
                    except ValueError:
                        temp[key].append(nullValue)
                else:
                    temp[key].append(record[key])

    output = {}
    for key in temp.keys():
        output[mappings[key]] = array(temp[key])

    return output

"""
Example configuration file textread.ini:
[File]
Filename=data.csv
# How the columns are arranged and named in the input file - this can change:
Columns=ind,unit,year,month,day,hour,minute,intensity,lat,lon,pressure,maxspeed,surfcode

[Data]
# What the data should look like - this doesn't change:
Fields=indicator,unit,year,month,day,hour,minute,intensity,lat,lon,pressure,max_speed,csurf_code

[Mappings]
# This sets out how 'Headings' map to 'Fields'
ind=indicator
maxspeed=max_speed
surfcode=csurf_code

"""

def savegrid(fname, X, lon, lat, delta, delimiter=' ', nodata=-9999,
             fmt='%.10e', coords='latlon'):
    """
    savegrid(filename, X, lon,lat,delta,delimiter=' ',nodata=-9999,fmt='%.10e'):
    Save formatted data to an ascii grid format file.
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    savegrid(filename, X, lon, lat, delta, delimiter=' ', nodata=-9999,
             fmt='%.10e, coords='latlon')
    """
    if pylab.is_string_like(fname):
        if fname.endswith('.gz'):
            import gzip
            fh = gzip.open(fname, 'wb')
        else:
            fh = file(fname, 'w')
    elif hasattr(fname, 'seek'):
        fh = fname
    else:
        raise ValueError('fname must be a string or file handle')

    if coords=='UTM':
        zone, xllcorner, yllcorner = LLtoUTM(lat.min(), lon.min())
        delta = convert(delta, "deg", "m")
    else:
        # Assume geographic coordinates
        xllcorner = lon.min()
        yllcorner = lat.min()

    fh.write('ncols         '+str(len(lon))+'\n')
    fh.write('nrows         '+str(len(lat))+'\n')
    fh.write('xllcorner     '+str(xllcorner)+'\n')
    fh.write('yllcorner     '+str(yllcorner)+'\n')
    fh.write('cellsize      '+str(delta)+'\n')
    fh.write('NODATA_value  '+str(nodata)+'\n')
    X = asarray(X)
    origShape = None
    if len(X.shape)==1:
        origShape = X.shape
        X.shape = len(X), 1
    for row in X:
        fh.write(delimiter.join([fmt%val for val in row]) + '\n')

    if origShape is not None:
        X.shape = origShape

def readgrid(filename, delimiter=None):
    """readgrid(filename, delimiter=None):
    Read formatted data from an ascii grid format file.
    Returns the longitude and latitude of the grid and the data values
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    longitude, latitude, data = readgrid(filename, [delimiter])
    """

    try:
        f = open(filename, 'r')
    except:
        raise IOError, "Cannot open %s"%filename
        return

    metadata = {}
    metadata["ncols"] = []
    metadata["nrows"] = []
    metadata["xllcorner"] = []
    metadata["yllcorner"] = []
    metadata["cellsize"] = []
    metadata["NODATA_value"] = []

    for i in range(0,6):
        line = f.readline()
        contents = line.split()
        label = contents[0]
        metadata[label] = float(contents[1])

    lon0 = metadata["xllcorner"]
    lon = numpy.array(range(int(metadata["ncols"])), dtype=float)
    lon = lon*metadata["cellsize"]+lon0
    lat0 = metadata["yllcorner"]
    lat = numpy.array(range(int(metadata["nrows"])), dtype=float)
    lat = lat*metadata["cellsize"]+lat0
    lat = numpy.flipud(lat)

    data = numpy.zeros([metadata["nrows"], metadata["ncols"]], dtype=float)

    for i in range(int(metadata["nrows"])):
        row = numpy.zeros([metadata["ncols"]], dtype=float)
        line = f.readline()
        for j, val in enumerate(line.split(delimiter)):
            value = float(val)
            if value == metadata["NODATA_value"]:
                value = Nan
            row[j] = value
        data[i,:] = row

    return lon, lat, data

class SampleGrid:
    """SampleGrid:
    Description: Sample data from an ascii grid file
    The files have 6 header lines describing the data, followed by the
    data in a gridded format.

    Headers:
    ncols
    nrows
    xllcorner
    yllcorner
    cellsize
    NODATA_value

    Usage:
    grid = SampleGrid(filename)
    value = grid.sampleGrid(lon, lat)

    Parameters: file (filename) containing gridded data in ascii grd
                format

    Members:
    sampleGrid(cLon, cLat): sample a value from the grid at the given
                            cLon, cLat
        At this time, does not interplolate from teh input grid to the
        given location

    Methods:
    sampleGrid(cLon, cLat): sample a value from the grid at the given
                            cLon, cLat
        At this time, does not interplolate from the input grid to the
        given location.

    Internal Methods:

    """

    def __init__(self, file):
        """
        Read in the data and ensure it's the right way around.
        """
        self.lon, self.lat, self.grid = readgrid(file)
        self.lat = numpy.flipud(self.lat)  # This is to ensure searchsorted works as expected.
        self.grid = numpy.flipud(self.grid)

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Sample data from an ascii grid file \
    The files have 6 header lines describing the data, followed by the \
    data in a gridded format.\
    \
    Headers:\
    ncols\
    nrows\
    xllcorner\
    yllcorner\
    cellsize\
    NODATA_value\
    \
    Usage:\
    grid = SampleGrid(filename)\
    value = grid.sampleGrid(lon, lat)'

    def sampleGrid(self, cLon, cLat):
        """sampleGrid(self, cLon, cLat):
        Grab nearest value to given location.
        """
        indi = self.lon.searchsorted(cLon)-1
        indj = self.lat.searchsorted(cLat)-1

        return self.grid[indj, indi]

def getUser():
    """getUser():
    Determine user name for incorporation into metadata strings
    """
    import getpass
    return getpass.getuser()

def getOS():
    """getOS():
    Determine operating system
    """
    import sys
    return sys.platform

def getDirSep():
    """getDirSep:
    Determine the directory separator for file manipulation
    """
    opsys = getOS()
    if opsys == "win32":
        dirsep = "\\"
    else:
        dirsep = "/"
    return dirsep

def checkPath(path):
    """checkPath(path):
    Check the path given exists and has the correct
    ending - i.e. path.endswith(dirsep) where dirsep is defined in
    the getDirSep function.
    """
    if path.endswith(getDirSep()):
        pass
    else:
        path += getDirSep()

    if os.path.isdir(path):
        pass
    else:
        try:
            os.mkdirs(path)
        except OSError:
            raise OSError, ("Cannot create %s"%path)

    return path

def getPath():
    """getPath():
    Prompt user for path to store data (if not already defined in the
    configuration file)
    """
