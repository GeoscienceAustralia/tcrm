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

 Title: files.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 04/23/08 2:38:PM
 Description: Helper functions for working with various files.

 Version: 75
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2009-04-29
 Modification: Added flStartLog

 Version: $Rev: 685 $
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2009-11-27 3:19:PM
 Modification: All exceptions are logged using the logging.exception() level
                (equivalent to logging.error)

 ModifiedBy: Nicholas Summons, nicholas.summons@ga.gov.au
 ModifiedDate: 2011-03-17
 Modification: 'logging.exception' replaced with 'logging.error' unless in try/except clause.
               This avoids error since logging.exception produces traceback that requires try/except error catching.

 ModifiedBy: Nicholas Summons, nicholas.summons@ga.gov.au
 ModifiedDate: 2011-06-22
 Modification: Replaced obsolete md5 package with hashlib

$Id: files.py 685 2012-03-29 04:22:32Z carthur $
"""
import os, sys, pdb, logging, traceback

try:
    import hashlib
    md5_constructor = hashlib.md5
except ImportError:
    import md5
    md5_constructor = md5.new
from time import time, ctime, localtime, strftime
import datetime
import ConfigParser
import numpy
import csv
import inspect

__version__ = '$Id: files.py 685 2012-03-29 04:22:32Z carthur $'

logger = logging.getLogger( )


def flModulePath(level=1):
    """
    Get the path of the module <level> levels above this function

    Input: level - level in the stack of the module calling this function
           (default = 1, function calling flModulePath)
    Output: path, basename and extension of the file containing the module
    Example: path, base, ext = flModulePath( )
             Calling flModulePath() from "/foo/bar/baz.py" produces the result
             "/foo/bar", "baz", ".py"
    """
    filename = os.path.realpath( sys._getframe( level ).f_code.co_filename )
    path, fname = os.path.split( filename )
    base, ext = os.path.splitext( fname )
    path = path.replace( os.path.sep, '/' )
    return path, base, ext

def flModuleName(level=1):
    """
    Get the name of the module <level> levels above this function

    Input: level - level in the stack of the module calling this function
           (default = 1, function calling flModuleName)
    Output: module name (as str)
    Example: mymodule = flModuleName( )
    """
    package = sys._getframe( level ).f_code.co_name
    return package

def flProgramVersion(level=len( inspect.stack( ) ) - 1 ):
    """
    Function to return the __version__ string from the parent
    program, where it is defined.
    If it is not defined, return an empty string.
    
    Input: level - level in the stack of the main script 
           (default = maximum level in the stack)
    Output: version string (defined as the __version__ global variable)

    Example: my_program_version = flProgramVersion( )
    """
    f = sys._getframe( level )
    if f.f_globals.has_key( '__version__' ):
        return f.f_globals[ '__version__' ]
    else:
        return ''

def flLoadFile(fname, comments='%', delimiter=None, skiprows=0):
    """
    Load ASCII data from fname into an array and return the array. The
    data must be regular, same number of values in every row.
    Code originated from pylab.py in matplotlib
    by John D. Hunter <jdhunter@ace.bsd.uhicago.edu>
    modified by Geoff Xu, 2005

    Input: filename, comment indicator, delimiter, skiprows
    Output: array of data, each column represents the columns in the
            input file
    Example: data = flLoadFile('/home/user/data.csv', comments='#',
                               delimiter=',', skiprows=0)
    """
    if os.path.isfile( fname ):
        pass
    else:
        logger.error( 'Filename %s is not a file'%( fname ) )
        raise IOError( 'Filename %s is not a file'%( fname ) )

    logger.debug( 'Loading %s'%( fname ) )
    if os.path.isfile( fname ):
        if fname.endswith( '.gz' ):
            import gzip
            try:
                fh = gzip.open( fname )
            except IOError:
                logger.exception( 'Cannot open %s '%( fname ) )
                flLogFatalError( traceback.format_exc( ).splitlines( ) )
        else:
            try:
                fh = open( fname )
            except IOError:
                logger.exception( 'Cannot open %s'%( fname ) )
                flLogFatalError( traceback.format_exc( ).splitlines( ) )
    elif hasattr( fname, 'seek' ):
        fh = fname
    else:
        logger.error( 'Filename must be a string or file handle' )
        raise IOError( 'Filename must be a string or file handle' ) 
        flLogFatalError( traceback.format_exc( ).splitlines( ) )
    data = []
    i = -1
    while True:
        line = fh.readline( )
        if len(line) == 0:
            break
        line = line.rstrip( '\n' )
        i = i + 1
        commentIdx = line.find( comments )
        if (commentIdx != -1):
            line = line[:commentIdx]
        line = line.strip( )
        if i < skiprows: continue
        if not len( line ): continue
        row = []
        for val in line.split( delimiter ):
            if val == 'NaN': val = Nan
            try:
                row.append( float( val ) )
            except ValueError:
                row.append( str( val ) )
        thisLen = len( row )
        data.append( row )
    data = numpy.array( data )
    r,c = data.shape
    if r == 1 or c == 1:
        data.shape = max( [r,c] ),
    fh.close( )
    return data

def flSaveFile(filename, data, header=None, delimiter=' ', fmt='%.18e'):
    """
    Save the data in X to file fname using fmt string to convert the
    data to strings, originated from pylab.py in matplotlib
    by John D. Hunter <jdhunter@ace.bsd.uhicago.edu>
    modified by Geoff Xu, 2006

    Input: filename, array of data, header line (default=None),
           delimiter (default=' '), format string (default='%.18e')
    Output: None
    Example: flSaveFile(file, data, header='List of fields',
                        delimiter=',', fmt='%5.2f')
    """
    try:
        directory, fname = os.path.split( filename )
    except AttributeError:
        logger.exception( 'Input filename is not a string' )
        flLogFatalError( traceback.format_exc( ).splitlines( ) )
    if not os.path.isdir( directory ):
        try:
            os.makedirs( directory )
        except:
            logger.exception( 'Cannot build path: %s'%( directory ) )
            flLogFatalError( traceback.format_exc( ).splitlines( ) )
    logger.debug( 'Saving data to %s'%( filename ) )

    if type( filename ) == str:
        if fname.endswith( '.gz' ):
            import gzip
            try:
                fh = gzip.open( filename, 'wb' )
            except IOError: 
                logger.exception( 'Cannot open %s'%( filename ) )
                flLogFatalError( traceback.format_exc( ).splitlines( ) )
        else:
            try:
                fh = open( filename, 'w' )
            except IOError:
                logger.exception( 'Cannot open %s'%( filename ) )
                flLogFatalError( traceback.format_exc( ).splitlines( ) )

    elif hasattr( filename, 'seek' ):
        fh = filename
    else:
        logger.error( 'Filename must be a string or file handle' )
        raise IOError( 'Filename must be a string or file handle' )

    if header:
        fh.write( '%' + header + '\n' )

    X = numpy.asarray( data )
    origShape = None
    if len( X.shape ) == 1:
        origShape = X.shape
        X.shape = len( X ), 1

    for row in X:
        try:
            if type(fmt) == list:
                fh.write( delimiter.join( [f%v for f,v in zip(fmt,row)] ) + '\n' )
            elif type(fmt) == str:
                fh.write( delimiter.join( [fmt%val for val in row] ) + '\n' )
            else:
                logger.exception( "Mismatch between format string and values in _write" )
                raise TypeError, "Mismatch between format string and values in _write"
        except ValueError:
            logger.exception( "Cannont write data to file" )
            flLogFatalError( traceback.format_exc( ).splitlines( ) )
            raise ValueError( "Cannont write data to file" )

    fh.close( )
    if origShape is not None:
        X.shape = origShape


def flGetStat(filename, CHUNK=2**16):
    """
    Get basic statistics of filename - namely directory, name (excluding
    base path), md5sum and the last modified date. Useful for checking
    if a file has previously been processed.

    Input: filename, chunk size (for md5sum calculation)
    Output: path, name, md5sum, modification date
    Example: dir, name, md5sum, moddate = flGetStat(filename)
    """
    try:
        fh = open( filename )
        fh.close( )
    except:
        logger.exception("Cannot open %s"%( filename ) )
        raise IOError("Cannot open %s"%( filename ) )

    try:
        directory, fname = os.path.split( filename )
    except:
        logger.exception( 'Input file is not a string' )
        raise TypeError( 'Input file is not a string' )

    try:
        si = os.stat( filename )
    except IOError:
        logger.exception( 'Input file is not a valid file: %s'%( filename ) )
        raise IOError( 'Input file is not a valid file: %s'%( filename ) )

    moddate = ctime( si.st_mtime )
    m = md5_constructor( )
    f = open( filename, 'rb' )

    while 1:
        chunk = f.read( CHUNK )
        if not chunk:
            break
        m.update( chunk )
    md5sum = m.hexdigest( )

    return directory, fname, md5sum, moddate

def flConfigFile(extension='.ini', prefix='', level=len(inspect.stack())):
    """
    Build a configuration filename (default extension .ini) based on the
    name and path of the function/module calling this function. Can also
    be useful for setting log file names automatically.
    If prefix is passed, this is preprended to the filename.

    Input: extension (default=.ini), prefix (default is empty)
    Output: Full path of calling function/module, with the source file's
    extension replaced with extension, and optionally prefix inserted
    after the last path separator
    Example: configFile = flConfigFile('.ini')
             Calling flConfigFile from /foo/bar/baz.py should 
             return /foo/bar/baz.ini
    """
    if level is None:
        level = len(inspect.stack())

    path, base, ext = flModulePath( level )
    config_file = os.path.join( path, prefix + base + extension )
    config_file = config_file.replace( os.path.sep, '/' )
    return config_file

def flStartLog(logFile, logLevel, verbose=False, datestamp=False, newlog=True):
    """
    Start logging to logFile all messages of logLevel and higher.
    Setting verbose=True will report all messages to STDOUT as well

    Input: logFile - full path to log file
           logLevel - string specifiying one of the standard Python logging
               levels ('NOTSET','DEBUG','INFO','WARNING','ERROR','CRITICAL')
           verbose - boolean: True will echo all logging calls to STDOUT
           datestamp - boolean: True will include a timestamp of the creation
                       time in the filename
           newlog - boolean: True will create a new log file each time this
                    function is called. False will append to the existing file.
    Output: None
    Example: flStartLog('/home/user/log/app.log','INFO',verbose=True)
    """
    if datestamp:
        b, e = os.path.splitext(logFile)
        curdate = datetime.datetime.now()
        curdatestr = curdate.strftime( '%Y%m%d%H%M' )
        # The lstrip on the extension is required as splitext leaves it on.
        logFile = "%s.%s.%s"%( b, curdatestr, e.lstrip( '.' ) )

    logDir = os.path.dirname( os.path.realpath( logFile ) )
    if not os.path.isdir( logDir ):
        try:
            os.makedirs( logDir )
        except OSError:
            # Unable to create the directory, so stick it in the 
            # current working directory:
            path, fname = os.path.split( logFile )
            logFile = os.path.join( os.getcwd( ), fname )

    if newlog:
        mode = 'w'
    else:
        mode = 'a'

    logging.basicConfig(level=getattr( logging, logLevel ),
                        format='%(asctime)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=logFile,
                        filemode=mode)
    logger = logging.getLogger( )

    if len( logger.handlers ) < 2:
        # Assume that the second handler is a StreamHandler for verbose
        # logging. This ensures we do not create multiple StreamHandler
        # instances that will *each* print to STDOUT
        if verbose and sys.stdout.isatty():
            # If set to true, all logging calls will also be printed to the
            # console (i.e. STDOUT)
            console = logging.StreamHandler( )
            console.setLevel( getattr( logging, logLevel ) )
            formatter = logging.Formatter( '%(asctime)s: %(message)s',
                                           '%Y-%m-%d %H:%M:%S', )
            console.setFormatter( formatter )
            logger.addHandler( console )

    logger.info( 'Started log file %s (detail level %s)'%( logFile, logLevel ) )
    logger.info( 'Running %s (pid %d)'%( sys.argv[0], os.getpid( ) ) )
    logger.info( 'Version %s'% ( flProgramVersion( ) ) )
    return logger

def flLogFatalError(tblines):
    """
    Log the error messages normally reported in a traceback
    so that all error messages can be caught.
    the input 'tblines' is created by calling
    traceback.format_exc().splitlines()
    """
    for line in tblines:
        logger.critical( line.lstrip( ) )
    sys.exit( )

def flModDate(filename, dateformat='%Y-%m-%d %H:%M:%S'):
    """
    Return the update date of the input file

    Input: filename - file name (full path)
           dateformat - (optional) format string for the date
    Output: File modification date/time as a string
    Example: modDate = flModDate( 'C:/foo/bar.csv' , dateformat='%Y-%m-%dT%H:%M:%S' )
    """
    try:
        si = os.stat( filename )
    except IOError:
        logger.exception( 'Input file is not a valid file: %s'%( filename ) )
        raise IOError( 'Input file is not a valid file: %s'%( filename ) )
    moddate = localtime( si.st_mtime )

    return strftime( dateformat, moddate )

def flSize(filename):
    """
    Return the size of the input file in bytes
    
    Input: filename - file name (full path)
    Output: file size in bytes
    Example: file_size = flSize( 'C:/foo/bar.csv' )
    """
    try:
        si = os.stat( filename )
    except WindowsError:
        logger.exception( 'Input file is not a valid file: %s'%( filename ) )
        raise IOError( 'Input file is not a valid file: %s'%( filename ) )
    except IOError:
        logger.exception( 'Input file is not a valid file: %s'%( filename ) )
        raise IOError( 'Input file is not a valid file: %s'%( filename ) )
    else:
        size = si.st_size

    return size


