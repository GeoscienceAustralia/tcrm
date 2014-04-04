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

 Title: process.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 11/14/07 4:30:PM
 Description: Provides functions to control processing of a large number
 of files. It records details of each file processed in a text file
 which can then be looked up at a later time to check if it has been
 processed.  It uses the MD5 checksum value of the file as the primary
 indicator of modifications to a file.
 Constraints: Python versions <= 2.4 should use the md5 module. Later
 versions should use the hashlib module.
 SeeAlso: files.py

 Version :$Rev: 642 $

 $Id: process.py 642 2012-02-21 07:54:04Z nsummons $
"""
import os, sys, time, pdb, logging

from files import flGetStat, flModDate
__version__ = '$Id: process.py 642 2012-02-21 07:54:04Z nsummons $'
global gDatFile
global gProcessedFiles
global g_archive_dir
global g_archive_timestamp
global g_archive_date_format
logger = logging.getLogger( )

gProcessedFiles = { }
g_archive_dir = ''
g_archive_date_format = '%Y%m%d%H%M'
g_archive_timestamp = True

def pSetProcessedEntry( directory, filename, attribute, value ):
    """
    Update the attribute of the given file with the given value.
    Input: directory, filename, attribute, value
    Output: None (gProcessedFiles is updated)
    Example: pSetProcessedEntry(directory, filename, attribute, value)
    """
    global gDatFile
    global gProcessedFiles
    if directory in gProcessedFiles:
        if filename in gProcessedFiles[directory]:
            if attribute in gProcessedFiles[directory][filename]:
                gProcessedFiles[directory][filename][attribute] = value
            else:
                gProcessedFiles[directory][filename].update( {attribute:value} )
        else:
            gProcessedFiles[directory].update( {filename:{attribute:value}} )
    else:
        gProcessedFiles.update( {directory:{filename:{attribute:value}}} )

def pGetProcessedEntry( directory, filename, attribute ):
    """
    Retrieve the value of an attribute for a file from the
    gProcessedFiles dictionary.
    Input: directory, name and attribute of the file to retrieve
    Output: Value (if it exists) of the attribute
    Example: value = pGetProcessedEntry(directory, filename, attribute)
    """
    global gDatFile
    global gProcessedFiles
    try:
        value = gProcessedFiles[directory][filename][attribute]
        rc = value
    except KeyError:
        rc = 0
    return rc


def pGetProcessedFiles( datFileName=None ):
    """
    Retrieve a list of processed files from a dat file. This will also
    set the global gDatFile.
    Input: datFileName
    Output: True if successfully read the dat file, False otherwise
    Example: rc = pGetProcessedFiles(gDatFile)
    """
    global gDatFile
    global gProcessedFiles
    rc = 0
    if datFileName:
        gDatFile = datFileName
        try:
            fh = open( datFileName )

        except IOError:
            logger.warn( "Couldn't open dat file %s"%( datFileName ) )
            return rc
        else:
            logger.debug( "Getting previously-processed files from %s"%( datFileName ) )
            for line in fh:
                line.rstrip( '\n' )
                directory, filename, moddate, md5sum = line.split('|')
                pSetProcessedEntry( directory, filename, 'moddate',
                                   moddate.rstrip( '\n' ) )
                pSetProcessedEntry( directory, filename, 'md5sum',
                                   md5sum.rstrip( '\n' ) )
            rc = 1
            fh.close()

    else:
        logger.info("No dat file name provided - all files will be processed")

        return rc
    return rc

def pWriteProcessedFile( filename ):
    """pWriteProcessedFile(filename):
    Write the various attributes of the given file to gDatFile
    Input: filename that has been processed
    Output: True if the attributes of the file are successfully stored
            in gProcessedFiles and written to gDatFile, False otherwise.
    Example: rc = pWriteProcessedFile(filename)
    """
    global gDatFile
    global gProcessedFiles
    rc = 0
    if gDatFile:
        directory,fname,md5sum,moddate = flGetStat( filename )
        try:
            fh = open( gDatFile,'a' )
        except:
            logger.info( "Cannot open %s"%( gDatFile ) )

        else:
            pSetProcessedEntry( directory, fname, 'md5sum', md5sum )
            pSetProcessedEntry( directory, fname, 'moddate', moddate )
            fh.write( '|'.join( [directory, fname, moddate, md5sum] ) + '\n' )
            fh.close()
            rc = 1
    else:
        logger.warn( "Dat file name not provided. Can't record %s as processed."%( filename ) )

    return rc

def pDeleteDatFile( ):
    """
    Delete the existing dat file - defined in the gDatFile variable
    (list of previously-processed files).
    Input: None
    Output: True if existing dat file successfully deleted,
            False otherwise
    Example: rc = pDeleteDatFile()
    """
    global gDatFile
    global gProcessedFiles
    rc = 0
    if os.unlink( gDatFile ):
        rc = 1
    else:
        logger.warn( "Cannot remove dat file %s"%( gDatFile ) )
    return rc

def pAlreadyProcessed( directory, filename, attribute, value ):
    """
    Determine if a file has already been processed (i.e. it is stored in
    gProcessedFiles)
    Input: directory, filename, attribute and corresponding value of a
           file.
    Output: True if the value matches that stored in gProcessedFiles,
            False otherwise
    Example: rc = pAlreadyProcessed(directory, filename, attribute,
                                    value)
    """
    global gDatFile
    global gProcessedFiles
    rc = False
    if pGetProcessedEntry( directory, filename, attribute ) == value:
        rc = True
    else:
        rc = False
    return rc

def pArchiveDir( archive_dir=None ):
    """
    pArchiveDir
    Set or get the archive directory.
    Input: archive directory (if setting it)
    Output: original archive directory
    """
    global g_archive_dir
    
    if archive_dir:
        g_archive_dir = archive_dir
        g_archive_dir.rstrip("/\\")
        if (not os.path.isdir(g_archive_dir) ):
            try:
                os.makedirs(g_archive_dir)
            except:
                logger.exception( "Cannot create %s"%( g_archive_dir ) )
    rc = g_archive_dir
    
    return rc


def pArchiveDateFormat( date_format=None ):
    """pArchiveDateFormat:
    Set or get archive date format.
    Input: archive date format (if setting it)
    Output: original archive date format
    """
    global g_archive_date_format
    if ( date_format ):
        g_archive_date_format = date_format
    rc = g_archive_date_format
    return rc


def pArchiveTimestamp( timestamp=False ):
    """pArchiveTimestamp
    Set or get archive timstamp flag.
    Input: true or false (if setting it)
    Output: original setting
    """
    global g_archive_timestamp
    if ( timestamp ):
        g_archive_timestamp = timestamp
    rc = g_archive_timestamp
    return rc

def pMoveFile( origin, destination ):
    """pMoveFile: move a single file"""
    try:
        os.rename( origin, destination )
    except:
        logger.warn( "Error moving %s to %s"%( origin, destination ) )
        raise
        rc = 0
    else:
        logger.debug( "%s moved to %s"%( origin, destination ) )
        rc = 1

    return rc

def pArchiveFile( filename ):
    """pArchiveFile
    Move the file to the archive directory (if specified), inserting a
    timestamp in the name.
    """
    path, ext = os.path.splitext( filename )
    path, base = os.path.split( path )
    archive_dir = pArchiveDir( )
    ext = ext.lstrip('.')
    if ( archive_dir ):
        if os.path.isdir( archive_dir ):
            pass
        else:
            try:
                os.makedirs( archive_dir )
            except:
                logger.critcal( "Cannot create %s"%( archive_dir ) )
                raise

    if ( pArchiveTimestamp( ) ):
        archive_date = flModDate( filename, g_archive_date_format )
        archive_file_name = os.path.join( archive_dir, "%s.%s.%s"%( base, archive_date, ext ) )
    else:
        archive_file_name = os.path.join( archive_dir, "%s.%s"%( base, ext ) )

    rc = pMoveFile( filename, archive_file_name )
    return rc
