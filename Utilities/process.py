"""
:mod:`process` -- control processing of files
=============================================

.. module:: process
    :synopsis: Provides functions to control processing of files. The
               module records details of each file processed in a text file,
               which can then be looked up at a later time to determin if a file
               has previously been processed. It uses the MD5 checksum value of
               the file as the primary indicator of modifications to a file.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import sys
import time
import logging

from files import flGetStat, flModDate

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

def pSetProcessedEntry(directory, filename, attribute, value):
    """
    Update the attribute of the given file with the given value.

    :param str directory: Base directory of the file.
    :param str filename: File name.
    :param str attribute: Name of the file attribute to be updated.
    :param str value: Attribute value.

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

    :param str directory: Path name of the file.
    :param str filename: File name to retrieve the details of.
    :param str attribute: Attribute to retrieve.

    :returns: Attribute value

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

    :param str datFileName: Name of a data file to read from.

    :returns: True if successfully read the data file, False
              otherwise.
    :rtype: bool

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
    """
    Write the various attributes of the given file to `gDatFile`

    :param str filename: Name of file that has been processed.

    :returns: True if the attributes of the file are successfully
              stored in gProcessedFiles and written to gDatFile, False
              otherwise.
    :rtype: bool

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
    Delete the existing data file - defined in the `gDatFile` variable
    (list of previously-processed files).

    :return: True if existing dat file successfully deleted,
             False otherwise
    :rtype: bool

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

    :param str directory: Base path of the file being checked.
    :param str filename: Name of the file.
    :param str attribute: Attribute name to be checked.
    :param str value: Value of the attribute to be tested.

    :return: True if the value matches that stored in gProcessedFiles,
             False otherwise.
    :rtype: boolean

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
    Set or get the archive directory. If setting the directory, its
    existence is checked and the directory is created.

    :param str archive_dir: Archive directory (if setting it).

    :return: The archive directory.
    :rtype: str

    :raises OSError: if the directory cannot be created

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
                raise OSError

    rc = g_archive_dir

    return rc


def pArchiveDateFormat( date_format=None ):
    """
    Set or get archive date format. Archived files can optionally have
    a date string inserted into the filename to retain all files with
    the same name, but different timestamps.

    :param str date_format: archive date format (if setting it)

    :return: archive date format
    :rtype: str

    """
    global g_archive_date_format
    if ( date_format ):
        g_archive_date_format = date_format
    rc = g_archive_date_format
    return rc


def pArchiveTimestamp( timestamp=False ):
    """
    Set or get archive timstamp flag. If the flag is `True`, then
    files that are to be archived will have a timestamp inserted into
    the file name.

    :param bool timestamp: `True` or `False` (if setting it)

    :return: The value of `g_archive_timestamp`
    :rtype: bool

    """
    global g_archive_timestamp
    if ( timestamp ):
        g_archive_timestamp = timestamp
    rc = g_archive_timestamp
    return rc

def pMoveFile( origin, destination ):
    """
    Move a single file to an archive directory.

    :param str origin: Full path of the file to be moved.
    :param str destination: Full path of the file destination.

    :return: `True` if the file is successfully moved, `False` otherwise.
    :rtype: bool

    """
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
    """
    Move the file to the archive directory (if specified), inserting a
    timestamp in the name.

    :param str filename: Full path of the file to be archived.

    :return: `True` if the file is successfully moved, `False` otherwise.
    :rtype: bool

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
