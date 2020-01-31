"""
:mod:`tsmultipliers` -- apply site-exposure multipliers to time series output
=============================================================================

.. module:: tsmultipliers
    :synopsis: Multiply the wind speed in a timeseries file by the
               appropriate multiplier values. Still a very rudimentary
               process.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import sys
import logging
import argparse
import traceback

from os.path import join as pjoin, dirname, realpath, isdir

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

import numpy as np
from Utilities.files import flStartLog
from Utilities.config import ConfigParser
from Utilities.dynarray import DynamicRecArray
from Utilities import shapefile
from Utilities import pathLocator


OUTPUTFMT = ['%s', '%s', '%9.5f', '%9.5f',
              '%6.2f', '%6.2f', '%6.2f', '%6.2f',
              '%7.2f']
INPUTFMT = ('|U16', '|U16', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')
INPUTNAMES = ('Station', 'Time', 'Longitude', 'Latitude', 'Speed', 'UU', 'VV',
              'Bearing', 'Pressure')

MINMAX_NAMES = ('Station', 'Time', 'Longitude', 'Latitude',
                'Speed', 'UU', 'VV', 'Bearing', 'Pressure')
MINMAX_TYPES = ['|U16', '|U16',  'f8', 'f8',  'f8', 'f8', 'f8', 'f8', 'f8']
MINMAX_FMT = ['%s', '%s', '%9.5f', '%9.5f',
              '%6.2f', '%6.2f', '%6.2f', '%6.2f',
              '%7.2f']

def tsmultiply(inputFile, multipliers, outputFile):
    """
    Apply multipliers to a single file. Values are combined then written
    back to the source file.

    :param str inputFile: Path to the input timeseries file. This will need
                          to contain the values of the three multipliers
                          (topography, terrain and shielding) for each of
                          eight directions.
    :param tuple multipliers: The eight combined multiplier values for the
                              location.
    :param str outputFile: Destination for the processed file.

    :returns: The records corresponding to the maximum (localised) wind speed
              and minimum pressure.

    """

    log.info("Processing {0}".format(inputFile))
    tsdata = np.genfromtxt(inputFile, dtype=INPUTFMT, names=INPUTNAMES,
                           delimiter=',', skip_header=1)
    tstep = tsdata['Time']
    lon = tsdata['Longitude']
    lat = tsdata['Latitude']
    gust = tsdata['Speed']
    uu = tsdata['UU']
    vv = tsdata['VV']
    bear = tsdata['Bearing']
    pressure = tsdata['Pressure']
    bear = 2*np.pi - (np.arctan2(-vv, -uu) - np.pi/2)
    bear = (180./np.pi) * np.mod(bear, 2.*np.pi)

    mcn, mcne, mce, mcse, mcs, mcsw, mcw, mcnw = multipliers

    # Apply multipliers:

    ii = np.where((bear >= 0.0) & (bear < 45.))
    gust[ii] *= (1./45.)*(mcn*(bear[ii] - 0.0) +
                          mcne*(45. - bear[ii]))
    ii = np.where((bear >= 45.0) & (bear < 90.))
    gust[ii] *= (1./45.)*(mcne*(bear[ii] - 45.0) +
                          mce*(90. - bear[ii]))
    ii = np.where((bear >= 90.0) & (bear < 135.))
    gust[ii] *= (1./45.)*(mce*(bear[ii] - 90.0) +
                          mcse*(135. - bear[ii]))
    ii = np.where((bear >= 135.0) & (bear < 180.))
    gust[ii] *= (1./45.)*(mcse*(bear[ii] - 135.0) +
                          mcs*(180. - bear[ii]))
    ii = np.where((bear >= 180.0) & (bear < 225.))
    gust[ii] *= (1./45.)*(mcs*(bear[ii] - 180.0) +
                          mcsw*(225. - bear[ii]))
    ii = np.where((bear >= 225.0) & (bear < 270.))
    gust[ii] *= (1./45.)*(mcsw*(bear[ii] - 225.0) +
                          mcw*(270. - bear[ii]))
    ii = np.where((bear >= 270.0) & (bear < 315.))
    gust[ii] *= (1./45.)*(mcw*(bear[ii] - 270.0) +
                          mcnw*(315. - bear[ii]))
    ii = np.where((bear >= 315.0) & (bear <= 360.))
    gust[ii] *= (1./45.)*(mcnw*(bear[ii] - 315.0) +
                          mcn*(360. - bear[ii]))


    ii = np.where(gust==0)
    bear[ii] = 0

    data = np.array([tsdata['Station'], tstep, lon, lat, gust, uu, vv, bear, pressure]).T

    maxidx = np.argmax(gust)
    minidx = np.argmin(pressure)
    maxdata = data[maxidx, :]
    mindata = data[minidx, :]

    header = 'Station,Time,Longitude,Latitude,Speed,UU,VV,Bearing,Pressure'
    np.savetxt(outputFile, data, fmt='%s', delimiter=',',
               header=header, comments='')

    return maxdata, mindata


def process_timeseries(config_file):
    """
    Process a set of timeseries files to include the multiplier values.

    The combined multiplier values are stored in a shape file as fields,
    and records are keyed by the same code that is used to select
    stations for sampling.

    :param str config_file: Path to a configuration file.

    """

    config = ConfigParser()
    config.read(config_file)

    stnFile = config.get('Input', 'LocationFile')
    key_name = config.get('Timeseries', 'StationID')
    inputPath = pjoin(config.get('Output', 'Path'),
                                  'process', 'timeseries')
    outputPath = pjoin(inputPath, 'local')

    if not isdir(outputPath):
        try:
            os.makedirs(outputPath)
        except OSError:
            raise

    log.info("Loading stations from %s"%stnFile)
    log.info("Timeseries data will be written into %s"%outputPath)

    directions = ['n', 'ne', 'e', 'se', 's', 'sw', 'w', 'nw']

    sf = shapefile.Reader(stnFile)
    field_names = [sf.fields[i][0] for i in range(1, len(sf.fields))]
    try:
        key_index = field_names.index(key_name)
    except ValueError:
        log.exception("{0} not a field in {1}".format(key_name, stnFile))
        raise

    min_data = DynamicRecArray(dtype={'names': MINMAX_NAMES,
                                      'formats': MINMAX_TYPES})
    max_data = DynamicRecArray(dtype={'names': MINMAX_NAMES,
                                      'formats': MINMAX_TYPES})

    records = sf.records()
    indexes = []
    for dir in directions:
        fieldname = 'm4_%s' % dir
        indexes.append(field_names.index(fieldname))

    for record in records:
        stnId = record[key_index]
        inputFile = pjoin(inputPath, 'ts.{0}.csv'.format(stnId))
        outputFile = pjoin(outputPath, 'ts.{0}.csv'.format(stnId))
        if os.path.isfile(inputFile):
            # Load multipliers for this location:
            mvals = [float(record[i]) for i in indexes]
            maxdata, mindata = tsmultiply(inputFile, tuple(mvals), outputFile)
            min_data.append(tuple(mindata))
            max_data.append(tuple(maxdata))

        else:
            log.debug("No timeseries file for {0}".format(stnId))

    # Save local minima/maxima
    maxfile = pjoin(outputPath, 'local_maxima.csv')
    minfile = pjoin(outputPath, 'local_minima.csv')
    maxheader = ('Station,Time,Longitude,Latitude,Speed,'
                        'UU,VV,Bearing,Pressure')
    np.savetxt(maxfile, max_data.data, fmt=MINMAX_FMT, delimiter=',',
               header=maxheader, comments='')
    np.savetxt(minfile, min_data.data, fmt=MINMAX_FMT, delimiter=',',
               header=maxheader, comments='')


def startup():
    """
    Parse command line arguments and call the :func:`main` function.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file',
                        help='Path to configuration file')
    parser.add_argument('-v', '--verbose', help='Verbose output',
                        action='store_true')
    parser.add_argument('-d', '--debug', help='Allow pdb traces',
                        action='store_true')
    args = parser.parse_args()
    config_file = args.config_file
    config = ConfigParser()
    config.read(config_file)

    rootdir = pathLocator.getRootDirectory()
    os.chdir(rootdir)

    logfile = config.get('Logging', 'LogFile')
    logdir = dirname(realpath(logfile))
    # If log file directory does not exist, create it
    if not isdir(logdir):
        try:
            os.makedirs(logdir)
        except OSError:
            logfile = pjoin(os.getcwd(), 'tsmultipliers.log')

    logLevel = config.get('Logging', 'LogLevel')
    verbose = config.getboolean('Logging', 'Verbose')
    datestamp = config.getboolean('Logging', 'Datestamp')
    debug = False

    if args.verbose:
        verbose = True

    if args.debug:
        debug = True

    flStartLog(logfile, logLevel, verbose, datestamp)

    if debug:
        process_timeseries(config_file)
    else:
        try:
            process_timeseries(config_file)
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())

    log.info("Completed {0}".format(sys.argv[0]))
if __name__ == "__main__":
    startup()
