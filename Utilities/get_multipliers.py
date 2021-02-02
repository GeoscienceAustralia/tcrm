## get_multipliers.py
'''
This script collects wind multipliers for the chosen tiles, converts them to 
geotiffs, merges them into a single tile for each wind direction, then combines 
the shielding, terrain and topography for each wind direction.

Make sure the following modules are loaded into the environment prior to running:
module load openmpi/1.8.4
module load python/2.7.6
module load python/2.7.6-matplotlib
module load geos
module load gdal/1.11.1-python
'''

from shutil import copyfile
import glob
import os
import argparse
import traceback
from os.path import join as pjoin, realpath, isdir, dirname
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from osgeo.gdal_array import BandReadAsArray, CopyDatasetInfo, BandWriteArray

import numpy as np
import numpy.ma as ma
import logging as log
from Utilities.config import ConfigParser
from Utilities import pathLocator
from Utilities.files import flStartLog

#class generateWindMultipliers():
'''
    Collect components of the wind multipliers, and generate m4 fields for the whole connected
    area
'''

def startup():
    """
    Parse the command line arguments and call the :func:`main`
    function.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file',
                        help='Path to configuration file')
    parser.add_argument('-v', '--verbose', help='Verbose output',
                        action='store_true')
    parser.add_argument('-d', '--debug', help='Allow pdb traces',
                        action='store_true')
    args = parser.parse_args()

    configFile = args.config_file
    config = ConfigParser()
    config.read(configFile)

    rootdir = pathLocator.getRootDirectory()
    os.chdir(rootdir)

    logfile = config.get('Logging','LogFile')
    logdir = dirname(realpath(logfile))

    # If log file directory does not exist, create it
    if not isdir(logdir):
        try:
            os.makedirs(logdir)
        except OSError:
            logfile = pjoin(os.getcwd(), 'tcrm.log')

    logLevel = config.get('Logging', 'LogLevel')
    verbose = config.getboolean('Logging', 'Verbose')
    datestamp = config.getboolean('Logging', 'Datestamp')
    debug = False

    if args.verbose:
        verbose = True

    if args.debug:
        debug = True

    flStartLog(logfile, logLevel, verbose, datestamp)
    # Switch off minor warning messages
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    if debug:
        main(configFile)
    else:
        try:
            main(configFile)
        except ImportError as e:
            log.critical("Missing module: {0}".format(e.strerror))
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())

def checkOutputFolders(OUTPUT_dir):
    '''
    Looks for the existance of the output folder/s specified in the config file,
    and if they don't exist, creates them.

    :param str OUTPUT_dir: path to the output directory
    '''
    dir_check = os.path.isdir(OUTPUT_dir)
    if dir_check == False:
        os.makedirs(OUTPUT_dir)
        log.info('Creating directories for outputs')

    # Assume if one is missing, they are all missing
    dir_check = os.path.isdir(OUTPUT_dir + '/shielding')
    if dir_check == False:
        os.makedirs(OUTPUT_dir + '/shielding')
        os.makedirs(OUTPUT_dir + '/terrain')           
        os.makedirs(OUTPUT_dir + '/topographic')

def copyTranslateMultipliers(configFile, type_mapping, output_path):
    '''  
    Copy wind multipliers from directory specified in the configuration file, to an
    output directory.
    Once the files have been copied, they are translated into Geotiffs

    :param str configFile: Path to configuration file
    :param str type_mapping: dict containing the three wind multiplier inputs
    '''
    config = ConfigParser()
    config.read(configFile)

    # Check for wind multiplier file path in config file
    if config.has_option('Input', 'Multipliers'):
        WMPath = config.get('Input', 'Multipliers')
        log.info('Using multiplier files from %s', WMPath)
    else:
        log.info('Using default multiplier files from /g/data/fj6/multipliers/')
        WMPath = '/g/data/fj6/multipliers/'

    tiles = config.get('Input', 'Tiles')
    tiles = [item.strip() for item in tiles.split(',')]
    log.info('Multipliers will be written out to %s', output_path)
    checkOutputFolders(output_path)
    
    for tile in tiles:
        for wm in type_mapping:
            var = type_mapping[wm]
            pathn = WMPath + wm + '/'
            log.debug('Copying files from %s', pathn)
            log.debug('Beginning to translate tiles into Geotiff')
            for file in glob.glob(pathn + tile + '*'):
                file_break = file.split('/')
                output = file_break[-1]
                output_name = wm + '/' + output
                copyfile(file, output_path + output_name)
                os.system('gdal_translate -of GTiff NETCDF:{0}{1}/{2}:{3} {4}{5}.tif'
                          .format(output_path, wm, output, var, output_path, 
                                  output_name[:-3])) # -3 to drop '.nc'
                log.debug('%s translated to Geotif', output_name[:-3])

def mergeWindMultipliers(type_mapping, dirns, output_path):
    '''
    Merge the Geotiff tiles together for each wind direction.

    :param str type_mapping: dict containing the three wind multiplier inputs
    :param str dirns: list of eight ordinal directions for wind
    :param str output_path: path to the output directory
    '''
    
    for wm in type_mapping:
        pathn = output_path + wm + '/'
        for dirn in dirns:
            common_dir = glob.glob(pathn + '*_' + dirn + '.tif')
            filelist = " ".join(common_dir)
            log.info('Merging %s %s files', wm, dirn)
            os.system('gdal_merge.py -of GTiff -ot float32 -n -9999 -a_nodata -9999 -o  {0}{1}_{2}.tif {3}'
                      .format(pathn, dirn, wm, filelist))

def combineDirections(dirns, output_path):
    '''
    Multiply geotiffs for terrain, topographic and shielding into a single geotiff.
    Output files are named "m4_" direction.

    :param str dirns: list of eight ordinal directions for wind
    :param str output_path: path to the output directory
    '''
    log.info('Multipliers will be written to %s', output_path)
    for dirn in dirns:
        log.info('working on %s', dirn)
        ds1 = gdal.Open('{0}terrain/{1}_terrain.tif'.format(output_path, dirn))
        ds2 = gdal.Open('{0}topographic/{1}_topographic.tif'.format(output_path, dirn))
        ds3 = gdal.Open('{0}shielding/{1}_shielding.tif'.format(output_path, dirn))
        band1 = ds1.GetRasterBand(1)
        band2 = ds2.GetRasterBand(1)
        band3 = ds3.GetRasterBand(1)
        data1 = BandReadAsArray(band1)
        data2 = BandReadAsArray(band2)
        data3 = BandReadAsArray(band3)
        m_data1 = ma.masked_values(data1, -9999)
        m_data2 = ma.masked_values(data2, -9999)
        m_data3 = ma.masked_values(data3, -9999)

        dataOut = m_data1 * m_data2 * m_data3

        driver = gdal.GetDriverByName("GTiff")
        log.info('Writing m4_%s.tiff', dirn)
        dsOut = driver.Create('{0}m4_{1}.tif'.format(output_path, dirn), ds1.RasterXSize, ds1.RasterYSize, 1, band1.DataType)
        CopyDatasetInfo(ds1,dsOut)
        bandOut=dsOut.GetRasterBand(1)
        bandOut.SetNoDataValue(-9999)
        BandWriteArray(bandOut, dataOut.data)

def main(configFile):
    """
    Main function for collecting and processing qind multipliers
    
    :param str configFile: Path to configuration file.
    """
    
    config = ConfigParser()
    config.read(configFile)
    output_path = config.get('Output', 'Path')
    
    type_mapping = {'shielding': 'Ms', 'terrain': 'Mz', 'topographic': 'Mt'}
    dirns = ['e', 'n', 'ne', 'nw', 's', 'se', 'sw', 'w']

    copyTranslateMultipliers(configFile, type_mapping, output_path)
    mergeWindMultipliers(type_mapping, dirns, output_path)
    combineDirections(dirns, output_path)

if __name__ == "__main__":
    startup()
