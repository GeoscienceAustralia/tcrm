"""
:mod:`processMultipliers` -- combine multipliers with wind speed data
=====================================================================

.. module:: processMultipliers
    :synopsis: Process a wind field file and combine with directional
               site-exposure multipliers to evaluate local wind
               speed.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

Combine the regional wind speed data with the site-exposure
multipliers, taking care to select the directional multiplier that
corresponds to the direction of the maximum wind speed.

This version assumes that the site-exposure multipliers are given as a
combined value (i.e. ``Ms * Mz * Mh``), and the files are ERDAS
Imagine-format files ('*.img'). Further, the files are assumed to have
the file name ``m4_<dir>.img``, where <dir> is the direction (n, ne, e,
se, s, sw, w or nw).

Requires the Python GDAL bindings, Numpy, netCDF4 and the :mod:`files`
and :mod:`config` modules from TCRM. It assumes :mod:`Utilities` can
be found in the ``PYTHONPATH`` directory.
    Make sure the following modules are loaded into the environment prior to running:
    module load openmpi/1.8.4
    module load python/2.7.6
    module load python/2.7.6-matplotlib
    module load pypar/26Feb15-2.7.6-1.8.4
    module load geos
    module load gdal/1.11.1-python

"""

from shutil import copyfile, rmtree
import glob
import os
from os.path import join as pjoin, dirname, realpath, isdir, splitext
import time
import logging as log
import argparse
import traceback
from functools import wraps, reduce

from Utilities.files import flStartLog
from Utilities.config import ConfigParser
from Utilities import pathLocator

import numpy as np
import numpy.ma as ma

from osgeo import osr, gdal, gdalconst
from osgeo.gdal_array import BandReadAsArray, CopyDatasetInfo, BandWriteArray
from gdal import *

from netCDF4 import Dataset

gdal.UseExceptions()

syn_indices = {
    0: {'dir': 'n', 'min': 0., 'max': 22.5, 'fill': 0},
    1: {'dir': 'ne', 'min': 22.5, 'max': 67.5, 'fill': 1},
    2: {'dir': 'e', 'min': 67.5, 'max': 112.5, 'fill': 2},
    3: {'dir': 'se', 'min': 112.5, 'max': 157.5, 'fill': 3},
    4: {'dir': 's', 'min': 157.5, 'max': 202.5, 'fill': 4},
    5: {'dir': 'sw', 'min': 202.5, 'max': 247.5, 'fill': 5},
    6: {'dir': 'w', 'min': 247.5, 'max': 292.5, 'fill': 6},
    7: {'dir': 'nw', 'min': 292.5, 'max': 337.5, 'fill': 7},
    8: {'dir': 'n', 'min': 337.5, 'max': 360., 'fill': 0},
    9: {'dir': 'max', 'fill': 0}}

def timer(f):
    """
    Basic timing functions for entire process
    """
    @wraps(f)
    def wrap(*args, **kwargs):
        t1 = time.time()
        res = f(*args, **kwargs)

        tottime = time.time() - t1
        msg = "%02d:%02d:%02d " % \
            reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
                   [(tottime,), 60, 60])

        log.info("Time for {0}: {1}".format(f.__name__, msg))
        return res

    return wrap

class getMultipliers():
    '''
    This script collects wind multipliers for the chosen tiles, converts them to
    geotiffs, merges them into a single tile for each wind direction, then combines
    the shielding, terrain and topography for each wind direction.
    '''

    def __init__(self, configFile):

        config = ConfigParser()
        config.read(configFile)

        # Check for wind multiplier file path in config file
        if config.has_option('Input', 'RawMultipliers'):
            self.WMPath = config.get('Input', 'RawMultipliers')
            log.info('Using multiplier files from {0}'.format(self.WMPath))
        else:
            log.info('Using default multiplier files from /g/data/fj6/multipliers/')
            self.WMPath = '/g/data/fj6/multipliers/'

    def checkOutputFolders(self, working_dir, type_mapping):
        '''
        Looks for the existance of the output folder/s specified in the config file,
        and if they don't exist, creates them.

        :param str working_dir: path to the output directory
        :param dict type_mapping: dict of shielding, terrain, topographic
        '''
        dir_check = os.path.isdir(working_dir)
        if dir_check == False:
            os.makedirs(working_dir)
            log.info('Creating directories for outputs')
        else:
            log.info('Using existing output directories')

        # Assume if one is missing, they are all missing
        dir_check = os.path.isdir(working_dir + '/shielding')
        if dir_check is False:
            for wm in type_mapping:
                os.makedirs('{0}/{1}'.format(working_dir, wm))
        else:
            #Else clean up the folders to make sure we don't mix up old and new
            log.info('Deleting the contents of the shielding, terrain and topographic '
                     'folders in {0} to avoid using old files'.format(working_dir))
            for wm in type_mapping:
                files = glob.glob('{0}/{1}/*'.format(working_dir, wm))
                for f in files:
                    os.remove(f)

    def copyTranslateMultipliers(self, tiles, configFile, type_mapping, working_dir):
        '''
        Copy wind multipliers from directory specified in the configuration file, to an
        output directory.
        Once the files have been copied, they are translated into Geotiffs

        :param str configFile: Path to configuration file
        :param str type_mapping: dict containing the three wind multiplier inputs
        :param str working_dir: path to the output directory
        '''
        for tile in tiles:
            log.info('tile = {0}'.format(tile))
            for wm in type_mapping:
                log.debug('type mapping = {0}'.format(wm))
                var = type_mapping[wm]
                pathn = self.WMPath + wm + '/'
                log.debug('Copying files from {0}'.format(pathn))
                log.debug('Beginning to translate tiles into Geotiff')
                for file in glob.glob(pathn + tile + '*'):
                    file_break = file.split('/')
                    output = file_break[-1]
                    output_name = wm + '/' + output
                    copyfile(file, working_dir + output_name)
                    os.system('gdal_translate -a_srs EPSG:4326 -of GTiff '
                              'NETCDF:{0}{1}/{2}:{3} {4}{5}.tif'
                              .format(working_dir, wm, output, var, working_dir,
                                      output_name[:-3])) # -3 to drop '.nc'
                    
                    log.info('%s translated to Geotif', output_name[:-3])

    def mergeWindMultipliers(self, type_mapping, dirns, working_dir):
        '''
        Merge the Geotiff tiles together for each wind direction.

        :param str type_mapping: dict containing the three wind multiplier inputs
        :param str dirns: list of eight ordinal directions for wind
        :param str working_dir: path to the output directory

        '''

        for wm in type_mapping:
            pathn = working_dir + wm + '/'
            for dirn in dirns:
                common_dir = glob.glob(pathn + '*_' + dirn + '.tif')
                filelist = " ".join(common_dir)
                log.info('Merging {0} {1} files'.format(wm, dirn))
                os.system('gdal_merge.py -of GTiff -ot float32 -n -9999 '
                          '-a_nodata -9999 -o  {0}{1}_{2}.tif {3}'
                          .format(pathn, dirn, wm, filelist))
        log.debug('Finished merging the wind multiplier tiles')

    def combineDirections(self, dirns, working_dir):
        '''
        Multiply geotiffs for terrain, topographic and shielding into a single geotiff.
        Output files are named "m4_" direction.

        :param str dirns: list of eight ordinal directions for wind
        :param str output_path: path to the output directory
        '''
        log.info('Multipliers will be written to {0}'.format(working_dir))
        for dirn in dirns:
            log.info('working on %s', dirn)
            ds1 = gdal.Open('{0}terrain/{1}_terrain.tif'.format(working_dir, dirn))
            ds2 = gdal.Open('{0}topographic/{1}_topographic.tif'.format(working_dir, dirn))
            ds3 = gdal.Open('{0}shielding/{1}_shielding.tif'.format(working_dir, dirn))
            band1 = ds1.GetRasterBand(1)
            band2 = ds2.GetRasterBand(1)
            band3 = ds3.GetRasterBand(1)
            data1 = BandReadAsArray(band1)
            data2 = BandReadAsArray(band2)
            data3 = BandReadAsArray(band3)
            m_data1 = ma.masked_values(data1, -9999)
            m_data2 = ma.masked_values(data2, -9999)
            m_data3 = ma.masked_values(data3, -9999)
            
            log.debug('Size of data arrays: terrain = {0}, topographic = {1}, shielding = {2}'
                      .format(m_data1.shape, m_data2.shape, m_data3.shape))

            dataOut = m_data1 * m_data2 * m_data3

            driver = gdal.GetDriverByName("GTiff")
            log.info('Writing m4_{0}.tiff'.format(dirn))
            dsOut = driver.Create('{0}m4_{1}.tif'.format(working_dir, dirn), ds1.RasterXSize,
                                  ds1.RasterYSize, 1, band1.DataType)
            CopyDatasetInfo(ds1, dsOut)
            bandOut = dsOut.GetRasterBand(1)
            bandOut.SetNoDataValue(-9999)
            BandWriteArray(bandOut, dataOut.data)

def generate_syn_mult_img(tl_x, tl_y, delta, dir_path, shape,
                          indices=syn_indices,
                          every_fill=None):
    """

    :param x_tl:  top left x
    :param y_tl: top left y
    :param shape: The array shape
    :return:

    """
    tl_y = np.asarray([tl_y]) # top left y
    tl_x = np.asarray([tl_x]) # top left x

    multiplier_values = np.zeros(shape)

    for value in indices.items():
        if every_fill is None:
            fill = value[1]['fill']
        else:
            fill = every_fill
        multiplier_values.fill(fill)
        img_name = 'm4_' + value[1]['dir'] + '.tif'
        file_path = pjoin(dir_path, img_name)
        createRaster(multiplier_values, tl_x, tl_y,
                     delta, -delta,
                     filename=file_path)

def createRaster(array, x, y, dx, dy, epsg = 4326, filename=None, nodata=-9999):
    """
    Create an in-memory raster for processing. By default, we assume
    the input array is in geographic coordinates, using WGS84 spatial
    reference system.

    :param array: Data array to be stored.
    :type  array: :class:`numpy.ndarray`
    :param x: x-coordinate of the array values.
    :type  x: :class:`numpy.ndarray`
    :param y: y-coordinate of the array values - should be a negative
              value.
    :type  y: :class:`numpy.ndarray`
    :param float dx: Pixel size in x-direction.
    :param float dy: Pixel size in y-direction.
    :param int epsg: EPSG code of the spatial reference system
                     of the input array (default=4326, WGS84)
    :param filename: Optional path
     to store the data in.
    :type  filename: str or None

    """
    if filename:
        log.debug("Creating raster: {0}".format(filename))
    else:
        log.debug("Creating in-memory raster")
    rows, cols = array.shape
    originX, originY = x[0], y[-1]
    if filename:
        _, ext = splitext(filename)
    if filename and ext == '.tif':
        driver = gdal.GetDriverByName('GTiff')
        tempRaster = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    elif filename and ext == '.img':
        driver = gdal.GetDriverByName('HFA')
        tempRaster = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    else:
        driver = gdal.GetDriverByName('MEM')
        tempRaster = driver.Create('', cols, rows, 1, gdal.GDT_Float32)

    tempRaster.SetGeoTransform((originX, dx, 0,
                                originY, 0, dy))
    tempBand = tempRaster.GetRasterBand(1)
    tempBand.WriteArray(array[::int(np.sign(dy) * 1)])
    tempBand.SetNoDataValue(nodata)
    tempRasterSRS = osr.SpatialReference()
    tempRasterSRS.ImportFromEPSG(epsg)
    tempRaster.SetProjection(tempRasterSRS.ExportToWkt())

    log.debug("Spatial reference system is:")
    log.debug(tempRasterSRS.ExportToWkt())
    tempBand.FlushCache()
    return tempRaster


def loadRasterFile(raster_file, fill_value=1):
    """
    Load a raster file and return the data as a :class:`numpy.ndarray`.
    No prorjection information is returned, just the actual data as an
    array.

    :param str raster_file: Path to the raster file to load.
    :param fill_value: Value to replace `nodata` values with (default=1).
    :returns: 2-d array of the data values.
    :rtype: :class:`numpy.ndarray`

    """

    log.debug("Loading raster data from {0} into array".format(raster_file))
    ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()

    nodata = band.GetNoDataValue()

    if nodata is not None:
        np.putmask(data, data == nodata, fill_value)

    del ds
    return data


def loadRasterFileBandLonLat(raster_file, fill_value=1):
    """
    Load a raster file and return the data as a :class:`numpy.ndarray`.
    No prorjection information is returned, just the actual data as an
    array.

    :param str raster_file: Path to the raster file to load.
    :param fill_value: Value to replace `nodata` values with (default=1).
    :returns: 2-d array of the data values.
    :rtype: :class:`numpy.ndarray`

    """

    log.debug("Loading raster data from {0} into array".format(raster_file))
    ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    nodata = band.GetNoDataValue()

    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*gt[5]
    maxx = gt[0] + width*gt[1] + height*gt[2]
    maxy = gt[3]


    if nodata is not None:
        np.putmask(data, data == nodata, fill_value)

    del ds
    return minx, miny, maxx, maxy, data


def calculateBearing(uu, vv):
    """
    Calculate the wind direction from the u (eastward) and v
    (northward) components of the wind speed.

    :param uu: :class:`numpy.ndarray` of eastward values
    :param vv: :class:`numpy.ndarray` of northward values.

    :returns: Direction the wind is coming from, zero northwards, positive
              clockwise, in degrees.  This is the direction the wind is
              blwing from.  So if thw wind is speeding north, the
              bearing is 180 deg.
    :rtype: :class:`numpy.ndarray`

    """
    bearing = 2 * np.pi - (np.arctan2(-vv, -uu) - np.pi / 2)
    bearing = (180. / np.pi) * np.mod(bearing, 2. * np.pi)
    return bearing


@timer
def reprojectDataset(src_file, match_filename, dst_filename,
                     resampling_method = gdalconst.GRA_Bilinear, 
                     match_projection = None):
    """
    Reproject a source dataset to match the projection of another
    dataset and save the projected dataset to a new file.

    :param src_filename: Filename of the source raster dataset, or an
                         open :class:`gdal.Dataset`
    :param match_filename: Filename of the dataset to match to, or an
                           open :class:`gdal.Dataset`
    :param str dst_filename: Destination filename.
    :param resampling_method: Resampling method. Default is bilinear
                              interpolation.

    """

    log.debug("Reprojecting {0}".format(repr(src_file)))
    log.debug("Match raster: {0}".format(repr(match_filename)))
    log.debug("Output raster: {0}".format(dst_filename))

    if isinstance(src_file, str):
        src = gdal.Open(src_file, gdal.GA_ReadOnly)
    else:
        src = src_file
    srcBand = src.GetRasterBand(1)
    srcBand.SetNoDataValue(-9999)
    src_proj = src.GetProjection()

    # We want a section of source that matches this:
    if isinstance(match_filename, str):
        match_ds = gdal.Open(match_filename, gdal.GA_ReadOnly)
    else:
        match_ds = match_filename
    matchBand = match_ds.GetRasterBand(1)
    matchBand.SetNoDataValue(-9999)

    if match_projection:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(match_projection)
        match_proj = srs.ExportToWkt()
    else:
        match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize

    # Output / destination
    drv = gdal.GetDriverByName('GTiff')
    dst = drv.Create(dst_filename, wide, high, 1, gdal.GDT_Float32)
    dst.SetGeoTransform(match_geotrans)
    dst.SetProjection(match_proj)
    dstBand = dst.GetRasterBand(1)
    dstBand.SetNoDataValue(-9999)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, resampling_method)

    del dst  # Flush
    if isinstance(match_filename, str):
        del match_ds
    if isinstance(src_file, str):
        del src

    return

@timer
def processMult(wspd, uu, vv, lon, lat, working_dir, m4_max_file = 'm4_ne.tif'):
    """

    The lat and lon values are the top left corners of the cells
    The speed arrays are in bottom to top format

    :param wspd: The gust speed
    :param uu: x component of the wind speed
    :param vv: y component of the wind speed
    :param lon: list of raster longitude values
    :param lat:  list of raster latitude values
    :param m4_max_file: Multiplier file used for reprojection.
    :param working_dir: The working output directory
    :return:
    """

    # This gives different bearing values
    # thank the bearings in the result tuple
    bearing = calculateBearing(uu, vv)

    delta = lon[1] - lon[0]

    log.debug('Create rasters from the netcdf gust file variables')
    wind_raster_file = pjoin(working_dir, 'region_wind.tif')
    wind_raster = createRaster(np.flipud(wspd), lon, lat, delta, delta,
                               filename = wind_raster_file)
    bear_raster = createRaster(np.flipud(bearing), lon, lat, delta, delta)
    uu_raster = createRaster(np.flipud(uu), lon, lat, delta, delta)
    vv_raster = createRaster(np.flipud(vv), lon, lat, delta, delta)

    log.info("Reprojecting regional wind data")
    wind_prj_file = pjoin(working_dir, 'gust_prj.tif')
    bear_prj_file = pjoin(working_dir, 'bear_prj.tif')
    uu_prj_file = pjoin(working_dir, 'uu_prj.tif')
    vv_prj_file = pjoin(working_dir, 'vv_prj.tif')

    log.info('Reproject the wind speed and bearing data to match '
             'the input wind multipliers')
    m4_max_file = pjoin(working_dir, m4_max_file)

    reprojectDataset(wind_raster, m4_max_file, wind_prj_file)
    reprojectDataset(bear_raster, m4_max_file, bear_prj_file,
                     resampling_method = gdalconst.GRA_NearestNeighbour)
    reprojectDataset(uu_raster, m4_max_file, uu_prj_file,
                     resampling_method = gdalconst.GRA_NearestNeighbour)
    reprojectDataset(vv_raster, m4_max_file, vv_prj_file,
                     resampling_method = gdalconst.GRA_NearestNeighbour)

    wind_prj_ds = gdal.Open(wind_prj_file, gdal.GA_ReadOnly)
    wind_prj = wind_prj_ds.GetRasterBand(1)
    bear_prj_ds = gdal.Open(bear_prj_file, gdal.GA_ReadOnly)
    bear_prj = bear_prj_ds.GetRasterBand(1)
    uu_prj_ds = gdal.Open(uu_prj_file, gdal.GA_ReadOnly)
    uu_prj = uu_prj_ds.GetRasterBand(1)
    vv_prj_ds = gdal.Open(vv_prj_file, gdal.GA_ReadOnly)
    vv_prj = vv_prj_ds.GetRasterBand(1)
    wind_proj = wind_prj_ds.GetProjection()
    wind_geot = wind_prj_ds.GetGeoTransform()

    wind_data = wind_prj.ReadAsArray()
    bear_data = bear_prj.ReadAsArray()
    uu_data = uu_prj.ReadAsArray()
    vv_data = vv_prj.ReadAsArray()
    bearing = calculateBearing(uu_data, vv_data)

    # The local wind speed array:
    local = np.zeros(wind_data.shape, dtype = 'float32')

    indices = {
        0: {'dir': 'n', 'min': 0., 'max': 22.5},
        1: {'dir': 'ne', 'min': 22.5, 'max': 67.5},
        2: {'dir': 'e', 'min': 67.5, 'max': 112.5},
        3: {'dir': 'se', 'min': 112.5, 'max': 157.5},
        4: {'dir': 's', 'min': 157.5, 'max': 202.5},
        5: {'dir': 'sw', 'min': 202.5, 'max': 247.5},
        6: {'dir': 'w', 'min': 247.5, 'max': 292.5},
        7: {'dir': 'nw', 'min': 292.5, 'max': 337.5},
        8: {'dir': 'n', 'min': 337.5, 'max': 360.}
    }
    log.info("Processing all directions")
    for i in list(indices.keys()):
        dn = indices[i]['dir']
        log.info("Processing {0}".format(dn))
        m4_file = pjoin(working_dir, 'm4_{0}.tif'.format(dn.lower()))
        m4 = loadRasterFile(m4_file)
        idx = np.where((bear_data >= indices[i]['min']) &
                       (bear_data < indices[i]['max']))

        local[idx] = wind_data[idx] * m4[idx]
    rows, cols = local.shape
    output_file = pjoin(working_dir, 'local_wind.tif')
    log.info("Creating output file: {0}".format(output_file))
    # Save the local wind field to a raster file with the SRS of the
    # multipliers
    drv = gdal.GetDriverByName("GTiff")
    dst_ds = drv.Create(output_file, cols, rows, 1,
                        gdal.GDT_Float32, ['BIGTIFF=NO', 'SPARSE_OK=TRUE'])
    dst_ds.SetGeoTransform(wind_geot)
    dst_ds.SetProjection(wind_proj)
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue(-9999)
    dst_band.WriteArray(local)

    # dst_band.FlushCache()

    del dst_ds
    log.info("Completed")

    return output_file

class run():
    
    def __init__(self):

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

        self.configFile = args.config_file
        config = ConfigParser()
        config.read(self.configFile)


        logfile = config.get('Logging', 'LogFile')
        logdir = dirname(realpath(logfile))

        # If log file directory does not exist, create it
        if not isdir(logdir):
            try:
                os.makedirs(logdir)
            except OSError:
                logfile = pjoin(os.getcwd(), 'processMultipliers.log')

        logLevel = config.get('Logging', 'LogLevel')
        verbose = config.getboolean('Logging', 'Verbose')
        datestamp = config.getboolean('Logging', 'Datestamp')

        if args.verbose:
            verbose = True

        flStartLog(logfile, logLevel, verbose, datestamp)
        # Switch off minor warning messages
        import warnings
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
        warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
        warnings.filterwarnings("ignore", category=UserWarning,
                                module="matplotlib")

        warnings.filterwarnings("ignore", category=RuntimeWarning)

        self.working_dir = config.get('Output', 'Working_dir')
        self.gust_file = config.get('Input', 'Gust_file')

        tiles = config.get('Input', 'Tiles')
        self.tiles = [item.strip() for item in tiles.split(',')]
        log.debug('List of tiles to be processed: {0}'.format(self.tiles))
        log.info('Multipliers will be written out to %s', self.working_dir)

        # Get the multipliers, and process them if need be
        self.type_mapping = {'shielding': 'Ms', 'terrain': 'Mz', 'topographic': 'Mt'}
        self.dirns = ['e', 'n', 'ne', 'nw', 's', 'se', 'sw', 'w']

        rootdir = pathLocator.getRootDirectory()
        os.chdir(rootdir)

        try:
            self.main()
        except ImportError as e:
            log.critical("Missing module: {0}".format(e.strerror))
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())

    @timer
    def main(self):
        """
        Main function to combine the multipliers with the regional wind
        speed data.

        :param str configFile: Path to configuration file.

        """
        log.debug('Instantiating getMultipliers class')
        gM = getMultipliers(self.configFile)
        log.debug('Running checkOutputFolders')
        gM.checkOutputFolders(self.working_dir, self.type_mapping)
        log.debug('Running Translate Multipliers')
        gM.copyTranslateMultipliers(self.tiles, self.configFile,
                                    self.type_mapping, self.working_dir)
        log.debug('Running mergeWindMultipliers')
        gM.mergeWindMultipliers(self.type_mapping, self.dirns, self.working_dir)
        log.debug('Running combineDirections')
        gM.combineDirections(self.dirns, self.working_dir)

        # Load the wind data:
        log.info("Loading regional wind data from {0}".format(self.gust_file))
        ncobj = Dataset(self.gust_file, 'r')

        lat = ncobj.variables['lat'][:]
        lon = ncobj.variables['lon'][:]

        delta = lon[1] - lon[0]
        lon = lon - delta / 2.
        lat = lat - delta / 2.

        # Wind speed:
        wspd = ncobj.variables['vmax'][:]

        # Components:
        uu = ncobj.variables['ua'][:]
        vv = ncobj.variables['va'][:]

        processMult(wspd, uu, vv, lon, lat, self.working_dir)

if __name__ == "__main__":
    run()
