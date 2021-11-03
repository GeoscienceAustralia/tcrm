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
se, s, sw, w or nw). Alternatively, computed site-exposure multipliers
can also be provided with all 8 direction data in 'e', 'n', 'ne', 'nw',
's', 'se', 'sw' and 'w' sequence bands.

AWS S3 location can also be specified for input and output of this
module following GDAL Virtual File Systems convention
(https://gdal.org/user/virtual_file_systems.html#vsis3). The location
should be specified in /vsis3/bucket/key format. AWS authentication
should be provided in ~/.aws/credentials property file and AWS region
should be set in ~/.aws/config property file under 'default' profile.
Aleternative ways of AWS authentication can also be found in
https://gdal.org/user/virtual_file_systems.html#vsis3.
AWS user access keys should be correctly configured in "Identity and
Access Management (IAM)" => "Users" => "Security credentials" =>
"Access keys" with correct pemissions to access S3.
Example of ~/.aws/credentials:
    [default]
    aws_access_key_id = <access_key>
    aws_secret_access_key = <secret_access_key>
Example of ~/.aws/config:
    [default]
    region = ap-southeast-2
    output = json

Requires the Python GDAL bindings, Numpy, netCDF4 and the :mod:`files`
and :mod:`config` modules from TCRM. It assumes :mod:`Utilities` can
be found in the ``PYTHONPATH`` directory.
    Make sure the following modules are loaded into the environment prior to running:
    module load openmpi/1.8.4
    module load python/2.7.6
    module load python/2.7.6-matplotlib
    module load geos
    module load gdal/1.11.1-python

"""

import glob
import logging as log
import math
import os
import queue
import tempfile
import threading
import time
import traceback
from concurrent import futures
from functools import wraps, reduce
from os.path import join as pjoin, dirname, realpath, isdir, splitext
from shutil import copyfile

import argparse
import boto3
import numpy as np
import numpy.ma as ma
from botocore.exceptions import ClientError
from netCDF4 import Dataset
from osgeo import osr, gdal, gdalconst
from osgeo.gdal_array import BandReadAsArray, CopyDatasetInfo, BandWriteArray

from Utilities import pathLocator
from Utilities.AsyncRun import AsyncRun
from Utilities.config import ConfigParser
from Utilities.files import flStartLog

threadLock_gust = threading.Lock()
threadLock_bear = threading.Lock()
threadLock_m4 = threading.Lock()
threadLock_out = threading.Lock()

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

        # Initialise variables to None that will be used later
        self.temporary_raw_multipliers_directory = None
        self.temporary_working_directory = None
        self.temporary_working_directory_name = None
        self.s3_upload_path = None
        self._s3_client = None
        self.computed_wm_path = None
        # Initialise to default values
        self.extent_applied = False
        self.process_multi_version = "2"
        self.max_working_threads = 4
        self.processing_segment_size = 256
        self.warp_memory_limit = 500

        # Check for computed wind multiplier file path in config file
        if config.has_option('Input', 'Multipliers'):
            self.computed_wm_path = config.get('Input', 'Multipliers')
            log.info('Using computed wind multiplier files from {0}'.format(self.computed_wm_path))
            if config.has_option('Input', 'Extent'):
                self.extent = config.geteval('Input', 'Extent')
                log.info('Provided extent {0}'.format(self.extent))
            else:
                self.extent = dict()
        # Check for wind multiplier file path in config file
        elif config.has_option('Input', 'RawMultipliers'):
            self.WMPath = config.get('Input', 'RawMultipliers')
            log.info('Using multiplier files from {0}'.format(self.WMPath))
        else:
            log.info('Using default multiplier files from /g/data/fj6/multipliers/')
            self.WMPath = '/g/data/fj6/multipliers/'

        if config.has_option('Input', 'ExtentApplied'):
            self.extent_applied = (config.get('Input', 'ExtentApplied').lower() == 'yes')
        if config.has_option('ProcessMultipliers', 'ProcessMultiVersion'):
            self.process_multi_version = config.get('ProcessMultipliers', 'ProcessMultiVersion')
            log.info('Using ProcessMultiVersion {0}'.format(self.process_multi_version))
        if config.has_option('ProcessMultipliers', 'MaxWorkingThreads'):
            self.max_working_threads = config.getint('ProcessMultipliers', 'MaxWorkingThreads')
            log.info('Using MaxWorkingThreads {0}'.format(self.max_working_threads))
        if config.has_option('ProcessMultipliers', 'ProcessingSegmentSize'):
            self.processing_segment_size = config.getint('ProcessMultipliers', 'ProcessingSegmentSize')
            log.info('Using ProcessingSegmentSize {0}'.format(self.processing_segment_size))
        if config.has_option('ProcessMultipliers', 'WarpMemoryLimit'):
            self.warp_memory_limit = config.getint('ProcessMultipliers', 'WarpMemoryLimit')
            log.info('Using WarpMemoryLimit {0}'.format(self.warp_memory_limit))

        if config.has_option('Output', 'Temporary_working_dir'):
            self.temporary_working_directory_name = config.get('Output', 'Temporary_working_dir')


    def get_s3_client(self):
        """
        Returns service client for S3. It eliminates initialising service client if AWS
        path is not used.
        """
        if self._s3_client is None:
            self._s3_client = boto3.client('s3')
        return self._s3_client

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

        # For computed_wm_path creating sub-directories is not necessary
        if self.computed_wm_path is not None:
            return

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

    @timer
    def build_combined_multipliers_for_all_directions(self, output_path):
        """
        Create Geotiffs for wind multiplier (terrain, topographic and shielding combined) into
        a single Geotiff from 8-band source file by applying specified extent.
        Output files are named "m4_" direction.

        :param str output_path: path to the output directory
        :returns combined direction file location
        """
        if self.computed_wm_path.endswith('.tif') and self.extent_applied:
            # If S3 download
            if self.computed_wm_path.startswith('/vsis3/'):
                downloaded_computed_wm_path = self.download_from_s3(self.computed_wm_path,output_path)
                log.info('Wind multipliers with applied extent is downloded to {0} from {1}'
                         .format(downloaded_computed_wm_path, self.computed_wm_path))
                return downloaded_computed_wm_path
            else:
                log.info('Wind multipliers with applied extent located at {0}'.format(self.computed_wm_path))
                return self.computed_wm_path
        # Copy source data to local directory limited by extent
        combined_direction_file = pjoin(output_path,'m4_source.tif')
        log.info('Wind multipliers within extent will be written to {0}'.format(combined_direction_file))
        translate_options = gdal.TranslateOptions(format='GTiff',
                                                  outputType=gdal.GDT_Float32,
                                                  outputSRS='EPSG:4326',
                                                  bandList=[1, 2, 3, 4, 5, 6, 7, 8],
                                                  projWin=[self.extent['xMin'], self.extent['yMax'], self.extent['xMax'], self.extent['yMin']])
        gdal.Translate(combined_direction_file, self.computed_wm_path, options=translate_options)
        return combined_direction_file

    def check_if_wm_source_from_s3(self, type_mapping, tiles, dirns):
        """
        Function to check if wind multiplier location (WMPath) is referring to s3 location. If S3,
        this function downloads relevant input files from 'shielding', 'terrain' and 'topographic'
        sub-directories in 3 async call, place them in local temporary directory and return
        that temporary directory path. As opening NETCDF file from /vsis3 is buggy, source
        .nc files are downloaded and processed.

        :param Dict[str,str] type_mapping: Dict of multipliers to extract sub-directory and suffix
        :param List[str] tiles: Array of tiles specified in config to get the file name
        :param List[str] dirns: Array of string representing 8 directions to be used in file name suffix
        """
        if not self.WMPath.startswith('/vsis3/'):
            return
        if not self.WMPath.endswith('/'):
            self.WMPath = self.WMPath + '/'
        # Keep the reference with self otherwise the temporary directory will be deleted.
        self.temporary_raw_multipliers_directory = tempfile.TemporaryDirectory(prefix='RawMulti-')
        log.info('Copying to temporary directory %s as RawMultiplier is specified in S3 %s',
                 self.temporary_raw_multipliers_directory.name, self.WMPath)
        self.get_s3_client() # Dummy call so that AsyncRun can't create multiple service client
        async_runs = []
        for type_name in type_mapping:
            files_to_download = []
            os.mkdir(pjoin(self.temporary_raw_multipliers_directory.name, type_name))
            for timage in tiles:
                for dir in dirns:
                    files_to_download.append('{0}_{1}_{2}.nc'.format(timage, type_mapping[type_name].lower(), dir))
            # Download files from S3 with async call
            args = {
                's3_source_directory_path': self.WMPath + type_name + '/',
                'destination_directory': pjoin(self.temporary_raw_multipliers_directory.name, type_name),
                'files_to_download': files_to_download
            }
            thread = AsyncRun(self.download_files_from_s3, args)
            thread.start()
            async_runs.append(thread)
        # Wait for threads to complete
        for thread in async_runs:
            thread.join()
        log.info('RawMultiplier copy completed to temporary directory %s', self.temporary_raw_multipliers_directory.name)
        self.WMPath = self.temporary_raw_multipliers_directory.name + os.path.sep

    def check_if_gust_file_from_s3(self, gust_file, working_dir):
        """
        Function to check if Gust_file location is referring to s3 location. If S3 location specified
        .nc file will be downloaded. Cause: opening NETCDF file from /vsis3 is buggy with Dataset()
        and gdal.open() so source must be downloaded. Most likely NetCDF4 library don't support
        correctly although HDF5 claimed to support correctly.

        :param str gust_file: Location of Gust_file specified in config
        :param str working_dir: Working_dir property value provided in config to write output
        """
        if gust_file.startswith('/vsis3/'):
            return self.download_from_s3(gust_file, working_dir)
        return gust_file

    def check_if_working_dir_in_s3(self, working_dir):
        """
        Function to check if output directory is on S3 create a temporary working
        directory for intermidiate outputs

        :param str working_dir: Working_dir property value provided in config to write output
        """
        if working_dir.startswith('/vsis3/'):
            # Keep the reference with self otherwise the temporary directory will be deleted.
            if self.temporary_working_directory_name is None:
                self.temporary_working_directory = tempfile.TemporaryDirectory(prefix='Multipliers-')
                log.info('Creating temporary working_dir %s as S3 working_dir specified %s',
                         self.temporary_working_directory.name, working_dir)
                self.temporary_working_directory_name = self.temporary_working_directory.name
            self.s3_upload_path = working_dir
            if self.temporary_working_directory_name.endswith(os.path.sep):
                return self.temporary_working_directory_name
            else:
                return self.temporary_working_directory_name + os.path.sep
        return working_dir

    @timer
    def upload_files_to_s3(self, local_directory, s3_destination_directory_path, files_to_upload):
        """
        Function to upload files from local directory to s3.

        :param str local_directory: Local directory path containing files to upload.
        :param str s3_destination_directory_path: S3 path of the destination directory.
        :param List[str] files_to_upload: List of file names to upload.
        """
        if not s3_destination_directory_path.endswith('/'):
            s3_destination_directory_path = s3_destination_directory_path + '/'
        [bucket_name, bucket_key, file_name] = self.s3_path_segments_from_vsis3(s3_destination_directory_path)
        try:
            s3_client = self.get_s3_client()
            for file_name in files_to_upload:
                local_path = pjoin(local_directory, file_name)
                log.info("Uploading file to S3 bucket: {0}, key: {1}, local path: {2}"
                         .format(bucket_name, bucket_key + file_name, local_path))
                s3_client.upload_file(local_path, bucket_name, bucket_key + file_name)
        except ClientError as e:
            log.exception("S3 write error: {0}".format(file_name))
            raise e

    def download_from_s3(self, s3_source_path, destination_directory):
        """
        Function to download a S3 file into local directory.

        :param str s3_source_path: S3 path of the file.
        :param str destination_directory: Local directory location to
        """
        [bucket_name, bucket_key, file_name] = self.s3_path_segments_from_vsis3(s3_source_path)
        file_path = pjoin(destination_directory, file_name)
        log.info("Downloading file from S3 bucket: {0}, key: {1}, local path: {2}"
                 .format(bucket_name, bucket_key, file_path))
        try:
            s3_client = self.get_s3_client()
            s3_client.download_file(bucket_name, bucket_key, file_path)
        except ClientError as e:
            log.exception("S3 read error: {0}".format(file_name))
            raise e
        return file_path

    def download_files_from_s3(self, s3_source_directory_path, destination_directory, files_to_download):
        """
        Function to download a S3 file into local directory.

        :param str s3_source_directory_path: S3 path of the directory contining files.
        :param str destination_directory: Local directory location to store the downloaded files
        :param List[str] files_to_download: List of files to download from S3
        """
        if not s3_source_directory_path.endswith('/'):
            s3_source_directory_path = s3_source_directory_path + '/'
        [bucket_name, bucket_key, file_name] = self.s3_path_segments_from_vsis3(s3_source_directory_path)
        try:
            s3_client = self.get_s3_client()
            for file_name in files_to_download:
                file_path = pjoin(destination_directory, file_name)
                log.info("Downloading file from S3 bucket: {0}, key: {1}, local path: {2}"
                         .format(bucket_name, bucket_key + file_name, file_path))
                s3_client.download_file(bucket_name, bucket_key + file_name, file_path)
        except ClientError as e:
            log.exception("S3 read error: {0}".format(file_name))
            raise e
        return file_path

    def s3_path_segments_from_vsis3(self, s3_path):
        """
        Function to extract bucket name, key and filename from path specified using
        GDAL Virtual File Systems conventions

        :param str s3_path: Path to S3 location starting with /vsis3/.
        """
        s3_path_segments = s3_path.split('/')
        if s3_path_segments[0] != '' or s3_path_segments[1] != 'vsis3':
            raise ValueError('Invalid path: ', [s3_path, s3_path_segments])
        file_name = s3_path_segments[-1]
        bucket_name = s3_path_segments[2]
        bucket_key = '/'.join(s3_path_segments[3:])
        return bucket_name, bucket_key, file_name


def computeOutputExtentIfInvalid(extent, gust_file, computed_wm_path):
    """
    Check validity of extent. If extent is not valid, output image extent is computed from
    regional wind data (gust file) and wind multiplier image extents.

    :param dict extent: extent provided in config file
    :param str gust_file: file path of regional wind data / gust file
    :param str computed_wm_path: file path of wind multiplier image
    """
    if 'xMin' in extent and 'xMax' in extent and 'yMin' in extent and 'yMax' in extent:
        return extent

    log.info('Invalid extent. Calculating default extent.')
    ncobj = Dataset(gust_file, 'r')
    lat = ncobj.variables['lat'][:]
    lon = ncobj.variables['lon'][:]
    xMinGust = min(lon)
    xMaxGust = max(lon)
    yMinGust = min(lat)
    yMaxGust = max(lat)
    log.info('Extent from regional wind data '
             '{{\'xMin\': {0}, \'xMax\': {1}, \'yMin\': {2}, \'yMax\': {3} }}'
             .format(xMinGust, xMaxGust, yMinGust, yMaxGust))
    del ncobj

    # Calculate extent of wind multiplier image
    ds = gdal.Open(computed_wm_path, gdal.GA_ReadOnly)
    widthWM = ds.RasterXSize
    heightWM = ds.RasterYSize
    gt = ds.GetGeoTransform()
    xMinWM = gt[0]
    yMinWM = gt[3] + widthWM * gt[4] + heightWM * gt[5]
    xMaxWM = gt[0] + widthWM * gt[1] + heightWM * gt[2]
    yMaxWM = gt[3]
    log.info('Extent from wind multiplier image '
             '{{\'xMin\': {0}, \'xMax\': {1}, \'yMin\': {2}, \'yMax\': {3} }}'
             .format(xMinWM, xMaxWM, yMinWM, yMaxWM))
    del ds

    # Take only intersecting extent of provided extent and wind multiplier image extent
    extent = dict(xMin=max(xMinGust, xMinWM),
                  xMax=min(xMaxGust, xMaxWM),
                  yMin=max(yMinGust, yMinWM),
                  yMax=min(yMaxGust, yMaxWM))
    log.info('Applying effective extent {0}'.format(extent))
    return extent

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


def loadRasterFile(raster_file, fill_value=1, band_number=1):
    """
    Load a raster file and return the data as a :class:`numpy.ndarray`.
    No prorjection information is returned, just the actual data as an
    array.

    :param str raster_file: Path to the raster file to load.
    :param fill_value: Value to replace `nodata` values with (default=1).
    :param int band_number: Band number in raster file (default=1).
    :returns: 2-d array of the data values.
    :rtype: :class:`numpy.ndarray`

    """
    log.debug("Loading raster data from {0} band {1} into array".format(raster_file, band_number))
    ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(band_number)
    data = band.ReadAsArray()

    nodata = band.GetNoDataValue()

    if nodata is not None:
        np.putmask(data, data == nodata, fill_value)

    del ds
    return data


def loadAllBandArrayData(band_sources, fill_value=1, segment_info=None):
    """
    Load all 8 band data within specified segment and return the data
    as a :class:`numpy.ndarray`. No prorjection information is
    returned, just the actual data as an array.

    :param List[object] band_sources: reference to 8 bands in array
    :param fill_value: Value to replace `nodata` values with (default=1).
    :param List[int] segment_info: location and size of segment
    :returns: 2-d array of the data values.
    :rtype: :class:`numpy.ndarray`
    """
    data_array = []
    with threadLock_m4:
        for band in band_sources:
            if segment_info is None:
                data = band.ReadAsArray()
            else:
                data = band.ReadAsArray(segment_info[0], segment_info[1], segment_info[2], segment_info[3])
            nodata = band.GetNoDataValue()
            if nodata is not None:
                np.putmask(data, data == nodata, fill_value)
            data_array.append(data)
    return data_array


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
                     resampling_method=gdalconst.GRA_Bilinear,
                     match_projection=None, warp_memory_limit=0.0):
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
    gdal.ReprojectImage(src, dst, src_proj, match_proj, resampling_method, WarpMemoryLimit=warp_memory_limit)

    del dst  # Flush
    if isinstance(match_filename, str):
        del match_ds
    if isinstance(src_file, str):
        del src

    return

@timer
def processMult(wspd, uu, vv, lon, lat, working_dir, m4_max_file = 'm4_ne.tif', combined_mulipliers_file=None):
    """
    The lat and lon values are the top left corners of the cells
    The speed arrays are in bottom to top format

    :param wspd: The gust speed
    :param uu: x component of the wind speed
    :param vv: y component of the wind speed
    :param lon: list of raster longitude values
    :param lat:  list of raster latitude values
    :param working_dir: The working output directory
    :param m4_max_file: Multiplier file used for reprojection.
    :param string combined_mulipliers_file: Path of raster file containing wind
        multiplier data of all 8 directions
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
    if combined_mulipliers_file is not None:
        m4_max_file = combined_mulipliers_file
    else:
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
    band_numbers_for_indices_in_geotiff = [2, 3, 1, 6, 5, 7, 8, 4, 2]
    log.info("Processing all directions")
    for i in list(indices.keys()):
        dn = indices[i]['dir']
        log.info("Processing {0}".format(dn))
        if combined_mulipliers_file is None:
            m4_file = pjoin(working_dir, 'm4_{0}.tif'.format(dn.lower()))
            m4 = loadRasterFile(m4_file)
        else:
            m4 = loadRasterFile(combined_mulipliers_file, band_number=band_numbers_for_indices_in_geotiff[i])
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

@timer
def processMultV2(wspd, uu, vv, lon, lat, working_dir, dirns,
                  max_working_threads, processing_segment_size, warp_memory_limit,
                  m4_max_file = 'm4_ne.tif', combined_mulipliers_file=None):
    """
    The lat and lon values are the top left corners of the cells
    The speed arrays are in bottom to top format.
    V2 for processing by image segment instead of loading
    whole image into memory.

    :param wspd: The gust speed
    :param uu: x component of the wind speed
    :param vv: y component of the wind speed
    :param lon: list of raster longitude values
    :param lat:  list of raster latitude values
    :param working_dir: The working output directory
    :param m4_max_file: Multiplier file used for reprojection.
    :param string combined_mulipliers_file: Path of raster file containing wind
        multiplier data of all 8 directions
    :return string: path of local wind file
    """
    # This gives different bearing values
    # thank the bearings in the result tuple
    bearing = calculateBearing(uu, vv)
    delta = lon[1] - lon[0]

    log.debug('Create rasters from the netcdf gust file variables')
    wind_raster_file = pjoin(working_dir, 'region_wind.tif')
    wind_raster = createRaster(np.flipud(wspd), lon, lat, delta, delta,
                               filename=wind_raster_file)
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
    if combined_mulipliers_file is not None:
        m4_max_file = combined_mulipliers_file
    else:
        m4_max_file = pjoin(working_dir, m4_max_file)

    future_requests = []
    with futures.ThreadPoolExecutor(max_workers=max_working_threads) as e:
        m4_max_file_obj = gdal.Open(m4_max_file, gdal.GA_ReadOnly)
        reprojectDataset(wind_raster, m4_max_file_obj, wind_prj_file,
                         warp_memory_limit=warp_memory_limit)
        reprojectDataset(bear_raster, m4_max_file_obj, bear_prj_file,
                         warp_memory_limit=warp_memory_limit,
                         resampling_method=gdalconst.GRA_NearestNeighbour)
        future_requests.append(e.submit(reprojectDataset, uu_raster, m4_max_file_obj, uu_prj_file,
                                        warp_memory_limit=warp_memory_limit,
                                        resampling_method=gdalconst.GRA_NearestNeighbour))
        future_requests.append(e.submit(reprojectDataset, vv_raster, m4_max_file_obj, vv_prj_file,
                                        warp_memory_limit=warp_memory_limit,
                                        resampling_method=gdalconst.GRA_NearestNeighbour))

        wind_prj_ds = gdal.Open(wind_prj_file, gdal.GA_ReadOnly)
        wind_prj = wind_prj_ds.GetRasterBand(1)
        bear_prj_ds = gdal.Open(bear_prj_file, gdal.GA_ReadOnly)
        bear_prj = bear_prj_ds.GetRasterBand(1)
        wind_proj = wind_prj_ds.GetProjection()
        wind_geot = wind_prj_ds.GetGeoTransform()

        cols = wind_prj.XSize
        rows = wind_prj.YSize

        output_file = pjoin(working_dir, 'local_wind.tif')
        log.info("Creating output file: {0}".format(output_file))
        # Save the local wind field to a raster file with the SRS of the
        # multipliers
        drv = gdal.GetDriverByName("GTiff")
        dst_ds = drv.Create(output_file, cols, rows, 1,
                            gdal.GDT_Float32, ['BIGTIFF=YES', 'SPARSE_OK=TRUE'])
        dst_ds.SetGeoTransform(wind_geot)
        dst_ds.SetProjection(wind_proj)
        dst_band = dst_ds.GetRasterBand(1)
        dst_band.SetNoDataValue(-9999)
        print('processMultV2', dst_ds.GetProjection())

        log.info("Reading bands")
        source_dir_bands = []
        m4_ds_arr = []
        if combined_mulipliers_file is None:
            for dn in dirns:
                m4_file = pjoin(working_dir, 'm4_{0}.tif'.format(dn.lower()))
                m4_ds = gdal.Open(m4_file, gdal.GA_ReadOnly)
                source_dir_bands.append(m4_ds.GetRasterBand(1))
                m4_ds_arr.append(m4_ds)
        else:
            dn_ix = 0
            m4_ds = gdal.Open(combined_mulipliers_file, gdal.GA_ReadOnly)
            for dn in dirns:
                dn_ix = dn_ix + 1
                source_dir_bands.append(m4_ds.GetRasterBand(dn_ix))

        total_segments = int(math.ceil(1.0 * cols / processing_segment_size)
                             * math.ceil(1.0 * rows / processing_segment_size))
        segment_count = 0
        segment_queue = queue.Queue(total_segments)
        for y_offset in range(0, rows, processing_segment_size):
            height = rows - y_offset if y_offset + processing_segment_size > rows else processing_segment_size
            for x_offset in range(0, cols, processing_segment_size):
                segment_count = segment_count + 1
                width = cols - x_offset if x_offset + processing_segment_size > cols else processing_segment_size
                segment_queue.put([x_offset, y_offset, width, height, segment_count, total_segments])

        log.info("Lunching {0} segmented task in {1} worker threads".format(total_segments, max_working_threads))
        for _ in range(max_working_threads):
            future_requests.append(e.submit(call_process_multiplier_segment, segment_queue, source_dir_bands, wind_prj, bear_prj, dst_band))

        futures.wait(future_requests, return_when='FIRST_EXCEPTION')
        for task in future_requests:
            task.result()  # Called to obtain exception information if any
        dst_ds.FlushCache()
        del dst_ds

    log.info("Completed")
    return output_file


def call_process_multiplier_segment(segment_queue, source_dir_band, wind_prj, bear_prj, dst_band):
    while not segment_queue.empty():
        processMultiplierSegment(segment_queue.get(), source_dir_band, wind_prj, bear_prj, dst_band)
    dst_band.FlushCache()

def processMultiplierSegment(segment, source_dir_band, wind_prj, bear_prj, dst_band):
    """
    Calculates local wind multiplier data by image segments
    and writes to corresponding segment of output file

    :param segment: image segment specified by [x_offset, y_offset,
                    width, height, segment_count, total_segments]
    :param source_dir_band: 8 band array representing wind mulitpliers
                    data in 8 directions
    :param wind_prj: band representing gust data
    :param bear_prj: band representing bear data
    :param dst_band: band of output file
    """
    band_numbers_for_indices_in_geotiff = [2, 3, 1, 6, 5, 7, 8, 4, 2]
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
    [x_offset, y_offset, width, height, segment_id, total_segments] = segment
    with threadLock_gust:
        wind_data = wind_prj.ReadAsArray(x_offset, y_offset, width, height)
    with threadLock_bear:
        bear_data = bear_prj.ReadAsArray(x_offset, y_offset, width, height)
    m4_all = loadAllBandArrayData(source_dir_band, segment_info=segment)
    local = np.zeros([height, width], dtype='float32')
    for i in list(indices.keys()):
        m4 = m4_all[band_numbers_for_indices_in_geotiff[i] - 1]
        idx = np.where((bear_data >= indices[i]['min']) &
                       (bear_data < indices[i]['max']))
        local[idx] = wind_data[idx] * m4[idx]
    with threadLock_out:
        dst_band.WriteArray(local, x_offset, y_offset)
    if segment_id % int(math.ceil(total_segments / 100.0)) == 0:
        dst_band.FlushCache()
        log.info('Progress: {0:.2f} %'.format((segment_id * 100.0) / total_segments))

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

        if config.has_option('Input', 'Tiles'):
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
        self.working_dir = gM.check_if_working_dir_in_s3(self.working_dir)
        self.gust_file = gM.check_if_gust_file_from_s3(self.gust_file, self.working_dir)
        log.debug('Running checkOutputFolders')
        gM.checkOutputFolders(self.working_dir, self.type_mapping)
        combined_mulipliers_file = None

        if gM.computed_wm_path is None:
            gM.check_if_wm_source_from_s3(self.type_mapping, self.tiles, self.dirns)
            log.debug('Running Translate Multipliers')
            gM.copyTranslateMultipliers(self.tiles, self.configFile,
                                        self.type_mapping, self.working_dir)
            log.debug('Running mergeWindMultipliers')
            gM.mergeWindMultipliers(self.type_mapping, self.dirns, self.working_dir)
            log.debug('Running combineDirections')
            gM.combineDirections(self.dirns, self.working_dir)
        else:
            gM.extent = computeOutputExtentIfInvalid(gM.extent, self.gust_file, gM.computed_wm_path)
            log.debug('Running build_combined_multipliers_for_all_directions')
            combined_mulipliers_file = gM.build_combined_multipliers_for_all_directions(self.working_dir)

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

        del ncobj

        if gM.process_multi_version == "1":
            processMult(wspd, uu, vv, lon, lat, self.working_dir,
                        combined_mulipliers_file=combined_mulipliers_file)
        else:
            processMultV2(wspd, uu, vv, lon, lat, self.working_dir, self.dirns,
                          gM.max_working_threads, gM.processing_segment_size,
                          gM.warp_memory_limit,
                          combined_mulipliers_file=combined_mulipliers_file)

        if gM.s3_upload_path is not None:
            gM.upload_files_to_s3(self.working_dir, gM.s3_upload_path,
                                  ['local_wind.tif', 'region_wind.tif', 'vv_prj.tif',
                                   'uu_prj.tif', 'bear_prj.tif', 'gust_prj.tif'])

if __name__ == "__main__":
    run()
