from osgeo import osr, gdal, gdalconst
from osgeo.gdal_array import BandReadAsArray, CopyDatasetInfo, BandWriteArray
from netCDF4 import Dataset
import os
import xarray as xr
import numpy as np
from tqdm import tqdm
from Utilities.config import ConfigParser
from Utilities.files import flStartLog
import argparse
from os.path import join as pjoin, dirname, realpath, isdir, splitext
import traceback
import logging as log


def downscale_multipliers(src_file, match_file, dst_file):
    """
    Downscales and clips GDAL compatable file and saves it as a GEOTIFF.

    Params:
        - src_file: filepath of the input file to be downscaled
        - match_file: filepath of the file that the input is transformed to match
        - dst_file: output filepath
    """
    # load src
    src = gdal.Open(src_file, gdal.GA_ReadOnly)

    # load match info
    ncobj = Dataset(match_file, 'r')
    lat = ncobj.variables['lat'][:]
    lon = ncobj.variables['lon'][:]
    delta = lon[1] - lon[0]
    lon = lon - delta / 2.
    lat = lat - delta / 2.

    dx = lon[1] - lon[0]
    dy = lat[1] - lat[0]
    originX, originY = lon[0], lat[0]
    epsg = 4326

    wide = len(lon)
    high = len(lat)

    # Output / destination
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)

    # warp
    drv = gdal.GetDriverByName('GTiff')
    dst = drv.Create(dst_file, wide, high, 8, gdal.GDT_Float32)
    dst.SetGeoTransform((originX, dx, 0, originY, 0, dy))
    dst.SetProjection(srs.ExportToWkt())
    dstBand = dst.GetRasterBand(1)
    dstBand.SetNoDataValue(-9999)

    gdal.ReprojectImage(src, dst, src.GetProjection(), dst.GetProjection(), gdalconst.GRA_Bilinear)


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
        warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

        warnings.filterwarnings("ignore", category=RuntimeWarning)

        self.working_dir = config.get('Output', 'Working_dir')
        self.gust_dir = config.get('Input', 'Gust_dir')
        self.mult_file = config.get('Input', 'Multipliers')

        try:
            self.main()
        except ImportError as e:
            log.critical("Missing module: {0}".format(e.strerror))
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())
    
    def main(self):
        low_res_mult_file = os.path.join(self.working_dir, "low_res_m4.tif")
        gust_files = [os.path.join(self.gust_dir, fn) for fn in os.listdir(self.gust_dir) if fn.startswith("gust")]

        log.info("Downscaling multipliers")
        gdal.SetConfigOption('GDAL_NUM_THREADS', "16")
        downscale_multipliers(self.mult_file, gust_files[0], low_res_mult_file)
        gdal.SetConfigOption('GDAL_NUM_THREADS', "1")

        # indices, band numbers, and directions for using wind multipliers
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

        # load in wind multiplier data
        # reducing the resolution seems to ignore the -9999 nodata causing some -9999 to be averaged with real data
        # these points (along with other nodata points) are set to 1
        ds = gdal.Open(low_res_mult_file, gdal.GA_ReadOnly)
        bands = []
        for i in range(1, 9):
            band = ds.GetRasterBand(i).ReadAsArray()
            band[band < 0] = 1
            bands.append(band)

        log.info("Applying multipliers")
        # loop through the gust files and apply the wm
        for gust_file in tqdm(gust_files):
            gust = xr.load_dataset(gust_file)
            wind_data = gust.vmax
            local = wind_data.copy()

            bearing = 2 * np.pi - (np.arctan2(-gust.va, -gust.ua) - np.pi / 2)
            bearing = (180. / np.pi) * np.mod(bearing, 2. * np.pi)

            for i in list(indices.keys()):
                idx = np.where((bearing >= indices[i]['min']) & (bearing < indices[i]['max']))
                m4 = bands[band_numbers_for_indices_in_geotiff[i] - 1]
                local.data[idx] = wind_data.data[idx] * m4[idx]

            outds = xr.Dataset()
            outds["vmax"] = local
            outds.to_netcdf(gust_file.replace("gust", "wm_gust"))


if __name__ == "__main__":
    run()