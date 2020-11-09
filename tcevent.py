"""
:mod:`tcevent` -- run the windfield module for a single TC track
================================================================

.. module:: tcevent
    :synopsis: Run the wind field module for a single TC track.

Run the :mod:`wind.windmodels` module to calculate the wind field for a
single TC event. The track of the TC is interpolated to a fine
temporal resolution, then the maximum wind field evaluated.

Data at selected points within the model domain can be extracted
at each time step, giving a time history of wind speed, direction
and estimated sea level pressure for the location(s).

See the :ref:`Scenario modelling <scenariomodelling>` section of
the TCRM User Guide for details on running this script.

"""

import logging as log
log.getLogger('matplotlib').setLevel(log.WARNING)
from functools import reduce

import os
import time
import argparse
import traceback

from functools import wraps
from os.path import join as pjoin, realpath, isdir, dirname

from Utilities import pathLocator
from Utilities.config import ConfigParser
from Utilities.files import flStartLog
from Utilities.version import version
from Utilities.progressbar import SimpleProgressBar as ProgressBar
from Evaluate import interpolateTracks

__version__ = version()

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
          reduce(lambda ll, b : divmod(ll[0], b) + ll[1:],
                        [(tottime,), 60, 60])

        log.info("Time for {0}: {1}".format(f.__name__, msg) )
        return res

    return wrap

def doOutputDirectoryCreation(configFile):
    """
    Create all the necessary output folders.

    :param str configFile: Name of configuration file.
    :raises OSError: If the directory tree cannot be created.

    """

    config = ConfigParser()
    config.read(configFile)

    outputPath = config.get('Output', 'Path')

    log.info('Output will be stored under %s', outputPath)

    subdirs = ['tracks', 'windfield', 'plots', 'plots/timeseries',
               'log', 'process', 'process/timeseries']

    if not isdir(outputPath):
        try:
            os.makedirs(outputPath)
        except OSError:
            raise
    for subdir in subdirs:
        if not isdir(realpath(pjoin(outputPath, subdir))):
            try:
                os.makedirs(realpath(pjoin(outputPath, subdir)))
            except OSError:
                raise

def doTimeseriesPlotting(configFile):
    """
    Run functions to plot time series output.

    :param str configFile: Path to configuration file.

    """

    config = ConfigParser()
    config.read(configFile)

    outputPath = config.get('Output', 'Path')
    timeseriesPath = pjoin(outputPath, 'process', 'timeseries')
    plotPath = pjoin(outputPath, 'plots', 'timeseries')
    log.info("Plotting time series data to {0}".format(plotPath))
    from PlotInterface.plotTimeseries import plotTimeseries
    plotTimeseries(timeseriesPath, plotPath)

def doWindfieldPlotting(configFile):
    """
    Plot the wind field on a map.

    :param str configFile: Path to the configuration file.

    :Note: the file name is assumed to be 'gust.interp.nc'

    """
    from netCDF4 import Dataset
    import numpy as np
    from PlotInterface.maps import saveWindfieldMap

    config = ConfigParser()
    config.read(configFile)
    outputPath = config.get('Output', 'Path')
    windfieldPath = pjoin(outputPath, 'windfield')

    inputFile = config.get('DataProcess', 'InputFile')
    if inputFile.endswith(".nc"):
        # We have a netcdf track file. Work under the assumption it is
        # drawn directly from TCRM.
        trackFile = os.path.basename(inputFile)
        trackId = trackFile.split('.')[1]
        gustFile = 'gust.{0}.nc'.format(trackId)
        outputWindFile = pjoin(windfieldPath, gustFile)
    else:
        # Note the assumption about the file name!
        outputWindFile = pjoin(windfieldPath, 'gust.001-00001.nc')
    plotPath = pjoin(outputPath, 'plots', 'maxwind.png')

    f = Dataset(outputWindFile, 'r')

    xdata = f.variables['lon'][:]
    ydata = f.variables['lat'][:]

    vdata = f.variables['vmax'][:]

    gridLimit = None
    if config.has_option('Region','gridLimit'):
        gridLimit = config.geteval('Region', 'gridLimit')
        ii = np.where((xdata >= gridLimit['xMin']) &
                      (xdata <= gridLimit['xMax']))
        jj = np.where((ydata >= gridLimit['yMin']) &
                      (ydata <= gridLimit['yMax']))
        [xgrid, ygrid] = np.meshgrid(xdata[ii], ydata[jj])
        ig, jg = np.meshgrid(ii, jj)
        vdata = vdata[jg, ig]
    else:
        [xgrid, ygrid] = np.meshgrid(xdata, ydata)
    map_kwargs = dict(llcrnrlon=xgrid.min(),
                      llcrnrlat=ygrid.min(),
                      urcrnrlon=xgrid.max(),
                      urcrnrlat=ygrid.max(),
                      projection='merc',
                      resolution='i')
    title = "Maximum wind speed"
    cbarlabel = "Wind speed ({0})".format(f.variables['vmax'].units)
    levels = np.arange(30, 101., 5.)
    saveWindfieldMap(vdata, xgrid, ygrid, title, levels,
                     cbarlabel, map_kwargs, plotPath)

@timer
def main(configFile):
    """
    Main function to execute the :mod:`wind`.

    :param str configFile: Path to configuration file.

    """
    config = ConfigParser()
    config.read(configFile)
    doOutputDirectoryCreation(configFile)

    trackFile = config.get('DataProcess', 'InputFile')
    source = config.get('DataProcess', 'Source')
    delta = 1/12.
    outputPath = pjoin(config.get('Output','Path'), 'tracks')
    outputTrackFile = pjoin(outputPath, "tracks.interp.nc")

    # This will save interpolated track data in TCRM format:
    interpTrack = interpolateTracks.parseTracks(configFile, trackFile,
                                                source, delta,
                                                outputTrackFile,
                                                interpolation_type='akima')

    showProgressBar = config.get('Logging', 'ProgressBar')

    pbar = ProgressBar('Calculating wind fields: ', showProgressBar)

    def status(done, total):
        pbar.update(float(done)/total)

    import wind
    wind.run(configFile, status)

    import impact
    impact.run_optional(config)

    if config.getboolean('WindfieldInterface', 'PlotOutput'):
        doWindfieldPlotting(configFile)

    if config.getboolean('Timeseries', 'Extract'):
        doTimeseriesPlotting(configFile)

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
            log.critical("Missing module: {0}".format(e))
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())


if __name__ == "__main__":
    startup()


