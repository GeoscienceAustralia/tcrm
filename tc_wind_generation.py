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
from Utilities.parallel import attemptParallel, disableOnWorkers

from Utilities.progressbar import SimpleProgressBar as ProgressBar

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

@disableOnWorkers
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

@disableOnWorkers
def plotFirstWindfields(configFile, limit=10):
    """Plot windfield maps for the first ``limit`` generated gust files."""
    from netCDF4 import Dataset
    import numpy as np
    from PlotInterface.maps import saveWindfieldMap

    config = ConfigParser()
    config.read(configFile)
    outputPath = config.get('Output', 'Path')
    windfieldPath = pjoin(outputPath, 'windfield')
    plotPath = pjoin(outputPath, 'plots', 'windfield')

    if not isdir(plotPath):
        os.makedirs(plotPath)

    gustFiles = sorted(
        pjoin(windfieldPath, f) for f in os.listdir(windfieldPath)
        if f.startswith('gust.') and f.endswith('.nc'))

    for gustFile in gustFiles[:limit]:
        with Dataset(gustFile, 'r') as f:
            xdata = f.variables['lon'][:]
            ydata = f.variables['lat'][:]
            vdata = f.variables['vmax'][:]

            gridLimit = None
            if config.has_option('Region', 'gridLimit'):
                gridLimit = config.geteval('Region', 'gridLimit')
                ii = np.where((xdata >= gridLimit['xMin']) &
                              (xdata <= gridLimit['xMax']))
                jj = np.where((ydata >= gridLimit['yMin']) &
                              (ydata <= gridLimit['yMax']))
                xgrid, ygrid = np.meshgrid(xdata[ii], ydata[jj])
                ig, jg = np.meshgrid(ii, jj)
                vdata = vdata[jg, ig]
            else:
                xgrid, ygrid = np.meshgrid(xdata, ydata)

            map_kwargs = dict(llcrnrlon=xgrid.min(),
                              llcrnrlat=ygrid.min(),
                              urcrnrlon=xgrid.max(),
                              urcrnrlat=ygrid.max(),
                              projection='merc',
                              resolution='i')
            title = os.path.basename(gustFile).replace('.nc', '')
            cbarlabel = "Wind speed ({0})".format(f.variables['vmax'].units)
            maxwind = float(np.nanmax(vdata))
            upper = max(35., np.ceil(maxwind / 5.) * 5.)
            levels = np.arange(0., upper + 5., 5.)
            outfile = pjoin(plotPath, os.path.basename(gustFile).replace('.nc', '.png'))
            saveWindfieldMap(vdata, xgrid, ygrid, title, levels,
                             cbarlabel, map_kwargs, outfile)

    log.info("Plotted %d windfield files to %s", min(limit, len(gustFiles)),
             plotPath)

@timer
def main(configFile):
    """
    Main function to execute the :mod:`wind`.

    :param str configFile: Path to configuration file.

    """
    config = ConfigParser()
    config.read(configFile)
    doOutputDirectoryCreation(configFile)

    # Use the existing 1-hour interpolated tracks under Output/Path/tracks.
    # This wind-generation entry point intentionally does not create or
    # overwrite track files.
    comm.barrier()

    showProgressBar = config.get('Logging', 'ProgressBar')

    pbar = ProgressBar('Calculating wind fields: ', showProgressBar)

    def status(done, total):
        pbar.update(float(done)/total)

    import wind
    wind.run(configFile, status)

    # Generate PNG maps for only the first 10 gust files. Do not plot
    # timeseries from this windfield-only driver.
    plotFirstWindfields(configFile, limit=10)

def startup():
    """
    Parse the command line arguments and call the :func:`main`
    function.

    """
    print("starting up")
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

    global MPI, comm
    MPI = attemptParallel()
    import atexit
    atexit.register(MPI.Finalize)
    comm = MPI.COMM_WORLD

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


