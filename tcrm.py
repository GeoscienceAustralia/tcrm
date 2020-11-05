"""
:mod:`tcrm` -- the main TCRM interface
======================================

.. module:: tcrm
    :synopsis: The main interface to TCRM. Determines which components
               of the model system to execute, based on the configuration
               settings specified.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""
# This package needs patching to run on python 2.6
import logging as log
log.getLogger('matplotlib').setLevel(log.WARNING)
log.getLogger('shapely').setLevel(log.WARNING)
from functools import reduce

import Utilities.datasets as datasets
import traceback
import argparse
import time
import sys
import os

from os.path import join as pjoin, realpath, isdir, dirname
from functools import wraps
from Utilities.progressbar import SimpleProgressBar as ProgressBar
from Utilities.files import flStartLog, flLoadFile
from Utilities.config import ConfigParser
from Utilities.parallel import attemptParallel, disableOnWorkers
from Utilities.version import version
from Utilities import pathLocator

import matplotlib
matplotlib.use('Agg')  # Use matplotlib backend

# Set Basemap data path if compiled with py2exe
if pathLocator.is_frozen():
    os.environ['BASEMAPDATA'] = pjoin(
        pathLocator.getRootDirectory(), 'mpl-data', 'data')

def checkModules():
    """
    Check that compiled extensions are available

    :param list moduleList: list of modules to check for
    
    :return:
    """
    try:
        import Utilities.akima
        log.debug("Compiled modules were found.")
    except ImportError:
        log.critical("Unable to import compiled modules. "
                     "Make sure to compile the extensions using "
                     "`python installer/setup.py build_ext -i` "
                     "or the `compile.cmd` script")
        sys.exit(1) # or raise

def timer(f):
    """
    A simple timing decorator for the entire process.

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


# Set global version string (for output metadata purposes):
__version__  = version() # pylint: disable-msg=C0103



@disableOnWorkers
def doDataDownload(configFile):
    """
    Check and download the data files listed in the configuration file.
    Datasets are listed in the `Input` section of the configuration
    file, with the option `Datasets`. There must also be a corresponding
    section in the configuration file that inlcudes the url, path where
    the dataset will be stored and the filename that will be stored, e.g.::

        [Input]
        Datasets=IBTRACS

        [IBTRACS]
        URL=ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r05/wmo/csv/Allstorms.ibtracs_wmo.v03r05.csv.gz
        filename=Allstorms.ibtracs_wmo.v03r05.csv
        path=input

    This will attempt to download the gzipped csv file from the given URL
    and save it to the given filename, in the 'input' folder under the current
    directory. Gzipped files are automatically unzipped.


    :param str configFile: Name of configuration file.
    :raises IOError: If the data cannot be downloaded.


    """

    log.info('Checking availability of input data sets')

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')

    datasets.loadDatasets(configFile)
    for dataset in datasets.DATASETS:
        if not dataset.isDownloaded():
            log.info('Input file %s is not available', dataset.filename)
            try:
                log.info('Attempting to download %s', dataset.filename)

                pbar = ProgressBar('Downloading file %s: ' % dataset.filename,
                                   showProgressBar)

                def status(fn, done, size):
                    pbar.update(float(done)/size)

                dataset.download(status)
                log.info('Download successful')
            except IOError:
                log.error('Unable to download %s. Maybe a proxy problem?',
                          dataset.filename)
                sys.exit(1)


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

    subdirs = ['tracks', 'hazard', 'windfield', 'plots', 'plots/hazard',
               'plots/stats', 'log', 'process', 'process/timeseries',
               'process/dat']

    if not isdir(outputPath):
        try:
            os.makedirs(outputPath)
        except (OSError, FileExistsError):
            raise
    for subdir in subdirs:
        if not isdir(realpath(pjoin(outputPath, subdir))):
            try:
                os.makedirs(realpath(pjoin(outputPath, subdir)))
            except (OSError, FileExistsError):
                raise


def doTrackGeneration(configFile):
    """
    Do the tropical cyclone track generation in :mod:`TrackGenerator`.

    The track generation settings are read from *configFile*.

    :param str configFile: Name of configuration file.

    """

    log.info('Starting track generation')

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')

    pbar = ProgressBar('Simulating cyclone tracks: ', showProgressBar)

    def status(done, total):
        pbar.update(float(done)/total)

    import TrackGenerator
    TrackGenerator.run(configFile, status)

    pbar.update(1.0)
    log.info('Completed track generation')


def doWindfieldCalculations(configFile):
    """
    Do the wind field calculations, using :mod:`wind`. The wind
    field settings are read from *configFile*.

    :param str configFile: Name of configuration file.

    """

    log.info('Starting wind field calculations')

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')

    pbar = ProgressBar('Calculating wind fields: ', showProgressBar)

    def status(done, total):
        pbar.update(float(done)/total)

    import wind
    wind.run(configFile, status)

    pbar.update(1.0)
    log.info('Completed wind field calculations')


@disableOnWorkers
def doDataProcessing(configFile):
    """
    Parse the input data and turn it into the necessary format
    for the model calibration step, using the :mod:`DataProcess` module.

    :param str configFile: Name of configuration file.

    """

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')

    pbar = ProgressBar('Processing data files: ', showProgressBar)

    log.info('Running Data Processing')

    from DataProcess.DataProcess import DataProcess
    dataProcess = DataProcess(configFile, progressbar=pbar)
    dataProcess.processData()

    log.info('Completed Data Processing')
    pbar.update(1.0)


@disableOnWorkers
def doDataPlotting(configFile):
    """
    Plot the pre-processed input data. Requires the Data Processing step to have
    been executed first (``Actions -- DataProcess = True``)

    :param str configFile: Name of configuration file.

    """

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')
    pbar = ProgressBar('Plotting results: ', showProgressBar)

    outputPath = config.get('Output', 'Path')

    statsPlotPath = pjoin(outputPath, 'plots', 'stats')
    processPath = pjoin(outputPath, 'process')

    pRateData = flLoadFile(pjoin(processPath, 'pressure_rate'))
    pAllData = flLoadFile(pjoin(processPath, 'all_pressure'))
    bRateData = flLoadFile(pjoin(processPath, 'bearing_rate'))
    bAllData = flLoadFile(pjoin(processPath, 'all_bearing'))
    sRateData = flLoadFile(pjoin(processPath, 'speed_rate'))
    sAllData = flLoadFile(pjoin(processPath, 'all_speed'))
    freq = flLoadFile(pjoin(processPath, 'frequency'))


    indLonLat = flLoadFile(pjoin(processPath, 'cyclone_tracks'),
                           delimiter=',')
    indicator = indLonLat[:, 0]
    lonData = indLonLat[:, 1]
    latData = indLonLat[:, 2]

    jdayobs = flLoadFile(pjoin(processPath, 'jday_obs'), delimiter=',')
    jdaygenesis = flLoadFile(pjoin(processPath, 'jday_genesis'), delimiter=',')


    from PlotInterface.plotStats import PlotPressure, PlotBearing, \
        PlotSpeed, PlotFrequency, PlotDays, PlotLonLat

    log.info('Plotting pressure data')
    pbar.update(0.05)
    PrsPlot = PlotPressure(statsPlotPath, "png")
    PrsPlot.plotPressure(pAllData)
    PrsPlot.plotPressureRate(pRateData)
    PrsPlot.plotMinPressure(indicator, pAllData)

    #FIXME: To be moved into `PlotPressure` class.
    #plotting.minPressureLat(pAllData, latData)

    log.info('Plotting bearing data')
    pbar.update(0.15)
    BearPlot = PlotBearing(statsPlotPath, "png")
    BearPlot.plotBearing(bAllData)
    BearPlot.plotBearingRate(bRateData)

    log.info('Plotting speed data')
    pbar.update(0.25)
    SpeedPlot = PlotSpeed(statsPlotPath, "png")
    SpeedPlot.plotSpeed(sAllData)
    SpeedPlot.plotSpeedRate(sRateData)

    log.info('Plotting longitude and latitude data')
    pbar.update(0.45)

    # FIXME: To be moved to it's own class in PlotStats
    LLPlot = PlotLonLat(statsPlotPath, "png")
    LLPlot.plotLonLat(lonData, latData, indicator)

    pbar.update(0.65)

    log.info('Plotting frequency data')
    pbar.update(0.85)
    FreqPlot = PlotFrequency(statsPlotPath, "png")
    FreqPlot.plotFrequency(freq[:, 0], freq[:, 1])

    DayPlot = PlotDays(statsPlotPath, "png")
    DayPlot.plotJulianDays(jdayobs, jdaygenesis)

    pbar.update(1.0)

@disableOnWorkers
def doStatistics(configFile):
    """
    Calibrate the model with the :mod:`StatInterface` module.

    :param str configFile: Name of configuration file.

    """
    from DataProcess.CalcTrackDomain import CalcTrackDomain

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')
    getRMWDistFromInputData = config.getboolean('RMW',
                                                'GetRMWDistFromInputData')

    log.info('Running StatInterface')
    pbar = ProgressBar('Calibrating model: ', showProgressBar)

    # Auto-calculate track generator domain
    CalcTD = CalcTrackDomain(configFile)
    domain = CalcTD.calcDomainFromFile()

    pbar.update(0.05)

    from StatInterface import StatInterface
    statInterface = StatInterface.StatInterface(configFile,
                                                autoCalc_gridLimit=domain)
    statInterface.kdeGenesisDate()
    pbar.update(0.4)

    statInterface.kdeOrigin()
    pbar.update(0.5)

    statInterface.cdfCellBearing()
    pbar.update(0.6)

    statInterface.cdfCellSpeed()
    pbar.update(0.7)

    statInterface.cdfCellPressure()
    pbar.update(0.8)

    statInterface.calcCellStatistics()

    if getRMWDistFromInputData:
        statInterface.cdfCellSize()

    pbar.update(1.0)
    log.info('Completed StatInterface')


def doHazard(configFile):
    """
    Do the hazard calculations (extreme value distribution fitting)
    using the :mod:`hazard` module.

    :param str configFile: Name of configuration file.

    """

    log.info('Running HazardInterface')

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')
    pbar = ProgressBar('Performing hazard calculations: ', showProgressBar)

    def status(done, total):
        pbar.update(float(done)/total)

    import hazard
    hazard.run(configFile)

    log.info('Completed HazardInterface')
    pbar.update(1.0)

@disableOnWorkers
def doHazardPlotting(configFile):
    """
    Do the hazard plots (hazard maps and curves for all locations within
    the model domain). Plotting is performed by the
    :mod:`PlotInterface.AutoPlotHazard` module.

    :param str configFile: Name of configuration file.

    """

    config = ConfigParser()
    config.read(configFile)

    log.info('Plotting Hazard Maps')

    showProgressBar = config.get('Logging', 'ProgressBar')
    pbar = ProgressBar('Plotting hazard maps: ', showProgressBar)
    pbar.update(0.0)

    from PlotInterface.AutoPlotHazard import AutoPlotHazard
    plotter = AutoPlotHazard(configFile, progressbar=pbar)
    plotter.plotMap()
    plotter.plotCurves()

    pbar.update(1.0)

def doDatabaseUpdate(configFile):
    """
    Build a database containing info on the events, locations, return
    period wind speeds and tracks.

    :param str configFile: Name of the configuration file.

    """

    log.info("Creating hazard database")
    import database
    database.run(configFile)

    log.info("Created and populated database")


def doEvaluation(configFile):
    """
    Do the track model evaluation processing, using :mod:`Evaluate`.

    To perform this stage, it is recommended to generate an event set
    with several simulated years (e.g. 10, 20, 30 years) in each
    simulation.

    A good setting might be::

        [Actions]
        ExecuteTrackGenerator=True
        ExecuteEvaluate=True

        [TrackGenerator]
        NumSimulations=1000
        YearsPerSimulation=50

    This will generate 1000 simulations each with 50 years of simulated
    TC activity. :mod:`Evaluate` will then compare pressure distributions,
    track density, landfall rates and longitude crossing rates for the
    input dataset and the full 1000 simulations.

    :param str configFile: Name of the configuration file.

    """

    log.info("Running Evaluation")

    import Evaluate
    Evaluate.run(configFile)


@timer
def main(configFile='main.ini'):
    """
    Main interface of TCRM that allows control and interaction with the
    5 interfaces: DataProcess, StatInterface, TrackGenerator,
    WindfieldInterface and HazardInterface

    :param str configFile: Name of file containing configuration settings
                           for running TCRM

    """

    log.info("Starting TCRM")
    log.info("Configuration file: %s", configFile)

    doOutputDirectoryCreation(configFile)

    config = ConfigParser()
    config.read(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'DownloadData'):
        doDataDownload(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'DataProcess'):
        doDataProcessing(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'ExecuteStat'):
        doStatistics(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'ExecuteTrackGenerator'):
        doTrackGeneration(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'ExecuteWindfield'):
        doWindfieldCalculations(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'ExecuteHazard'):
        doHazard(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'PlotData'):
        doDataPlotting(configFile)

    comm.barrier()

    if config.getboolean('Actions', 'CreateDatabase'):
        doDatabaseUpdate(configFile)

    comm.barrier()
    if config.getboolean('Actions', 'ExecuteEvaluate'):
        doEvaluation(config)

    comm.barrier()

    if config.getboolean('Actions', 'PlotHazard'):
        doHazardPlotting(configFile)

    comm.barrier()

    log.info('Completed TCRM')


def startup():
    """
    Parse command line arguments, set up logging and attempt
    to execute the main TCRM functions.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', help='The configuration file')
    parser.add_argument('-v', '--verbose', help='Verbose output',
                        action='store_true')
    parser.add_argument('-d', '--debug', help='Allow pdb traces',
                        action='store_true')
    args = parser.parse_args()

    configFile = args.config_file

    rootdir = pathLocator.getRootDirectory()
    os.chdir(rootdir)

    config = ConfigParser()
    config.read(configFile)

    logfile = config.get('Logging', 'LogFile')
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
    if comm.size > 1 and comm.rank > 0:
        logfile += '-' + str(comm.rank)
        verbose = False  # to stop output to console
    else:
        pass
        #codeStatus = status()
        #print __doc__ + codeStatus

    flStartLog(logfile, logLevel, verbose, datestamp)

    # Switch off minor warning messages
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    checkModules()

    if debug:
        main(configFile)
    else:
        try:
            main(configFile)
        except ImportError as e:
            log.critical("Missing module: {0}".format(e))
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())
        except Exception:  # pylint: disable=W0703
            # Catch any exceptions that occur and log them (nicely):
            tblines = traceback.format_exc().splitlines()
            for line in tblines:
                log.critical(line.lstrip())


if __name__ == "__main__":
    startup()
