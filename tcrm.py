#!/usr/bin/env python
"""
Tropical Cyclone Risk Model
Copyright (c) 2013 Commonwealth of Australia (Geoscience Australia)
"""

import Utilities.datasets as datasets
import logging as log
import subprocess
import traceback
import argparse
import time
import sys
import os

from os.path import join as pjoin, realpath, isdir, dirname
from Utilities.progressbar import SimpleProgressBar as ProgressBar
from Utilities.files import flStartLog, flLoadFile
from Utilities.config import ConfigParser
from Utilities import pathLocator

# Set Basemap data path if compiled with py2exe
if pathLocator.is_frozen():
    os.environ['BASEMAPDATA'] = pjoin(
        pathLocator.getRootDirectory(), 'mpl-data', 'data')


UPDATE_MSG = """
----------------------------------------------------------
Your TCRM version is not up-to-date. The last 3 things that
have been fixed are:

%s

To update your version of TCRM, try running:

git pull

----------------------------------------------------------
"""


def timer(func):
    """
    Basic timing functions for entire process
    """
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        log.debug("%s took %0.0f s", func.func_name, t2 - t1)
        return res

    return wrapper


def status():
    """
    Check status TCRM of code.
    """
    msg = ''

    def git(command):
        with open(os.devnull) as devnull:
            return subprocess.check_output(['git ' + command],
                                           shell=True,
                                           stderr=devnull)

    try:
        git('fetch')  # update status of remote
        behind = int(git('rev-list HEAD...origin/master --count'))
        recent = git('log --pretty=format:" - %s (%ad)" --date=relative ' +
                     'origin/master HEAD~3..HEAD')
        if behind != 0:
            msg = UPDATE_MSG % recent
    except subprocess.CalledProcessError:
        pass

    return msg


def attemptParallel():
    """
    Attempt to load Pypar globally as `pp`. If Pypar loads
    successfully, then a call to `pypar.finalize` is registered to be
    called at exit of the Python interpreter. This is to ensure that
    MPI exits cleanly.

    If pypar cannot be loaded then a dummy `pp` is created.
    """
    # pylint: disable=W0621,W0601,C0111,C0103,C0321,R0201
    global pp

    try:
        # load pypar for everyone

        import pypar as pp

        # success! now ensure a clean MPI exit

        import atexit
        atexit.register(pp.finalize)

    except ImportError:

        # just in case user is trying to run with mpirun anyhow
        # NOTE: this only works with OpenMPI

        if os.getenv('OMPI_COMM_WORLD_RANK'):
            raise

        # no pypar and not run with mpirun, use the dummy one

        class DummyPypar(object):
            def size(self): return 1
            def rank(self): return 0
            def barrier(self): pass

        pp = DummyPypar()


def disableOnWorkers(f):
    """
    Disable function calculation on workers. Function will
    only be evaluated on the master.
    """
    def wrap(*args, **kwargs):
        if pp.size() > 1 and pp.rank() > 0:
            return
        else:
            return f(*args, **kwargs)
    return wrap


@disableOnWorkers
def doDataDownload(configFile):
    """
    Check and download the data files.
    """
    log.info('Checking availability of input data sets')

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')

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
        except OSError:
            raise
    for subdir in subdirs:
        if not isdir(realpath(pjoin(outputPath, subdir))):
            try:
                os.makedirs(realpath(pjoin(outputPath, subdir)))
            except OSError:
                raise


def doTrackGeneration(configFile):
    """
    Do the tropical cyclone track generation.

    The track generation settings are read from *configFile*.
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
    Do the wind field calculations.

    The wind field settings are read from *configFile*.
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
    for the model calibration step.
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
    Plot the data.
    """
    import matplotlib
    matplotlib.use('Agg')  # Use matplotlib backend

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

    indLonLat = flLoadFile(pjoin(processPath, 'cyclone_tracks'),
                           delimiter=',')
    indicator = indLonLat[:, 0]
    lonData = indLonLat[:, 1]
    latData = indLonLat[:, 2]

    from PlotInterface.plotStats import PlotData
    plotting = PlotData(statsPlotPath, "png")

    log.info('Plotting pressure data')
    pbar.update(0.05)

    plotting.plotPressure(pAllData, pRateData)
    plotting.scatterHistogram(
        pAllData[1:], pAllData[:-1], 'prs_scatterHist', allpos=True)
    plotting.scatterHistogram(
        pRateData[1:], pRateData[:-1], 'prsRate_scatterHist')
    plotting.minPressureHist(indicator, pAllData)
    plotting.minPressureLat(pAllData, latData)


    log.info('Plotting bearing data')
    pbar.update(0.15)

    plotting.plotBearing(bAllData, bRateData)

    log.info('Plotting speed data')
    pbar.update(0.25)

    plotting.plotSpeed(sAllData, sRateData)

    log.info('Plotting longitude and lattitude data')
    pbar.update(0.45)

    plotting.plotLonLat(lonData, latData, indicator)

    log.info('Plotting quantiles for pressure, bearing, and speed')
    pbar.update(0.65)

    plotting.quantile(pRateData, "Pressure")
    plotting.quantile(bRateData, "Bearing")
    plotting.quantile(sRateData, "Speed")

    log.info('Plotting frequency data')
    pbar.update(0.85)

    try:
        freq = flLoadFile(pjoin(processPath, 'frequency'))
        years = freq[:, 0]
        frequency = freq[:, 1]
        plotting.plotFrequency(years, frequency)
    except IOError:
        log.warning("No frequency file available - skipping this stage")

    pbar.update(1.0)

@disableOnWorkers
def doStatistics(configFile):
    """
    Calibrate the model.
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


@disableOnWorkers
def doHazard(configFile):
    """
    Do the hazard calculations.
    """
    log.info('Running HazardInterface')

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')
    pbar = ProgressBar('Performing hazard calculations: ', showProgressBar)

    from HazardInterface.HazardInterface import HazardInterface
    hzdinterface = HazardInterface(configFile)
    hzdinterface.calculateWindHazard()

    log.info('Completed HazardInterface')
    pbar.update(1.0)

@disableOnWorkers
def doHazardPlotting(configFile):
    """
    Do the hazard plots.
    """
    import matplotlib
    matplotlib.use('Agg')  # Use matplotlib backend

    config = ConfigParser()
    config.read(configFile)

    log.info('Plotting Hazard Maps')

    showProgressBar = config.get('Logging', 'ProgressBar')
    pbar = ProgressBar('Plotting hazard maps: ', showProgressBar)
    pbar.update(0.0)

    from PlotInterface.AutoPlotHazard import AutoPlotHazard
    plotter = AutoPlotHazard(configFile, progressbar=pbar)
    plotter.plot()

    pbar.update(1.0)


@timer
def main(configFile='main.ini'):
    """
    Main interface of TCRM that allows control and interaction with the
    5 interfaces: DataProcess, StatInterface, TrackGenerator,
    WindfieldInterface and HazardInterface

    Input: configFile  - file containing configuration
                         settings for running TCRM
    Output: None
    Example: main('main.ini')
    """

    log.info("Starting TCRM")
    log.info("Configuration file: %s", configFile)

    doOutputDirectoryCreation(configFile)

    config = ConfigParser()
    config.read(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'DownloadData'):
        doDataDownload(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'DataProcess'):
        doDataProcessing(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'ExecuteStat'):
        doStatistics(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'ExecuteTrackGenerator'):
        doTrackGeneration(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'ExecuteWindfield'):
        doWindfieldCalculations(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'ExecuteHazard'):
        doHazard(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'PlotData'):
        doDataPlotting(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'PlotHazard'):
        doHazardPlotting(configFile)

    pp.barrier()

    log.info('Completed TCRM')


def startup():
    parser = argparse.ArgumentParser()
    parser.add_argument('configFile', help='the configuration file')
    parser.add_argument('-v', '--verbose', help='verbose output',
                        action='store_true')
    args = parser.parse_args()

    configFile = args.configFile

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
            logfile = pjoin(os.getcwd(), 'main.log')

    logLevel = config.get('Logging', 'LogLevel')
    verbose = config.getboolean('Logging', 'Verbose')

    if args.verbose:
        verbose = True

    if not verbose:
        logLevel = 'ERROR'
        verbose = True

    attemptParallel()

    if pp.size() > 1 and pp.rank() > 0:
        logfile += '-' + str(pp.rank())
        verbose = False  # to stop output to console
    else:
        codeStatus = status()
        print __doc__ + codeStatus

    flStartLog(logfile, logLevel, verbose)

    # Switch off minor warning messages
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")

    try:
        main(configFile)
    except Exception:  # pylint: disable=W0703
        # Catch any exceptions that occur and log them (nicely):
        tblines = traceback.format_exc().splitlines()
        for line in tblines:
            log.critical(line.lstrip())


if __name__ == "__main__":
    startup()
