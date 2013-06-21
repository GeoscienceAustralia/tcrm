#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Title: main.py - main executable for TCRM
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2008-07-21 2:25:PM
Description: A text-based interface executable in Python for the
Tropical Cyclone Risk Model. Allows users to interact with all available
methods from DataProcess, StatInterface, TrackGenerator,
WindfieldInterface and HazardInterface.

Usage:


Version :$Rev: 826 $

$Id: main.py 826 2012-03-26 02:06:55Z nsummons $

Copyright - Geoscience Australia, 2008
"""

__version__ = "$Id: main.py 826 2012-03-26 02:06:55Z nsummons $"

import gc
gc.set_debug(gc.DEBUG_SAVEALL)

import os
import io
import sys
import logging
import traceback
import matplotlib

from os.path import join as pjoin, dirname

from Utilities.files import flConfigFile, flStartLog, flLoadFile
from Utilities import pathLocator
from Utilities.progressbar import ProgressBar
from config import ConfigParser

matplotlib.use('Agg')  # Use matplotlib backend

# Set Basemap data path if compiled with py2exe
if pathLocator.is_frozen():
    os.environ['BASEMAPDATA'] = pjoin(
        pathLocator.getRootDirectory(), 'mpl-data', 'data')

log = logging.getLogger('main')

def attemptParallel():
    """
    Attempt to load Pypar globally as `pp`. If Pypar loads successfully, then a
    call to `pypar.finalize` is registered to be called at exit of the Python
    interpreter. This is to ensure that MPI exits cleanly.

    If pypar cannot be loaded then a dummy `pp` is created.
    """
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
            def size(self):
                return 1

            def rank(self):
                return 0

            def barrier(self):
                pass

        pp = DummyPypar()


def disableOnWorkers(f):
    def wrap(*args, **kwargs):
        if pp.size() > 1 and pp.rank() > 0:
            return
        else:
            return f(*args, **kwargs)
    return wrap


@disableOnWorkers
def doOutputDirectoryCreation(configFile):
    # Set up output folders - at this time, we create *all* output
    # folders that may be required:

    config = ConfigParser()
    config.read(configFile)

    outputPath = config.get('Output', 'Path')

    log.info('Output will be stored under %s' % outputPath)

    subdirs = ['tracks', 'hazard', 'windfield', 'plots', 'plots/hazard',
               'plots/stats', 'log', 'process', 'process/timeseries',
               'process/dat']

    if not os.path.isdir(outputPath):
        try:
            os.makedirs(outputPath)
        except OSError:
            raise
    for subdir in subdirs:
        if not os.path.isdir(os.path.realpath(pjoin(outputPath, subdir))):
            try:
                os.makedirs(os.path.realpath(pjoin(outputPath, subdir)))
            except OSError:
                raise


def doTrackGeneration(configFile):
    """
    Do the tropical cyclone track generation.

    The track generation settings are read from *configFile*.
    """

    log.info('Starting track generation')

    import TrackGenerator
    TrackGenerator.run(configFile)

    log.info('Completed track generation')


def doWindfieldCalculations(configFile):
    """
    Do the wind field calculations.

    The wind field settings are read from *configFile*.
    """

    log.info('Starting wind field calculations')

    import WindfieldInterface
    WindfieldInterface.run(configFile)

    log.info('Completed wind field calculations')


@disableOnWorkers
def doDataProcessing(configFile):

    config = ConfigParser()
    config.read(configFile)

    outputPath      = config.get('Output', 'Path')
    showProgressBar = config.get('Logging', 'ProgressBar')
 
    pbar = ProgressBar('(1/6) Processing data files: ', showProgressBar)

    log.info('Running Data Processing')

    from DataProcess.DataProcess import DataProcess
    data_process = DataProcess(configFile, progressbar=pbar)
    data_process.processData()

    log.info('Completed Data Processing')
    pbar.update(1.0)


@disableOnWorkers
def doDataPlotting(configFile):
    config = ConfigParser()
    config.read(configFile)

    outputPath      = config.get('Output', 'Path')
    showProgressBar = config.get('Logging', 'ProgressBar')
 
    pbar = ProgressBar('Plotting data files: ', showProgressBar)

    statsPlotPath = pjoin(outputPath, 'plots', 'stats')
    processPath = pjoin(outputPath, 'process')

    pRateData = flLoadFile(pjoin(processPath, 'pressure_rate'))
    pAllData = flLoadFile(pjoin(processPath, 'all_pressure'))
    bRateData = flLoadFile(pjoin(processPath, 'bearing_rate'))
    bAllData = flLoadFile(pjoin(processPath, 'all_bearing'))
    sRateData = flLoadFile(pjoin(processPath, 'speed_rate'))
    sAllData = flLoadFile(pjoin(processPath, 'all_speed'))
    pbar.update(0.625)

    indLonLat = flLoadFile(pjoin(processPath, 'cyclone_tracks'),
                           delimiter=',')
    indicator = indLonLat[:, 0]
    lonData = indLonLat[:, 1]
    latData = indLonLat[:, 2]

    from PlotInterface.plotStats import PlotData
    plotting = PlotData(statsPlotPath, "png")
    plotting.plotPressure(pAllData, pRateData)
    plotting.minPressureHist(indicator, pAllData)
    plotting.minPressureLat(pAllData, latData)
    pbar.update(0.75)

    plotting.plotBearing(bAllData, bRateData)
    plotting.plotSpeed(sAllData, sRateData)
    # plotStats.plotSpeedBear(sAllData, bAllData, statsPlotPath)
    plotting.plotLonLat(lonData, latData, indicator)
    pbar.update(0.875)

    plotting.quantile(pRateData, "Pressure")
    plotting.quantile(bRateData, "Bearing")
    plotting.quantile(sRateData, "Speed")
    try:
        freq = flLoadFile(pjoin(processPath, 'frequency'))
    except IOError:
        log.warning("No frequency file available - skipping this stage")
    else:
        years = freq[:, 0]
        frequency = freq[:, 1]
        plotting.plotFrequency(years, frequency)

    pbar.update(1.0)

@disableOnWorkers
def doStatistics(configFile):
    from DataProcess.CalcTrackDomain import CalcTrackDomain

    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')
    getRMWDistFromInputData = config.getboolean('RMW', 'GetRMWDistFromInputData')

    log.info('Running StatInterface')
    # Auto-calculate track generator domain
    pbar = ProgressBar('(2/6) Calculating statistics:', showProgressBar)
    CalcTD = CalcTrackDomain(configFile)
    TG_domain = CalcTD.calcDomainFromFile()

    from StatInterface import StatInterface
    statInterface = StatInterface.StatInterface(configFile,
                                                autoCalc_gridLimit=TG_domain,
                                                progressbar=pbar)
    statInterface.kdeGenesisDate()
    pbar.update(0.6)
    statInterface.kdeOrigin()
    pbar.update(0.7)
    statInterface.cdfCellBearing()
    pbar.update(0.8)
    statInterface.cdfCellSpeed()
    pbar.update(0.9)

    statInterface.cdfCellPressure()

    if getRMWDistFromInputData: 
        statInterface.cdfCellSize()

    pbar.update(1.0)
    log.info('Completed StatInterface')


@disableOnWorkers
def doHazard(configFile):
    log.info('Running HazardInterface')
    from HazardInterface.HazardInterface import HazardInterface
    hzdinterface = HazardInterface(configFile)
    hzdinterface.calculateWindHazard()
    log.info('Completed HazardInterface')


@disableOnWorkers
def doHazardPlotting(configFile):
    config = ConfigParser()
    config.read(configFile)

    showProgressBar = config.get('Logging', 'ProgressBar')

    log.info('Plotting Hazard Maps')
    pbar = ProgressBar('(6/6) Plotting results:      ', showProgressBar)
    from PlotInterface.AutoPlotHazard import AutoPlotHazard
    plot_hazard = AutoPlotHazard(configFile, progressbar=pbar)
    plot_hazard.plot()


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
    log.info("Configuration file: %s" % configFile)

    doOutputDirectoryCreation(configFile)

    config = ConfigParser()
    config.read(configFile)

    pp.barrier()

    if config.getboolean('Actions', 'DataProcess'):
        doDataProcessing(configFile)
        #doDataPlotting(configFile)

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

    if config.getboolean('Actions', 'PlotHazard'):
        doHazardPlotting(configFile)

    pp.barrier()

    log.info('Completed TCRM')

if __name__ == "__main__":

    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = 'main.ini'
        # If no filename is specified and default filename doesn't exist =>
        # raise error
        if not os.path.exists(configFile):
            ERROR_MSG = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError(ERROR_MSG)
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        ERROR_MSG = "Configuration file '" + configFile + "' not found"
        raise IOError(ERROR_MSG)

    TCRM_DIR = pathLocator.getRootDirectory()
    os.chdir(TCRM_DIR)

    config = ConfigParser()
    config.read(configFile)

    LOG_FILE = config.get('Logging', 'LogFile')
    LOG_FILE_DIR = os.path.dirname(os.path.realpath(LOG_FILE))
    # If log file directory does not exist, create it
    if not os.path.isdir(LOG_FILE_DIR):
        try:
            os.makedirs(LOG_FILE_DIR)
        except OSError:
            LOG_FILE = pjoin(os.getcwd(), 'main.log')

    logLevel = config.get('Logging', 'LogLevel')
    verbose = config.getboolean('Logging', 'Verbose')

    attemptParallel()

    if pp.size() > 1 and pp.rank() > 0:
        LOG_FILE += '-' + str(pp.rank())
        verbose = False  # to stop output to console

    flStartLog(LOG_FILE, logLevel, verbose)

    # Switch off minor warning messages
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    warnings.filterwarnings("ignore", category=UserWarning, module="numpy")

    try:
        main(configFile)
    except:
        # Catch any exceptions that occur and log them (nicely):
        TB_LINES = traceback.format_exc().splitlines()
        for line in TB_LINES:
            log.critical(line.lstrip())
        # sys.exit(1)

    #for o in gc.garbage:
    #    print('%s %s' % (type(o), str(o)))
