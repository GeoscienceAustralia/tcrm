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

import os
import io
import sys
import logging as log
import traceback
import matplotlib

from os.path import join as pjoin, dirname
from ConfigParser import RawConfigParser

from Utilities.config import cnfGetIniValue
from Utilities.files import flConfigFile, flStartLog, flLoadFile
from Utilities import pathLocator
from DataProcess.CalcTrackDomain import CalcTrackDomain
from DataProcess.CalcFrequency import CalcFrequency
from Utilities.progressbar import ProgressBar

matplotlib.use('Agg') # Use matplotlib backend

# Set Basemap data path if compiled with py2exe
if pathLocator.is_frozen():
    os.environ['BASEMAPDATA'] = os.path.join(pathLocator.getRootDirectory(), 'mpl-data', 'data')

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

    from WindfieldInterface import WindfieldInterface
    nfiles = cnfGetIniValue(configFile, 'WindfieldInterface', 'NumberofFiles', 1)
    wfinterface = WindfieldInterface.WindfieldInterface(configFile, nfiles)

    log.info('Completed wind field calculations')


def main(config_file='main.ini'):
    """
    Main interface of TCRM that allows control and interaction with the
    5 interfaces: DataProcess, StatInterface, TrackGenerator,
    WindfieldInterface and HazardInterface
    Input: configFile  - file containing configuration
                         settings for running TCRM
    Output: None
    Example: main('main.ini')
    """

    logger = log.getLogger()
    # Temporarily define a level so that we can write to the log file
    # some general info about the program...
    log.addLevelName(100, '')
    logger.log(100, "Starting TCRM")

    logger.info("Configuration file: %s"%config_file)

    # Set up output folders - at this time, we create *all* output
    # folders that may be required:
    output_path = cnfGetIniValue(config_file, 'Output', 'Path',
                                os.path.join(os.getcwd(), 'output'))
    logger.info("Output will be stored under %s"%output_path)
    subdirs = ['tracks', 'hazard', 'windfield', 'plots', 'plots/hazard',
               'plots/stats', 'log', 'process', 'process/timeseries',
               'process/dat']

    if not os.path.isdir(output_path):
        try:
            os.makedirs(output_path)
        except OSError:
            raise
    for subdir in subdirs:
        if not os.path.isdir(os.path.realpath(os.path.join(output_path, subdir))):
            try:
                os.makedirs(os.path.realpath(os.path.join(output_path, subdir)))
            except OSError:
                raise

    show_progress_bar = cnfGetIniValue( config_file, 'Logging', 'ProgressBar', True )
    # Execute Data Process:
    if cnfGetIniValue(config_file, 'Actions', 'DataProcess', False):
        pbar = ProgressBar('(1/6) Processing data files: ', show_progress_bar )

        from DataProcess import DataProcess
        logger.info('Running Data Processing')
        data_process = DataProcess.DataProcess(config_file, progressbar=pbar)
        data_process.processData()

        statsPlotPath = os.path.join(output_path, 'plots', 'stats')
        processPath = os.path.join(output_path, 'process')

        pRateData = flLoadFile(os.path.join(processPath, 'pressure_rate'))
        pAllData = flLoadFile(os.path.join(processPath, 'all_pressure'))
        bRateData = flLoadFile(os.path.join(processPath, 'bearing_rate'))
        bAllData = flLoadFile(os.path.join(processPath, 'all_bearing'))
        sRateData = flLoadFile(os.path.join(processPath, 'speed_rate'))
        sAllData = flLoadFile(os.path.join(processPath, 'all_speed'))
        pbar.update(0.625)

        indLonLat = flLoadFile(os.path.join(processPath, 'cyclone_tracks'),
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
        #plotStats.plotSpeedBear(sAllData, bAllData, statsPlotPath)
        plotting.plotLonLat(lonData, latData, indicator)
        pbar.update(0.875)

        plotting.quantile(pRateData, "Pressure")
        plotting.quantile(bRateData, "Bearing")
        plotting.quantile(sRateData, "Speed")
        try:
            freq = flLoadFile(os.path.join(processPath, 'frequency'))
        except IOError:
            logger.warning("No frequency file available - skipping this stage")
        else:
            years = freq[:, 0]
            frequency = freq[:, 1]
            plotting.plotFrequency(years, frequency)
        logger.info('Completed Data Processing')
        pbar.update(1.0)

    # Execute StatInterface:
    if cnfGetIniValue(config_file, 'Actions', 'ExecuteStat', False):
        logger.info('Running StatInterface')
        # Auto-calculate track generator domain
        pbar = ProgressBar('(2/6) Calculating statistics:', show_progress_bar )
        CalcTD = CalcTrackDomain(config_file)
        TG_domain = CalcTD.calcDomainFromFile()

        from StatInterface import StatInterface
        statInterface = StatInterface.StatInterface(config_file,
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
        if cnfGetIniValue(config_file, 'RMW', 'GetRMWDistFromInputData', False):
            statInterface.cdfCellSize()
        pbar.update(1.0)
        logger.info('Completed StatInterface')

    # Execute TrackGenerator:
    if cnfGetIniValue(config_file, 'Actions', 'ExecuteTrackGenerator', False):
        doTrackGeneration(config_file)

    # Execute Windfield:
    if cnfGetIniValue(config_file, 'Actions', 'ExecuteWindfield', False):
        doWindfieldCalculations(config_file)

    # Execute Hazard:
    if cnfGetIniValue(config_file, 'Actions', 'ExecuteHazard', False):
        logger.info('Running HazardInterface')
        from HazardInterface.HazardInterface import HazardInterface
        hzdinterface = HazardInterface(config_file)
        hzdinterface.calculateWindHazard()
        logger.info('Completed HazardInterface')

    # Plot Hazard:
    if cnfGetIniValue(config_file, 'Actions', 'PlotHazard', False):
        logger.info('Plotting Hazard Maps')
        pbar = ProgressBar('(6/6) Plotting results:      ', show_progress_bar )
        from PlotInterface.AutoPlotHazard import AutoPlotHazard
        plot_hazard = AutoPlotHazard(config_file, progressbar=pbar)
        plot_hazard.plot()

    logger.info('Completed TCRM')
    #exits program if not do matplotlib plots

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

         # no pypar, create a dummy one

         class DummyPypar(object):
             def size(self): return 1
             def rank(self): return 0
             def barrier(self): pass

         pp = DummyPypar()


if __name__ == "__main__":

    attemptParallel()

    try:
        CONFIG_FILE = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        CONFIG_FILE = 'main.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(CONFIG_FILE):
            ERROR_MSG = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, ERROR_MSG
    # If config file doesn't exist => raise error
    if not os.path.exists(CONFIG_FILE):
        ERROR_MSG = "Configuration file '" + CONFIG_FILE +"' not found"
        raise IOError, ERROR_MSG

    TCRM_DIR = pathLocator.getRootDirectory()
    os.chdir(TCRM_DIR)

    LOG_FILE = cnfGetIniValue(CONFIG_FILE, 'Logging', 'LogFile', 'main.log')
    LOG_FILE_DIR = os.path.dirname(os.path.realpath(LOG_FILE))
    # If log file directory does not exist, create it
    if not os.path.isdir(LOG_FILE_DIR):
        try:
            os.makedirs(LOG_FILE_DIR)
        except OSError:
            LOG_FILE = os.path.join(os.getcwd(), 'main.log')

    logLevel = cnfGetIniValue(CONFIG_FILE, 'Logging', 'LogLevel', 'INFO')
    verbose = cnfGetIniValue(CONFIG_FILE, 'Logging', 'Verbose', False)

    if pp.rank() > 0:
        LOG_FILE += '-' + str(pp.rank())
        verbose = False # to stop output to console

    flStartLog(LOG_FILE, logLevel, verbose)

    # Switch off minor warning messages
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    warnings.filterwarnings("ignore", category=UserWarning, module="numpy")

    try:
        main(CONFIG_FILE)
    except:
        # Catch any exceptions that occur and log them (nicely):
        TB_LINES = traceback.format_exc().splitlines()
        LOGGER = log.getLogger()
        for line in TB_LINES:
            LOGGER.critical(line.lstrip())
        #sys.exit(1)
