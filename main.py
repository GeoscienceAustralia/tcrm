#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

Version :$Rev: 652 $

$Id: main.py 652 2011-10-31 05:35:50Z nsummons $

Copyright - Geoscience Australia, 2008
"""

import os, sys, pdb, logging, traceback

from Utilities.config import cnfGetIniValue
from Utilities.files import flConfigFile, flStartLog, flLoadFile
from Utilities import pathLocator
from DataProcess.CalcTrackDomain import CalcTrackDomain
from DataProcess.CalcFrequency import CalcFrequency
from Utilities.progressbar import ProgressBar

# Set Basemap data path if compiled with py2exe
if pathLocator.is_frozen():
    os.environ['BASEMAPDATA'] = os.path.join(pathLocator.getRootDirectory(), 'mpl-data', 'data')

__version__ = "$Id: main.py 652 2011-10-31 05:35:50Z nsummons $"

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

    logger = logging.getLogger()
    # Temporarily define a level so that we can write to the log file
    # some general info about the program...
    logging.addLevelName(100, '')
    logger.log(100, "Starting TCRM")

    logger.info("Configuration file: %s"%configFile)

    # Set up output folders - at this time, we create *all* output
    # folders that may be required:
    outputPath = cnfGetIniValue(configFile, 'Output', 'Path',
                                os.path.join(os.getcwd(), 'output'))
    logger.info("Output will be stored under %s"%outputPath)
    subdirs = ['tracks', 'hazard', 'windfield', 'plots', 'plots/hazard', 'plots/stats',
               'log', 'process', 'process/timeseries', 'process/dat']

    if not os.path.isdir(outputPath):
        try:
            os.makedirs(outputPath)
        except OSError:
            raise
    for subdir in subdirs:
        if not os.path.isdir(os.path.realpath(os.path.join(outputPath, subdir))):
            try:
                os.makedirs(os.path.realpath(os.path.join(outputPath, subdir)))
            except OSError:
                raise

    # Execute Data Process:
    if cnfGetIniValue(configFile, 'Actions', 'DataProcess', False):
        pbar = ProgressBar('(1/6) Processing data files: ')
        
        from DataProcess import DataProcess
        logger.info('Running Data Processing')
        dP = DataProcess.DataProcess(configFile, progressbar=pbar)
        dP.processData()

        statsPlotPath = os.path.join(outputPath, 'plots', 'stats')
        processPath = os.path.join(outputPath, 'process')
        
        pRateData = flLoadFile(os.path.join(processPath,'pressure_rate'))
        pAllData = flLoadFile(os.path.join(processPath,'all_pressure'))
        bRateData = flLoadFile(os.path.join(processPath,'bearing_rate'))
        bAllData = flLoadFile(os.path.join(processPath,'all_bearing'))
        sRateData = flLoadFile(os.path.join(processPath,'speed_rate'))
        sAllData = flLoadFile(os.path.join(processPath,'all_speed'))
        pbar.update(0.625)
        
        indLonLat = flLoadFile(os.path.join(processPath,'cyclone_tracks'),delimiter=',')
        indicator = indLonLat[:,0]
        lonData = indLonLat[:,1]
        latData = indLonLat[:,2]

        import PlotInterface.plotStats as plotStats
        plotStats.plotPressure(pAllData,pRateData,statsPlotPath)
        plotStats.minPressureHist(indicator,pAllData,statsPlotPath)
        plotStats.minPressureLat(pAllData,latData,statsPlotPath)
        pbar.update(0.75)
        
        plotStats.plotBearing(bAllData,bRateData,statsPlotPath)
        plotStats.plotSpeed(sAllData,sRateData,statsPlotPath)
        #plotStats.plotSpeedBear(sAllData,bAllData,statsPlotPath)
        plotStats.plotLonLat(lonData,latData,indicator,statsPlotPath)
        pbar.update(0.875)
        
        plotStats.quantile(pRateData,statsPlotPath,"Pressure")
        plotStats.quantile(bRateData,statsPlotPath,"Bearing")
        plotStats.quantile(sRateData,statsPlotPath,"Speed")
        try:
            freq = flLoadFile(os.path.join(processPath,'frequency'))
        except IOError:
            logger.warning("No frequency file available - skipping this stage")
        else:
            years = freq[:,0]
            frequency = freq[:,1]
            plotStats.plotFrequency(years,frequency,statsPlotPath)
        logger.info('Completed Data Processing')
        pbar.update(1.0)

    # Execute StatInterface:
    if cnfGetIniValue(configFile, 'Actions', 'ExecuteStat', False):
        logger.info('Running StatInterface')
        # Auto-calculate track generator domain
        pbar = ProgressBar('(2/6) Calculating statistics:')
        CalcTD = CalcTrackDomain(configFile)
        autoCalc_TG_gridLimit = CalcTD.calc()

        from StatInterface import StatInterface
        statInterface = StatInterface.StatInterface(configFile, autoCalc_gridLimit=autoCalc_TG_gridLimit, progressbar=pbar)
        statInterface.kdeOrigin()
        pbar.update(0.7)
        statInterface.cdfCellBearing()
        pbar.update(0.8)
        statInterface.cdfCellSpeed()
        pbar.update(0.9)
        statInterface.cdfCellPressure()
        if cnfGetIniValue(configFile, 'RMW', 'GetRMWDistFromInputData', False):
            statInterface.cdfCellSize()
        pbar.update(1.0)
        logger.info('Completed StatInterface')

    # Execute TrackGenerator:
    if cnfGetIniValue(configFile, 'Actions', 'ExecuteTrackGenerator', False):
        logger.info('Running TrackGenerator')
        # Auto-calculate track generator domain if not calculated previously
        try:
            autoCalc_TG_gridLimit
        except NameError:
            CalcTD = CalcTrackDomain(configFile)
            autoCalc_TG_gridLimit = CalcTD.calc()        
        from TrackGenerator.trackSimulation import trackSimulation
        numSimulations = cnfGetIniValue(configFile, 'TrackGenerator',
                                        'NumSimulations', 50)
        dt = cnfGetIniValue(configFile, 'TrackGenerator', 'TimeStep', 1)
        tsteps = cnfGetIniValue(configFile, 'TrackGenerator',
                                'NumTimeSteps', 360)
        trackPath = os.path.join(outputPath, 'tracks')
        yrsPerSim = cnfGetIniValue(configFile, 'TrackGenerator',
                                   'YearsPerSimulation', 10)
        frequency = cnfGetIniValue(configFile, 'TrackGenerator', 'Frequency', '')
        if frequency == '':
            logger.info('No genesis frequency specified -> auto-calculating')
            CalcF = CalcFrequency(configFile, autoCalc_TG_gridLimit)
            frequency = CalcF.calc()
            logger.info("Estimated annual genesis frequency for Track Generator domain: %s"%frequency)
        format = cnfGetIniValue(configFile, 'TrackGenerator',
                                'Format', 'csv')
        trackSimulation(configFile, numSimulations, frequency, yrsPerSim, trackPath,
                        format, dt=dt, tsteps=tsteps, autoCalc_gridLimit=autoCalc_TG_gridLimit)
        logger.info('Completed TrackGenerator')

    # Execute Windfield:
    if cnfGetIniValue(configFile, 'Actions', 'ExecuteWindfield', False):
        logger.info('Running WindfieldInterface')
        from WindfieldInterface import WindfieldInterface
        nfiles = cnfGetIniValue(configFile, 'WindfieldInterface',
                                'NumberofFiles', 1)
        wfinterface = WindfieldInterface.WindfieldInterface(configFile, nfiles)
        logger.info('Completed WindfieldInterface')

    # Execute Hazard:
    if cnfGetIniValue(configFile, 'Actions', 'ExecuteHazard', False):
        logger.info('Running HazardInterface')
        from HazardInterface.HazardInterface import HazardInterface
        hzdinterface = HazardInterface(configFile)
        hzdinterface.calculateWindHazard()
        logger.info('Completed HazardInterface')

    # Plot Hazard:
    if cnfGetIniValue(configFile, 'Actions', 'PlotHazard', False):
        logger.info('Plotting Hazard Maps')
        from PlotInterface.AutoPlotHazard import AutoPlotHazard
        plotHazard = AutoPlotHazard(configFile)
        plotHazard.plot()

    logger.info('Completed TCRM')
    #exits program if not do matplotlib plots

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = 'main.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg
    
    tcrm_dir = pathLocator.getRootDirectory()
    os.chdir(tcrm_dir)

    logFile = cnfGetIniValue(configFile, 'Logging', 'LogFile', 'main.log')
    logFileDir = os.path.dirname(os.path.realpath(logFile))
    # If log file directory does not exist, create it
    if not os.path.isdir(logFileDir):
        try:
            os.makedirs(logFileDir)
        except OSError:
            logFile = os.path.join(os.getcwd(), 'main.log')
    
    flStartLog(logFile,
               cnfGetIniValue(configFile, 'Logging', 'LogLevel', 'INFO'),
               cnfGetIniValue(configFile, 'Logging', 'Verbose', False))
    
    # Switch off minor warning messages
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="pytz")
    
    try:
        main(configFile)
    except:
        # Catch any exceptions that occur and log them (nicely):
        tblines = traceback.format_exc().splitlines()
        logger = logging.getLogger()
        for line in tblines:
            logger.critical(line.lstrip())
        #sys.exit(1)
