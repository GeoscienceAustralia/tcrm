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


Title: WindfieldInterface.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2008-12-18 08:48:14
Description: Interface between main.py and generateWindfield module for
             TCRM

SeeAlso: generateWindfield.py, main.py
Constraints:

Version: $Rev: 819 $
ModifiedDate: $Date: 2012-03-15 15:03:34 +1100 (Thu, 15 Mar 2012) $
ModifiedBy: $Author: carthur $

$Id: WindfieldInterface.py 819 2012-03-15 04:03:34Z carthur $
"""

import os, sys, pdb, logging

import generateWindfield as gW

import Utilities.metutils as metutils
import Utilities.maputils as maputils

from numpy import transpose, array, meshgrid, pi
from PlotInterface.plotWindfield import plotWindfield, plotPressurefield
from PlotInterface.plotTimeseries import plotTimeseries

from Utilities.grid import grdSave
from Utilities.files import flConfigFile, flModDate, flProgramVersion, flGetStat
from Utilities.config import cnfGetIniValue
from Utilities.columns import colReadCSV
from Utilities.process import pGetProcessedFiles, pAlreadyProcessed, pWriteProcessedFile
import Utilities.nctools as nctools
from Utilities.progressbar import ProgressBar


class WindfieldInterface:


    def __init__(self, configFile='WindfieldInterface.ini', nfiles=1, show_progress_bar=False):
        """
        Initiate and run generateWindfield for nfiles simulations. If
        there is only a single simulation, the module will automatically
        plot the resulting windfield, and save the U and V components.
        Otherwise only the gust wind speed will be stored and the rest
        of the data discarded (this may change in future versions).
        """

        self.configFile = configFile
        show_progress_bar = cnfGetIniValue( self.configFile, 'Logging', 'ProgressBar', True )
        self.pbar = ProgressBar("(4/6) Generating windfields: ", show_progress_bar)
        self.logger = logging.getLogger()
        self.logger.info("Initialising WindfieldInterface")
        self.datFile = cnfGetIniValue(self.configFile, 'Process', 'DatFile', configFile.split('.')[0] + '.dat')
        self.excludePastProcessed = cnfGetIniValue(self.configFile, 'Process', 'ExcludePastProcessed', True)
        self.source = cnfGetIniValue(self.configFile, 'WindfieldInterface', 'Source')
        self.outputPath = os.path.join(cnfGetIniValue(self.configFile, 'Output', 'Path'), 'windfield')
        self.nfiles = nfiles

        self.WF = gW.generateWindfield(configFile)
        if self.nfiles == 1:
            trackFile = cnfGetIniValue(self.configFile, 'WindfieldInterface',
                                       'TrackFile')

            saveAllTimes = cnfGetIniValue(self.configFile,'WindfieldInterface','SaveAllTimes',False)
            plotOutput = cnfGetIniValue(self.configFile,'WindfieldInterface','PlotOutput',False)

            self.generateSingleWindfield(trackFile, saveAllTimes, plotOutput)
        else:
            trackPath = cnfGetIniValue(self.configFile, 'WindfieldInterface', 'TrackPath', os.path.join(cnfGetIniValue(self.configFile, 'Output', 'Path'), 'tracks'))
            self.logger.info("Processing %d files in %s"%(self.nfiles, trackPath))
            trackFiles = os.listdir(trackPath)
            # For convenience, put them in some sort of order:
            trackFiles.sort()

            if len(trackFiles) < nfiles:
                # For the case where someone is unable to count:
                self.logger.warning("Number of available track files is less than requested number")
                self.logger.warning("Resetting requested number to %d"%len(trackFiles))
                # Reset number of files to the available number of track files in the directory:
                self.nfiles = len(trackFiles)
            if self.nfiles == 0:
                self.logger.error("No track files available - not proceeding with WindfieldInterface")
                sys.exit(1)
            if self.nfiles == 1:
                # Handle the case of someone only putting a single track file in the directory:
                trackFile = os.path.join(trackPath, trackFiles[0])
                self.generateSingleWindfield(trackFile)
            else:
                self.generateMultipleWindfields(trackPath, trackFiles)


    def __doc__(self):
        return "An interface between main and generateWindfield modules"

    def _loadFile(self,inputFile):
        """
        Load the data from a track file into an array suitable for
        digestion by generateWindfield.  The data source is defined in
        the initialisation function (above). At this time, we only pass
        a bare minimum of defining variables to generateWindfield. As
        other components become available, we may include options to
        vary parameters such as the beta value or thetaMax.
        """
        self.logger.debug("Loading data from %s"%inputFile)

        trackData = colReadCSV(self.configFile, inputFile,self.source)

        index = array(trackData['index'], int)
        age = array(trackData['age'], float)
        lon = array(trackData['lon'], float)
        lat = array(trackData['lat'], float)
        speed = array(trackData['speed'], float)
        penv = array(trackData['penv'], float)
        pressure = array(trackData['pressure'], float)

        try:
            speedUnits = cnfGetIniValue(self.configFile, self.source, 'SpeedUnits', 'mps')
        except NameError:
            pass
        else:
            speed = metutils.convert(speed, speedUnits, 'mps')

        try:
            prUnits = cnfGetIniValue(self.configFile, self.source, 'PressureUnits', 'hPa')
        except NameError:
            pass
        else:
            pressure = metutils.convert(pressure, prUnits, 'Pa')
            penv = metutils.convert(penv, prUnits, 'Pa')
        # At this time, assume the storm motion direction is given in
        # degrees true.  Thus the storm motion direction must be
        # converted to radians.
        bearing = array(trackData['bearing'])*pi/180.
        bearing = maputils.bearing2theta(bearing)

        rmax=array(trackData['rmax'],float)

        return transpose(array([index, age, lon, lat, speed, bearing,
                                pressure, penv, rmax]))

    def generateSingleWindfield(self,trackFile,saveAllTimes=False,plotOutput=True):
        """
        Run generatWindfield for a single track file
        Automatically saves gust, U and V data, and plots resulting
        gust speeds
        """

        data = self._loadFile(trackFile)

        speed, UU, VV, pressure, lon, lat = self.WF.calcwind(data, saveAllTimes)

        # Save final data (i.e. maximum wind speeds, minimum pressures)
        # to netCDF file
        dimensions = {0:{'name':'lat','values':lat,'dtype':'f','atts':{'long_name':'Latitude','units':'degrees_north'} },
                      1:{'name':'lon','values':lon,'dtype':'f','atts':{'long_name':'Longitude','units':'degrees_east'} } }
        variables = {0:{'name':'vmax','dims':('lat','lon'),
                        'values':array(speed),'dtype':'f',
                        'atts':{'long_name':'Maximum 3-second gust wind speed',
                                'units':'m/s'} },
                     1:{'name':'ua','dims':('lat','lon'),
                        'values':array(UU),'dtype':'f',
                        'atts':{'long_name':'Maximum eastward wind',
                                'units':'m/s'} },
                     2:{'name':'va','dims':('lat','lon'),
                        'values':array(VV),'dtype':'f',
                        'atts':{'long_name':'Maximum northward wind',
                                'units':'m/s'} },
                     3:{'name':'psl','dims':('lat','lon'),
                        'values':array(pressure),'dtype':'f',
                        'atts':{'long_name':'Minimum air pressure at sea level',
                                'units':'hPa'} } }
        nctools.ncSaveGrid(os.path.join(self.outputPath,'tc.nc'), dimensions, variables)

        # Plot up the final results for verifying the output:
        if plotOutput:
            [gridX,gridY] = meshgrid(lon,lat)
            plotWindfield(gridX, gridY, speed,title="Windfield", fileName=os.path.join(self.outputPath,'windfield.png'))
            plotPressurefield(gridX, gridY, pressure,title="Pressure field", fileName=os.path.join(self.outputPath,'pressure.png'))

        # If timeseries data was collected, save that as well:
        if cnfGetIniValue(self.configFile, 'Timeseries', 'Extract', False):
            tsPath = os.path.join(cnfGetIniValue(self.configFile, 'Output', 'Path'), 'process', 'timeseries')
            plotPath =os.path.join(cnfGetIniValue(self.configFile, 'Output', 'Path'), 'plots')
            plotTimeseries(tsPath, plotPath)
        return

    def generateMultipleWindfields(self, directory, trackFiles):
        """
        Run generateWindfield multiple times for multiple track files
        Only saves the maximum gust speed, and only plots output if
        requested in the configuration settings.
        After completing each file, write it's details to a dat file to
        record that it has been processed.
        """

        rc = pGetProcessedFiles(self.datFile)

        for n in range(self.nfiles):
            inputFile = os.path.join(directory, trackFiles[n])

            directory, fname, md5sum, moddate = flGetStat(inputFile)
            if pAlreadyProcessed(directory,fname, 'md5sum', md5sum) & self.excludePastProcessed:
                self.logger.info("Already processed %s"%inputFile)
                pass
            else:
                self.logger.info("Processing %s"%inputFile)
                data = self._loadFile(inputFile)
                outputFile = os.path.join(self.outputPath, 'gust.%04d.nc'%n)
                speed, UU, VV, pressure, lon, lat = self.WF.calcwind(data, False)
                #nctools.ncSaveGrid(outputFile, lon, lat, speed, 'vmax', 'm/s',longname='Maximum 3-second gust wind speed')

                inputFileDate = flModDate( inputFile )
                gatts = {'history':'TCRM hazard simulation - synthetic event wind field',
                         'version':flProgramVersion( ),
                         'Python_ver':sys.version,
                         'track_file':'%s (modified %s)'%( inputFile, inputFileDate ),
                         'radial_profile':self.WF.profileType,
                         'boundary_layer':self.WF.windFieldType,
                         'beta':self.WF.beta}
                dimensions = {0:{'name':'lat','values':lat,'dtype':'f',
                                 'atts':{'long_name':'Latitude',
                                         'units':'degrees_north'} },
                              1:{'name':'lon','values':lon,'dtype':'f',
                                 'atts':{'long_name':'Longitude',
                                         'units':'degrees_east'} } }
                variables = {0:{'name':'vmax','dims':('lat','lon'),
                                'values':numpy.array(speed),'dtype':'f',
                                'atts':{'long_name':'Maximum 3-second gust wind speed',
                                        'units':'m/s'} } }
                nctools.ncSaveGrid(os.path.join(self.outputPath,'gust.%04d.nc'%n),
                                    dimensions, variables, gatts=gatts)

                if cnfGetIniValue(self.configFile, 'WindfieldInterface',
                                  'PlotOutput', False):
                    plotWindfield(gridX, gridY, speed, title="Windfield",
                                  fileName=os.path.join(self.outputPath, 'windfield.%04d.png'%(n+1)))
                pWriteProcessedFile(inputFile)
            self.pbar.update((n+1)/float(self.nfiles))

