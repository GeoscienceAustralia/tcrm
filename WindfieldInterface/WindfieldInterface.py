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

from Utilities.files import flModDate, flProgramVersion, flGetStat
from Utilities.config import cnfGetIniValue
from Utilities.columns import colReadCSV
from Utilities.process import pGetProcessedFiles, pAlreadyProcessed, pWriteProcessedFile
import Utilities.nctools as nctools
from Utilities.progressbar import ProgressBar


class WindfieldInterface:


    def __init__(self, config_file='WindfieldInterface.ini', nfiles=1, 
                 show_progress_bar=False):
        """
        Initiate and run generateWindfield for nfiles simulations. If
        there is only a single simulation, the module will automatically
        plot the resulting windfield, and save the U and V components.
        Otherwise only the gust wind speed will be stored and the rest
        of the data discarded (this may change in future versions).
        """

        self.config_file = config_file
        show_progress_bar = cnfGetIniValue( self.config_file, 'Logging', 
                                           'ProgressBar', True )
        self.pbar = ProgressBar("(4/6) Generating windfields: ", 
                                show_progress_bar)
        self.logger = logging.getLogger()
        self.logger.info("Initialising WindfieldInterface")
        self.datFile = cnfGetIniValue(self.config_file, 'Process', 'DatFile',
                                      config_file.split('.')[0] + '.dat')
        
        self.excludePastProcessed = cnfGetIniValue(self.config_file, 'Process',
                                                   'ExcludePastProcessed', True)
        
        self.source = cnfGetIniValue(self.config_file, 'WindfieldInterface', 'Source')
        self.output_path = os.path.join(cnfGetIniValue(self.config_file,
                                                       'Output', 'Path'),
                                                       'windfield')
        self.nfiles = nfiles

        self.WF = gW.generateWindfield(config_file)
        if self.nfiles == 1:
            track_file = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                       'TrackFile')

            saveAllTimes = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                          'SaveAllTimes',False)
            plot_output = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                         'PlotOutput',False)

            self.generateSingleWindfield(track_file, saveAllTimes, plot_output)
        else:
            track_path = cnfGetIniValue(self.config_file, 'WindfieldInterface', 'TrackPath',
                                        os.path.join(cnfGetIniValue(self.config_file,
                                                                    'Output', 'Path'),
                                                                    'tracks'))
                                       
            self.logger.info("Processing %d files in %s"%(self.nfiles, track_path))
            track_files = os.listdir(track_path)
            # For convenience, put them in some sort of order:
            track_files.sort()

            if len(track_files) < nfiles:
                # For the case where someone is unable to count:
                self.logger.warning("Number of available track files is less than requested number")
                self.logger.warning("Resetting requested number to %d"%len(track_files))
                # Reset number of files to the available number of track files in the directory:
                self.nfiles = len(track_files)
            if self.nfiles == 0:
                self.logger.error("No track files available - not proceeding with WindfieldInterface")
                sys.exit(1)
            if self.nfiles == 1:
                # Handle the case of someone only putting a single track file in the directory:
                track_file = os.path.join(track_path, track_files[0])
                self.generateSingleWindfield(track_file)
            else:
                self.generateMultipleWindfields(track_path, track_files)


    def __doc__(self):
        return "An interface between main and generateWindfield modules"

    def load_file(self, track_file):
        """
        Load the data from a track file into an array suitable for
        digestion by generateWindfield.  The data source is defined in
        the initialisation function (above). At this time, we only pass
        a bare minimum of defining variables to generateWindfield. As
        other components become available, we may include options to
        vary parameters such as the beta value or thetaMax.
        """
        self.logger.debug("Loading data from %s"%track_file)

        track_data = colReadCSV(self.config_file, track_file, self.source)

        index = array(track_data['index'], int)
        age = array(track_data['age'], float)
        lon = array(track_data['lon'], float)
        lat = array(track_data['lat'], float)
        speed = array(track_data['speed'], float)
        penv = array(track_data['penv'], float)
        pressure = array(track_data['pressure'], float)

        try:
            speed_units = cnfGetIniValue(self.config_file, self.source, 
                                        'SpeedUnits', 'mps')
        except NameError:
            pass
        else:
            speed = metutils.convert(speed, speed_units, 'mps')

        try:
            pr_units = cnfGetIniValue(self.config_file, self.source, 
                                      'PressureUnits', 'hPa')
        except NameError:
            pass
        else:
            pressure = metutils.convert(pressure, pr_units, 'Pa')
            penv = metutils.convert(penv, pr_units, 'Pa')
        # At this time, assume the storm motion direction is given in
        # degrees true.  Thus the storm motion direction must be
        # converted to radians.
        bearing = array(track_data['bearing'])*pi/180.
        bearing = maputils.bearing2theta(bearing)

        rmax = array(track_data['rmax'], float)

        return transpose(array([index, age, lon, lat, speed, bearing,
                                pressure, penv, rmax]))

    def generateSingleWindfield(self, track_file, saveAllTimes=False,
                                plot_output=True):
        """
        Run generatWindfield for a single track file
        Automatically saves gust, U and V data, and plots resulting
        gust speeds
        """

        data = self.load_file(track_file)

        speed, UU, VV, pressure, lon, lat = self.WF.calcwind(data, saveAllTimes)

        # Save final data (i.e. maximum wind speeds, minimum pressures)
        # to netCDF file
        dimensions = {0:{'name':'lat', 'values':lat, 'dtype':'f',
                             'atts':{'long_name':'Latitude',
                                         'units':'degrees_north'} },
                      1:{'name':'lon', 'values':lon, 'dtype':'f',
                             'atts':{'long_name':'Longitude',
                                         'units':'degrees_east'} } }
        variables = {0:{'name':'vmax', 'dims':('lat', 'lon'),
                        'values':array(speed), 'dtype':'f',
                        'atts':{'long_name':'Maximum 3-second gust wind speed',
                                'units':'m/s'} },
                     1:{'name':'ua', 'dims':('lat', 'lon'),
                        'values':array(UU), 'dtype':'f',
                        'atts':{'long_name':'Maximum eastward wind',
                                'units':'m/s'} },
                     2:{'name':'va', 'dims':('lat', 'lon'),
                        'values':array(VV), 'dtype':'f',
                        'atts':{'long_name':'Maximum northward wind',
                                'units':'m/s'} },
                     3:{'name':'psl', 'dims':('lat', 'lon'),
                        'values':array(pressure), 'dtype':'f',
                        'atts':{'long_name':'Minimum air pressure at sea level',
                                'units':'hPa'} } }
        nctools.ncSaveGrid(os.path.join(self.output_path, 'tc.nc'), 
                               dimensions, variables)

        # Plot up the final results for verifying the output:
        if plot_output:
            [gridX, gridY] = meshgrid(lon, lat)
            plotWindfield(gridX, gridY, speed, title="Windfield", 
                          fileName=os.path.join(self.output_path, 
                                                'windfield.png'))
                                                
            plotPressurefield(gridX, gridY, pressure, title="Pressure field", 
                              fileName=os.path.join(self.output_path, 
                                                    'pressure.png'))

        # If timeseries data was collected, save that as well:
        if cnfGetIniValue(self.config_file, 'Timeseries', 'Extract', False):
            ts_path = os.path.join(cnfGetIniValue(self.config_file, 
                                                  'Output', 'Path'), 
                                  'process', 'timeseries')
            plot_path = os.path.join(cnfGetIniValue(self.config_file, 
                                                    'Output', 'Path'), 'plots')
            plotTimeseries(ts_path, plot_path)
        return

    def generateMultipleWindfields(self, directory, track_files):
        """
        Run generateWindfield multiple times for multiple track files
        Only saves the maximum gust speed, and only plots output if
        requested in the configuration settings.
        After completing each file, write it's details to a dat file to
        record that it has been processed.
        """

        rc = pGetProcessedFiles(self.datFile)

        for n in range(self.nfiles):
            input_file = os.path.join(directory, track_files[n])

            directory, fname, md5sum, moddate = flGetStat(input_file)
            if pAlreadyProcessed(directory, fname, 'md5sum', md5sum) & self.excludePastProcessed:
                self.logger.info("Already processed %s"%input_file)
            else:
                self.logger.info("Processing %s"%input_file)
                data = self.load_file(input_file)
                output_file = os.path.join(self.output_path, 'gust.%04d.nc'%n)
                speed, UU, VV, pressure, lon, lat = self.WF.calcwind(data, False)

                inputFileDate = flModDate( input_file )
                gatts = {'history':'TCRM hazard simulation - synthetic event wind field',
                         'version':flProgramVersion( ),
                         'Python_ver':sys.version,
                         'track_file':'%s (modified %s)'%( input_file, inputFileDate ),
                         'radial_profile':self.WF.profileType,
                         'boundary_layer':self.WF.windFieldType,
                         'beta':self.WF.beta}
                dimensions = {0:{'name':'lat', 'values':lat, 'dtype':'f',
                                 'atts':{'long_name':'Latitude',
                                         'units':'degrees_north'} },
                              1:{'name':'lon', 'values':lon, 'dtype':'f',
                                 'atts':{'long_name':'Longitude',
                                         'units':'degrees_east'} } }
                variables = {0:{'name':'vmax', 'dims':('lat', 'lon'),
                                'values':array(speed), 'dtype':'f',
                                'atts':{'long_name':'Maximum 3-second gust wind speed',
                                        'units':'m/s'} } }
                nctools.ncSaveGrid(output_file, dimensions, variables, 
                                   gatts=gatts)

                if cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                  'PlotOutput', False):
                    [gridX, gridY] = meshgrid(lon, lat)                  
                    plotWindfield(gridX, gridY, speed, title="Windfield",
                                  fileName=os.path.join(self.output_path, 
                                                        'windfield.%04d.png'%(n)))
                pWriteProcessedFile(input_file)
            self.pbar.update((n+1)/float(self.nfiles))

