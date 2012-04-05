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


Title: StatInterface.py

Author: Geoff Xu
CreationDate: 2005-12-14
Description:  Defines the class for the interface of all n-by-n
              statistical classes for representing historical cyclone
              data. This is a wrapper class that allows the user to
              interact with all the classes related to generating new
              samples of cyclone data from historical data.

ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2006-10-31
Modification: Added descriptive header and metadata

Version: $Rev: 499$
ModifiedBy: Nariman Habili, nariman.habili@ga.gov.au
ModifiedDate: 2006-11-29
Modification: KDEParameters and SamplingParameters instantiated directly.
              File names Stat parameters passed by dictionaries instead
              of individual file names.
              Upgrade to ndarray
              Conformance with style guide

Version: 530
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-07-04
Modification: Removed deprecated parameterType variable. Only seen in
              initialisation.

Version: 626
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-10-02 11:23:AM
Modification:

Version: 285
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-04-10
Modification: Changed logging method

Version: $Rev: 810 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-08-19 10:03:AM
Modification: Removed generation of unneccesary CDF files (namely
              all_cell_cdf_<parameter> and all_cell_cdf_<parameter>_rate
              The only required CDF's are for initial values for sampling
              at the initial time step. All other values are calculated
              using autoregressive model.
SeeAlso: (related programs)
Constraints:
References:

$Id: StatInterface.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

import KDEOrigin
import KDEParameters
#import SamplingOrigin
import GenerateDistributions

from Utilities.config import cnfGetIniValue

class StatInterface:
    """
    Parameters
    ----------
    None

    Members
    -------
    lon_lat : string (file name including path)
        latitude and longitude data for cyclone origins
    init_bearing : string (file name including path)
        initial bearings of cyclones
    bearing_no_init : string (file name including path)
        bearings of cyclones with no initial bearing
    all_bearing : string (file name including path)
        all bearings of cyclones
    init_speed : string (file name including path)
        initial speeds of cyclones
    speed_no_init : string (file name including path)
        speeds of cyclones with no initial speed
    all_speed : string (file name including path)
        all speeds of cyclones
    init_pressure : string (file name including path)
        initial pressures of cyclones
    pressure_no_init : string (file name including path)
        pressures of cyclones with no initial pressure
    all_pressure : string (file name including path)
        all pressures of cyclones
    kdeType : string
        kernel density estimation type
    ns : int
        number of samples

    Methods
    -------
    kdeBearing()
        generate KDEs relating to the bearings of cyclones
    kdeSpeed()
        generate KDEs relating to the speeds of cyclones
    kdePressure()
        generate KDEs relating to the pressures of cyclones
    kdeCell(cellNo)
        generate cyclone parameter KDEs relating to a particular cell

    """
    def __init__(self, configFile, autoCalc_gridLimit=None, progressbar=None):
        """
        Initialize the data and variables required for the interface
        """
        self.configFile = configFile
        self.logger = logging.getLogger()
        self.progressbar = progressbar
        self.logger.info("Initialising StatInterface")

        self.kdeType = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'kdeType', 'Biweight')
        self.kde2DType = cnfGetIniValue(self.configFile, 'StatInterface',
                                        'kde2DType', 'Gaussian')
        self.ns = cnfGetIniValue(self.configFile, 'StatInterface',
                                 'Samples', 50000)
        minSamplesCell = cnfGetIniValue(self.configFile, 'StatInterface',
                                        'minSamplesCell', 100)
        self.kdeStep = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'kdeStep', 0.2)
        self.outputPath = cnfGetIniValue(self.configFile, 'Output', 'Path')
        self.processPath = os.path.join(self.outputPath, 'process')
        missingValue = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'MissingValue', sys.maxint)

        gridLimitStr = cnfGetIniValue(self.configFile, 'StatInterface', 'gridLimit', '')
        if gridLimitStr is not '':
            try:
                self.gridLimit = eval(gridLimitStr)
            except SyntaxError:
                self.logger.exception('Error! gridLimit is not a dictionary' )
        else:
            self.gridLimit = autoCalc_gridLimit
            self.logger.info('No gridLimit specified - using automatic selection: ' + str(self.gridLimit))
        
        try:
            gridSpace = eval(cnfGetIniValue(self.configFile, 'StatInterface',
                                            'gridSpace'))
            gridInc = eval(cnfGetIniValue(self.configFile, 'StatInterface',
                                          'gridInc'))
        except SyntaxError:
            self.logger.exception('Error! gridSpace or gridInc not dictionaries' )
            raise

        self.generateDist = GenerateDistributions.GenerateDistributions( \
                               self.configFile, self.gridLimit, gridSpace, gridInc,
                               self.kdeType, minSamplesCell, missingValue)

    def __doc__(self):
        """
        Documentation on what this class does
        """
        return 'A wrapper class that allows users to interact with all the\
               classes relating with generating new samples of cyclone data from\
               historical data as well as statistics from the historical data'

    def kdeOrigin(self):
        """
        Generate 2D KDEs relating to the origin of cyclones
        """
        self.logger.info('Generating 2D PDF of TC origins')
        self.logger.debug('Reading data from %s'%os.path.join(self.processPath,'origin_lon_lat'))
        
        kde = KDEOrigin.KDEOrigin(self.configFile, self.kde2DType, self.gridLimit, 0.1, progressbar=self.progressbar)
        kde.generateKDE(None, save=True)
        kde.generateCdf()

    # def samplingOrigin(self):
        # """
        # sampling of initial cyclone origns
        # """
        # self.logger.info('Sampling TC origins')
        # self.logger.debug('Reading data from %s', os.path.join(self.processPath, 'originPDF.txt'))
        # self.logger.debug('Outputting data into %s', os.path.join(self.processPath, 'sample_origin_lonlat'))

        # samplingOrigin = SamplingOrigin.SamplingOrigin(kdeOrigin=os.path.join(self.processPath, 'originPDF.txt'))
        # samplingOrigin.generateSamples(self.ns,os.path.join(self.processPath, 'sample_origin_lonlat'))
        
    def kdeGenesisDate(self):
        """
        Generate CDFs relating to the genesis day-of-year of cyclones
        """
        self.logger.info('Generating CDFs for TC genesis day')
        self.logger.debug('Reading data from %s', os.path.join( self.processPath, 'jdays' ) )
        jdays = os.path.join( self.processPath, 'jdays' )
        kde = KDEParameters.KDEParameters( self.kdeType )
        kde.generateGenesisDateCDF( jdays, bw=14, genesisKDE=os.path.join( self.processPath, 'cdfGenesisDays' ) )

    def cdfCellBearing(self):
        """
        generate CDFs relating to the bearing of cyclones
        """
        self.logger.info('Generating CDFs for TC bearing')
        self.logger.debug('Reading data from %s', os.path.join(self.processPath, 'init_lon_lat'))
        self.logger.debug('Reading data from %s', os.path.join(self.processPath, 'init_bearing'))
        self.logger.debug('Outputting data into %s', os.path.join(self.processPath, 'all_cell_cdf_init_bearing'))

        lonLat = os.path.join(self.processPath, 'init_lon_lat')
        pList = os.path.join(self.processPath, 'init_bearing')
        self.generateDist.allDistributions(lonLat, pList, 'init_bearing', 1, True)

    def cdfCellSpeed(self):
        """
        generate CDFs relating to the pressures of cyclones
        """
        self.logger.info('Generating CDFs for TC speed')
        self.logger.debug('Reading data from %s',
                           os.path.join(self.processPath, 'init_lon_lat'))
        self.logger.debug('Reading data from %s',
                           os.path.join(self.processPath, 'init_speed'))
        self.logger.debug('Outputting data into %s',
                           os.path.join(self.processPath,
                                        'all_cell_cdf_init_speed'))

        lonLat = os.path.join(self.processPath, 'init_lon_lat')
        pList = os.path.join(self.processPath, 'init_speed')
        self.generateDist.allDistributions(lonLat, pList, 'init_speed',
                                           self.kdeStep)

    def cdfCellPressure(self):
        """
        generate CDFs relating to the pressures of cyclones
        """
        self.logger.info('Generating CDFs for TC pressure')
        self.logger.debug('Reading data from %s',
                           os.path.join(self.processPath, 'origin_lon_lat'))
        self.logger.debug('Reading data from %s',
                           os.path.join(self.processPath, 'init_pressure'))
        self.logger.debug('Outputting data into %s',
                           os.path.join(self.processPath,
                                        'all_cell_cdf_init_pressure'))

        lonLat = os.path.join(self.processPath, 'origin_lon_lat')
        pList = os.path.join(self.processPath, 'init_pressure')
        self.generateDist.allDistributions(lonLat, pList, 'init_pressure',
                                           self.kdeStep)

    def cdfCellSize(self):
        """
        generate CDFs relating to the pressures of cyclones
        """
        self.logger.info('Generating CDFs for TC size')
        self.logger.debug('Reading data from %s',
                           os.path.join(self.processPath, 'origin_lon_lat'))
        self.logger.debug('Reading data from %s',
                           os.path.join(self.processPath, 'init_rmax'))
        self.logger.debug('Outputting data into %s',
                           os.path.join(self.processPath,
                                        'all_cell_cdf_init_rmax'))

        lonLat = os.path.join(self.processPath, 'origin_lon_lat')
        pList = os.path.join(self.processPath, 'init_rmax')
        self.generateDist.allDistributions(lonLat, pList, 'init_rmax',
                                           self.kdeStep)


if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg

    logging.basicConfig(level=logging.getattr(cnfGetIniValue(configFile, 'Logging','LogLevel','DEBUG')),
                        format='%(asctime)s %(name)-15s: %(levelname)-8s %(message)s',
                        filename=cnfGetIniValue(configFile,'Logging','LogFile',__file__.rstrip('.py') + '.log'),
                        filemode='w')
    si = StatInterface(configFile)

    si.kdeOrigin()
    #si.samplingOrigin()
    si.cdfCellBearing()
    si.cdfCellSpeed()
    si.cdfCellPressure()
    si.cdfCellSize()
