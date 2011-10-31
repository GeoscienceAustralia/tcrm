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


Title: GenerateDistributions.py - generate the distributions of
       parameters.

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-12-07
Description: Generate the cumulative distribution functions (CDF's) for
             a given parameter for each cell in the lat-lon grid
             (defined by gridLimit and gridSpace).  This uses the method
             of kernel density estimators to determine the
             distributions.

Version: $Rev: 644 $

Version: $Rev: 644 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2006-12-14
Modification: Added _expandCell and _checkGridLimits to permit
              increasing the population of cells that initially have
              insufficient observations to generate a quality PDF.

ModifiedBy: Nariman Habili, nariman.habili@ga.gov.au
ModifiedDate: 2007-03-22
Modification: lonLat, parameterList, parameterName are now pass passed
              in at allDistributions.

Version: $Rev: 644 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 10/04/08 11:42:AM
Modification: Changed logging method

SeeAlso: (related programs)
Constraints:

$Id: GenerateDistributions.py 644 2011-10-31 05:32:50Z nsummons $
"""

import os, sys, pdb, logging

import Utilities.stats as stats
import KDEParameters
import pylab
from scipy import array, arange, where, size, transpose, concatenate
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile, flSaveFile, flStartLog
from Utilities.AsyncRun import AsyncRun


class GenerateDistributions:
    """
    Description: Generate the cumulative distribution functions (CDF's)
    for a given parameter for each cell in the lat-lon grid (defined by
    gridLimit and gridSpace).  This uses the method of kernel density
    estimators to determine the distributions.

    Parameters:
    files : dictionary
        Contains input and output files read from config files
    lonLat : string or 2d-array
        If string, a filename containing the lon/lat of all cyclone
        observations
        If array, contains lon/lat pairs of all cyclone observations
    parameterList: string or array
        If string, a filename containing the values of a parameter for
        all cyclone observations.
        If array, contains all values of a parameter for all cyclone
        observations
    paramterName : string
        Name of the parameter being examined
    gridLimit : dictionary
        lon/lat bounds of the grid over which the parameter
        distributions will be determined.
    gridSpace : dictionary
        lon/lat grid spacing (in degrees)
    kdeType : string
        Name of the kernel density estimator to use

    Members:
    files : dictionary
        Contains input and output files read from config files
    lonLat : string or 2d-array
        If string, a filename containing the lon/lat of all cyclone
        observations
        If array, contains lon/lat pairs of all cyclone observations
    parameterList: string or array
        If string, a filename containing the values of a parameter for
        all cyclone observations.
        If array, contains all values of a parameter for all cyclone
        observations
    paramterName : string
        Name of the parameter being examined
    gridLimit : dictionary
        lon/lat bounds of the grid over which the parameter
        distributions will be determined.
    gridSpace : dictionary
        lon/lat grid spacing (in degrees)
    kdeType : string
        Name of the kernel density estimator to use

    Methods:
    allDistributions : calculate and save the CDF for each cell for a
                       given parameter.
    extractParameter : Identify the data in the full dataset
                       corresponding to the cell under investigation and
                       return the values of the given parameter, with
                       null values removed.

    Internal Methods:
    _expandCell : Obtain the indices of observations from adjacent
                  cells.  This is called when there are insufficient
                  observations in a cell to generate a PDF.
    _checkCellLimits : Check the bounds of the cell are within the
                       limits of the region under investigation - if not
                       then set to the limits
    """
    
    def __init__(self, configFile, gridLimit, gridSpace, gridInc, kdeType,
                 minSamplesCell, missingValue=sys.maxint):
        """
        Initialise required fields
        """
        self.logger = logging.getLogger()
        
        self.logger.info("Initialising GenerateDistributions")
        self.gridSpace = gridSpace
        self.gridLimit = gridLimit
        self.gridInc = gridInc
        self.kdeType = kdeType
        self.minSamplesCell = minSamplesCell
        self.outputPath = cnfGetIniValue(configFile, 'Output', 'Path')
        self.kdeParameter = KDEParameters.KDEParameters(kdeType)

        self.missingValue = missingValue

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Generate the cumulative distribution functions for a parameter in\
        each cell of the grid over the area of investigation.'

    def allDistributions(self, lonLat, parameterList, parameterName=None,
                         kdeStep=0.1, angular=False, plotParam=False):
        """
        Calculate a distribution for each individual cell and store in a
        file or return
        """
        if parameterName:
            self.logger.debug("Running allDistributions for %s"%parameterName)
        else:
            self.logger.debug("Running allDistributions")

        if type(lonLat) is str:
            self.logger.debug("Loading lat/lon data from file")
            self.lonLat = array(flLoadFile(lonLat, delimiter=','))
        else:
            self.lonLat = lonLat

        if type(parameterList) is str:
            self.logger.debug("Loading parameter data from file: %s" %
                          parameterList)
            self.pList = array(flLoadFile(parameterList))
        else:
            self.pList = parameterList

        self.pName = parameterName

        maxCellNum = stats.maxCellNum(self.gridLimit, self.gridSpace)

        # Writing CDF dataset for all individual cell number into files
        self.logger.debug("Writing CDF dataset for all individual cells into files")

        for cellNum in xrange(0, maxCellNum + 1):
            self.logger.debug("Processing cell number %i"%cellNum)

            # Generate cyclone parameter data for the cell number
            self.extractParameter(cellNum)

            # Estimate cyclone parameter data using KDE
            # The returned array contains the grid, the PDF and the CDF
            cdf = self.kdeParameter.generateKDE(self.parameter, kdeStep,
                                                angular=angular)
            if plotParam:
                self._plotParameter(cellNum, kdeStep)
            self.logger.debug('size of parameter array = %d: size of cdf array = %d'
                          % (self.parameter.size,cdf.size))

            cellNumlist = []
            for i in range(len(cdf)):
                cellNumlist.append(cellNum)
            if cellNum == 0:
                results = transpose(array([cellNumlist, cdf[:,0], cdf[:,2]]))
            else:
                self.logger.debug('size of results array = %s'%str(results.size))
                results = concatenate((results, transpose(array([cellNumlist,
                                                                 cdf[:,0],
                                                                 cdf[:,2]]))))

        if parameterName == None:
            self.logger.debug("Returning CDF dataset for all individual cell numbers")
            return results
        else:
            cdfHeader = "Cell_Number, CDF_" + self.pName + "_x, CDF_" + \
                        self.pName + "_y"
            allCellCdfOutput = os.path.join(self.outputPath,
                                            'process',
                                            'all_cell_cdf_' + self.pName)
            args = {"filename":allCellCdfOutput, "data":results,
                    "header":cdfHeader, "delimiter":",", "fmt":"%f"}
            fl = AsyncRun(flSaveFile, args)
            fl.start()
            self.logger.debug("Writing CDF dataset for all individual cell numbers into files")

    def extractParameter(self, cellNum):
        """extractParameter(cellNum):
        Extracts the cyclone parameter data for the given cell.
        If the population of a cell is insufficient for generating a
        PDF, the bounds of the cell are expanded until the population is
        sufficient.

        Null/missing values are removed.
        """
        if not stats.validCellNum(cellNum, self.gridLimit, self.gridSpace):
            self.logger.critical("Invalid input on cellNum: cell number %i is out of range"%cellNum)
            raise InvalidArguments, 'Invalid input on cellNum: cell number is out of range'
        lon = self.lonLat[:,0]
        lat = self.lonLat[:,1]
        cellLon, cellLat = stats.getCellLonLat(cellNum, self.gridLimit,
                                               self.gridSpace)

        wLon = cellLon
        eLon = cellLon + self.gridSpace['x']
        nLat = cellLat
        sLat = cellLat - self.gridSpace['y']

        indij = where(((lat >= sLat) & (lat < nLat)) &
                      (lon >= wLon) & (lon < eLon))
        parameter_ = self.pList[indij]
        self.parameter = stats.statRemoveNum(array(parameter_),
                                             self.missingValue)

        while size(self.parameter ) <= self.minSamplesCell:
            self.logger.debug("Insufficient samples. Increasing the size of the cell")
            wLon_last = wLon
            eLon_last = eLon
            nLat_last = nLat
            sLat_last = sLat
            wLon, eLon, nLat, sLat = self._expandCell(lon, lat, wLon, eLon,
                                                      nLat, sLat)
            if (wLon == wLon_last) & (eLon == eLon_last) & (nLat == nLat_last) & (sLat == sLat_last):
                errMsg = "Insufficient grid points in selected domain to estimate storm statistics - please select a larger domain."
                self.logger.critical(errMsg)
                raise StopIteration, errMsg
            indij = where(((lat >= sLat) & (lat < nLat)) &
                          ((lon >= wLon) & (lon < eLon)))
            parameter_ = self.pList[indij]
            self.parameter = stats.statRemoveNum(array(parameter_),
                                                 self.missingValue)

        # Check to see if all values in the array are the same. If the
        # values are the same, bandwidth would be 0, and therefore KDE
        # cannot proceed
        while self.parameter.max() == self.parameter.min():
            self.logger.debug("Parameter values appear to be the same. Increasing the size of the cell")
            wLon_last = wLon
            eLon_last = eLon
            nLat_last = nLat
            sLat_last = sLat
            wLon, eLon, nLat, sLat = self._expandCell(lon, lat, wLon,
                                                      eLon, nLat, sLat)
            if (wLon == wLon_last) & (eLon == eLon_last) & (nLat == nLat_last) & (sLat == sLat_last):
                errMsg = "Insufficient grid points in selected domain to estimate storm statistics - please select a larger domain."
                self.logger.critical(errMsg)
                raise StopIteration, errMsg
            indij = where(((lat >= sLat) & (lat < nLat)) &
                          ((lon >= wLon) & (lon < eLon)))
            parameter_ = self.pList[indij]
            self.parameter = stats.statRemoveNum(array(parameter_),
                                                 self.missingValue)
        self.logger.debug("Number of valid observations in cell %s : %s" %
                      (str(cellNum), str(size(self.parameter))))


    def _plotParameter(self, cellNum, kdeStep):
        self.logger.debug("Plotting %s"%self.pName)
        pMin = self.parameter.min()
        pMax = self.parameter.max()
        rng = arange(pMin, pMax, kdeStep)
        if len(rng) < 10:
            rng = 10
        pylab.clf()
        pylab.hist(self.parameter,rng)
        pylab.title('Parameter: '+self.pName+' cell: '+str(cellNum))
        pylab.savefig(self.outputPath+self.pName+'.'+str(cellNum)+'.png')


    def _expandCell(self, lon, lat, wLon, eLon, nLat, sLat):
        """_expandCell(lon, lat, wLon, eLon, nLat, sLat):
        Obtain the indices of observations from adjacent cells.
        This is called when there are insufficient observations
        in a cell to generate a PDF.
        """
        self.logger.debug("Increasing cell size")
        wLon -= self.gridInc['x']
        eLon += self.gridInc['x']
        nLat += self.gridInc['y']
        sLat -= self.gridInc['y']

        wLon, eLon, nLat, sLat = self._checkGridLimits(wLon, eLon, nLat, sLat)

        return wLon, eLon, nLat, sLat

    def _checkGridLimits(self, wLon, eLon, nLat, sLat):
        """_checkGridLimits(wLon, eLon, nLat, sLat):
        Check the bounds of the cell are within the limits of the
        region under investigation - if not then set to the limits
        """
        self.logger.debug("Checking cell does not extend outside domain")
        if wLon < self.gridLimit['xMin']:
            wLon = self.gridLimit['xMin']
        if eLon > self.gridLimit['xMax']:
            eLon = self.gridLimit['xMax']
        if nLat > self.gridLimit['yMax']:
            nLat = self.gridLimit['yMax']
        if sLat < self.gridLimit['yMin']:
            sLat = self.gridLimit['yMin']
        return wLon, eLon, nLat, sLat

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

    flStartLog(cnfGetIniValue(configFile, 'Logging', 'LogFile', __file__.rstrip('.py') + '.log'),
               cnfGetIniValue(configFile, 'Logging', 'LogLevel', 'DEBUG'),
               cnfGetIniValue(configFile, 'Logging', 'Verbose', True))
    path = cnfGetIniValue(configFile, 'Output', 'Path')

    gridLimit = eval(cnfGetIniValue(configFile, 'Settings', 'gridLimit'))
    gridSpace = eval(cnfGetIniValue(configFile, 'Settings', 'gridSpace'))
    gridInc = eval(cnfGetIniValue(configFile, 'Settings', 'gridInc'))
    kdeType = cnfGetIniValue(configFile, 'Parameters', 'kdeType')
    kdeStep = cnfGetIniValue(configFile, 'Parameters', 'kdeStep', 0.1)
    minSamplesCell = cnfGetIniValue(configFile, 'Parameters',
                                    'minSamplesCell', 100)
    missingValue = cnfGetIniValue(configFile, 'StatInterface',
                                  'MissingValue', sys.maxint)
    gDist = GenerateDistributions(configFile, gridLimit, gridSpace, gridInc, kdeType,
                                  minSamplesCell, missingValue)


#    gDist.allDistributions(path+'all_lon_lat', path+'bearing_rate', 'bearing_rate', 5.)
#    gDist.allDistributions(path+'all_lon_lat', path+'pressure_rate', 'pressure_rate',kdeStep)
#    gDist.allDistributions(path+'all_lon_lat', path+'speed_rate', 'speed_rate',kdeStep)
#    gDist.allDistributions(path+'all_lon_lat', path+'rmax_rate', 'rmax_rate',kdeStep)

    gDist.allDistributions(os.path.join(path, 'init_lon_lat'),
                           os.path.join(path, 'init_bearing'),
                           'init_bearing', 5.)
    gDist.allDistributions(os.path.join(path, 'origin_lon_lat'),
                           os.path.join(path, 'init_pressure'),
                           'init_pressure', kdeStep)
    gDist.allDistributions(os.path.join(path, 'init_lon_lat'),
                           os.path.join(path, 'init_speed'),
                           'init_speed', kdeStep)

    gDist.allDistributions(os.path.join(path, 'origin_lon_lat'),
                           os.path.join(path, 'init_rmax'),
                           'init_rmax', kdeStep)
    gDist.allDistributions(os.path.join(path, 'all_lon_lat'),
                           os.path.join(path, 'all_bearing'),
                           'bearing', 5.)
    gDist.allDistributions(os.path.join(path, 'all_lon_lat'),
                           os.path.join(path, 'all_pressure'),
                           'pressure', kdeStep)
    gDist.allDistributions(os.path.join(path, 'all_lon_lat'),
                           os.path.join(path, 'all_speed'),
                           'speed', kdeStep)

    gDist.allDistributions(os.path.join(path, 'all_lon_lat'),
                           os.path.join(path, 'all_rmax'),
                           'rmax', kdeStep)
    logging.shutdown()
