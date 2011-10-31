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


Title: generateStats.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 11/27/07 2:58:PM
Description: Calculates the mean, variance and autocorrelation
statistics for all grid cells in the domain, for a given parameter
(speed, bearing, etc).  The statistics calculated herein provide
information used in generating synthetic TC tracks. The method is based
on the approached described by Hall & Jewson (2007)

References:
Hall, T.M. and Jewson, S. (2007): Statistical modelling of North
Atlantic tropical cyclone tracks.
Tellus A, 59(4), doi:10.1111/j.1600-0870.2007.00240.x, pp486-498

Constraints:
SeeAlso:

Version: $Rev: 644 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 10/04/08 11:42:AM
Modification: Changed logging method

$Id: generateStats.py 644 2011-10-31 05:32:50Z nsummons $
"""

import os, pdb, logging, sys, math

from numpy import *
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile
import Utilities.stats as stats


__version__ = '$Id: generateStats.py 644 2011-10-31 05:32:50Z nsummons $'

class parameters:
    """parameters:

    Description: Create an object that holds arrays of the statistical
    properties of each grid cell. There are arrays for both land and sea cells.

    Parameters:
    None

    Members:
    mu/lmu : Array of mean values of parameter for each grid cell
    sig/lsig : Array of variance values of parameter for each grid cell
    alpha/lalpha : Array of autoregression coefficients of parameter for
    each grid cell
    phi/lphi : Array of normalisation values for random variations

    Methods:
    None

    Internal methods:
    None
    """

    def __init__(self,numCells):
        self.mu = zeros(numCells)
        self.sig = zeros(numCells)
        self.alpha = zeros(numCells)
        self.phi = zeros(numCells)
        self.min = zeros(numCells)
        self.lmu = zeros(numCells)
        self.lsig = zeros(numCells)
        self.lalpha = zeros(numCells)
        self.lphi = zeros(numCells)
        self.lmin = zeros(numCells)

class GenerateStats:
    """GenerateStats:

    Description:

    Parameters:
    parameter: Array or String (representing filename)
               Contains the data which the statistical values will be
               based
    lonLat:    Array or String (representing filename)
               Contains the longitude and latitude of each of the
               observations in the array parameter
    gridLimit: dictionary containing limits of regional grid
               Contains keys 'xMin', 'xMax', 'yMin', 'yMax'
    gridSpace: Dictionary containing spacing of grid
               Contains keys 'x' and 'y'
    gridInc:   Dictionary containing increments for grid to be used when
               insufficient observations exist in the cell being
               analysed.
               Contain keys 'x' and 'y'
    minSample: Integer
               The minimum required number of observations in a given
               cell to calculate meaningful statistics.
    angular:   Boolean
               True = the data in array 'parameter' is angular in nature
               (e.g. bearings)
               False = the data in array 'parameter' is not angular.
    Members:
    Methods:
    Internal methods:
    """
    def __init__(self, configFile, parameter, lonLat, gridLimit,
                 gridSpace, gridInc, minSample=100, angular=False, missingValue=sys.maxint, 
                 progressbar=None, prgStartValue=0, prgEndValue=1):
        self.logger = logging.getLogger()
        self.logger.info('Initialising GenerateStats')

        if type(lonLat) is str:
            self.lonLat = array(flLoadFile(lonLat, delimiter=','))
        else:
            self.lonLat = lonLat
        if type(parameter) is str:
            self.param = array(flLoadFile(parameter))
        else:
            self.param = parameter

        self.gridLimit = gridLimit
        self.gridSpace = gridSpace
        self.gridInc = gridInc
        self.maxCell = stats.maxCellNum(self.gridLimit, self.gridSpace)
        #landmask = cnfGetIniValue(configFile, 'Input', 'LandMask')
        #lsfraction = stats.statCellFraction(self.gridLimit, self.gridSpace,
        #                                    landmask)
        self.minSample = minSample
        self.coeffs = parameters(self.maxCell+1)
        self.angular = angular
        self.missingValue = missingValue

        self.domain_warning_raised = False

        for i in range(self.maxCell + 1):
            self.coeffs.mu[i], self.coeffs.sig[i], self.coeffs.alpha[i], \
            self.coeffs.phi[i],self.coeffs.min[i] = self.calculate(i,False)
            self.coeffs.lmu[i], self.coeffs.lsig[i], self.coeffs.lalpha[i], \
            self.coeffs.lphi[i], self.coeffs.lmin[i] = self.calculate(i,True)
            if mod(i, 10) == 0:  # Periodically update progress bar
                if progressbar is not None:
                    progressbar.update((i+1)/float(self.maxCell+1), prgStartValue, prgEndValue)
        if progressbar is not None:
            progressbar.update(1.0, prgStartValue, prgEndValue)

    def calculate(self, cellNum, onLand):
        """calculate(cellNum):
        Calculate the required statistics (mean, variance,
        autocorrelation and regularized anomaly coefficient) for the
        given cell.
        """
        p = self.extractParameter(cellNum, onLand)
        if self.angular:
            mu = stats.circmean(radians(p))
            sig = stats.circstd(radians(p))
        else:
            mu = mean(p)
            sig = std(p)
        # Calculate the autocorrelations:
        alphas = correlate(p, p, 'full')
        n = len(p)
        # Grab only the lag-one autocorrelation coeff.
        alpha = alphas[n]/alphas.max()
        phi = sqrt(1 - alpha**2)
        mn = min(p)
        return mu, sig, alpha, phi, mn

    def extractParameter(self, cellNum, onLand):
        """extractParameter(cellNum):
        Extracts the cyclone parameter data for the given cell.
        If the population of a cell is insufficient for generating a
        PDF, the bounds of the cell are expanded until the population is
        sufficient.

        Null/missing values are removed.
        """
        if not stats.validCellNum(cellNum, self.gridLimit, self.gridSpace):
            self.logger.critical("Invalid input on cellNum: cell number %i is out of range"%cellNum)
            raise InvalidArguments, 'Invalid input on cellNum: cell number %i is out of range'%cellNum
        cellLon, cellLat = stats.getCellLonLat(cellNum, self.gridLimit,
                                               self.gridSpace)
        wLon = cellLon
        eLon = cellLon + self.gridSpace['x']
        nLat = cellLat
        sLat = cellLat - self.gridSpace['y']

        lon = self.lonLat[:,0]
        lat = self.lonLat[:,1]
        lsflag = self.lonLat[:,2]
        if onLand:
            ij = where(((lat >= sLat) & (lat < nLat)) &
                       (lon >= wLon) & (lon < eLon) & (lsflag>0))
        else:
            ij = where(((lat >= sLat) & (lat < nLat)) &
                       (lon >= wLon) & (lon < eLon) & (lsflag==0))
        p_ = self.param[ij]
        p = stats.statRemoveNum(array(p_), self.missingValue)

        while size(p) <= self.minSample:
            wLon_last = wLon
            eLon_last = eLon
            nLat_last = nLat
            sLat_last = sLat
            wLon, eLon, nLat, sLat = self._expandCell(lon, lat, wLon, eLon,
                                                      nLat, sLat)
            # Check if grid has reached maximum extent
            if (wLon == wLon_last) & (eLon == eLon_last) & (nLat == nLat_last) & (sLat == sLat_last):
                if onLand:
                    if not self.domain_warning_raised:
                        self.domain_warning_raised = True
                        self.logger.warning("Insufficient grid points over land in selected domain to estimate storm statistics - reverting to statistics for open ocean.")
                    return self.extractParameter(cellNum, False)
                else:
                    errMsg = "Insufficient grid points in selected domain to estimate storm statistics - please select a larger domain."
                    self.logger.critical(errMsg)
                    raise StopIteration, errMsg
            if onLand:
                ij = where(((lat >= sLat) & (lat < nLat)) & (lon >= wLon) &
                           (lon < eLon) & (lsflag>0))
            else:
                ij = where(((lat >= sLat) & (lat < nLat)) & (lon >= wLon) &
                           (lon < eLon) & (lsflag==0))
            p_ = self.param[ij]
            p = stats.statRemoveNum(array(p_), self.missingValue)

        # Check to see if all values in the array are the same. If the values
        # are the same, bandwidth would be 0, and therefore KDE cannot be generated
        while p.max() == p.min():
            wLon_last = wLon
            eLon_last = eLon
            nLat_last = nLat
            sLat_last = sLat
            wLon, eLon, nLat, sLat = self._expandCell(lon, lat, wLon, eLon,
                                                      nLat, sLat)
            # Check if grid has reached maximum extent
            if (wLon == wLon_last) & (eLon == eLon_last) & (nLat == nLat_last) & (sLat == sLat_last):
                if onLand:
                    if not self.domain_warning_raised:
                        self.domain_warning_raised = True
                        self.logger.warning("Insufficient grid points over land in selected domain to estimate storm statistics - reverting to statistics for open ocean.")
                    return self.extractParameter(cellNum, False)
                else:
                    errMsg = "Insufficient grid points in selected domain to estimate storm statistics - please select a larger domain."
                    self.logger.critical(errMsg)
                    raise StopIteration, errMsg
            if onLand:
                ij = where(((lat >= sLat) & (lat < nLat)) &
                           (lon >= wLon) & (lon < eLon) & (lsflag>0))
            else:
                ij = where(((lat >= sLat) & (lat < nLat)) &
                           (lon >= wLon) & (lon < eLon) & (lsflag==0))

            p_ = self.param[ij]
            p = stats.statRemoveNum(array(p_), self.missingValue)
        return p

    def _expandCell(self, lon, lat, wLon, eLon, nLat, sLat):
        """_expandCell(lon, lat, wLon, eLon, nLat, sLat):
        Obtain the indices of observations from adjacent cells.
        This is called when there are insufficient observations
        in a cell to generate a PDF.
        """
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
        if wLon < self.gridLimit['xMin']:
            wLon = self.gridLimit['xMin']
        if eLon > self.gridLimit['xMax']:
            eLon = self.gridLimit['xMax']
        if nLat > self.gridLimit['yMax']:
            nLat = self.gridLimit['yMax']
        if sLat < self.gridLimit['yMin']:
            sLat = self.gridLimit['yMin']
        return wLon, eLon, nLat, sLat

if __name__=="__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename does not exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file does not exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s',
                        filename=cnfGetIniValue(configFile, 'Logging',
                                                'LogFile',
                                                __file__.rstrip('.py') + '.ini'),
                        filemode='w')

    gridLimit = eval(cnfGetIniValue(configFile, 'Parameters', 'gridLimit'))
    gridSpace = eval(cnfGetIniValue(configFile, 'Parameters', 'gridSpace'))
    gridInc = eval(cnfGetIniValue(configFile, 'Parameters', 'gridInc'))
    path = cnfGetIniValue(configFile, 'Output', 'Path')
    latLon = os.path.join(path, 'all_lon_lat')
    pS = GenerateStats(configFile, os.path.join(path, 'all_pressure'), latLon,
                       gridLimit, gridSpace, gridInc)
    vS = GenerateStats(configFile, os.path.join(path, 'all_speed'), latLon,
                       gridLimit, gridSpace, gridInc)
    bS = GenerateStats(configFile, os.path.join(path, 'all_bearing'), latLon,
                       gridLimit, gridSpace, gridInc, angular=True)
    sS = GenerateStats(configFile, os.path.join(path, 'all_rmax'), latLon,
                       gridLimit, gridSpace, gridInc)
    dvS = GenerateStats(configFile, os.path.join(path, 'speed_rate'), latLon,
                        gridLimit, gridSpace, gridInc)
    dpS = GenerateStats(configFile, os.path.join(path, 'pressure_rate'), latLon,
                        gridLimit, gridSpace, gridInc)
    dbS = GenerateStats(configFile, os.path.join(path, 'bearing_rate'), latLon,
                        gridLimit, gridSpace, gridInc, angular=True)
    dsS = GenerateStats(configFile, os.path.join(path, 'rmax_rate'), latLon,
                        gridLimit, gridSpace, gridInc)
