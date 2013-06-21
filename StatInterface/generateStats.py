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

Version: $Rev: 810 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 10/04/08 11:42:AM
Modification: Changed logging method

$Id: generateStats.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, pdb, logging, sys, math

import numpy as np
from Utilities.files import flLoadFile
import Utilities.stats as stats
from config import ConfigParser


__version__ = '$Id: generateStats.py 810 2012-02-21 07:52:50Z nsummons $'

class parameters(object):
    """parameters:

    Description: Create an object that holds np.arrays of the statistical
    properties of each grid cell. There are np.arrays for both land and sea cells.

    Parameters:
    None

    Members:
    mu/lmu : np.array of mean values of parameter for each grid cell
    sig/lsig : np.array of variance values of parameter for each grid cell
    alpha/lalpha : np.array of autoregression coefficients of parameter for
    each grid cell
    phi/lphi : np.array of normalisation values for random variations

    Methods:
    None

    Internal methods:
    None
    """

    def __init__(self,numCells):
        self.mu = np.zeros(numCells)
        self.sig = np.zeros(numCells)
        self.alpha = np.zeros(numCells)
        self.phi = np.zeros(numCells)
        self.min = np.zeros(numCells)
        self.lmu = np.zeros(numCells)
        self.lsig = np.zeros(numCells)
        self.lalpha = np.zeros(numCells)
        self.lphi = np.zeros(numCells)
        self.lmin = np.zeros(numCells)

class GenerateStats:
    """GenerateStats:

    Description:

    Parameters:
    parameter: np.array or String (representing filename)
               Contains the data which the statistical values will be
               based
    lonLat:    np.array or String (representing filename)
               Contains the longitude and latitude of each of the
               observations in the np.array parameter
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
               True = the data in np.array 'parameter' is angular in nature
               (e.g. bearings)
               False = the data in np.array 'parameter' is not angular.
    Members:
    Methods:
    Internal methods:
    """
    def __init__(self, parameter, lonLat, gridLimit,
                 gridSpace, gridInc, minSample=100, angular=False, missingValue=sys.maxint, 
                 progressbar=None, prgStartValue=0, prgEndValue=1, calculateLater=False):

        self.logger = logging.getLogger()
        self.logger.debug('Initialising GenerateStats')

        self.gridLimit = gridLimit
        self.gridSpace = gridSpace
        self.gridInc = gridInc
        self.maxCell = stats.maxCellNum(self.gridLimit, self.gridSpace)
        self.minSample = minSample
        self.coeffs = parameters(self.maxCell+1)
        self.angular = angular
        self.missingValue = missingValue

        self.domain_warning_raised = False

        self.progressbar = progressbar
        self.prgStartValue = prgStartValue
        self.prgEndValue = prgEndValue

        if not calculateLater:
            if type(lonLat) is str:
                self.lonLat = np.array(flLoadFile(lonLat, delimiter=','))
            else:
                self.lonLat = lonLat
            if type(parameter) is str:
                self.param = np.array(flLoadFile(parameter))
            else:
                self.param = parameter

            self.calculateStatistics()

    def calculateStatistics(self):
        progressbar = self.progressbar
        prgStartValue = self.prgStartValue
        prgEndValue = self.prgEndValue

        self.logger.debug('Calculating statistics for %i cells' % self.maxCell)
        for i in range(self.maxCell + 1):
            self.coeffs.mu[i], self.coeffs.sig[i], self.coeffs.alpha[i], \
            self.coeffs.phi[i],self.coeffs.min[i] = self.calculate(i,False)
            self.coeffs.lmu[i], self.coeffs.lsig[i], self.coeffs.lalpha[i], \
            self.coeffs.lphi[i], self.coeffs.lmin[i] = self.calculate(i,True)
            if np.mod(i, 10) == 0:  # Periodically update progress bar
                if progressbar is not None:
                    progressbar.update((i+1)/float(self.maxCell+1), prgStartValue, prgEndValue)
        if progressbar is not None:
            progressbar.update(1.0, prgStartValue, prgEndValue)
        self.logger.debug('Finished calculating statistics')

    def calculate(self, cellNum, onLand):
        """calculate(cellNum):
        Calculate the required statistics (mean, variance,
        autocorrelation and regularized anomaly coefficient) for the
        given cell.
        """
        p = self.extractParameter(cellNum, onLand)
        #self.logger.debug('Calculating statistics for cell %i from %i samples' % (cellNum, len(p)))
        if self.angular:
            mu = stats.circmean(np.radians(p))
            sig = stats.circstd(np.radians(p))
        else:
            mu = np.mean(p)
            sig = np.std(p)
        # Calculate the autocorrelations:
        alphas = np.correlate(p, p, 'full')
        n = len(p)
        # Grab only the lag-one autocorrelation coeff.
        alpha = alphas[n]/alphas.max()
        phi = np.sqrt(1 - alpha**2)
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
            ij = np.where(((lat >= sLat) & (lat < nLat)) &
                       (lon >= wLon) & (lon < eLon) & (lsflag>0))
        else:
            ij = np.where(((lat >= sLat) & (lat < nLat)) &
                       (lon >= wLon) & (lon < eLon) & (lsflag==0))
        p_ = self.param[ij]
        p = stats.statRemoveNum(np.array(p_), self.missingValue)

        while np.size(p) <= self.minSample:
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
                ij = np.where(((lat >= sLat) & (lat < nLat)) & (lon >= wLon) &
                           (lon < eLon) & (lsflag>0))
            else:
                ij = np.where(((lat >= sLat) & (lat < nLat)) & (lon >= wLon) &
                           (lon < eLon) & (lsflag==0))
            p_ = self.param[ij]
            p = stats.statRemoveNum(np.array(p_), self.missingValue)

        # Check to see if all values in the np.array are the same. If the values
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
                ij = np.where(((lat >= sLat) & (lat < nLat)) &
                           (lon >= wLon) & (lon < eLon) & (lsflag>0))
            else:
                ij = np.where(((lat >= sLat) & (lat < nLat)) &
                           (lon >= wLon) & (lon < eLon) & (lsflag==0))

            p_ = self.param[ij]
            p = stats.statRemoveNum(np.array(p_), self.missingValue)
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

    def load(self, filename):
        self.logger.debug('Loading statistics from %s' % filename)
        from scipy.io.netcdf import netcdf_file
        ncdf = netcdf_file(filename, 'r')
        for var in ncdf.variables.keys():
            setattr(self.coeffs, var, ncdf.variables[var][:].flatten())
        #TODO: maybe save and check grid settings?

    def save(self, filename, description=''):
        description = ' ' + description.strip()
        self.logger.debug('Saving' + description + ' statistics to %s' % filename)

        lon = np.arange(self.gridLimit['xMin'],self.gridLimit['xMax'],self.gridSpace['x'])
        lat = np.arange(self.gridLimit['yMax'],self.gridLimit['yMin'],-1*self.gridSpace['y'])

        nx = len(lon)
        ny = len(lat)

        dimensions = {0:{'name':'lat','values':lat,'dtype':'f',
                         'atts':{'long_name':'Latitude','units':'degrees_north'} },
                      1:{'name':'lon','values':lon,'dtype':'f',
                         'atts':{'long_name':'Longitude','units':'degrees_east'} } }

        variables = {0:{'name':'mu','dims':('lat','lon'),
                        'values':self.coeffs.mu.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Mean' + description,
                                'units':'m/s'} },
                     1:{'name':'alpha','dims':('lat','lon'),
                        'values':self.coeffs.alpha.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Lag-1 autocorrelation of' + description,
                                'units':''} },
                     2:{'name':'sig','dims':('lat','lon'),
                        'values':self.coeffs.sig.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Standard deviation' + description,
                                'units':'m/s'} },
                     3:{'name':'min','dims':('lat','lon'),
                        'values':self.coeffs.min.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Minimum' + description,
                                'units':'m/s'} },
                     4:{'name':'lmu','dims':('lat','lon'),
                        'values':self.coeffs.lmu.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Mean' + description +' (over land)',
                                'units':'m/s'} },
                     5:{'name':'lalpha','dims':('lat','lon'),
                        'values':self.coeffs.lalpha.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Lag-1 autocorrelation of' + description + ' (over land)',
                                'units':''} },
                     6:{'name':'lsig','dims':('lat','lon'),
                        'values':self.coeffs.lsig.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Standard deviation of' + description + ' (over land)',
                                'units':'m/s'} },
                     7:{'name':'lmin','dims':('lat','lon'),
                        'values':self.coeffs.lmin.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'Minimum' + description + ' (over land)',
                                'units':'m/s'} },
                     8:{'name':'cell', 'dims':('lat','lon'),
                        'values':np.arange(self.maxCell+1).reshape((ny,nx)),
                        'dtype':'i',
                        'atts':{'long_name':'Cell', 'units':''} },
                     9:{'name':'phi', 'dims':('lat','lon'),
                        'values':self.coeffs.phi.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'phi', 'units':''} },
                    10:{'name':'lphi', 'dims':('lat','lon'),
                        'values':self.coeffs.lphi.reshape((ny,nx)),
                        'dtype':'f',
                        'atts':{'long_name':'land phi', 'units':''} }
                     }

        import Utilities.nctools as nctools

        nctools.ncSaveGrid(filename, dimensions, variables,
                           nodata=self.missingValue,datatitle=None,dtype='f',
                           writedata=True, keepfileopen=False)

