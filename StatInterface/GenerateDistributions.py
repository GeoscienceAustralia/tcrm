"""
:mod:`GenerateDistributions` -- generate distributions of parameters
====================================================================

.. module:: GenerateDistributions
    :synopsis: Generate CDFs for parameters for each grid cell in the
               model domain.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>


Generate the cumulative distribution functions (CDF's) for a given
parameter for each cell in the lat-lon grid (defined by gridLimit and
gridSpace). This uses the method of kernel density estimators to
determine the distributions.

"""

import os
import sys
import logging
from os.path import join as pjoin

import Utilities.stats as stats
from . import KDEParameters
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile, flSaveFile, flStartLog

from netCDF4 import Dataset
import numpy as np

class GenerateDistributions(object):
    """
    Generate the cumulative distribution functions (CDF's)
    for a given parameter for each cell in the lat-lon grid (defined by
    gridLimit and gridSpace).  This uses the method of kernel density
    estimators to determine the distributions. The methods allow for
    extraction of paramter values from the parameter files created in
    :class:`DataProcess.DataProcess`, and calculation (and saving)
    distributions.

    :param str configFile: Path to configuration file.
    :param dict gridLimit: The bounds of the model domain. The
                           :class:`dict` should contain the keys
                           :attr:`xMin`, :attr:`xMax`, :attr:`yMin`
                           and :attr:`yMax`. The *x* variable bounds
                           the longitude and the *y* variable
                           bounds the latitude.
    :param dict gridSpace: The default grid cell size. The :class:`dict`
                           should contain keys of :attr:`x` and
                           :attr:`y`. The *x* variable defines the
                           longitudinal grid size, and the *y*
                           variable defines the latitudinal size.
    :param dict gridInc: The increment in grid size, for those cells
                         that do not contain sufficient observations
                         for generating distributions. The :class:`dict`
                         should contain the keys :attr:`x` and
                         :attr:`y`. The *x* variable defines the
                         longitudinal grid increment, and the *y*
                         variable defines the latitudinal increment.
    :param str kdeType: Name of the (univariate) kernel estimator to use
                        when generating the distribution. Must be one of
                        ``Epanechnikov``, ``Gaussian``, ``Biweight`` or
                        ``Triangular``.
    :param int minSamplesCell: Minimum number of valid observations
                               required to generate distributions.
                               If insufficient observations are found in
                               a grid cell, then it is incrementally
                               expanded until ``minSamplesCell`` is
                               reached.
    :param missingValue: Missing values have this value (default
                         :attr:`sys.maxint`).

    """

    def __init__(self, configFile, gridLimit, gridSpace, gridInc, kdeType,
                 minSamplesCell=40, missingValue=sys.maxsize):
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

    def allDistributions(self, lonLat, parameterList, parameterName=None,
                         kdeStep=0.1, angular=False, periodic=False,
                         plotParam=False):
        """
        Calculate a distribution for each individual cell and store in a
        file or return the distribution.

        :param lonLat: The longitude/latitude of all observations in
                       the model domain. If a string is given, then
                       it is the path to a file containing the
                       longitude/latitude information. If an array
                       is given, then it should be a 2-d array
                       containing the data values.
        :type  lonLat: str or :class:`numpy.ndarray`
        :param parameterList: Parameter values. If a string is given,
                              then it is the path to a file containing
                              the values. If an array is passed, then it
                              should hold the parameter values.
        :type  parameterList: str or :class:`numpy.ndarray`
        :param str parameterName: Optional. If given, then the
                                  cell distributions will be saved to a
                                  file with this name. If absent,
                                  the distribution values are returned.
        :param kdeStep: Increment of the ordinate values at which
                        the distributions will be calculated.
        :type  kdeStep: float, default=`0.1`
        :param angular: Does the data represent an angular measure
                        (e.g. bearing).
        :type  angular: boolean, default=``False``
        :param periodic: Does the data represent some form of periodic
                         data (e.g. day of year). If given, it should
                         be the period of the data (e.g. for annual data,
                         ``periodic=365``).
        :type  periodic: boolean or float, default=``False``
        :param boolean plotParam: Plot the parameters. Default is ``False``.

        :returns: If no ``parameterName`` is given returns ``None``
                  (data are saved to file), otherwise
                  :class:`numpy.ndarray`.


        """
        if parameterName:
            self.logger.debug("Running allDistributions for %s",
                              parameterName)
        else:
            self.logger.debug("Running allDistributions")

        if isinstance(lonLat, str):
            self.logger.debug("Loading lat/lon data from file")
            self.lonLat = np.array(flLoadFile(lonLat, delimiter=','))
        else:
            self.lonLat = lonLat

        if isinstance(parameterList, str):
            self.logger.debug("Loading parameter data from file: %s",
                              parameterList)
            self.pList = np.array(flLoadFile(parameterList))
        else:
            self.pList = parameterList

        self.pName = parameterName

        if len(self.pList) != len(self.lonLat):
            errmsg = ("Parameter data and "
                      "Lon/Lat data are not the same length "
                      "for {}.".format(parameterName))
            self.logger.critical(errmsg)
            raise IndexError(errmsg)

        maxCellNum = stats.maxCellNum(self.gridLimit, self.gridSpace)

        # Writing CDF dataset for all individual cell number into files
        self.logger.debug(("Writing CDF dataset for all individual "
                           "cells into files"))

        for cellNum in range(0, maxCellNum + 1):
            self.logger.debug("Processing cell number %i", cellNum)

            # Generate cyclone parameter data for the cell number
            self.extractParameter(cellNum)

            # Estimate cyclone parameter data using KDE
            # The returned array contains the grid, the PDF and the CDF
            cdf = self.kdeParameter.generateKDE(self.parameter, kdeStep,
                                                angular=angular,
                                                periodic=periodic)
            if plotParam:
                self._plotParameter(cellNum, kdeStep)
            self.logger.debug(('size of parameter array = %d: '
                               'size of cdf array = %d'), 
                              self.parameter.size, cdf.size)

            cellNumlist = []
            for i in range(len(cdf)):
                cellNumlist.append(cellNum)
            if cellNum == 0:
                results = np.transpose(np.array([cellNumlist,
                                                 cdf[:, 0], cdf[:, 2]]))
            else:
                self.logger.debug('size of results = %s', str(results.size))
                results = np.concatenate((results,
                                          np.transpose(np.array([cellNumlist,
                                                                 cdf[:, 0],
                                                                 cdf[:, 2]]))))

        if parameterName == None:
            self.logger.debug(("Returning CDF dataset for all "
                               "individual cell numbers"))
            return results
        else:
            cdfHeader = "Cell_Number, CDF_" + self.pName + "_x, CDF_" + \
                        self.pName + "_y"
            allCellCdfOutput = pjoin(self.outputPath, 'process',
                                     'all_cell_cdf_' + self.pName)

            args = {"filename":allCellCdfOutput, "data":results,
                    "header":cdfHeader, "delimiter":",", "fmt":"%f"}

            self.logger.debug(("Writing CDF dataset for all individual "
                               "cell numbers into files"))
            flSaveFile(**args)

            # Save to netcdf too

            filename = allCellCdfOutput + '.nc'

            ncdf = Dataset(filename, 'w')

            ncdf.createDimension('cell', len(results[:, 0]))
            cell = ncdf.createVariable('cell', 'i', ('cell',))
            cell[:] = results[:, 0]

            x = ncdf.createVariable('x', 'f', ('cell',))
            x[:] = results[:, 1]

            y = ncdf.createVariable('CDF', 'f', ('cell',))
            y[:] = results[:, 2]

            ncdf.close()

    def extractParameter(self, cellNum):
        """
        Extracts the cyclone parameter data for the given cell.
        If the population of a cell is insufficient for generating a
        PDF, the bounds of the cell are expanded until the population is
        sufficient.

        Null/missing values are removed.

        :param int cellNum: The cell number to process.
        :returns: None. The :attr:`parameter` attribute is updated.
        :raises IndexError: if the cell number is not valid
                            (i.e. if it is outside the possible
                            range of cell numbers).
        """
        if not stats.validCellNum(cellNum, self.gridLimit, self.gridSpace):
            self.logger.critical(("Invalid input on cellNum: "
                                  "cell number %i is out of range")%cellNum)
            raise IndexError('Invalid input on cellNum: '
                               'cell number is out of range')
        lon = self.lonLat[:, 0]
        lat = self.lonLat[:, 1]
        cellLon, cellLat = stats.getCellLonLat(cellNum, self.gridLimit,
                                               self.gridSpace)

        wLon = cellLon
        eLon = cellLon + self.gridSpace['x']
        nLat = cellLat
        sLat = cellLat - self.gridSpace['y']

        indij = np.where(((lat >= sLat) & (lat < nLat)) &
                         ((lon >= wLon) & (lon < eLon)))
        parameter_ = self.pList[indij]
        self.parameter = stats.statRemoveNum(np.array(parameter_),
                                             self.missingValue)

        while np.size(self.parameter) <= self.minSamplesCell:
            self.logger.debug(("Insufficient samples. Increasing the "
                               "size of the cell"))
            wLon_last = wLon
            eLon_last = eLon
            nLat_last = nLat
            sLat_last = sLat
            wLon, eLon, nLat, sLat = self._expandCell(lon, lat, wLon, eLon,
                                                      nLat, sLat)
            if ((wLon == wLon_last) & (eLon == eLon_last) &
                    (nLat == nLat_last) & (sLat == sLat_last)):
                errMsg = ("Insufficient grid points in selected domain to "
                          "estimate storm statistics - please select a larger "
                          "domain. Samples = %i / %i")%(np.size(self.parameter),
                                                        self.minSamplesCell)
                self.logger.critical(errMsg)
                raise StopIteration(errMsg)
            indij = np.where(((lat >= sLat) & (lat < nLat)) &
                             ((lon >= wLon) & (lon < eLon)))
            parameter_ = self.pList[indij]
            self.parameter = stats.statRemoveNum(np.array(parameter_),
                                                 self.missingValue)

        # Check to see if all values in the array are the same. If the
        # values are the same, bandwidth would be 0, and therefore KDE
        # cannot proceed
        while self.parameter.max() == self.parameter.min():
            self.logger.debug(("Parameter values appear to be the same. "
                               "Increasing the size of the cell"))
            wLon_last = wLon
            eLon_last = eLon
            nLat_last = nLat
            sLat_last = sLat
            wLon, eLon, nLat, sLat = self._expandCell(lon, lat, wLon,
                                                      eLon, nLat, sLat)
            if ((wLon == wLon_last) & (eLon == eLon_last) &
                    (nLat == nLat_last) & (sLat == sLat_last)):
                errMsg = ("Insufficient grid points in selected domain "
                          "to estimate storm statistics - "
                          "please select a larger domain.")
                self.logger.critical(errMsg)
                raise StopIteration(errMsg)
            indij = np.where(((lat >= sLat) & (lat < nLat)) &
                             ((lon >= wLon) & (lon < eLon)))
            parameter_ = self.pList[indij]
            self.parameter = stats.statRemoveNum(np.array(parameter_),
                                                 self.missingValue)
        self.logger.debug("Number of valid observations in cell %s : %s",
                          str(cellNum), str(np.size(self.parameter)))


    def _plotParameter(self, cellNum, kdeStep):
        import pylab
        self.logger.debug("Plotting %s"%self.pName)
        pMin = self.parameter.min()
        pMax = self.parameter.max()
        rng = np.arange(pMin, pMax, kdeStep)
        if len(rng) < 10:
            rng = 10
        pylab.clf()
        pylab.hist(self.parameter, rng)
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
        configfile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = ("No configuration file specified. "
                         "please type: python main.py {config filename}.ini")
            raise IOError(error_msg)
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError(error_msg)

    flStartLog(cnfGetIniValue(configFile, 'Logging',
                              'LogFile', __file__.rstrip('.py') + '.log'),
               cnfGetIniValue(configFile, 'Logging', 'LogLevel', 'DEBUG'),
               cnfGetIniValue(configFile, 'Logging', 'Verbose', True))
    path = cnfGetIniValue(configFile, 'Output', 'Path')

    gridLim = eval(cnfGetIniValue(configFile, 'Region', 'gridLimit'))
    gridSp = eval(cnfGetIniValue(configFile, 'Region', 'gridSpace'))
    gridinc = eval(cnfGetIniValue(configFile, 'Region', 'gridInc'))
    kdetype = cnfGetIniValue(configFile, 'StatInterface', 'kdeType')
    kdestep = cnfGetIniValue(configFile, 'StatInterface', 'kdeStep', 0.1)
    minSamples = cnfGetIniValue(configFile, 'StatInterface',
                                'minSamplesCell', 100)
    mv = cnfGetIniValue(configFile, 'StatInterface',
                        'MissingValue', sys.maxsize)
    gDist = GenerateDistributions(configFile, gridLim, gridSp,
                                  gridinc, kdetype,
                                  minSamples, mv)


    gDist.allDistributions(pjoin(path, 'init_lon_lat'),
                           pjoin(path, 'init_bearing'),
                           'init_bearing', 5.)
    gDist.allDistributions(pjoin(path, 'origin_lon_lat'),
                           pjoin(path, 'init_pressure'),
                           'init_pressure', kdestep)
    gDist.allDistributions(pjoin(path, 'init_lon_lat'),
                           pjoin(path, 'init_speed'),
                           'init_speed', kdestep)

    gDist.allDistributions(pjoin(path, 'origin_lon_lat'),
                           pjoin(path, 'init_rmax'),
                           'init_rmax', kdestep)
    gDist.allDistributions(pjoin(path, 'all_lon_lat'),
                           pjoin(path, 'all_bearing'),
                           'bearing', 5.)
    gDist.allDistributions(pjoin(path, 'all_lon_lat'),
                           pjoin(path, 'all_pressure'),
                           'pressure', kdestep)
    gDist.allDistributions(pjoin(path, 'all_lon_lat'),
                           pjoin(path, 'all_speed'),
                           'speed', kdestep)

    gDist.allDistributions(pjoin(path, 'all_lon_lat'),
                           pjoin(path, 'all_rmax'),
                           'rmax', kdestep)
    logging.shutdown()
