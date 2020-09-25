"""
:mod:`generateStats` -- calculation of statistical values
=========================================================

.. module:: generateStats
    :synopsis: Calculates the mean, variance and
               autocorrelation statistics for all grid
               cells in the domain, for a given parameter.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>


"""

import logging
import sys

import numpy as np

import Utilities.stats as stats

from Utilities.files import flLoadFile

from scipy.stats import scoreatpercentile as percentile, circmean, circstd
from PlotInterface.curves import RangeCurve, saveFigure

def acf(p, nlags=1):
    """
    Autocorrelation coefficient

    :param p: array of values to calculate autocorrelation coefficient
    :type p: 1-d :class:`numpy.ndarray`

    """
    ar = np.array([1]+[np.corrcoef(p[:-i], p[i:])[0,1] for i in range(1, nlags + 1)])
    #ar = np.correlate(p, p, 'full')
    #n = len(p)
    ## Grab only the lag-one autocorrelation coeff.
    #ar = ar[n-1:(n+nlags)]/ar.max()
    return ar


class parameters(object):
    """
    Description: Create an object that holds :class:`numpy.ndarray`s
    of the statistical properties of each grid cell. There are
    :class:`numpy.ndarray` for both land and sea cells.

    mu/lmu : array of mean values of parameter for each grid cell
    sig/lsig : array of variance values of parameter for each grid cell
    alpha/lalpha : array of autoregression coefficients of parameter for
    each grid cell
    phi/lphi : array of normalisation values for random variations

    """

    def __init__(self, numCells):
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
    """
    Generate the main statistical distributions across the grid domain.

    :type  parameter: :class:`numpy.ndarray` or str
    :param parameter: contains the data on which the statistical
                      values will be based. If `str`, then
                      represents the name of a file that contains
                      the data
    :type  lonLat: :class:`numpy.ndarray` or str
    :param lonLat: Contains the longitude and latitude of each
                   of the observations in the :class:`numpy.ndarray`
                   parameter
    :param dict gridLimit: dictionary containing limits of regional grid
               Contains keys 'xMin', 'xMax', 'yMin', 'yMax'
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

    :param int minSample: Minimum number of valid observations
                          required to generate distributions. If
                          insufficient observations are found in a
                          grid cell, then it is incrementally expanded
                          until ``minSample`` is reached.
    :param boolean angular: If ``True`` the data represents an angular
                            variable (e.g. bearings). Default is ``False``.

    """

    def __init__(self, parameter, lonLat, gridLimit,
                 gridSpace, gridInc, minSample=100, angular=False,
                 missingValue=sys.maxsize, progressbar=None,
                 prgStartValue=0, prgEndValue=1, calculateLater=False):

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
        """
        Cycle through the cells and calculate the statistics
        for the variable.

        """
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

    def plotStatistics(self, output_file):

        p = stats.statRemoveNum(np.array(self.param), self.missingValue)
        a = p - np.mean(p)
        pmin = p.min()
        pmax = p.max()
        amin = a.min()
        amax = a.max()
        abins = np.linspace(amin, amax, 50)
        bins = np.linspace(pmin, pmax, 50)
        hist = np.empty((len(bins) - 1, self.maxCell))
        ahist = np.empty((len(abins) - 1, self.maxCell))
        x = np.arange(11)
        alpha = np.empty((11, self.maxCell))
        aalpha = np.empty((11, self.maxCell))

        for i in range(self.maxCell + 1):
            p = self.extractParameter(i, 0)
            a = p - np.mean(p)
            hist[:, i - 1], b = np.histogram(p, bins, density=True)
            ahist[:, i - 1], b = np.histogram(a, abins, density=True)
            alpha[:, i - 1] = acf(p, 10)
            aalpha[:, i - 1] = acf(a, 10)

        mhist = np.mean(hist, axis=1)
        uhist = percentile(hist, per=95, axis=1)
        lhist = percentile(hist, per=5, axis=1)

        mahist = np.mean(ahist, axis=1)
        uahist = percentile(ahist, per=95, axis=1)
        lahist = percentile(ahist, per=5, axis=1)

        malpha = np.mean(alpha, axis=1)
        ualpha = percentile(alpha, per=95, axis=1)
        lalpha = percentile(alpha, per=5, axis=1)

        maalpha = np.mean(aalpha, axis=1)
        uaalpha = percentile(aalpha, per=95, axis=1)
        laalpha = percentile(aalpha, per=5, axis=1)

        fig = RangeCurve()
        fig.add(bins[:-1], mhist, uhist, lhist, "Values", "Probability", "")
        fig.add(abins[:-1], mahist, uahist, lahist, "Anomalies", "Probability", "")
        fig.add(x, malpha, ualpha, lalpha, "Lag", "Autocorrelation", "ACF of values")
        fig.add(x, maalpha, uaalpha, laalpha, "Lag", "Autocorrelation", "ACF of anomalies")
        fig.plot()

        saveFigure(fig, output_file + '.png')

    def calculate(self, cellNum, onLand):
        """
        Calculate the required statistics (mean, variance,
        autocorrelation and regularized anomaly coefficient) for the
        given cell.

        :param int cellNum: The cell number to process.
        :param boolean onLand: If ``True``, then the cell is (mostly
                               or entirely) over land. If ``False``,
                               the cell is over water.

        :returns: mean, standard deviation, autocorrelation, residual
                  correlation and the minimum parameter value.
        """

        p = self.extractParameter(cellNum, onLand)

        if self.angular:
            mu = circmean(np.radians(p))
            sig = circstd(np.radians(p))
        else:
            mu = np.mean(p)
            sig = np.std(p)

        # Calculate the autocorrelations:
        alphas = np.correlate(p, p, 'full')
        n = len(p)

        # Grab only the lag-one autocorrelation coeff.
        alpha = alphas[n]/alphas.max()
        alpha = acf(p)[-1]
        phi = np.sqrt(1 - alpha**2)
        mn = min(p)

        return mu, sig, alpha, phi, mn

    def extractParameter(self, cellNum, onLand):
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
            self.logger.critical("Invalid input on cellNum: cell number %i is out of range"%cellNum)
            raise IndexError('Invalid input on cellNum: cell number %i is out of range'%cellNum)
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
                    errMsg = ("Insufficient grid points in selected "
                              "domain to estimate storm statistics - "
                              "please select a larger domain.")
                    self.logger.critical(errMsg)
                    raise StopIteration(errMsg)

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
                    raise StopIteration(errMsg)
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
        """
        Load pre-calculated statistics from a netcdf file.

        :param str filename: Path to the netcdf-format file containing
                             the statistics.

        """

        self.logger.debug('Loading statistics from %s' % filename)
        from netCDF4 import Dataset
        ncdf = Dataset(filename, 'r')
        for var in list(ncdf.variables.keys()):
            setattr(self.coeffs, var, ncdf.variables[var][:].flatten())
        #TODO: maybe save and check grid settings?

    def save(self, filename, description=''):
        """
        Save parameters to a netcdf file for later access.

        :param str filename: Path to the netcdf file to be created.
        :param str description: Name of the parameter.

        """

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
                        'dtype':'i8',
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
                           nodata=self.missingValue,
                           datatitle=None, writedata=True,
                           keepfileopen=False)

