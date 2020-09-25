"""
:mod:`KDEParameters` -- generate KDE of cyclone parameters
==========================================================

.. module:: KDEParameters
    :synopsis: Generate KDEs of TC parameters.

.. moduleauthor:: Geoff Xu <geoff.xu@ga.gov.au>

Generates the probability density functions (using kernel density
estimation) of given cyclone parameters (speed, pressure, bearing, etc).
Each of these PDF's is converted to a cumulative density function for
use in other sections.

.. note:: In changing from the previous KPDF module to statsmodels, the 
          bandwidth calculation gives substantially different values for
          univariate data. For test data, the updated functions give a 
          smaller bandwidth value compared to KPDF.

"""

import os
import sys
import logging

import numpy as np
import Utilities.stats as stats

from statsmodels.nonparametric.kde import kernel_switch
import statsmodels.nonparametric.bandwidths as smbw
import statsmodels.api as sm

from Utilities.files import flLoadFile, flSaveFile
from Utilities.config import cnfGetIniValue

LOG = logging.getLogger()

class KDEParameters(object):
    """
    Generates the probability density functions (using kernel density
    estimation) of given cyclone parameters (speed, pressure, bearing,
    etc).  Each of these PDF's is converted to a cumulative density
    function for use in other sections.

    :param str kdeType: Name of the (univariate) kernel estimator to use
                        when generating the distribution. Must be one of
                        ``Epanechnikov``, ``Gaussian``, ``Biweight`` or
                        ``Triangular``.

    """


    def __init__(self, kdeType):
        """
        Initialize the logger and ensure the requested KDE type exists.
        
        This uses the `statsmodels.nonparametric.kernel_density` library
        
        """
        LOG.info("Initialising KDEParameters")
        kernels = kernel_switch.keys()
        if kdeType in kernels:
            LOG.debug(f"Using {kdeType} to generate distribution")
            self.kdeType = kdeType
        else:
            msg = (f"Invalid kernel type: {kdeType} \n"
                   f"Valid kernels are {repr(kernels)}")
            LOG.error(msg)
            raise NotImplementedError(msg)

    def generateKDE(self, parameters, kdeStep, kdeParameters=None,
                    cdfParameters=None, angular=False, periodic=False,
                    missingValue=sys.maxsize):
        """
        Generate a PDF and CDF for a given parameter set using the
        method of kernel density estimators. Optionally return the PDF
        and CDF as an array, or write both to separate files.

        :param parameters: Parameter values. If a string is given,
                           then it is the path to a file containing
                           the values. If an array is passed, then it
                           should hold the parameter values.

        :param kdeStep: Increment of the ordinate values at which
                        the distributions will be calculated.
        :type  kdeStep: float, default=`0.1`
        :param str kdeParameters: Optional. If given, then the
                                  cell distributions will be saved to a
                                  file with this name. If absent,
                                  the distribution values are returned.
        :param str cdfParameters: Optional. If given, then the
                                  cell distributions will be saved to a
                                  file with this name. If absent,
                                  the distribution values are returned.
        :param angular: Does the data represent an angular measure
                        (e.g. bearing).
        :type  angular: boolean, default=``False``
        :param periodic: Does the data represent some form of periodic
                         data (e.g. day of year). If given, it should
                         be the period of the data (e.g. for annual data,
                         ``periodic=365``).
        :type  periodic: boolean or int, default=``False``
        :param missingValue: Missing values have this value (default
                         :attr:`sys.maxint`).

        returns: If ``kdeParameters`` is given, returns ``None``
                  (data are saved to file), otherwise
                  :class:`numpy.ndarray` of the parameter grid, the PDF and CDF.

        """

        LOG.debug("Running generateKDE")
        if type(parameters) is str:
            self.parameters = stats.statRemoveNum(flLoadFile(parameters,
                                                             '%', ','),
                                                  missingValue)
        else:
            if parameters.size <= 1:
                LOG.error("Insufficient members in parameter list")
                raise IndexError("Insufficient members in parameter list")

            self.parameters = stats.statRemoveNum(parameters, missingValue)

        if angular:
            xmin = 0.0
            xmax = 360.0
        elif periodic:
            xmin = 0.0
            xmax = periodic
        else:
            xmin = self.parameters.min()
            xmax = self.parameters.max()

        LOG.debug("xmin=%7.3f, xmax=%7.3f, kdeStep=%7.3f" %
                  (xmin, xmax, kdeStep))
        if periodic:
            x = np.arange(1, periodic + 1, kdeStep)
            self.grid = np.concatenate([x - periodic, x, x + periodic])
            self.parameters = np.concatenate([self.parameters - periodic,
                                              self.parameters,
                                              self.parameters + periodic])
        else:
            self.grid = np.arange(xmin, xmax, kdeStep)

        if self.grid.size < 2:
            LOG.critical("Grid for CDF generation is a single value")
            LOG.critical("xmin=%7.3f, xmax=%7.3f, kdeStep=%7.3f", 
                         xmin, xmax, kdeStep)
            raise ValueError

        #bw = KPDF.UPDFOptimumBandwidth(self.parameters)
        bw = stats.bandwidth(self.parameters)
        self.pdf = self._generatePDF(self.grid, bw, self.parameters)

        if periodic:
            idx = int(periodic/kdeStep)
            self.pdf = 3.0*self.pdf[idx:2*idx]
            self.grid = self.grid[idx:2*idx]

        self.cy = stats.cdf(self.grid, self.pdf)
        if kdeParameters is None:
            return np.transpose(np.array([self.grid, self.pdf, self.cy]))
        else:
            # Assume both kdeParameters and cdfParameters are defined as files:
            LOG.debug("Saving KDE and CDF data to files")
            flSaveFile(kdeParameters,
                       np.transpose(np.array([self.grid, self.pdf])))
            flSaveFile(cdfParameters,
                       np.transpose(np.array([self.grid, self.cy])))

    def generateGenesisDateCDF(self, genDays, lonLat, bw=None, genesisKDE=None):
        """
        Calculate the PDF of genesis day using KDEs.
        Since the data is periodic, we use a simple method to include
        the periodicity in estimating the PDF. We prepend and append
        the data to itself, then use the central third of the PDF and
        multiply by three to obtain the required PDF. Probably not
        quite exact, but it should be sufficient for our purposes.

        :param str genDays: Name of file containing genesis days
                            (as day of year).
        :param lonLat: Array of genesis longitudes and latitudes.
        :param float bw: Optional. Bandwidth of the KDE to use.
        :param str genesisKDE: Optional. File name to save resulting CDF to.
        :type  lonLat: :class:`numpy.ndarray`

        :returns: :class:`numpy.ndarray` containing the days, the PDF and CDF
                  of the genesis days.
        """

        data = flLoadFile(genDays)
        days = np.arange(1, 366)
        ndays = np.concatenate([days - 365, days, days + 365])
        ndata = np.concatenate([data - 365, data, data + 365])

        if bw is None:
            bw = stats.bandwidth(self.parameters)

        kde = sm.nonparametric.KDEUnivariate(self.parameters)
        kde.fit(kernel=self.kdeType, bw=bw, fft=False, 
                gridsize=len(grid), clip=(min(grid), max(grid)), cut=0)
        #try:
        #    kdeMethod = getattr(KPDF, "UPDF%s" % self.kdeType)
        #except AttributeError:
        #    LOG.exception(("Invalid input on option: "
        #                   "KDE method UPDF%s does not exist"),
        #                  self.kdeType)
        #    raise
            
        veceval = np.vectorize(kde.evaluate)
        pdf = np.nan_to_num(veceval(grid))
        
        # Actual PDF to return
        apdf = 3.0*pdf[365:730]
        cy = stats.cdf(days, apdf)
        if genesisKDE is None:
            return np.transpose(np.array(np.concatenate([days, apdf, cy])))
        else:
            # Assume both kdeParameters and cdfParameters are defined as files:
            LOG.debug("Saving KDE and CDF data to files")
            flSaveFile(genesisKDE, np.transpose(np.array([days, cy])))


    def _generatePDF(self, grid, bw, dataset):
        """
        Sub-function that generates the PDFs of kernel
        density estimation from raw dataset
        """
        LOG.debug("Generating PDF")
        if bw <= 0:
            LOG.critical("bw = %d. Bandwidth cannot be negative or zero", bw)
            raise ValueError('bw = %d. Bandwidth cannot be negative or zero'%bw)

        kde = sm.nonparametric.KDEUnivariate(dataset)

        kde.fit(kernel=self.kdeType, fft=False, 
                gridsize=len(grid), clip=(min(grid), max(grid)), cut=0)
        veceval = np.vectorize(kde.evaluate)
        pdf = veceval(grid)
        return np.nan_to_num(pdf)
        
        #try:
        #    
        #    #kdeMethod = getattr(KPDF, "UPDF%s" %self.kdeType)
        #except AttributeError:
        #    LOG.exception(("Invalid input on option: "
        #                   "KDE method UPDF%s does not exist"),
        #                  self.kdeType)
        #    raise

        #return kdeMethod(dataset, grid, bw)

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = ("No configuration file specified, please type: "
                         "python main.py {config filename}.ini")
            raise IOError(error_msg)
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError(error_msg)

    logging.basicConfig(
        level=logging.getattr(cnfGetIniValue(configFile, 'Logging',
                                             'LogLevel', 'DEBUG')),
        format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s',
        filename=cnfGetIniValue(configFile, 'Logging', 'LogFile',
                                __file__.rstrip('.py') + '.log'),
        filemode='w')

    path = cnfGetIniValue(configFile, 'Output', 'Path')

#    init_bearing = path+'init_bearing'
#    kde_init_bearing = path+'kde_init_bearing'
#    cdf_init_bearing = cpath+'cdf_init_bearing'
#    all_bearing = path+'all_bearing'
#    kde_all_bearing = path+'kde_all_bearing'
#    cdf_all_bearing = path+'cdf_all_bearing'
#    bearing_no_init = path+'bearing_no_init'
#    kde_no_init_bearing = path+'kde_no_init_bearing'
#    cdf_no_init_bearing = path+'cdf_no_init_bearing'
    pressure_rate = os.path.join(path, 'pressure_rate')
#    bearing_rate = path+'bearing_rate'
#    speed_rate = path+'speed_rate'
    kde_pressure_rate = os.path.join(path, 'kde_pressure_rate')
    cdf_pressure_rate = os.path.join(path, 'cdf_pressure_rate')
#    kde_bearing_rate = path+'kde_bearing_rate'
#    cdf_bearing_rate = path+'cdf_bearing_rate'
#    kde_speed_rate = path+'kde_speed_rate'
#    cdf_speed_rate = path+'cdf_speed_rate'

    kdeStep = 0.1

    k = KDEParameters(cnfGetIniValue(configFile, 'Parameters',
                                     'kdeType', 'Biweight'), kdeStep)

    k.generateKDE(pressure_rate, kdeStep, kde_pressure_rate, cdf_pressure_rate)
    #k.generateKDE(bearing_rate, kde_bearing_rate, cdf_bearing_rate)
    #k.generateKDE(speed_rate, kde_speed_rate, cdf_speed_rate)
    #k.plotKdeInit()
    #k.plotKdeAll()
    #k.plotKdeNoInit()
    #k.plotCdf()
    #pylab.show()
