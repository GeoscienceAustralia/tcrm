"""
:mod:`KDEOrigin` -- kernel density estimation for genesis probability
=====================================================================

.. module:: KDEOrigin
    :synopsis: Kernel density estimation for genesis probability.

.. moduleauthor:: Geoff Xu <geoff.xu@ga.gov.au>

Calculate a genesis probability distribution, based on the observed
genesis locations and applying a 2-d kernel density estimation method.

"""

from os.path import join as pjoin
import logging
import numpy as np

from Utilities.files import flLoadFile
from Utilities.nctools import ncSaveGrid
from Utilities.config import ConfigParser

from statsmodels.nonparametric.kernel_density import KDEMultivariate

LOGGER = logging.getLogger(__name__)

def getOriginBandwidth(data):
    """
    Calculate the optimal bandwidth for kernel density estimation
    from data.

    :param data: :class:`numpy.ndarray` of data points for training data

    :returns: Bandwidth parameter.
    """
    dens = KDEMultivariate(data=data, var_type='cc', bw='cv_ml')
    return dens.bw

class KDEOrigin(object):
    """
    Initialise the class for generating the genesis probability distribution.
    Initialisation will load the required data (genesis locations) and
    calculate the optimum bandwidth for the kernel density method.

    :param str configFile: Path to the configuration file.
     :param dict gridLimit: The bounds of the model domain. The
                           :class:`dict` should contain the keys
                           :attr:`xMin`, :attr:`xMax`, :attr:`yMin`
                           and :attr:`yMax`. The *x* variable bounds
                           the longitude and the *y* variable
                           bounds the latitude.
    :param float kdeStep: Increment of the ordinate values at which
                          the distributions will be calculated.
                          Default=`0.1`
    :param lonLat: If given, a 2-d array of the longitude and latitude
                   of genesis locations. If not given, attempt to load
                   an ``init_lon_lat`` file from the processed files.
    :param progressbar: A :meth:`SimpleProgressBar` object to print
                        progress to STDOUT.
    :type  lonLat: :class:`numpy.ndarray`
    :type  progressbar: :class:`Utilities.progressbar` object.


    """

    def __init__(self, configFile, gridLimit, kdeStep, lonLat=None,
                 progressbar=None):
        """

        """
        self.progressbar = progressbar
        LOGGER.info("Initialising KDEOrigin")
        self.x = np.arange(gridLimit['xMin'], gridLimit['xMax'], kdeStep)
        self.y = np.arange(gridLimit['yMax'], gridLimit['yMin'], -kdeStep)

        self.kdeStep = kdeStep
        self.kde = None
        self.pdf = None
        self.cz = None

        self.configFile = configFile
        self.config = ConfigParser()
        self.config.read(configFile)

        if lonLat is None:
            # Load the data from file:
            self.outputPath = self.config.get('Output', 'Path')
            self.processPath = pjoin(self.outputPath, 'process')
            LOGGER.debug("Loading " + pjoin(self.processPath, 'init_lon_lat'))
            ll = flLoadFile(pjoin(self.processPath, 'init_lon_lat'), '%', ',')
            self.lonLat = ll[:, 0:2]
        else:
            self.lonLat = lonLat[:, 0:2]

        ii = np.where((self.lonLat[:, 0] >= gridLimit['xMin']) &
                      (self.lonLat[:, 0] <= gridLimit['xMax']) &
                      (self.lonLat[:, 1] >= gridLimit['yMin']) &
                      (self.lonLat[:, 1] <= gridLimit['yMax']))

        self.lonLat = self.lonLat[ii]

        self.bw = getOriginBandwidth(self.lonLat)
        LOGGER.info("Bandwidth: %s", repr(self.bw))


    def generateKDE(self, save=False, plot=False):
        """
        Generate the PDF for cyclone origins using kernel density
        estimation technique then save it to a file path provided by
        user.

        :param float bw: Optional, bandwidth to use for generating the PDF.
                         If not specified, use the :attr:`bw` attribute.
        :param boolean save: If ``True``, save the resulting PDF to a
                             netCDF file called 'originPDF.nc'.
        :param boolean plot: If ``True``, plot the resulting PDF.

        :returns: ``x`` and ``y`` grid and the PDF values.

        """

        self.kde = KDEMultivariate(self.lonLat, bw=self.bw, var_type='cc')
        xx, yy = np.meshgrid(self.x, self.y)
        xy = np.vstack([xx.ravel(), yy.ravel()])
        pdf = self.kde.pdf(data_predict=xy)
        pdf = pdf.reshape(xx.shape)

        self.pdf = pdf.transpose()

        if save:
            dimensions = {
                0: {
                    'name': 'lat',
                    'values': self.y,
                    'dtype': 'f',
                    'atts': {
                        'long_name':' Latitude',
                        'units': 'degrees_north'
                    }
                },
                1: {
                    'name': 'lon',
                    'values': self.x,
                    'dtype': 'f',
                    'atts': {
                        'long_name': 'Longitude',
                        'units': 'degrees_east'
                    }
                }
            }

            variables = {
                0: {
                    'name': 'gpdf',
                    'dims': ('lat', 'lon'),
                    'values': np.array(pdf),
                    'dtype': 'f',
                    'atts': {
                        'long_name': 'TC Genesis probability distribution',
                        'units': ''
                        }
                    }
                }

            ncSaveGrid(pjoin(self.processPath, 'originPDF.nc'),
                       dimensions, variables)

        if plot:
            from PlotInterface.maps import FilledContourMapFigure, \
                saveFigure, levels

            lvls, exponent = levels(pdf.max())

            [gx, gy] = np.meshgrid(self.x, self.y)

            map_kwargs = dict(llcrnrlon=self.x.min(),
                              llcrnrlat=self.y.min(),
                              urcrnrlon=self.x.max(),
                              urcrnrlat=self.y.max(),
                              projection='merc',
                              resolution='i')

            cbarlabel = r'Genesis probability ($\times 10^{' + \
                        str(exponent) + '}$)'
            figure = FilledContourMapFigure()
            figure.add(pdf*(10**-exponent), gx, gy, 'TC Genesis probability',
                       lvls*(10**-exponent), cbarlabel, map_kwargs)
            figure.plot()

            outputFile = pjoin(self.outputPath, 'plots',
                               'stats', 'originPDF.png')
            saveFigure(figure, outputFile)

        return self.x, self.y, self.pdf

    def generateCdf(self, save=False):
        """
        Generate the CDFs corresponding to PDFs of cyclone origins,
        then save it on a file path provided by user

        :param boolean save: If ``True``, save the CDF to a netcdf file
                             called 'originCDF.nc'. If ``False``, return
                             the CDF.

        """
        xx, yy = np.meshgrid(self.x, self.y)
        xy = np.vstack([xx.ravel(), yy.ravel()])
        self.cz = self.kde.cdf(data_predict=xy)

        if save:
            outputFile = pjoin(self.processPath, 'originCDF.nc')
            dimensions = {
                0: {
                    'name': 'lat',
                    'values': self.y,
                    'dtype': 'f',
                    'atts': {
                        'long_name': 'Latitude',
                        'units': 'degrees_north'
                    }
                },
                1: {
                    'name': 'lon',
                    'values': self.x,
                    'dtype': 'f',
                    'atts': {
                        'long_name': 'Longitude',
                        'units':'degrees_east'
                    }
                }
            }

            variables = {
                0: {
                    'name': 'gcdf',
                    'dims': ('lat', 'lon'),
                    'values': np.array(self.cz),
                    'dtype': 'f',
                    'atts': {
                        'long_name': ('TC Genesis cumulative '
                                      'distribution'),
                        'units': ''
                        }
                    }
                }

            ncSaveGrid(outputFile, dimensions, variables)
        else:
            return self.cz

    def updateProgressBar(self, step, stepMax):
        """
        Callback function to update progress bar from C code

        :param int n: Current step.
        :param int nMax: Maximum step.

        """
        if self.progressbar:
            self.progressbar.update(step/float(stepMax), 0.0, 0.7)

