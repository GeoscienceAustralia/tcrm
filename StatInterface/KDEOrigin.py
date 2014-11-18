"""
:mod:`KDEOrigin` -- kernel density estimation for genesis probability
=====================================================================

.. module:: KDEOrigin
    :synopsis: Kernel density estimation for genesis probability.

.. moduleauthor:: Geoff Xu <geoff.xu@ga.gov.au>

Calculate a genesis probability distribution, based on the observed
genesis locations and applying a 2-d kernel density estimation method.

"""

import os
import logging
import numpy as np

import Utilities.stats as stats
import Utilities.KPDF as KPDF

from Utilities.files import flLoadFile, flStartLog
from Utilities.grid import grdSave
from Utilities.nctools import ncSaveGrid
from Utilities.config import ConfigParser

from statsmodels.nonparametric.kernel_density import KDEMultivariate

def getOriginBandwidth(data):
    """
    Calculate the optimal bandwidth for kernel density estimation
    from data.

    :param data: :class:`numpy.ndarray` of data points for training data

    :returns: Bandwidth parameter.
    """
    dens = KDEMultivariate(data=data, var_type='cc', bw='cv_ml')
    return dens.bw

class KDEOrigin:
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

    def __init__(self, configFile, gridLimit, kdeStep, lonLat=None, progressbar=None):
        """
        
        """
        self.logger = logging.getLogger()
        self.progressbar = progressbar
        if self.progressbar:
            KPDF.set_callback(self.updateProgressBar)
        self.logger.info("Initialising KDEOrigins")
        self.configFile = configFile
        self.x = np.arange(gridLimit['xMin'], gridLimit['xMax'], kdeStep)
        self.y = np.arange(gridLimit['yMax'], gridLimit['yMin'], -kdeStep)

        self.kdeStep = kdeStep

        config = ConfigParser()
        config.read(configFile)

        if lonLat is None:
            self.outputPath = config.get('Output', 'Path')
            self.processPath = os.path.join(self.outputPath, 'process')
            self.logger.debug("Loading "+os.path.join(self.processPath,
                                                  'init_lon_lat'))
            ll = flLoadFile(os.path.join(self.processPath, 'init_lon_lat'),
                            '%', ',')
            self.lonLat = ll[:,0:2]
        else:
            self.lonLat = lonLat[:,0:2]

        ii = np.where((self.lonLat[:,0] >= gridLimit['xMin']) &
                      (self.lonLat[:,0] <= gridLimit['xMax']) &
                      (self.lonLat[:,1] >= gridLimit['yMin']) &
                      (self.lonLat[:,1] <= gridLimit['yMax']))
        
        self.lonLat = self.lonLat[ii]
        
        self.bw = getOriginBandwidth(self.lonLat)
        self.logger.debug("Bandwidth: {0}".format(self.bw))


    def generateKDE(self, bw=None, save=False, plot=False):
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

        kde = KDEMultivariate(self.lonLat, bw=self.bw, var_type='cc')
        xx, yy = np.meshgrid(self.x, self.y)
        xy = np.vstack([xx.ravel(), yy.ravel()])
        pdf = kde.pdf(data_predict=xy)
        pdf = pdf.reshape(xx.shape)

        self.pdf = pdf.transpose()

        if save:
            outputFile = os.path.join(self.processPath, 'originPDF.nc')
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

            ncSaveGrid(outputFile, dimensions, variables)

        if plot:

            from PlotInterface.maps import FilledContourMapFigure, saveFigure
            # Automatically determine appropriate contour levels
            min_lvls = 6.0
            lvls_options = np.array([1.0, 0.5, 0.25, 0.2, 0.1])
            pdfMax = pdf.max()
            exponent = int(np.floor(np.log10(pdfMax)))
            significand = pdfMax * 10**-exponent
            lvl_step = lvls_options[np.where((significand/lvls_options) > min_lvls)[0][0]]
            lvls = np.arange(0, pdf.max(), lvl_step*(10.0**exponent))
            [gx,gy] = np.meshgrid(self.x,self.y)
            map_kwargs = dict(llcrnrlon=self.x.min(),
                              llcrnrlat=self.y.min(),
                              urcrnrlon=self.x.max(),
                              urcrnrlat=self.y.max(),
                              projection='merc',
                              resolution='i')
            
            cbarlabel = r'Genesis probability ($\times 10^{' + str(exponent) + '}$)'
            title = 'TC Genesis probability'
            figure = FilledContourMapFigure()
            figure.add(pdf, gx, gy, title, lvls, cbarlabel, map_kwargs)
            figure.plot()

            outputFile = os.path.join(self.outputPath, 'plots', 'stats', 'originPDF_fill.png')
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
        self.cz = stats.cdf2d(self.x, self.y, self.pdf)
        if save:
            self.logger.debug("Saving origin CDF to file")
            grdSave(self.processPath+'originCDF.txt', self.cz, self.x,
                    self.y, self.kdeStep)

        if save:
            outputFile = os.path.join(self.processPath, 'originCDF.nc')
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

            variables =  {
                0: {
                    'name': 'gcdf',
                    'dims': ('lat','lon'),
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

    def updateProgressBar(self, n, nMax):
        """
        Callback function to update progress bar from C code

        :param int n: Current step.
        :param int nMax: Maximum step. 
        
        """
        if self.progressbar:
            self.progressbar.update(n/float(nMax), 0.0, 0.7)

