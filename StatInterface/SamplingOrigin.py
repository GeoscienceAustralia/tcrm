"""
:mod:`SamplingOrigin` -- Generate random TC origins
===================================================

.. module:: SamplingOrigin
    :synopsis: Generate a sample of random TC origins from
               CDFs derived using kernel density estimation.

.. moduleauthor:: Geoff Xu <geoff.xu@ga.gov.au>

Define the class for sampling tropical cyclone origins. 

"""

import os, sys, pdb, logging

import time
from Utilities.files import flLoadFile, flSaveFile
from Utilities.grid import grdRead, grdReadFromNetcdf
import numpy as np
import scipy
import Utilities.stats as stats
from Utilities.config import ConfigParser

class SamplingOrigin:
    """
    Class for generating samples of TC origins.

    :param kdeOrigin: Name of a file containing TC genesis PDF
                      data, or a 2-d array containing the PDF. 
    :type  kdeOrigin: str or :class:`numpy.ndarray`
    :param x: Longitude coordinates of the grid on which the PDF
              is defined.
    :param y: Latitude coordinates of the grid on which the PDF
              is defined.
    :type  x: :class:`numpy.ndarray`
    :type  y: :class:`numpy.ndarray`
    

    """

    def __init__(self, kdeOrigin=None, x=None, y=None):
        """
        Initialise the array of probabilities of genesis, plus the
        lon/lat arrays.
        """
        self.logger=logging.getLogger()

        if type(kdeOrigin) == str:
            self.logger.debug("Loading PDF from %s"%kdeOrigin)
            try:
                if kdeOrigin.endswith('nc'):
                    self.x, self.y, self.z = grdReadFromNetcdf(kdeOrigin)
                else:
                    self.x, self.y, self.z = grdRead(kdeOrigin)
            except IOError:
                self.logger.critical('Error! Files relating to cdf of cyclone parameters does not exist, please generate KDE of cyclone parameters first.')
                raise
            self._calculateCDF()  # calculate CDF of (x,Px) and (y,Py)
        elif type(kdeOrigin) == np.ndarray:
            self.x = x
            self.y = y
            self.z = kdeOrigin
            self._calculateCDF()  # calculate CDF of (x,Px) and (y,Py)
        else:
            self.x = np.array([])
            self.y = np.array([])
            self.z = np.array([])

    def setKDEOrigins(self, kdeOriginX=None, kdeOriginY=None, kdeOriginZ=None,
                      outputPath=None):
        """
        Set kernel density estimation origin parameters.

        :param kdeOriginX: x coordinates of kde result generated
                           from :class:`KDEOrigin`
        :param kdeOriginY: y coordinates of kde result generated
                           from :class:`KDEOrigin`
        :param kdeOriginZ: z coordinates of kde result generated
                           from :class:`KDEOrigin`
        :param outputPath: Path to output folder to load PDF file.

        :type  kdeOriginX: str or :class:`numpy.ndarray`
        :type  kdeOriginY: str or :class:`numpy.ndarray`
        :type  kdeOriginZ: str or :class:`numpy.ndarray`
        :type  outputPath: str
        
        """
        if outputPath:
            try:
                self.x, self.y. self.z = grdRead(os.path.join(self.outputPath,
                                                              'originPDF.txt'))
            except IOError:
                self.logger.critical('Error! Files relating to KDE of cyclone origins does not exist. Execute KDE of cyclone origins first.')
                raise
            self._calculateCDF() #calculate CDF of (x,Px) and (y,Py)
        elif type(kdeOriginZ) == str:
            try:
                self.x = flLoadFile(kdeOriginX)
                self.y = flLoadFile(kdeOriginY)
                self.z = flLoadFile(kdeOriginZ)
            except IOError:
                self.logger.critical('Error! Files relating to CDF of cyclone parameters do not exist. Generate KDE of cyclone parameters first.')
                raise
            self._calculateCDF()  # calculate CDF of (x,Px) and (y,Py)
        elif type(kdeOriginZ) == np.ndarray:
            self.x, self.y, self.z = grdRead(kdeOriginZ)
            self._calculateCDF()  # calculate CDF of (x,Px) and (y,Py)
        else:
            self.logger.error("No input arguments")
            raise

    def generateOneSample(self):
        """Generate a random cyclone origin."""
        # generate 2 uniform random variables
        unifX = scipy.rand()
        unifY = scipy.rand()

        xi = np.array(self.cdfX).searchsorted(unifX)
        yj = self.cdfY[xi, :].searchsorted(unifY)

        return self.x[xi], self.y[yj] #lon, lat

    def ppf(self, q1, q2):
        """
        Percent point function on 2-d grid (inverse of CDF).

        :param float q1: Quantile for the x-coordinate.
        :param float q2: Quantile for the y-coordinate.

        :returns: Longitude & latitude of the given quantile values.

        """
        xi = self.cdfX.searchsorted(q1)
        yj = self.cdfY[xi, :].searchsorted(q2)
        return self.x[xi], self.y[yj] #lon, lat

    def cdf(self, x, y):
        """
        Return CDF value at the given location.

        :param float x: x-ccordinate.
        :param float y: y-coordinate.

        :returns: CDF values for x & y at the given location.
        
        """
        
        # crude, this should be an interpolation
        xi = self.x.searchsorted(x) - 1
        yi = self.y.searchsorted(y) - 1
        return self.cdfX[xi], self.cdfY[xi, yi]

    def generateSamples(self, ns, outputFile=None):
        """
        Generate random samples of cyclone origins.

        :param int ns: Number of samples to generate.
        :param str outputFile: If given, save the samples to the file.

        :returns: :class:`numpy.ndarray` containing longitude and
                  latitude of a random sample of TC origins.

        :raises ValueError: If :attr:`ns` <= 0.
        :raises IndexError: If an invalid index is returned when generating
                            uniform random values. 
                            
        """
        if ns <= 0:
            self.logger.error('Invalid input on ns: number of sample cannot be zero or negative')
            raise ValueError

        # Generate 2 vectors of uniform random variables
        unifX = scipy.rand(ns)
        unifY = scipy.rand(ns)

        self.oLon = np.empty(ns, 'd')
        self.oLat = np.empty(ns, 'd')

        # For each random variable
        try:
            for i in xrange(ns):
                xi = self.cdfX.searchsorted(unifX[i])
                yj = self.cdfY[xi, :].searchsorted(unifY[i])
                if (i % (ns/100)) == 0 and i != 0:
                    self.logger.debug("Processing %ith element"%i)
                self.oLon[i] = self.x[xi]
                self.oLat[i] = self.y[yj]

        except IndexError:
            self.logger.debug("i = %s"%str(i))
            self.logger.debug("unifX = %s"%str(unifX[i]))
            self.logger.debug("unifY = %s"%str(unifY[i]))
            self.logger.debug("cdfY[xi,:] = %s"%str(self.cdfY[xi, :]))
            raise

        if outputFile:
            flSaveFile(outputFile, np.transpose([self.oLon, self.oLat]),
                       fmt='%2.3f', header='Origin Lon, Origin Lat',
                       delimiter=',')
        else:
            lonLat = np.empty([ns, 2], 'd')
            lonLat[:, 0] = self.oLon
            lonLat[:, 1] = self.oLat
            return lonLat

    def _calculateCDF(self):
        """Calculate Py and CDFy beforehand to remove the need of
        repeated calculation later
        """
        # sum along the column of z to get sum(z(i,:))
        # (check 'help sum' if need)
        px = self.z.sum(axis=1)
        # calculate CDF of (x,Px)
        cdfX = stats.cdf(self.x, px)
        # define Py & CDFy with nx by ny
        py = np.zeros(self.z.shape, 'd')
        cdfY = np.zeros(self.z.shape, 'd')
        # Py=conditional distribution,  CDFy = CDF of Y
        try:
            for i in xrange(len(self.x)):
                for j in xrange(len(self.z[i, :])):
                    if px[i] == 0:
                        py[i,j] = 0
                    else:
                        py[i,j] = self.z[i, j]/px[i]
                cdfTemp = stats.cdf(self.y, py[i, :])
                for j in xrange(len(cdfTemp)):
                    cdfY[i,j] = cdfTemp[j]
        except IndexError:
            self.logger.debug("i = %s"%str(i))
            self.logger.debug("j = %s"%str(j))
            self.logger.debug("p_y[%s, %s] = %s"%(str(i), str(j), str(py[i, j])))
            self.logger.debug("z[%s, %s] = %s"%(str(i), str(j), str(self.z[i, j])))
            self.logger.debug("p_x[%s] = %s"%(str(i), str(px[i])))
            self.logger.debug("cdfy dim = %s"%(str(cdfY.shape)))
            self.logger.debug("p_y dim = %s"%(str(py.shape)))
            self.logger.debug("cdfx dim = %s"%(str(cdfX.shape)))
            self.logger.debug("p_x dim = %s"%(str(px.shape)))

            raise

        self.cdfX = cdfX
        self.cdfY = cdfY
        return
