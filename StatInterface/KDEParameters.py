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


Title: KDEParameters.py - generate kernel density estimations of cyclone
parameters

Author: Geoff Xu, geoff.xu@ga.gov.au
CreationDate: 2006-01-05
Description:
Generates the probability density functions (using kernel density
estimation) of given cyclone parameters (speed, pressure, bearing, etc).
Each of these PDF's is converted to a cumulative density function for
use in other sections.

Version: $Rev: 810 $

ModifiedDate: 2006-11-29
ModifiedBy: N. Habili, nariman.habili@ga.gov.au
Modifications: Added accessor methods, upgraded to ndarray and conforms
               with style guide

ModifiedDate: 2006-12-07
ModifiedBy: C. Arthur, craig.arthur@ga.gov.au
Modifications: Added a generic generateKDE() method, tidied _generatePDF()

ModifiedDate: 2007-05-01
ModifiedBy: N. Habili, nariman.habili@ga.gov.au
Modifications: Removed plotting methods
               Removed generateKDEAll, generateKDEInit, and
               generateKDENoInit

$Id: KDEParameters.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

from scipy import array, arange, transpose
import numpy
import Utilities.stats as stats
import Utilities.KPDF as KPDF

from Utilities.files import flLoadFile, flSaveFile
from Utilities.config import cnfGetIniValue


class KDEParameters:
    """
    Parameters
    ----------
    kdeType : string
        kernel density estimation type

    Members
    -------
    parameters : 1D array of float
        parameters of cyclones
    kdeType : string
        kernel density estimation type
    kdeParameters : string
        file path for kde parameters
    cdfParameters : string
        file path for cdf parameters

    Methods
    -------
    generateKde(self, parameters, kdeParameters, cdfParameters)
        Generate PDF and cdf for a cyclone parameter using kernel
        density estimation technique. It is then either returned or
        saved in a file path provided by user

    Internal Methods
    ----------------
    generatePDF(option,grid,dataset) : 1D array of float
        a sub function that generates the PDFs of kernel density
        estimation from raw dataset
    """

    def __init__(self, kdeType):
        """
        Initialize the logger and ensure the requested KDE type exists.
        """
        self.logger = logging.getLogger()
        self.logger.info("Initialising KDEParameters")

        if hasattr(KPDF, "UPDF%s" %kdeType):
            self.logger.debug("Using %s to generate distribution"%kdeType)
            self.kdeType = kdeType
        else:
            self.logger.error("Invalid KDE type: %s" %kdeType)
            raise NotImplementedError, "Invalid KDE type: %s" %kdeType

    def __doc__(self):
        """
        Documentation on what this class does
        """
        return "Calculate distributions for the cyclone parameters using \
                kernel density estimation technique"

    def generateKDE(self, parameters, kdeStep, kdeParameters=None,
                    cdfParameters=None, angular=False,
                    missingValue=sys.maxint):
        """
        Generate a PDF and CDF for a given parameter set using the
        method of kernel density estimators.
        Optionally return the PDF and CDF as an array, or write both
        to separate files.
        """

        self.logger.debug("Running generateKDE")
        if type(parameters) is str:
            self.parameters = stats.statRemoveNum(flLoadFile(parameters, '%', ','), missingValue)
        else:
            if parameters.size <= 1:
                self.logger.error("Insufficient members in parameter list")
                raise IndexError, "Insufficient members in parameter list"

            self.parameters = stats.statRemoveNum(parameters, missingValue)

        if angular:
            xmin = 0.0
            xmax = 360.0
        else:
            xmin = self.parameters.min()
            xmax = self.parameters.max()
        self.logger.debug("xmin=%7.3f, xmax=%7.3f, kdeStep=%7.3f" %
                           (xmin, xmax,kdeStep))
        self.grid = arange(xmin, xmax, kdeStep)
        if self.grid.size<2:
            self.logger.critical("Grid for CDF generation is a single value")
            self.logger.critical("xmin=%7.3f, xmax=%7.3f, kdeStep=%7.3f" %
                                  (xmin, xmax,kdeStep))
            raise ValueError
        bw = KPDF.UPDFOptimumBandwidth(self.parameters)
        self.pdf = self._generatePDF(self.grid, bw, self.parameters)
        self.cy = stats.cdf(self.grid, self.pdf)
        if kdeParameters is None:
            return transpose(array([self.grid, self.pdf, self.cy]))
        else:
            # Assume both kdeParameters and cdfParameters are defined as files:
            self.logger.debug("Saving KDE and CDF data to files")
            flSaveFile(kdeParameters, transpose(array([self.grid, self.pdf])))
            flSaveFile(cdfParameters, transpose(array([self.grid, self.cy])))

    def generateGenesisDateCDF( self, genDays, bw=None, genesisKDE=None ):
        """
        Calculate the PDF of genesis day using KDEs.
        Since the data is periodic, we use a simple method to include the 
        periodicity in estimating the PDF. We prepend and append the data
        to itself, then use the central third of the PDF and multiply by three to
        obtain the required PDF. Probably notquite exact, but it should be
        sufficient for our purposes. 
        """

        data = flLoadFile( genDays )
        days = arange( 1, 366 )
        ndays = numpy.concatenate( [days - 365, days, days + 365] )
        ndata = numpy.concatenate( [data - 365, data, data + 365] )

        if bw is None:
            bw = KPDF.UPDFOptimumBandwidth( ndata ) 

        try:
            kdeMethod = getattr(KPDF, "UPDF%s" %self.kdeType)
        except AttributeError:
            self.logger.exception("Invalid input on option: KDE method UPDF%s does not exist"%self.kdeType)
            raise
        pdf = kdeMethod( ndata, ndays, bw )
        # Actual PDF to return
        apdf = 3.0*pdf[365:730]
        cy = stats.cdf(days, apdf)
        if genesisKDE is None:
            return transpose(array(numpy.concatenate( [days, apdf, cy] ) ))
        else:
            # Assume both kdeParameters and cdfParameters are defined as files:
            self.logger.debug("Saving KDE and CDF data to files")
            flSaveFile(genesisKDE, transpose(numpy.concatenate([days, pdf])))
            #flSaveFile(cdfParameters, transpose(array([days, cy])))
        

    def _generatePDF(self, grid, bw, dataset):
        """
        Sub-function that generates the PDFs of kernel
        density estimation from raw dataset
        """
        self.logger.debug("Generating PDF")
        if bw <= 0:
            self.logger.critical("bw = %d. Bandwidth cannot be negative or zero", bw)
            raise ValueError, 'bw = %d. Bandwidth cannot be negative or zero' %bw

        try:
            kdeMethod = getattr(KPDF, "UPDF%s" %self.kdeType)
        except AttributeError:
            self.logger.exception("Invalid input on option: KDE method UPDF%s does not exist" %self.kdeType)
            raise

        return kdeMethod(dataset, grid, bw)

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

    logging.basicConfig(level=logging.getattr(cnfGetIniValue(configFile, 'Logging', 'LogLevel', 'DEBUG')),
                        format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s',
                        filename=cnfGetIniValue(configFile, 'Logging', 'LogFile', __file__.rstrip('.py') + '.log'),
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
    pressure_rate = os.path.join(path,'pressure_rate')
#    bearing_rate = path+'bearing_rate'
#    speed_rate = path+'speed_rate'
    kde_pressure_rate = os.path.join(path,'kde_pressure_rate')
    cdf_pressure_rate = os.path.join(path,'cdf_pressure_rate')
#    kde_bearing_rate = path+'kde_bearing_rate'
#    cdf_bearing_rate = path+'cdf_bearing_rate'
#    kde_speed_rate = path+'kde_speed_rate'
#    cdf_speed_rate = path+'cdf_speed_rate'

    kdeStep = 0.1

    k = KDEParameters(cnfGetIniValue(configFile,'Parameters','kdeType','Biweight'), kdeStep)

    k.generateKDE(pressure_rate, kdeStep, kde_pressure_rate, cdf_pressure_rate)
    #k.generateKDE(bearing_rate, kde_bearing_rate, cdf_bearing_rate)
    #k.generateKDE(speed_rate, kde_speed_rate, cdf_speed_rate)
    #pdb.set_trace()
    #k.plotKdeInit()
    #k.plotKdeAll()
    #k.plotKdeNoInit()
    #k.plotCdf()
    #pylab.show()
