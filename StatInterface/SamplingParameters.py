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


Title: samplingParameters.py - generate samples of cyclone parameters
Author: Geoff Xu, geoff.xu@ga.gov.au
CreationDate: 2006-01-09
Description: Defines the class for sampling cyclone parameters. Can
             generate either a single sample, or an array of samples
             (for multiple cyclones)

ModifiedBy: Nariman Habili, nariman.habili@ga.gov.au
ModifiedDate: 2006-11-29
Modification: Added accessor methods, upgraded to ndarray and conformed
              with style guide

ModifiedBy: Nariman Habili, nariman.habili@ga.gov.au
ModifiedDate: 2007-04-02
Modification: Added setParameters

SeeAlso: (related programs)
Constraints:

$Id: SamplingParameters.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

from scipy import rand, ndarray, array
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile, flSaveFile

#-----------------------------------------------------------------------
# 'SamplingParameter' class:
#-----------------------------------------------------------------------
class SamplingParameters:
    """
    Parameters
    ----------
    cdfParameters : string (file name including path) or array
        CDF of cyclone parameters generated from kde classes
    ns : integer
        number of samples to generate

    Members
    -------
    xacy : 2D array of float
        xy coordinates of cdf of cyclone initial parameters
    ns : int
        number of samples to generate
    sample : 1D array of float
        random samples been generated from generateSample method

    Methods
    -------
    generateOneSample() : float
        Generate a random sample of cyclone parameters
    generateSamples(ns,sample_parameter_path)
        Generate random samples of cyclone initial parameters
    plotSamples(original_parameters,kde_parameters)
        plot the orignal historical data as well as the random samples
        generated for cyclone initial parameters

    Accessor Methods
    ----------------
    setParameters(cdfParameters) : string or array
        set parameters
    setName(n) : string
        set the name of plot
    name()
        get the name of plot
    setXlabel(xLabel) : string
        set the label name on x axis
    xLabel()
        get the label name on x axis
    setBarWidth(barWidth) : integer
        set the width of the bar on the plot
    barWidth()
        get the width of the bar on the plot
    """
    def __init__(self, cdfParameters=None):
        """Initialize the data needed for the plots including CDF of
        cyclone parameters
        """

        self.logger = logging.getLogger()

        if type(cdfParameters) == str:
            try:
                self.xacy = flLoadFile(cdfParameters)
            except IOError:
                self.logger.exception('Error! Files relating to cdf of cyclone parameters does not exist, please generate KDE of cyclone parameters first.')
                raise
        elif type(cdfParameters) == ndarray:
            self.xacy = cdfParameters
        else:
            self.xacy = array([])

    def __doc__(self):
        """Documentation on what this class does
        """

        return "Cyclone Parameters Sampling from the CDF"

    def setParameters(self, cdfParameters):
        """Set parameters
        """

        if type(cdfParameters) == str:
            try:
                self.xacy = flLoadFile(cdfParameters)
            except IOError:
                self.logger.exception('Error! Files relating to cdf of cyclone parameters does not exist, please generate KDE of cyclone parameters first.')
                raise
        else:
            self.xacy = cdfParameters

        return self

    def generateOneSample(self):
        """Generate a random sample of cyclone parameters
        """

        unif = rand()
        ind_kdf = self.xacy[:, 1].searchsorted(unif)
        return self.xacy[ind_kdf, 0]

    def generateSamples(self, ns, sample_parameter_path=None):
        """Generate random samples of cyclone initial parameters
        """

        if ns <= 0:
            raise ValueError, 'invalid input on ns: number of sample cannot be zero or negative'

        unif_s = rand(ns)

        ind_kdf = self.xacy[:, 1].searchsorted(unif_s)
        self.sample = self.xacy[ind_kdf, 0]

        if sample_parameter_path:
            flSaveFile(sample_parameter_path, self.sample)
        else:
            return self.sample


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

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s',
                        filename=cnfGetIniValue(configFile, 'Logging', 'LogFile', __file__.rstrip('.py') + '.log'),
                        filemode='w')

    path = cnfGetIniValue(configFile, 'Output','Path')

    cdfParameters = os.path.join(path,'all_cell_cdf_pressure_rate')

    sp = SamplingParameters(cdfParameters)

    #pdb.set_trace()
    os = sp.generateOneSample()
    s = sp.generateSamples(cnfGetIniValue(configFile, 'Parameters',
                                          'Samples'))
