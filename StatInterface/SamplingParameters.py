"""
:mod:`SamplingParameters` -- Sample TC parameters from distributions
====================================================================

.. module:: SamplingParamters
    :synopsis: Generate samples of TC parameters.

.. moduleauthor:: Geoff Xu <geoff.xu@ga.gov.au>

Defines the class for sampling cyclone parameters. Can generate either
a single sample, or an array of samples (for multiple cyclones).

"""

import os
import sys
import logging

import scipy
import numpy as np
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile, flSaveFile


class SamplingParameters:
    """
    Provides methods to sample one or many values from a CDF of
    parameter values.

    :param cdfParameters: Name of a file containing the CDF of
                          a parameter, or the actual CDF values.
    :type  cdfParamters: str, :class:`numpy.ndarray` or None


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
        elif type(cdfParameters) == np.ndarray:
            self.xacy = cdfParameters
        else:
            self.xacy = np.array([])

    def setParameters(self, cdfParameters):
        """
        Set parameters.

        :param cdfParameters: Name of a file containing the CDF of
                              a parameter, or the actual CDF values.
        :type  cdfParamters: str or :class:`numpy.ndarray`

        :raises IOError: If the CDF files do not exist.
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
        """Generate a single random sample of cyclone parameters."""

        unif = scipy.rand()
        ind_kdf = self.xacy[:, 1].searchsorted(unif)
        return self.xacy[ind_kdf, 0]

    def generateSamples(self, ns, sample_parameter_path=None):
        """
        Generate random samples of cyclone initial parameters.

        :param int ns: Number of samples to generate.
        :param str sample_parameter_path: Path to a file to save the
                                          sampled parameter values to.

        :returns: The sample values.
        :raises ValueError: If ns <= 0.

        """

        if ns <= 0:
            raise ValueError('invalid input on ns: number of sample cannot be zero or negative')

        unif_s = scipy.rand(ns)

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
            raise IOError(error_msg)
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError(error_msg)

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s',
                        filename=cnfGetIniValue(configFile, 'Logging', 'LogFile', __file__.rstrip('.py') + '.log'),
                        filemode='w')

    path = cnfGetIniValue(configFile, 'Output','Path')

    cdfParameters = os.path.join(path,'all_cell_cdf_pressure_rate')

    sp = SamplingParameters(cdfParameters)

    os = sp.generateOneSample()
    s = sp.generateSamples(cnfGetIniValue(configFile, 'Parameters',
                                          'Samples'))
