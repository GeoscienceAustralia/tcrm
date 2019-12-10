"""
:mod:`circularKDE` -- probability distributions for circular data
=================================================================

.. module:: circularKDE
    :synopsis: Probability distributions for circular data.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

Create a probability distribution function for radial-type
data, e.g. bearings. Returns the grid on which the PDF is defined, the
PDF itself and the corresponding CDF. By default, the grid is specified
on [0,2*pi).

"""

import numpy as np
from scipy.special import i0
import statsmodels.nonparametric.bandwidths as smbw

import stats

def circularKDE(parameters, kdeStep=np.pi/16.):
    """
    Create a probability distribution function for radial-type data,
    e.g. bearings. Returns the grid on which the PDF is defined, the
    PDF itself and the corresponding CDF.

    By default, the grid is specified on [0, 2*pi)

    :param parameters: :class:`numpy.ndarray` of parameter values.
    :param float kdeStep: Increment of the ordinate at which the
                          distribution will be estimated.

    :returns: :class:`numpy.ndarray` of the grid, the PDF and the CDF.

    """
    bw = stats.bandwidth(parameters)
    grid = np.arange(0, 2 * np.pi +kdeStep, kdeStep)
    pdf = np.empty(len(grid), 'float')
    chi = 1./(2 * np.pi * i0(bw))
    for k in parameters:
        kH = chi * np.exp(bw * np.cos(grid-k))
        pdf += kH/kH.sum()

    pdf = pdf/len(pdf)
    cy = stats.cdf(grid, pdf)

    return grid, pdf, cy

