"""
:mod:`GPD` -- Calculate extreme value distribution using GPD
============================================================

.. module: GPD
    :synopsis: Calculate parameters for the GPD distribution

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

Calculate the parameters for the GPD distribution. The function will 
not only fit the distribution parameters, but also calculate the return 
period values for specified return periods.

References:
Hosking, J. R. M., 1990: L-moments: Analysis and Estimation of
Distributions using Linear Combinations of Order Statistics. Journal
of the Royal Statistical Society, 52, 1, 105-124.

"""

import logging as log
import numpy as np
from scipy.stats import genpareto, scoreatpercentile

def gpdReturnLevel(intervals, mu, shape, scale, rate, npyr=365.25):
    """
    Calculate return levels for specified intervals for a distribution with
    the given threshold, scale and shape parameters.

    :param intervals: :class:`numpy.ndarray` or float of recurrence intervals
              to evaluate return levels for.
    :param float mu: Threshold parameter (also called location).
    :param float shape: Shape parameter.
    :param float scale: Scale parameter.
    :param float rate: Rate of exceedances (i.e. number of observations greater
                       than `mu`, divided by total number of observations).
    :param float npyr: Number of observations per year.

    :returns: return levels for the specified recurrence intervals.

    """

    rp = mu + (scale / shape) * (np.power(intervals * npyr * rate, shape) - 1.)
    return rp

def gpdfit(data, years, missingValue=-9999, minrecords=50, thresh=99.5):
    """
    Fit a Generalised Pareto Distribution to the data. For a quick evaluation,
    we use the 99.5th percentile as a threshold. 

    :param data: array of data values.
    :type data: :class:`numpy.ndarray`
    :param years: array of years for which to calculate return period values.
    :type years: :class:`numpy.ndarray`
    :param float missingValue: value to insert if fit does not converge.
    :param int minRecords: minimum number of valid observations required to
                           perform fitting.
    :param float thresh: Threshold for performing the fitting. Default is the
                     99.5th percentile

    Returns:
    --------

    :param Rp: `numpy.array` of return period wind speed values
    :param mu: location parameter
    :param params[2]: scale parameter
    :param params[0]: shape parameter
    """
    mu = scoreatpercentile(data, thresh)

    loc, scale, shp = [missingValue, missingValue, missingValue]
    Rp = missingValue * np.ones(len(years))

    if len(data[data > 0]) < minrecords:
        return Rp, loc, scale, shp

    rate = float(len(data[data > mu])) / float(len(data))

    try:
        shape, loc, scale = genpareto.fit(data[data > mu], floc = mu)
    except:
        return Rp, loc, scale, shp

    Rp = gpdReturnLevel(years, mu, shape, scale, rate)
    
    return Rp, mu, scale, shape

