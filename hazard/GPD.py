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

def gpdfit(data, years, numsim, missingValue=-9999,
           minrecords=50, threshold=99.5):
    """
    Fit a Generalised Pareto Distribution to the data. For a quick evaluation,
    we use the 99.5th percentile as a threshold. 

    :param data: array of data values.
    :type data: :class:`numpy.ndarray`
    :param years: array of years for which to calculate return period values.
    :param int numsim: number of simulations created.
    :type years: :class:`numpy.ndarray`
    :param float missingValue: value to insert if fit does not converge.
    :param int minrecords: minimum number of valid observations required to
                           perform fitting.
    :param float threshold: Threshold for performing the fitting. Default is 
                            the 99.5th percentile

    Returns:
    --------

    :param Rpeval: `numpy.array` of return period wind speed values
    :param location: location parameter
    :param scale: scale parameter
    :param shape: shape parameter
    """
    recs = data[data > 0]
    mu = scoreatpercentile(data, threshold)

    loc, scl, shp = [missingValue, missingValue, missingValue]
    Rp = missingValue * np.ones(len(years))

    log.debug("The length of the data currently is {0}".format(len(data)))

    if len(data) < minrecords:
        return Rp, loc, scl, shp

    # Fill each day that a cyclone isn't recorded with zero so we get 
    # the correct rate for the return periods
    datafilled = np.zeros(int(numsim * 365.25))
    datafilled[-len(data):] = data
    log.debug("The length of the filled data is {0}".format(len(datafilled)))

    rate = float(len(datafilled[datafilled > mu])) / float(len(datafilled))
    log.debug("The calculated rate is: {0}".format(rate))

    try:
        shape, location, scale = genpareto.fit(datafilled[datafilled > mu],
                                               floc = mu)
    except:
        return Rp, loc, scl, shp

    Rpeval = gpdReturnLevel(years, mu, shape, scale, rate)
    if shape > 0: # or Rpeval[0] < 0.0:
        return Rp, loc, scl, shp
    else:
        return Rpeval, location, scale, shape

