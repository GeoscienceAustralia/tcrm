"""
:mod:`evd` -- Calculate extreme value distributions
===================================================

.. module: evd
    :synopsis: Calculate parameters for the GEV distribution
               using L-moments

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

Calculate the parameters for the GEV distribution, using the method of
L-moments to estimate the parameters. The function will not only fit
the distribution parameters, but also calculate the return period
values for specified return periods.

References:
Hosking, J. R. M., 1990: L-moments: Analysis and Estimation of
Distributions using Linear Combinations of Order Statistics. Journal
of the Royal Statistical Society, 52, 1, 105-124.


"""

import logging as log
import numpy as np
from scipy.stats import genpareto, scoreatpercentile

try:
    import lmoments as lmom
except ImportError:
    # If not, load ported python version of required functions
    # (note: this runs about half the speed of the LMOMENTS fortran package)
    import Utilities.lmomentFit as lmom
    log.debug('LMOMENTS package not found - reverting to slower python' +
              'version of code')

def estimateEVD(v, years, missingValue=-9999., minRecords=50, yrspersim=1):
    """
    Calculate extreme value distribution parameters using the Lmoments module.
    Return period values are not calculated if the shape parameter is negative
    or zero.

    :param v: array of data values.
    :type v: :class:`numpy.ndarray`
    :param years: array of years for which to calculate return period values.
    :type years: :class:`numpy.ndarray`
    :param float missingValue: value to insert if fit does not converge.
    :param int minRecords: minimum number of valid observations required to
                           perform fitting.
    :param int yrspersim: data represent block maxima - this gives the length
                          of each block in years.

    :return: return period values
    :rtype: :class:`numpy.ndarray`
    :return: location, shape and scale parameters of the distribution
    :rtype: float

    """
    # Convert to float to prevent integer division & ensure consistent data
    # types for output variables
    yrspersim = float(yrspersim)
    missingValue = float(missingValue)
    years = np.array(years)

    # Initialise variables:
    loc, scale, shp = [missingValue, missingValue, missingValue]
    w = missingValue * np.ones(len(years))

    if (v.max() > 0.):
        ii = np.flatnonzero(v)
        # Only calculate l-moments for those grid points where the values are
        # not all equal, and where there are 50 or more valid (>0) values.
        if v[ii].min() != v[ii].max():
            if len(ii) >= minRecords:
                l1, l2, l3 = lmom.samlmu(v, 3)
                t3 = l3 / l2
                if (l2 <= 0.) or (np.abs(t3) >= 1.):
                    # Reject points where the second l-moment is negative
                    # or the ratio of the third to second is > 1.
                    log.debug("Invalid l-moments")
                else:
                    # Parameter estimation returns the location, scale and
                    # shape parameters
                    xmom = [l1, l2, t3]
                    loc, scale, shp = np.array(lmom.pelgev(xmom))
                    # We only store the values if the first parameter is
                    # finite (i.e. the location parameter is finite)
                    if not np.isfinite(loc):
                        loc, scale, shp = [missingValue,
                                           missingValue,
                                           missingValue]

    for i, t in enumerate(years):
        if shp <= 0.:
            w[i] = missingValue
        else:
            w[i] = (np.transpose(loc + (scale / shp) *
                    (1. - np.power(-1. * np.log(1. - (yrspersim / t)), shp))))

            # Replace any non-finite numbers with the missing value:
            if not np.isfinite(w[i]):
                w[i] = missingValue

    return w, loc, scale, shp

def gpdReturnLevel(intervals, mu, xi, sigma, rate, npyr=365.25):
    """
    Calculate return levels for specified intervals for a distribution with
    the given threshold, scale and shape parameters.

    :param intervals: :class:`numpy.ndarray` or float of recurrence intervals
              to evaluate return levels for.
    :param float mu: Threshold parameter (also called location).
    :param float xi: Shape parameter.
    :param float sigma: Scale parameter.
    :param float rate: Rate of exceedances (i.e. number of observations greater
                       than `mu`, divided by total number of observations).
    :param float npyr: Number of observations per year.

    :returns: return levels for the specified recurrence intervals.

    """

    rp = mu + (sigma / xi) * (np.power(intervals * npyr * rate, xi) - 1.)
    return rp

def gpdfit(data, rate, years, missingValue=-9999, minrecords=50, thresh=99.7):
    """
    Fit a Generalised Pareto Distribution to the data. For a quick evaluation,
    we use the 99.7th percentile as a threshold. 

    :param data: array of data values.
    :type data: :class:`numpy.ndarray`
    :param years: array of years for which to calculate return period values.
    :type years: :class:`numpy.ndarray`
    :param float missingValue: value to insert if fit does not converge.
    :param int minRecords: minimum number of valid observations required to
                           perform fitting.
    :param float thresh: Threshold for performing the fitting. Default is the
                     99.7th percentile

    :return: return period values
    :rtype: :class:`numpy.ndarray`
    :return: location, shape and scale parameters of the distribution,
             determined using MLE (as per :meth:`scipy.stats.rv_continuous`)
    :rtype: float
    """

    mu = scoreatpercentile(data, thresh)

    if len(data[data > 0]) < minrecords:
        return 
    loc, scale, shp = [missingValue, missingValue, missingValue]
    w = missingValue * np.ones(len(years))

    if len(data[data > 0]) < minrecords:
        return w, loc, scale, shp

    rate = float(len(data[data > mu])) / float(len(mu))

    try:
        params = genpareto.fit(data[data > mu], floc=mu)
    except:
        return w, loc, scale, shp

    w = gpdReturnLevel(years, mu, params[0], params[2], rate)
    
    return w, mu, params[2], params[0]
