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

# import logging as log
# import numpy as np
# from scipy.stats import genpareto, scoreatpercentile

# def gpdReturnLevel(intervals, mu, shape, scale, rate, npyr=365.25):
#     """
#     Calculate return levels for specified intervals for a distribution with
#     the given threshold, scale and shape parameters.

#     :param intervals: :class:`numpy.ndarray` or float of recurrence intervals
#               to evaluate return levels for.
#     :param float mu: Threshold parameter (also called location).
#     :param float shape: Shape parameter.
#     :param float scale: Scale parameter.
#     :param float rate: Rate of exceedances (i.e. number of observations greater
#                        than `mu`, divided by total number of observations).
#     :param float npyr: Number of observations per year.

#     :returns: return levels for the specified recurrence intervals.

#     """

#     rp = mu + (scale / shape) * (np.power(intervals * npyr * rate, shape) - 1.)
#     return rp


# def gpdRecurrenceIntervals(return_levels, mu, shape, scale, rate, npyr=365.25):
#     """
#     Calculate recurrence intervals for specified return levels for a distribution with
#     the given threshold, scale and shape parameters.

#     :param intervals: :class:`numpy.ndarray` or float of return levels
#               to evaluate recurrence intervals for.
#     :param float mu: Threshold parameter (also called location).
#     :param float shape: Shape parameter.
#     :param float scale: Scale parameter.
#     :param float rate: Rate of exceedances (i.e. number of observations greater
#                        than `mu`, divided by total number of observations).
#     :param float npyr: Number of observations per year.

#     :returns: recurrence intervals for the specified return levels.

#     """
#     ri = np.power((return_levels - mu) * (shape / scale) + 1, 1 / shape) / (npyr * rate)
#     return ri


# def gpdfit(data, years, numsim, missingValue=-9999,
#            minrecords=50, threshold=99.5):
#     """
#     Fit a Generalised Pareto Distribution to the data. For a quick evaluation,
#     we use the 99.5th percentile as a threshold. 

#     :param data: array of data values.
#     :type data: :class:`numpy.ndarray`
#     :param years: array of years for which to calculate return period values.
#     :param int numsim: number of simulations created.
#     :type years: :class:`numpy.ndarray`
#     :param float missingValue: value to insert if fit does not converge.
#     :param int minrecords: minimum number of valid observations required to
#                            perform fitting.
#     :param float threshold: Threshold for performing the fitting. Default is 
#                             the 99.5th percentile

#     Returns:
#     --------

#     :param Rpeval: `numpy.array` of return period wind speed values
#     :param location: location parameter
#     :param scale: scale parameter
#     :param shape: shape parameter
#     """
#     recs = data[data > 0]
#     mu = scoreatpercentile(data, threshold)

#     loc, scl, shp = [missingValue, missingValue, missingValue]
#     Rp = missingValue * np.ones(len(years))

#     log.debug("The length of the data currently is {0}".format(len(data)))
#     # breakpoint()
#     if len(data) < minrecords:
#         return Rp, loc, scl, shp

#     # Fill each day that a cyclone isn't recorded with zero so we get 
#     # the correct rate for the return periods
#     datafilled = np.zeros(int(numsim * 365.25))
#     datafilled[-len(data):] = data
#     log.debug("The length of the filled data is {0}".format(len(datafilled)))

#     rate = float(len(datafilled[datafilled > mu])) / float(len(datafilled))
#     log.debug("The calculated rate is: {0}".format(rate))

#     try:
#         shape, location, scale = genpareto.fit(datafilled[datafilled > mu],
#                                                floc = mu)
#     except:
#         return Rp, loc, scl, shp

#     Rpeval = gpdReturnLevel(years, mu, shape, scale, rate)

#     if shape > 0: # or Rpeval[0] < 0.0:
#         return Rp, loc, scl, shp
#     else:
#         return Rpeval, location, scale, shape

#     # return Rpeval, location, scale, shape


import numpy as np
from scipy.stats import genpareto
from scipy.stats.mstats import scoreatpercentile

# def gpdReturnLevel(intervals, mu, shape, scale, rate, npyr=365.25):
def gpdReturnLevel(intervals, mu, shape, scale, rate):
    """
    Your existing return-level formula.
    RL(T) = u + (σ/ξ) * ( (T * rate)^{ξ} - 1 )
    """
    return mu + (scale / shape) * ((intervals * rate)**shape - 1)


def gpdfit(data, intervals, numsim, missingValue=-9999,
           minrecords=50, threshold=99.5):

    recs = data[data > 0]
    mu = scoreatpercentile(data, threshold)

    loc, scl, shp = [missingValue, missingValue, missingValue]
    Rp = missingValue * np.ones(len(intervals))

    # Fill zeros for days without cyclones
    datafilled = np.zeros(int(numsim * 365.25))
    datafilled[-len(data):] = data

    rate = float(len(datafilled[datafilled > mu])) / float(len(datafilled))

    if len(data) < minrecords:
        return Rp, loc, scl, shp, rate, mu, datafilled

    try:
        shape, location, scale = genpareto.fit(
            datafilled[datafilled > mu], floc=mu
        )
    except:
        return Rp, loc, scl, shp, rate, mu, datafilled

    Rpeval = gpdReturnLevel(intervals, mu, shape, scale, rate)

    if shape > 0:
        return Rp, loc, scl, shp, rate, mu, datafilled
    else:
        return Rpeval, location, scale, shape, rate, mu, datafilled


# ---------------------------------------------------------------
# NEW: Bootstrap Confidence Interval for GPD Return Levels
# ---------------------------------------------------------------
def gpd_return_level_CI(datafilled, mu, shape, scale, rate,
                        intervals, B=1000, ci=95):
    """
    Compute bootstrap CI for GPD return levels.
    
    Parameters
    ----------
    datafilled : array
        Zero-filled data series
    mu, shape, scale : floats
        Fitted GPD parameters
    rate : float
        Exceedance rate (per day)
    years : array
        Return periods
    B : int
        Number of bootstrap samples
    ci : float
        Percent confidence interval (95, 99, etc.)
    """

    alpha = (100 - ci) / 2
    bootstrap_RL = np.zeros((B, len(intervals)))

    N = len(datafilled)

    for b in range(B):
        # ---- 1. simulate exceedances from fitted GPD ----
        nexceed = int(rate * N)
        y = genpareto.rvs(c=shape, loc=mu, scale=scale, size=nexceed)

        # ---- 2. reconstruct full time series ----
        sim_series = np.zeros(N)
        sim_series[-nexceed:] = y

        # ---- 3. refit GPD ----
        try:
            shp_b, loc_b, scl_b = genpareto.fit(
                sim_series[sim_series > mu], floc=mu
            )
        except:
            continue

        rate_b = float(np.sum(sim_series > mu)) / float(N)

        # ---- 4. compute bootstrap return levels ----
        bootstrap_RL[b, :] = gpdReturnLevel(intervals, mu, shp_b, scl_b, rate_b)

    # Remove possible rows of zeros (failed fits)
    bootstrap_RL = bootstrap_RL[np.any(bootstrap_RL != 0, axis=1)]

    # ---- 5. CI from percentiles ----
    lower = np.percentile(bootstrap_RL, alpha, axis=0)
    upper = np.percentile(bootstrap_RL, 100 - alpha, axis=0)

    return lower, upper