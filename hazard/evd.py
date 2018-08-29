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
from scipy.optimize import curve_fit

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

    :param v: array of data values. Values represent max events for each year 
              of simulation at a single grid box
    :type v: :class:`numpy.ndarray`
    :param years: array of years for which to calculate return period values.
    :type years: :class:`numpy.ndarray`
    :param float missingValue: value to insert if fit does not converge.
    :param int minRecords: minimum number of valid observations required to
                           perform fitting.
    :param int yrspersim: data represent block maxima - this gives the length
                          of each block in years.

    Returns:
    --------

    :param w: `numpy.array` of return period wind speed values
    :param loc: location parameter
    :param scale: scale parameter
    :param shp: shape parameter

    """
    # Convert to float to prevent integer division & ensure consistent data
    # types for output variables
    yrspersim = float(yrspersim)
    missingValue = float(missingValue)
    years = np.array(years)

    # Initialise variables:
    loc, scale, shp = [missingValue, missingValue, missingValue]
    w = missingValue * np.ones(len(years)) # Create empty return period array

    if (v.max() > 0.): # Check for valid data
        ii = np.flatnonzero(v) # Return indices for non-zero elements of the flattened aray
        # Only calculate l-moments for those grid points where the values are
        # not all equal, and where there are 50 or more valid (>0) values.
        if v[ii].min() != v[ii].max():
            if len(ii) >= minRecords:
                l1, l2, l3 = lmom.samlmu(v, 3) # find 3 l-moments
                # t3 = L-skewness (Hosking 1990)
                t3 = l3 / l2
                if (l2 <= 0.) or (np.abs(t3) >= 1.):
                    # Reject points where the second l-moment is negative
                    # or the ratio of the third to second is > 1, i.e. positive
                    # skew.
                    log.debug("Invalid l-moments")
                else:
                    # Parameter estimation returns the location, scale and
                    # shape parameters
                    xmom = [l1, l2, t3]
                    # Calculates the parameters of the distribution given its
                    # L-moments. GEV distribution calculated.
                    loc, scale, shp = np.array(lmom.pelgev(xmom))
                    # We only store the values if the first parameter is
                    # finite (i.e. the location parameter is finite)
                    if not np.isfinite(loc):
                        loc, scale, shp = [missingValue,
                                           missingValue,
                                           missingValue]

    # Calculate return period wind speeds
    for i, t in enumerate(years): 
        # if no valid fit was found, then there are no return period wind speeds
        if shp <= 0.:
            w[i] = missingValue
        # if a valid fit was found...
        else:
            # Calculate wind speed for each return period
            w[i] = (np.transpose(loc + (scale / shp) *
                    (1. - np.power(-1. * np.log(1. - (yrspersim / t)), shp))))

            # Replace any non-finite numbers with the missing value:
            if not np.isfinite(w[i]):
                w[i] = missingValue

    return w, loc, scale, shp

def powerfit(data, years, numsim, missingValue=-9999.,
             minrecords=50):

    """Fit a modified power law to empirical return period values

    w(t; a, b, c) = a - b * t^c
    
    where a > 0, 0 < b < a and c < 0.

    This function calculates the empirical return periods for the
    given set of wind speed values, then fits the function to the
    return period wind speeds. In this sense, the returned location,
    scale and shape parameters do not refer to the parameters of any 
    given distribution, but rather the coefficients of the fitted 
    function.

    
    :param data: :class:`numpy.ndarray` of data values
    :param years: :class:`numpy.ndarray` of years for which to calculate
                  return period values.
    :param int numsim: number of simulations created.
    :param float missingValue: value to insert if fit does not converge.
    :param int minrecords: minimum number of valid observations 
                           required to perform fitting.

    Returns:
    --------

    :param Rpeval: `numpy.array` of return period wind speed values
    :param location: location parameter
    :param scale: scale parameter
    :param shape: shape parameter

    """

    loc, scale, shp = [missingValue, missingValue, missingValue]
    w = missingValue * np.ones(len(years)) # Create empty return period array
    
    npyr = 365.25
    nobs = npyr * numsim
    
    wspd = np.zeros(int(nobs))
    wspd[-len(data):] = np.sort(data)

    emprp = 1./ (1. - np.arange(1, nobs + 1, 1)) / (nobs + 1) / npyr
    idx = np.where(wspd > np.mean(wspd[wspd > 0]))[0]

    def func(x, a, b, c):
        return a - b * np.power(x, c)

    try:
        pars, covar = curve_fit(func, emprp[idx], wspd[idx], 
                                p0=[np.max(data), np.max(data), -0.1],
                                maxfev=10000)
    except:
        return w, loc, scale, shp

    else:
        w = func(years, *pars)
        loc, scale, shp = pars
        return w, loc, scale, shp
    
