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
import sys
import types
import logging as log
import numpy as np
from scipy.stats import genpareto, scoreatpercentile
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

import inspect

try:
    import lmoments as lmom
except ImportError:
    # If not, load ported python version of required functions
    # (note: this runs about half the speed of the LMOMENTS fortran package)
    import Utilities.lmomentFit as lmom
    log.debug('LMOMENTS package not found - reverting to slower python' +
              'version of code')

# Average number of observations per year: daily observations
NPYR = 365.25 

class ExtremeValueDistribution(object):

    """
    Abstract extreme value distribution model
    """

    def __init__(self, intervals, numsim, nodata, minrecords):

        self.intervals = intervals
        self.numsim = numsim
        self.nodata = nodata
        self.minrecs = minrecords

    def calculate(self, data, *args, **kwargs):
        dims = data.shape
        Rp = np.zeros((len(self.intervals),) + dims.shape[1:], dtype='f')
        loc = np.zeros(dims[1:], dtype='f')
        scale = np.zeros(dims[1:], dtype='f')
        shp = np.zeros(dims[1:], dtype='f')


        for i in range(dims[1]):
            for j in range(dims[2]):
                if data[:, i, j].max() > 0:
                    w, l, sc, sh = self.fit(data[:, i, j],  *args, **kwargs)
                    Rp[:, i, j] = w
                    loc[i, j] = l
                    scale[i, j] = sc
                    shp[i, j] = sh

        return Rp, loc, scale, shp
                
    def fit(self, data, *args, **kwargs):
        """
        Calculate the average recurrence interval values, given a set of 
        data values
        
        :param data: :class:`numpy.ndarray` of data values to use to
        calculate return levels
        
        :returns: return levels at given intervals, along with
        parameters for the distribution (if calculated)

        """
        raise NotImplementedError

class GEVDistribution(ExtremeValueDistribution):
    """
    Generalised extreme value distribution, used for block maxima

    """
    def __init__(self, intervals, numsim, nodata, minrecords):
        ExtremeValueDistribution.__init__(self, intervals, numsim, nodata, minrecords)

    def fit(self, data):
        # Initialise variables:
        loc, scale, shp = [self.nodata] * 3
        w = nodata * np.ones(len(self.intervals)) # Create empty return period array

        if (data.max() > 0.): # Check for valid data
            ii = np.flatnonzero(data) # Return indices for non-zero elements of the flattened aray
            # Only calculate l-moments for those grid points where the values are
            # not all equal, and where there are 50 or more valid (>0) values.
            if data[ii].min() != data[ii].max():
                if len(ii) >= self.minrecords:
                    l1, l2, l3 = lmom.samlmu(data, 3) # find 3 l-moments
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
                            loc, scale, shp = [self.nodata] * 3

        # Calculate return period wind speeds
        for i, t in enumerate(self.intervals): 
            # if no valid fit was found, then there are no return period wind speeds
            if shp <= 0.:
                w[i] = self.nodata
            # if a valid fit was found...
            else:
                # Calculate wind speed for each return period
                w[i] = (np.transpose(loc + (scale / shp) *
                        (1. - np.power(-1. * np.log(1. - (1. / t)), shp))))

                # Replace any non-finite numbers with the missing value:
                if not np.isfinite(w[i]):
                    w[i] = self.nodata

        return w, loc, scale, shp

class GPDDistribution(ExtremeValueDistribution):
    """
    Generalised Pareto Distribution, fitted using peaks-over-threshold
    approach. The threshold is an arbitrary percentile of the data
    values (default =99.5)

    """

    def __init__(self, intervals, numsim, nodata, minrecords, threshold=99.5):
        ExtremeValueDistribution.__init__(self, intervals, numsim, nodata, minrecords)
        self.threshold = threshold

    def fit(self, data):
        recs = data[data > 0]
        mu = scoreatpercentile(recs, self.threshold)
        loc, scale, shp = [self.nodata]*3
        w = nodata * np.ones(len(self.intervals))

        if len(data) < self.minrecords:
            return Rp, loc, scale, shp

       # Fill each day that a cyclone isn't recorded with zero so we get 
       # the correct rate for the return periods
        datafilled = np.zeros(int(self.numsim * NPYR))
        datafilled[-len(data):] = data
        log.debug("The length of the filled data is {0}".format(len(datafilled)))

        rate = float(len(datafilled[datafilled > mu])) / float(len(datafilled))
        log.debug("The calculated rate is: {0}".format(rate))

        try:
            shape, location, scale = genpareto.fit(datafilled[datafilled > mu],
                                                   floc = mu)
        except:
            return w, loc, scale, shp

        w = gpdReturnLevel(self.intervals, mu, shape, scale, rate)
        if shape > 0: # or Rpeval[0] < 0.0:
            return w, loc, scl, shp
        else:
            return w, location, scale, shape


class EMPDistribution(ExtremeValueDistribution):

    """
    Empirical average recurrence intervals

    .. note:: Empirical ARI values cannot be calculated for ARIs
    greater than the number of simulated years. If using this
    distribution, the configuration option `Years` in the `Hazard`
    section must not contain values greater than the number of
    simulations (`TrackGenerator -- NumSimulation`).

    """
    def __init__(self, intervals, numsim, nodata, minrecords):
        ExtremeValueDistribution.__init__(self, intervals, numsim, nodata, minrecords)
    
    def fit(self, data):
        loc, scale, shp = [self.nodata] * 3
        w = self.nodata * np.ones(len(self.intervals))
        nobs = self.numsim * NPYR
        datafilled = np.zeros(int(nobs))
        datafilled[-len(data):] = np.sort(data)
        emprp = 1./ (1. - np.arange(1, int(nobs) + 1, 1) / (nobs + 1)) / NPYR

        rpfunc = interp1d(emprp, np.sort(datafilled), kind='linear', fill_value=self.nodata)

        w = rpfunc(intervals)

        return w, loc, scale, shp

class POWERDistribution(ExtremeValueDistribution):

    """
    Fit a function of the form w = A - B * ARI^C to the empirical 
    average recurrence intervals (ARI), where A > 0, 0 < B< A and C < 0.

    """

    def __init__(self, intervals, numsim, nodata, minrecords):
        ExtremeValueDistribution.__init__(self, intervals, numsim, nodata, minrecords)

    def fit(self, data):
        loc, scale, shp = [self.nodata] * 3
        w = self.nodata * np.ones(len(self.intervals)) 

        nobs = self.numsim * NPYR
        datafilled = np.zeros(int(nobs))
        datafilled[-len(data):] = np.sort(data)
        emprp = 1./ (1. - np.arange(1, int(nobs) + 1, 1) / (nobs + 1)) / NPYR

        idx = np.where(datafilled > np.mean(datafilled[datafilled > 0]))[0]
        
        def func(x, a, b, c):
            return a - b * np.power(x, c)

        try:
            pars, covar = curve_fit(func, emprp[idx], wspd[idx], 
                                p0=[np.max(data), np.max(data), -0.1],
                                maxfev=10000)
        except:
            return w, loc, scale, shp

        else:
            w = func(self.intervals, *pars)
            loc, scale, shp = pars
            return w, loc, scale, shp
        
def gevfit(data, intervals, nodata=-9999., minrecords=50, yrspersim=1):
    """
    Calculate extreme value distribution parameters using the Lmoments module.
    Return period values are not calculated if the shape parameter is negative
    or zero.

    :param data: array of data values. Values represent max events for each year 
              of simulation at a single grid box
    :type data: :class:`numpy.ndarray`
    :param intervals: array of years for which to calculate return period values.
    :type intervals: :class:`numpy.ndarray`
    :param float nodata: value to insert if fit does not converge.
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
    nodata = float(nodata)
    intervals = np.array(intervals)

    # Initialise variables:
    loc, scale, shp = [nodata, nodata, nodata]
    w = nodata * np.ones(len(intervals)) # Create empty return period array

    if (data.max() > 0.): # Check for valid data
        ii = np.flatnonzero(data) # Return indices for non-zero elements of the flattened aray
        # Only calculate l-moments for those grid points where the values are
        # not all equal, and where there are 50 or more valid (>0) values.
        if data[ii].min() != data[ii].max():
            if len(ii) >= minrecords:
                l1, l2, l3 = lmom.samlmu(data, 3) # find 3 l-moments
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
                        loc, scale, shp = [nodata, nodata, nodata]

    # Calculate return period wind speeds
    for i, t in enumerate(intervals): 
        # if no valid fit was found, then there are no return period wind speeds
        if shp <= 0.:
            w[i] = nodata
        # if a valid fit was found...
        else:
            # Calculate wind speed for each return period
            w[i] = (np.transpose(loc + (scale / shp) *
                    (1. - np.power(-1. * np.log(1. - (yrspersim / t)), shp))))

            # Replace any non-finite numbers with the missing value:
            if not np.isfinite(w[i]):
                w[i] = nodata

    return w, loc, scale, shp

def powerfit(data, intervals, numsim, nodata=-9999., minrecords=50):

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
    :param intervals: :class:`numpy.ndarray` of years for which to calculate
                  return period values.
    :param int numsim: number of simulations created.
    :param float nodata: value to insert if fit does not converge.
    :param int minrecords: minimum number of valid observations 
                           required to perform fitting.

    Returns:
    --------

    :param w: `numpy.array` of return period wind speed values
    :param location: location parameter
    :param scale: scale parameter
    :param shape: shape parameter

    """

    loc, scale, shp = [nodata, nodata, nodata]
    w = nodata * np.ones(len(intervals)) # Create empty return period array
    
    npyr = 365.25
    nobs = npyr * numsim
    
    wspd = np.zeros(int(nobs))
    wspd[-len(data):] = np.sort(data)

    emprp = 1./ (1. - np.arange(1, nobs + 1, 1) / (nobs + 1)) / npyr
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
        w = func(intervals, *pars)
        loc, scale, shp = pars
        return w, loc, scale, shp
    
def empfit(data, intervals, numsim, nodata=-9999., minrecords=50):
    """
    Calculate empirical ARI values for a collection 
    of wind speed records.
    
    

    :param data: :class:`numpy.ndarray` of data values
    :param intervals: :class:`numpy.ndarray` of years for which to calculate
                  return period values. The values will be determined 
                  empirically, then interpolated to these intervals.
    :param int numsim: number of simulations created.
    :param float nodata: value to insert if fit does not converge.
    :param int minrecords: minimum number of valid observations 
                           required to perform fitting.

    Returns:
    --------

    :param Rpeval: `numpy.array` of return period wind speed values
    :param location: location parameter
    :param scale: scale parameter
    :param shape: shape parameter

    """

    loc, scale, shp = [nodata, nodata, nodata]
    w = nodata * np.ones(len(intervals)) # Create empty return period array
    
    npyr = 365.25
    nobs = npyr * numsim
    
    wspd = np.zeros(int(nobs))
    wspd[-len(data):] = np.sort(data)

    emprp = 1./ (1. - np.arange(1, int(nobs) + 1, 1) / float((nobs + 1))) / npyr

    rpfunc = interp1d(emprp, wspd, kind='linear', fill_value=nodata)
    try:
        w = rpfunc(intervals)
    except ValueError:
        log.exception("Empirical ARI values cannot be greater than the length of simulation")
        log.exception("Highest empirical ARI: {0}".format(emprp[-1]))
        log.exception("ARIs requested: {0}".format(repr(intervals)))
        raise

    return w, loc, scale, shp

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

def gpdfit(data, intervals, numsim, nodata=-9999, minrecords=50, threshold=99.5):
    """
    Fit a Generalised Pareto Distribution to the data. For a quick evaluation,
    we use the 99.5th percentile as a threshold. 

    :param data: array of data values.
    :type data: :class:`numpy.ndarray`
    :param intervals: array of years for which to calculate return period values.
    :param int numsim: number of simulations created.
    :type intervals: :class:`numpy.ndarray`
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

    loc, scl, shp = [nodata, nodata, nodata]
    Rp = nodata * np.ones(len(intervals))

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

    Rpeval = gpdReturnLevel(intervals, mu, shape, scale, rate)
    if shape > 0: # or Rpeval[0] < 0.0:
        return Rp, loc, scl, shp
    else:
        return Rpeval, location, scale, shape

def islocal(func):
    return isinstance(func, types.FunctionType) \
        and func.__module__ == __name__ \
        and func.__name__.endswith('fit')

#EVFUNCS = dict(inspect.getmembers(sys.modules[__name__], \
#                                  predicate=islocal))

# Automatic discovery of models and required parameters


def allSubclasses(cls):
    """
    Recursively find all subclasses of a given class.
    """
    return cls.__subclasses__() + \
        [g for s in cls.__subclasses__() for g in allSubclasses(s)]

def evfunc(name):
    return EVFUNCS[name]

def evargs(name):
    from inspect import getargspec
    std = getargspec(ExtremeValueDistribution.__init__)[0]
    new = getargspec(evfunc(name).__init__)[0]
    params = [p for p in new if p not in std]
    return params

EVFUNCS = dict([(k.__name__.replace('Distribution', '').lower(), k)
               for k in allSubclasses(vars()['ExtremeValueDistribution'])])
