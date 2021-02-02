"""
:mod:`trackSize` -- determine radius to max winds
=================================================

.. module:: trackSize
    :synopsis: Determine the radius to maximum winds for each cyclone
               using either statistics from observations or a
               parametric model based on latitude and central pressure
               deficit.
.. moduleauthor:: Craig Arthur, <craig.arthur@ga.gov.au>

"""

import numpy as np
import logging
import scipy.stats as stats
import statsmodels.api as sm

LOG = logging.getLogger()

def rmax(dp, lat, eps, coeffs=[3.5843946536979779,-0.0045486143609339436,
                               0.78621467400844858, 0.0024030344245284741,
                               0.0015567629057007433]):
    """
    Calculate radius to maximum wind based on pressure deficit and
    latitude. This function allows for the random variate to be set
    when calling the function. Default coefficients for the functional
    form of ln(Rmw) are given, based on JTWC data for the southern hemisphere.
    
    ln(Rmw) = a + b*dp + c*exp(-d*dp^2) + f*|lat| + eps

    eps is not included in the coefficients (though that may be considered
    by some to be more logical), so that it can remain constant for a single
    TC event. 

    :param dp: Central pressure deficit (hPa)
    :param lat: Latitude of the storm (degrees)
    :param eps: random variate. This would normally be held constant
                for a single storm.
    :param coeffs: A list of coefficients for the functional form. Default
                   values are based on JTWC data from the southern hemisphere.

    :returns: radius to maximum wind value.

    """
    if len(coeffs) < 4:
        LOG.warn("Insufficient coefficients for rmw calculation!")
        LOG.warn("Using default values")
        coeffs = [3.5843946536979779,-0.0045486143609339436,
                  0.78621467400844858, 0.0024030344245284741,
                  0.0015567629057007433]

    if isinstance(dp, (np.ndarray, list)) and \
      isinstance(lat, (np.ndarray, list)):
        assert len(dp) == len(lat)
    yy = coeffs[0] + coeffs[1]*dp + coeffs[2] * np.exp(-coeffs[3] * dp * dp) +\
         coeffs[4] * np.abs(lat) + eps
    rm = np.exp(yy)
    return rm

def fitRmax(rmw, dp, lat):
    """
    Fit Rmw data to a function of pressure deficit and latitude.

    We fit a function of dp and latitude to ln(Rmw) values of the
    form:

    ln(Rmw) = a + b*dp + c*dp^2 + d*lat^2 + eps

    where eps is a random normal variate with zero mean and std. dev.
    describing the residual variance. 

    :param rmw: :class:`numpy.ndarray` of valid Rmw observations.
    :param dp: :class:`numpy.ndarray` of valid pressure deficit.
               observations, corresponding to the Rmw observations.
    :param lat: :class:`numpy.ndarray` of latitude observations for
                the Rmw observations.

    :returns: list of coefficients for the functional form and a
              std dev for the random variate.

    """

    assert len(dp) == len(lat)
    assert len(rmw) == len(dp)
    
    X = np.column_stack((dp, dp*dp, lat*lat))
    X = sm.add_constant(X)
    y = np.array(np.log(rmw))
    model = sm.OLS(y, X)
    results = model.fit()
    params = list(results.params)

    r = results.resid
    rf = stats.norm.fit(r, loc=np.mean(r), scale=np.std(r))

    params.append(rf[1])
    return params
