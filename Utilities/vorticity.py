"""
:mod:`vorticity` -- calculate vorticity of a vector field
=========================================================

.. module:: vorticity
    :synopsis: Calculates relative or absolute vorticity
               of a vector field.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import numpy as np
import metutils

__version__ = '$Id: vorticity.py 686 2012-03-29 04:24:59Z carthur $'

def relative(u, v, lon, lat):
    """
    Calculates the relative vorticity (curl(u,v)) of a wind field using
    a basic centred difference scheme. The centred differencing means
    the returned array is reduced in size by 2 elements in each
    dimension.

    :param u: 2-d array of eastward vector component.
    :param v: 2-d array of northward vector component.
    :param lon: 1-d array of longitudes of grid that defines the vector field.
    :param lat: 1-d array of latitudes of grid that defines the vector field.

    :type u: :class:`numpy.ndarray`
    :type v: :class:`numpy.ndarray`
    :type lon: :class:`numpy.ndarray`
    :type lat: :class:`numpy.ndarray`

    :return: 2-d :class:`numpy.ndarray` of relative vorticity values.

    """
    dx = np.zeros((len(lat), len(lon)-2))
    dy = np.zeros((len(lat)-2, len(lon)))
    du = np.zeros((len(lat)-2, len(lon)))
    dv = np.zeros((len(lat), len(lon)-2))
    zeta = np.zeros((len(lat)-2, len(lon)-2))

    for i in xrange(1, len(lon)-1):
        for j in xrange(0, len(lat)):
            dx[j, i-1] = metutils.convert((lon[i+1] - lon[i-1]) * \
                                          np.cos(np.pi*lat[j]/180.),
                                          "deg", "m")
            dv[j, i-1] = v[i+1, j] - v[i-1, j]

    for i in xrange(0, len(lon)):
        for j in xrange(1, len(lat)-1):
            dy[j-1, i] = metutils.convert((lat[j+1] - lat[j-1]), "deg", "m")
            du[j-1, i] = u[i, j+1] - u[i, j-1]

    for i in xrange(len(lat) - 2):
        for j in xrange(len(lon) - 2):
            zeta[i, j] = dv[i, j]/dx[i, j] - du[i, j]/dy[i, j]
    return zeta

def absolute(u, v, lon, lat):
    """
    Calculates the absolute vorticity (f + curl(u,v)) of a wind field using
    a basic centred difference scheme. The centred differencing means
    the returned array is reduced in size by 2 elements in each
    dimension.

    :param u: 2-d array of eastward vector component.
    :param v: 2-d array of northward vector component.
    :param lon: 1-d array of longitudes of grid that defines the vector field.
    :param lat: 1-d array of latitudes of grid that defines the vector field.

    :type u: :class:`numpy.ndarray`
    :type v: :class:`numpy.ndarray`
    :type lon: :class:`numpy.ndarray`
    :type lat: :class:`numpy.ndarray`

    :return: 2-d :class:`numpy.ndarray` of absolute vorticity values.

    """
    dx = np.zeros((len(lat), len(lon) - 2))
    dy = np.zeros((len(lat) - 2, len(lon)))
    du = np.zeros((len(lat) - 2, len(lon)))
    dv = np.zeros((len(lat), len(lon) - 2))
    zeta = np.zeros((len(lat) - 2, len(lon) - 2))

    for i in xrange(1, len(lon) - 1):
        for j in xrange(0, len(lat)):
            dx[j, i-1] = metutils.convert((lon[i+1] - lon[i-1]) * \
                                          np.cos(np.pi*lat[j]/180.),
                                          "deg", "m")
            dv[j, i-1] = v[i+1, j] - v[i-1, j]

    for i in xrange(0, len(lon)):
        for j in xrange(1, len(lat) - 1):
            dy[j-1, i] = metutils.convert((lat[j+1] - lat[j-1]), "deg", "m")
            du[j-1, i] = u[i, j+1] - u[i, j-1]

    for i in xrange(len(lat) - 2):
        for j in xrange(len(lon) - 2):
            zeta[i, j] = dv[i, j]/dx[i, j] - du[i, j]/dy[i, j] + \
               metutils.coriolis(lat[i+1])

    return zeta
