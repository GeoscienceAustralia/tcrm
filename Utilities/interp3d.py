"""
:mod:`interp3d` -- interpolate to a set of points in 3-dimensional space
========================================================================

.. module:: interp3d
    :synopsis: Interpolate data to a set of points in 3-d space.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

Use :func:`scipy.ndimage.interpolation.map-coordinates` to interpolate
data in three dimensions.

"""

import numpy as np
from scipy.ndimage.interpolation import map_coordinates

def interp3d(input_array, coords,
             scale=[360., 180., 365.],
             offset=[0., -90., 0.],
             prefilter=True):
    """
    Wrapper to :func:`scipy.ndimage.interpolation.map_coordinates`, which
    converts coordinates of points to indices that correspond to the
    array.  We assume that one is working with lon, lat, day data
    (i.e. initially designed to work with daily long term mean sea
    level pressure)

    :param input_array: A 3-d array of data at regular intervals,
                        representing the data to be evaluated.
    :param coords: A 3-by-n array of coordinates at which the data
                   in ``input_array`` will be interpolated to.
    :param list scale: A (list of) scale factor(s) that reduces the
                       ``coords`` values to the range of indices in
                       ``input_array``.
    :param list offset: A (list of) offset factor(s) that is
                        subtracted from the ``coords`` values before
                        adjusting the ``scale`` (above).
    :param boolean prefilter: If ``True`` (default), then apply a spline
                              filter before interpolation (necessary for
                              spline interpolation of order > 1).  If
                              ``False``, it is assumed that the input is
                              already filtered. Default is ``True``.

    :returns: 1-d array of values corresponding to the interpolated
              values at the points given in ``coords``.

    Example::

      >>> vals = interp3d(data, coords, scale=[360., 180., 365.],
      ...                 offset=[0., -90., 0.])

    """

    if input_array.ndim != 3:
        raise ValueError('Input array has incorrect shape')
    if coords.shape[0] != 3:
        raise ValueError('Coordinates of points must be 3-d')

    dims = input_array.shape
    indices = [d*(c - o) / s for d, c, o, s in
               zip(dims, coords, offset, scale)]

    values = map_coordinates(input_array, indices, mode='wrap',
                             prefilter=prefilter)
    dtype = input_array.dtype
    return np.array(values, dtype)

def _interp(data, coords, scale=[360., 180.], offset=[0., -90.]):
    """
    Wrapper to scipy.ndimage.interpolation.map_coordinates, which converts
    coordinates of points to indices that correspond to the array.
    We assume that one is working with lon, latdata (i.e. initially
    designed to work with daily long term mean sea level pressure)

    Input:
    input_array - a 2-d array of data at regular intervals, representing the
                  data to be evaluated
    coords - a 2xn array of coordinates at which the data in input_array
             will be interpolated
    scale - a scale factor that reduces the coords values to the range of
            indices in input_array
    offset - an offset factor that is subtracted from the coords values
             before adjusting the scale (above)

    Output:
    1-d array of values corresponding to the interpolated values
    at the points given in 'coords'

    Example: vals = interp2d( data, coords, scale=[360., 180.],
                              offset=[0.,-90.] )
    """

    if data.ndim != np.asarray(coords).ndim:
        raise ValueError('Input array and coordinates do not have matching dimensions')

    dims = np.array(data.shape)

    indices = [d*(c - o) / s for d, c, o, s in 
               zip(dims, coords, offset, scale)]

    values = map_coordinates(data, indices, mode='wrap')

    return values
