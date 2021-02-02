"""
:mod:`convolve` -- convolve two 2-d fields
==========================================

.. module:: convolve
    :synopsis: Convolve two 2-d fields to apply a smoother (e.g. directional
               filtering, gaussian smoother).

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import numpy as np
from scipy import signal

def getKernel(d, mtype, res=25., height=5.):
    """
    Define an appropriate kernel for smoothing the data

    :param str d: Wind direction - one of ['N','NE','E','SE','S','SW','W','NW']
    :param str m: Multiplier type - either 'terrain' or 'shield'
    :param float res: Resolution of the dataset (in metres) (default is 25)
    :param float height: Nominal height of buildings, for use in evaluating the
                    shielding multiplier (default is 5)

    :returns: kernel to apply to the raw multiplier values to generate
              directional multiplier values.
    :rtype: 2-d :class:`numpy.ndarray`

    """

    if dir not in ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']:
        raise ValueError("No valid direction provided")

    if mtype == "terrain":
        # We assume that the boundary layer develops over a 1 km distance:
        n = int(1000./res)
        g = np.zeros((2*n+1, 2*n+1))
        x = np.arange(-n, n+1)
        y = np.arange(n, -1*n-1, -1)
        xx, yy = np.meshgrid(x, np.transpose(y))
        rr = res*np.sqrt(xx**2+yy**2)
        bear = 90.-(180./np.pi)*np.arctan2(yy, xx)
        i = np.where(bear > 180.)
        bear[i] -= 360.

        if d == "S":
            ii = np.where((bear == 0.) & (rr <= 1000.))
        elif d == "SW":
            ii = np.where((bear == 45.) & (rr <= 1000.))
        elif d == "W":
            ii = np.where((bear == 90.) & (rr <= 1000.))
        elif d == "NW":
            ii = np.where((bear == 135.) & (rr <= 1000.))
        elif d == "N":
            ii = np.where(((bear == 180.) | (bear == -180.)) & (rr <= 1000.))
        elif d == "NE":
            ii = np.where((bear == -135.) & (rr <= 1000.))
        elif d == "E":
            ii = np.where((bear == 90.) & (rr <= 1000.))
        elif d == "SE":
            ii = np.where((bear == -45.) & (rr <= 1000.))

    elif mtype == "shield":
        # AS/NZS 1170.2 specifies the shielding distance at 20 times the
        # nominal height of the buildings - assume a 5 m height for
        # residential buildings, hence a 100 m radius.
        n = int(20. * height/res)
        g = np.zeros((2*n + 1, 2*n + 1))
        x = np.arange(-n, n + 1)
        y = np.arange(n, -1*n - 1, -1)
        xx, yy = np.meshgrid(x, np.transpose(y))
        rr = res * np.sqrt(xx**2 + yy**2)
        bear = 90. - (180./np.pi)*np.arctan2(yy, xx)
        i = np.where(bear > 180.)
        bear[i] -= 360.
        if d == "S":
            ii = np.where((bear >= -22.5) & (bear <= 22.5) &
                          (rr <= 20.*height))
        elif d == "SW":
            ii = np.where((bear >= 22.5) & (bear <= 67.5) &
                          (rr <= 20.*height))
        elif d == "W":
            ii = np.where((bear >= 67.5) & (bear <= 112.5) &
                          (rr <= 20.*height))
        elif d == "NW":
            ii = np.where((bear >= 112.5) & (bear <= 157.5) &
                          (rr <= 20.*height))
        elif d == "N":
            ii = np.where((bear >= 157.5) | (bear <= -157.5) &
                          (rr <= 20.*height))
        elif d == "NE":
            ii = np.where((bear >= -157.5) & (bear <= -112.5) &
                          (rr <= 20.*height))
        elif d == "E":
            ii = np.where((bear >= -112.5) & (bear <= -67.5) &
                          (rr <= 20.*height))
        elif d == "SE":
            ii = np.where((bear >= -67.5) & (bear <= -22.5) &
                          (rr <= 20.*height))

    g[ii] = 1.
    return g/g.sum()



def convolve(im, direction, mtype="terrain", res=25., height=5.):
    """
    Smooth a 2D array im by convolving with a kernel of size n.

    :param im: 2-d array of values to be smoothed.
    :param str dir: One of 'N','NE','E','SE','S','SW','W','NW' to define the
                    direction of the site multiplier to evaluate.
    :param str mtype: Model type = either "terrain" or "shield".
    :param float res: Resolution of the input grid dataset.
    :param float height: nominal height of the buildings to be used in
                         evaluating the shielding multiplier.
    :type  im: :class:`numpy.ndarray`

    :returns: 2-d array of convolved data (convolved with the appropriate
              kernel). The output array is the same size as the input array,
              with boundary values set to be filled to a value of 1.0
    :rtype: :class:`numpy.ndarray`

    """
    kernel = getKernel(direction, mtype, res, height)
    improc = signal.convolve2d(im, kernel, mode='same', boundary='fill',
                               fillvalue=1.0)
    return improc
