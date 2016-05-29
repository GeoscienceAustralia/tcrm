"""
:mod:`smooth` -- smooth a 2D field using a Gaussian filter
==========================================================

.. module: smooth
    :synopsis: Functions to smooth a 2D field using a Gaussian filter

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

"""

import numpy
from scipy import signal



def gaussKern(size):
    """
    Calculate a normalised Gaussian kernel to apply as a smoothing
    function.

    :param int size: the size of the kernel to use (how many points will be
                used in the smoothing operation).
    :returns: :class:`numpy.ndarray` normalised 2D kernel array for use in
               convolutions
    """
    size = int(size)
    x, y = numpy.mgrid[-size:size + 1, -size:size + 1]
    g = numpy.exp(-(x**2/float(size) + y**2/float(size)))
    return g / g.sum()

def smooth(im, n=15):
    """
    Smooth a 2D array `im` by convolving with a Gaussian kernel of size `n`.

    :param im: Array of values to be smoothed
    :type  im: :class:`numpy.ndarray`
    :param int n: Number of points to include in the smoothing.

    :returns: smoothed array (same dimensions as the input array)

    """
    g = gaussKern(n)
    improc = signal.convolve2d(im, g, mode='same', boundary='symm')
    return improc
