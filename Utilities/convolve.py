"""
:mod:`convolve` -- convolve two 2-d fields
==========================================

.. module:: convolve
    :synopsis: Convolve two 2-d fields to apply a smoother (e.g. directional
               filtering, gaussian smoother).

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import numpy
from scipy import signal

def getKernel(d,m,r=25.,h=5.):
    """
    Define an appropriate kernel for smoothing the data

    :param str d: Wind direction - one of ['N','NE','E','SE','S','SW','W','NW']
    :param str m: Multiplier type - either 'terrain' or 'shield'
    :param float r: Resolution of the dataset (in metres) (default is 25)
    :param float h: Nominal height of buildings, for use in evaluating the
                    shielding multiplier (default is 5)

    :returns: kernel to apply to the raw multiplier values to generate
              directional multiplier values.
    :rtype: 2-d :class:`numpy.ndarray`

    """

    if d not in ['N','NE','E','SE','S','SW','W','NW']:
        raise ValueError, "No valid direction provided"

    if m=="terrain":
        # We assume that the boundary layer develops over a 1 km distance:
        n = int(1000./r)
        g = numpy.zeros((2*n+1,2*n+1))
        x = numpy.arange(-n,n+1)
        y = numpy.arange(n,-1*n-1,-1)
        xx,yy = numpy.meshgrid(x,numpy.transpose(y))
        rr = r*numpy.sqrt(xx**2+yy**2)
        bear = 90.-(180./numpy.pi)*numpy.arctan2(yy,xx)
        i = numpy.where(bear>180.)
        bear[i] -= 360.

        if d=="S":
            ii = numpy.where((bear==0.) & (rr<=1000.))
        elif d=="SW":
            ii = numpy.where((bear==45.) & (rr<=1000.))
        elif d=="W":
            ii = numpy.where((bear==90.) & (rr<=1000.))
        elif d=="NW":
            ii = numpy.where((bear==135.) & (rr<=1000.))
        elif d=="N":
            ii = numpy.where(((bear==180.) | (bear==-180.)) & (rr<=1000.))
        elif d=="NE":
            ii = numpy.where((bear==-135.) & (rr<=1000.))
        elif d=="E":
            ii = numpy.where((bear==90.) & (rr<=1000.))
        elif d=="SE":
            ii = numpy.where((bear==-45.) & (rr<=1000.))

    elif m=="shield":
        # AS/NZS 1170.2 specifies the shielding distance at 20 times the
        # nominal height of the buildings - assume a 5 m height for
        # residential buildings, hence a 100 m radius.
        n = int(20.*h/r)
        g = numpy.zeros((2*n+1,2*n+1))
        x = numpy.arange(-n,n+1)
        y = numpy.arange(n,-1*n-1,-1)
        xx,yy = numpy.meshgrid(x,numpy.transpose(y))
        rr = r*numpy.sqrt(xx**2+yy**2)
        bear = 90.-(180./numpy.pi)*numpy.arctan2(yy,xx)
        i = numpy.where(bear>180.)
        bear[i] -= 360.
        if d=="S":
            ii = numpy.where((bear>=-22.5) & (bear<=22.5) & (rr<=20.*h))
        elif d=="SW":
            ii = numpy.where((bear>=22.5) & (bear<=67.5) & (rr<=20.*h))
        elif d=="W":
            ii = numpy.where((bear>=67.5) & (bear<=112.5) & (rr<=20.*h))
        elif d=="NW":
            ii = numpy.where((bear>=112.5) & (bear<=157.5) & (rr<=20.*h))
        elif d=="N":
            ii = numpy.where((bear>=157.5) | (bear<=-157.5) & (rr<=20.*h))
        elif d=="NE":
            ii = numpy.where((bear>=-157.5) & (bear<=-112.5) & (rr<=20.*h))
        elif d=="E":
            ii = numpy.where((bear>=-112.5) & (bear<=-67.5) & (rr<=20.*h))
        elif d=="SE":
            ii = numpy.where((bear>=-67.5) & (bear<=-22.5) & (rr<=20.*h))

    g[ii] = 1.
    return g/g.sum()



def convolve(im, direction, m="terrain", res=25.,height=5.):
    """
    Smooth a 2D array im by convolving with a kernel of size n.

    :param im: 2-d array of values to be smoothed.
    :param str dir: One of 'N','NE','E','SE','S','SW','W','NW' to define the
                    direction of the site multiplier to evaluate.
    :param str m: Model type = either "terrain" or "shield".
    :param float res: Resolution of the input grid dataset.
    :param float height: nominal height of the buildings to be used in evaluating
                         the shielding multiplier.
    :type  im: :class:`numpy.ndarray`

    :returns: 2-d array of convolved data (convolved with the appropriate kernel)
              The output array is the same size as the input array, with
              boundary values set to be filled to a value of 1.0
    :rtype: :class:`numpy.ndarray`

    """
    g = getKernel(direction,m,res,height)
    improc = signal.convolve2d(im, g, mode='same', boundary='fill', fillvalue=1.0)
    return(improc)
