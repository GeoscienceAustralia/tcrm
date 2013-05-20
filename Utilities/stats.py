#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Title: stats.py - helper functions for statistical methods
Author: Geoff Xu, geoff.xu@ga.gov.au
CreationDate: 2005-12-23
Description: Miscellaneous tools required for statistics-related classes.
SeeAlso:
Constraints:
Version: $Rev: 686 $

ModifiedBy: Geoff Xu, geoff.xu@ga.gov.au
ModifiedDate: 2006-01-25
Modifications:

ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2006-11-07
Modification: Changed to a module in utils

ModifiedBy: N. Habili, nariman.habili@ga.gov.au
ModifiedDate: 2006-11-30
Modification: Upgrade to numpy

ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2006-12-07
Modification: Added getCellNum, getMinRange, getMaxRange,
            getOccurence (last three deprecated)

ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2006-12-12
Modification: Added getCellLonLat, validCellNum

ModifiedBy: N. Habili, nariman.habili@ga.gov.au
ModifiedDate: 2006-12-22
Modification: Added maxCellNum
          Modified getCellNum. Index out of range will now assert.

Version: $Rev: 686 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-04-24 4:07:PM
Modification: Added statCellFraction

$Id: stats.py 686 2012-03-29 04:24:59Z carthur $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
import math
from scipy import array, arange, size, zeros
from numpy import *

from grid import grdRead

__version__ = '$Id: stats.py 686 2012-03-29 04:24:59Z carthur $'

logger = logging.getLogger()


"""
Functions:
-------
cdf(x,y) : 1D array of float
    Cumulative Density Function extracted from cdf_lin.m
cdf2d(x,y,z): 2D array of float
    2D Cumulative Density Function extracted from cdf2d.m
getCellNum(lon, lat, gridLimit, gridSpace): int
    Determine the cell number based on the lat/lon, the grid bounds
    and the grid spacing.
getCellLonLat(cellNum, gridLimit, gridSpace): 2D float
    Determine the lat/lon  of the northwestern corner of
    cellNum
validCellNum(cellNum, gridLimit, gridSpace): boolean
    Determine whether the given cell number can exist within
    the bounds of the region defined by gridSpace and gridLimit
maxCellNum(gridLimit, gridSpace): int
    Determine maximum cell number based on grid limits and spacing.
"""
def cdf(x, y):
    """
    Cumulative Density Function extracted from cdf_lin.m
    """

    h = abs(x[1] - x[0])
    y = h*y
    cy = y.cumsum()
    if cy[-1] == 0:
        return cy
    #normalize cy
    cy = cy/cy[-1]
    return cy

def cdf2d(x, y, z):
    """
    2D Cumulative Density Function extracted from cdf2d.m
    Assumes the grid is uniformly defined, i.e. dx and dy are
    constant.
    """
    if size(x) < 2 or size(y) < 2:
        logger.critical("X or Y grids are not arrays")
        raise TypeError, "X or Y grids are not arrays"
    grid_area = abs(x[1] - x[0])*abs(y[1] - y[0])
    grid_volume = grid_area*z

    cz = zeros([x.size, y.size], 'd')
    cz[:, 0] = (grid_volume[:,0]).cumsum()
    cz[0, :] = (grid_volume[0,:]).cumsum()

    for i in xrange(1, len(x)):
        for j in xrange(1, len(y)):
            cz[i,j] = cz[i-1,j] + cz[i,j-1] - cz[i-1,j-1] + grid_volume[i,j]

    if cz[-1,-1] == 0:
        return cz
    #normalize cy
    cz = cz/cz[-1,-1]
    return cz

def getCellNum(lon, lat, gridLimit, gridSpace):
    """
    Return the cell number given longitude and latititude
    """
    lon = int(math.floor(lon))
    lat = int(math.ceil(lat))

    if (lon < gridLimit['xMin'] or lon >= gridLimit['xMax'] or \
        lat <= gridLimit['yMin'] or lat > gridLimit['yMax']):
        raise ValueError, 'Invalid input on cellNum: cell number is out of range'

    j = abs((abs(lon) - abs(gridLimit['xMin'])))/abs(gridSpace['x'])
    i = abs((abs(lat) - abs(gridLimit['yMax'])))/abs(gridSpace['y'])

    return int(i*abs((gridLimit['xMax'] - gridLimit['xMin'])/gridSpace['x']) + j)

def getCellLonLat(cellNum, gridLimit, gridSpace):
    """
    Return the lon/lat of a given cell, based on gridLimit and gridSpace
    """
    if (cellNum < 0):
        raise IndexError, 'Index is negative'

    lat = arange(gridLimit['yMax'], gridLimit['yMin'], -gridSpace['y'])
    lon = arange(gridLimit['xMin'], gridLimit['xMax'], gridSpace['x'])
    indLat = cellNum/lon.size
    indLon = cellNum%lon.size
    return lon[indLon], lat[indLat]

def validCellNum(cellNum, gridLimit, gridSpace):
    """
    Checks whether the given cell number is valid
    """
    if (cellNum < 0):
        return False

    latCells = size(arange(gridLimit['yMax'], gridLimit['yMin'], -gridSpace['y']))
    lonCells = size(arange(gridLimit['xMin'], gridLimit['xMax'], gridSpace['x']))

    numCells = latCells*lonCells - 1
    if (cellNum > numCells):
        return False
    else:
        return True

def maxCellNum(gridLimit, gridSpace):
    """
    Get maximum cell number based on grid and grid space
    """
    latCells = size(arange(gridLimit['yMax'], gridLimit['yMin'], -gridSpace['y']))
    lonCells = size(arange(gridLimit['xMin'], gridLimit['xMax'], gridSpace['x']))

    return latCells*lonCells - 1

def getOccurence(occurList, indList):
    """
    Returns an array of indices corresponding to cyclone observations that
    """
    return array([occurList[ind] for ind in indList])

def statMaxRange(minval, maxval, step):
    """
    Returns the maximum value in a range
    used for arranging the cells to be evenly spaced
    """
    if maxval < minval:
        raise ValueError, 'Invalid minval maxval input: minval cannot be greater than maxval'
    if step <= 0:
        raise ValueError, 'Invalid step input: step cannot be 0 or negative number'

    ran = arange(minval, maxval + step, step)
    return ran[-1]

def statMinRange(minval, maxval, step):
    """
    Returns the minimum value in a range
    Used for arranging the cells to be evenly spaced
    """
    if maxval < minval:
        raise ValueError, 'Invalid minval maxval input: minval cannot be greater than maxval'
    if step <= 0:
        raise ValueError, 'Invalid step input: step cannot be 0 or negative number'
    ran = arange(maxval, minval - step, -step)
    return ran[-1]

def rMaxDist(mean, sig, maxrad=200.):
    """rMaxDist(mean, sig, maxrad=200.)
    Based on the logarithmic distribution reported by Willoughby & Rahn (2004)
    """
    x = arange(1., maxrad, 1.)
    mu = log(mean) - (0.5*sig**2)
    pdf = exp(-((log(x)-mu)**2)/(2*sig**2))/(x*sig*sqrt(2*pi))
    cy = cdf(x, pdf)
    return x, pdf, cy

def circmean(samples, high=2*pi, low=0):
    """
    Compute the circular mean for samples assumed to be in the range
    [low to high]
    """
    ang = (samples - low)*2*pi / (high-low)
    res = angle(mean(exp(1j*ang)))
    if (res < 0):
        res = res + 2*pi
    return res*(high-low)/2.0/pi + low

def circvar(samples, high=2*pi, low=0):
    """
    Compute the circular variance for samples assumed to be in the range
    [low to high]
    """
    ang = (samples - low)*2*pi / (high-low)
    res = mean(exp(1j*ang))
    V = 1-abs(res)
    return ((high-low)/2.0/pi)**2 * V

def circstd(samples, high=2*pi, low=0):
    """
    Compute the circular standard deviation for samples assumed to be
    in the range [low to high]
    """
    ang = (samples - low)*2*pi / (high-low)
    res = mean(exp(1j*ang))
    V = 1-abs(res)
    return ((high-low)/2.0/pi) * sqrt(V)

def statRemoveNum(a, Num=sys.maxint):
    """
    Remove all elements in an array for which value is Num
    """
    if shape(a) == ():
        raise ValueError, "Input array must be a 1-d array"
    tmp = a.compress(a <> Num)
    return tmp.compress(tmp<sys.maxint)

def statCellFraction(gridLimit, gridSpace, valueFile):
    """
    Calculate the fractional value of each grid cell, based on the
    values stored in valueFile.
    Input: gridLimit - dictionary of bounds of the grid
           gridSpace - resolution of the grid to calculate values.
           valueFile - path to the ascii grid file containing values to
                       sample
    Output: array of fractional values, with length equal to the number
            of cells
    Example:
    Notes: Still need to include bounds checking to ensure the valueFile
           data actually covers the gridLimits.
    """
    gLon, gLat, gData = grdRead(valueFile)
    nCells = maxCellNum(gridLimit, gridSpace) + 1
    output = zeros(nCells)
    for cellNum in xrange(nCells):
        cellLon, cellLat = getCellLonLat(cellNum, gridLimit, gridSpace)
        wLon = cellLon
        eLon = cellLon + gridSpace['x']
        nLat = cellLat
        sLat = cellLat - gridSpace['y']

        ii = where((gLon<=eLon) & (gLon>=wLon))
        jj = where((gLat<=nLat) & (gLat>=sLat))
        cellValues = gData[meshgrid(jj[0],ii[0])]

        if abs(cellValues).max() == 0:
            output[cellNum] = average(cellValues)
        else:
            output[cellNum] = average(cellValues)/abs(cellValues).max()
    return output

def probability(return_period):
    """Return an annual probability given a return period"""
    p = 1.0 - exp(-1.0/return_period)
    return p

def between(value, minval, maxval, fuzz=2, inclusive=True):
    """
    Test whether a value is within some range with some fuzziness at the edges
    to allow for floating point noise.

    The fuzziness is implemented by expanding the range at each end `fuzz` steps
    using the numpy.nextafter function. For example, with the inputs
    minval = 1, maxval = 2, and fuzz = 2; the range would be expanded to
    minval = 0.99999999999999978 and maxval = 2.0000000000000009 before doing
    comparisons.

    Parameters
    ----------
    val : float
        Value being tested.

    minval : float
        Lower bound of range. Must be lower than `maxval`.

    maxval : float
        Upper bound of range. Must be higher than `minval`.

    fuzz : int, optional
        Number of times to expand bounds using numpy.nextafter.

    inclusive : bool, optional
        Set whether endpoints are within the range.

    Returns
    -------
    between : bool
        True if `val` is between `minval` and `maxval`, False otherwise.

    From http://penandpants.com/category/python/numpy/
    """
    # expand bounds
    for _ in xrange(fuzz):
        minval = nextafter(minval, minval - 1e6)
        maxval = nextafter(maxval, maxval + 1e6)

    if inclusive:
        return minval <= value <= maxval

    else:
        return minval < value < maxval