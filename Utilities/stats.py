"""
:mod:`stats` -- helper functions for statistical methods
========================================================

.. module:: stats
    :synopsis: Miscellaneous tools required for statistics-related classes.

.. moduleauthor:: Geoff Xu <geoff.xu@ga.gov.au>

"""

import sys
import logging

import math
import numpy as np

from .grid import grdRead

logger = logging.getLogger()


"""
Functions:

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

    h = np.abs(x[1] - x[0])
    y = h * y
    cy = y.cumsum()
    if cy[-1] == 0:
        return cy
    # normalize cy
    cy = cy / cy[-1]
    return cy


def cdf2d(x, y, z):
    """
    2D Cumulative Density Function extracted from cdf2d.m
    Assumes the grid is uniformly defined, i.e. dx and dy are
    constant.
    """
    if np.size(x) < 2 or np.size(y) < 2:
        logger.critical("X or Y grids are not arrays")
        raise TypeError("X or Y grids are not arrays")
    grid_area = np.abs(x[1] - x[0]) * np.abs(y[1] - y[0])
    grid_volume = grid_area * z

    cz = np.zeros([x.size, y.size], 'd')
    cz[:, 0] = (grid_volume[:, 0]).cumsum()
    cz[0, :] = (grid_volume[0, :]).cumsum()

    for i in range(1, len(x)):
        for j in range(1, len(y)):
            cz[i, j] = cz[i - 1, j] + cz[i, j - 1] - \
                cz[i - 1, j - 1] + grid_volume[i, j]

    if cz[-1, -1] == 0:
        return cz
    # normalize cy
    cz = cz / cz[-1, -1]
    return cz


def getCellNum(lon, lat, gridLimit, gridSpace):
    """
    Return the cell number given longitude and latititude
    """
    lon = int(math.floor(lon))
    lat = int(math.ceil(lat))

    if (lon < gridLimit['xMin'] or lon >= gridLimit['xMax'] or
            lat <= gridLimit['yMin'] or lat > gridLimit['yMax']):
        raise ValueError('Invalid input on cellNum: cell number is out of range')

    j = abs((abs(lon) - abs(gridLimit['xMin']))) // abs(gridSpace['x'])
    i = abs((abs(lat) - abs(gridLimit['yMax']))) // abs(gridSpace['y'])

    return int(i * abs((gridLimit['xMax'] - gridLimit['xMin']) // gridSpace['x']) + j)


def getCellLonLat(cellNum, gridLimit, gridSpace):
    """
    Return the lon/lat of a given cell, based on gridLimit and gridSpace
    """
    if cellNum < 0:
        raise IndexError('Index is negative')

    lat = np.arange(gridLimit['yMax'], gridLimit['yMin'], -gridSpace['y'])
    lon = np.arange(gridLimit['xMin'], gridLimit['xMax'], gridSpace['x'])
    indLat = cellNum // lon.size
    indLon = cellNum % lon.size
    return lon[indLon], lat[indLat]


def validCellNum(cellNum, gridLimit, gridSpace):
    """
    Checks whether the given cell number is valid
    """
    if cellNum < 0:
        return False

    latCells = np.size(
        np.arange(gridLimit['yMax'], gridLimit['yMin'], -gridSpace['y']))
    lonCells = np.size(
        np.arange(gridLimit['xMin'], gridLimit['xMax'], gridSpace['x']))

    numCells = latCells * lonCells - 1
    if cellNum > numCells:
        return False
    else:
        return True


def maxCellNum(gridLimit, gridSpace):
    """
    Get maximum cell number based on grid and grid space
    """
    latCells = np.size(np.arange(gridLimit['yMax'],
                                 gridLimit['yMin'],
                                 -gridSpace['y']))
    lonCells = np.size(np.arange(gridLimit['xMin'],
                                 gridLimit['xMax'],
                                 gridSpace['x']))

    return latCells * lonCells - 1


def getOccurence(occurList, indList):
    """
    Returns an array of indices corresponding to cyclone observations that
    """
    return np.array([occurList[ind] for ind in indList])


def statMaxRange(minval, maxval, step):
    """
    Returns the maximum value in a range
    used for arranging the cells to be evenly spaced
    """
    if maxval < minval:
        raise ValueError('Invalid minval maxval input: minval cannot be greater than maxval')
    if step <= 0:
        raise ValueError('Invalid step input: step cannot be 0 or negative number')

    ran = np.arange(minval, maxval + step, step)
    return ran[-1]


def statMinRange(minval, maxval, step):
    """
    Returns the minimum value in a range
    Used for arranging the cells to be evenly spaced
    """
    if maxval < minval:
        raise ValueError('Invalid minval maxval input: minval cannot be greater than maxval')
    if step <= 0:
        raise ValueError('Invalid step input: step cannot be 0 or negative number')
    ran = np.arange(maxval, minval - step, -step)
    return ran[-1]


def rMaxDist(mean, sig, maxrad=200.):
    """rMaxDist(mean, sig, maxrad=200.)
    Based on the logarithmic distribution reported by Willoughby & Rahn (2004)
    """
    x = np.arange(1., maxrad, 1.)
    mu = np.log(mean) - (0.5 * sig ** 2)
    pdf = np.exp(-((np.log(x) - mu) ** 2) / (2 * sig ** 2)) / \
          (x * sig * np.sqrt(2 * np.pi))
    cy = cdf(x, pdf)
    return x, pdf, cy


def circmean(samples, high=2 * np.pi, low=0):
    """
    Compute the circular mean for samples assumed to be in the range
    [low to high]
    """
    ang = (samples - low) * 2 * np.pi / (high - low)
    res = np.angle(np.mean(np.exp(1j * ang)))
    if res < 0:
        res = res + 2 * np.pi
    return res * (high - low) / 2.0 / np.pi + low


def circvar(samples, high=2 * np.pi, low=0):
    """
    Compute the circular variance for samples assumed to be in the range
    [low to high]
    """
    ang = (samples - low) * 2 * np.pi / (high - low)
    res = np.mean(np.exp(1j * ang))
    V = 1 - np.abs(res)
    return ((high - low) / 2.0 / np.pi) ** 2 * V


def circstd(samples, high=2 * np.pi, low=0):
    """
    Compute the circular standard deviation for samples assumed to be
    in the range [low to high]
    """
    ang = (samples - low) * 2 * np.pi / (high - low)
    res = np.mean(np.exp(1j * ang))
    V = 1 - np.abs(res)
    return ((high - low) / 2.0 / np.pi) * np.sqrt(V)


def statRemoveNum(a, Num=sys.maxsize):
    """
    Remove all elements in an array for which value is Num
    """
    if np.shape(a) == ():
        raise ValueError("Input array must be a 1-d array")
    tmp = a.compress(a != Num)
    return tmp.compress(tmp < sys.maxsize)


def statCellFraction(gridLimit, gridSpace, valueFile):
    """
    Calculate the fractional value of each grid cell, based on the
    values stored in valueFile.
    :param dict gridLimit: Dictionary of bounds of the grid.
    :param dict gridSpace: Resolution of the grid to calculate values.
    :param str valueFile: Path to the ascii grid file containing values to sample.

    :returns: :class:`numpy.ndarray` of fractional values, with length equal to the number
              of cells

    Notes: Still need to include bounds checking to ensure the valueFile
    data actually covers the gridLimits.
    """
    gLon, gLat, gData = grdRead(valueFile)
    nCells = maxCellNum(gridLimit, gridSpace) + 1
    output = np.zeros(nCells)
    for cellNum in range(nCells):
        cellLon, cellLat = getCellLonLat(cellNum, gridLimit, gridSpace)
        wLon = cellLon
        eLon = cellLon + gridSpace['x']
        nLat = cellLat
        sLat = cellLat - gridSpace['y']

        ii = np.where((gLon <= eLon) & (gLon >= wLon))
        jj = np.where((gLat <= nLat) & (gLat >= sLat))
        cellValues = gData[np.meshgrid(jj[0], ii[0])]

        if abs(cellValues).max() == 0:
            output[cellNum] = np.average(cellValues)
        else:
            output[cellNum] = np.average(cellValues) / abs(cellValues).max()
    return output


def probability(return_period):
    """Return an annual probability given a return period"""
    p = 1.0 - np.exp(-1.0 / return_period)
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

    :param float val: Value being tested.
    :param float minval: Lower bound of range. Must be lower than `maxval`.
    :param float maxval: Upper bound of range. Must be higher than `minval`.
    :param int fuzz: Number of times to expand bounds using `numpy.nextafter`.
    :param boolean inclusive: Set whether endpoints are within the range.

    :returns: True if `val` is between `minval` and `maxval`, false otherwise.

    From http://penandpants.com/category/python/numpy/
    """
    # expand bounds
    for _ in range(fuzz):
        minval = np.nextafter(minval, minval - 1e6)
        maxval = np.nextafter(maxval, maxval + 1e6)

    if inclusive:
        return minval <= value <= maxval

    else:
        return minval < value < maxval

def bandwidth(data):
    """
    Calculate the bandwidth for a kernel density estimation, using the 
    normal reference method. 
    
    
    
    :param data: :class:`numpy.ndarray` of float values
    
    :returns: Float value of the "optimum" bandwidth for the kernel 
              density estimate
    
    """
    if not isinstance(data, np.ndarray):
        raise TypeError("Wrong input type to bandwidth()")
    if len(np.shape(data)) == 1:
        nobs = len(data)
        nvars = 1
    else:
        nobs, nvars = np.shape(data)
    X = np.std(data, axis=0)
    return 1.06 * X * nobs ** (- 1. / (4 + nvars))
    