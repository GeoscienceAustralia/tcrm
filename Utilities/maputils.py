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

Title: maputils.py

Author: Craig Arthur
Email: craig.arthur@ga.gov.au
CreationDate: 2006-11-03

Description: Contains mapping functions that operate on arrays
Supercedes some of the functions lying around in some unusual places
(like DataProcess)

Members:
- latLon2Azi(lat, lon, ieast) :
    Calculate the bearing (clockwise positive from true north) and
    distance between consecutive lat/lon pairs.
- bear2LonLat(bearing, distance, lon, lat) :
    Calculate the longitude and latitude of a new point given
    the bearing and distance from some origin lon/lat position
- latLon2XY(xr, yr,lat, lon) :
    Calculate the distance between consecutive lat/lon points - returns
    two vectors (x,y) containing the northward and eastward distances.
    Use xy2r to convert (x,y) to a distance.
- xy2r(x, y) :
    calculate the distance between consecutive members of the arrays
    (x,y)
- gridLatLonDist(clon,clat, lonarray, latarray) :
    Calculate the distance between a grid of lat/lon points and
    a given lat/lon point
- gridLatLonBear(clon,clat, lonarray, latarray) :
    Calculate the bearing of a grid of points described by the arrays
    (lonArray, latArray) from the point (cLon, cLat)
- bearing2theta(angle) :
    Convert from compass bearing to cartesian angle (uses radians)
- theta2bearing(angle) :
    Convert from cartesian angle to compass bearing (uses radians)

Still to be added: - error checking on function inputs
                   - xy2ll: inverse of ll2xy


SeeAlso: (related programs)
Constraints:

ModifiedBy: C. Arthur
ModifiedDate: 2006-11-29
Modification: Added gridLatLonBear(cLon, cLat, lonArray, latArray)
Version: 304

ModifiedBy: C. Arthur
ModifiedDate: 2007-01-18
Modification: Added dist2GC function, check if user requested degrees or
              radians to be returned from latlon2Azi
Version: $Rev: 686 $

$Id: maputils.py 686 2012-03-29 04:24:59Z carthur $
"""

import os, sys, pdb, logging

import numpy
import math
import metutils
import time

# C weave code disabled for now.  The code speeds up the windfield interface module by ~6% but
# does not appear to work on some systems.
#from scipy import weave
#from scipy.weave import converters

__version__ = '$Id: maputils.py 686 2012-03-29 04:24:59Z carthur $'
logger = logging.getLogger('maptools')
#class MapToolError(Exception): pass
#class ArrayMismatch(MapToolError): pass

def xy2r(x, y):
    """
    Given x and y arrays, returns the distance between consecutive
    elements.
    """

    #if len(x) != len(y):
    #   raise ArrayMismatch, "Input array sizes do not match"
    return numpy.sqrt(x**2 + y**2)

def latLon2Azi(lat, lon, ieast=1, azimuth=0, wantdeg=True):
    """
    latLon2Azi: returns the bearing and distance (in km) between
    consecutive members of the array pair (lat,lon)
    USE: [azi,length] = ll2azi(lat,lon,ieast)
          INPUT arrays:
            lat = arry of latitudes, decimal degrees
            lon = array of longitudes, decimal degrees
            ieast   = 1 for longitude increasing towards East
                      -1 for longitude increasing towards West
          OUTPUT arrays:
            azi    = azimuth (+ clockwise from North)
            length = distance (km) between 2 points (local Cartesian)
    NOTES: calls xy2r & ll2xy (hence will bomb at North and South poles)

    Matlab code by Andres Mendez
    Conversion to Python by C. Arthur 2006
    (following the Python code of G. Xu)
    """
    #if len(lat) != len(lon):
    #   raise ArrayMismatch, "Input array sizes do not match"

    xr = 0
    yr = 0

    yn, xe = latLon2XY(xr, yr, lat, lon, ieast, azimuth)
    length = xy2r(yn, xe)
    ### for azimuth calculation use atan2 which returns
    ### angle from -pi to pi. Rules for getting azimuth are:
    ### 1st quadrant (yn > 0 xe > 0):   0    <= angle <=  pi/2
    ###      Rule -1* + pi/2 maps to:  pi/2               0
    ### 2nd quadrant (yn < 0 xe > 0): -pi/2  <= angle <=  0
    ###      Rule -1* + pi/2 maps to:  pi                 pi/2
    ### 3rd quadrant (yn < 0 xe < 0): -pi    <= angle <= -pi/2
    ###      Rule -1* + pi/2 maps to:  3pi/2              pi
    ### 4th quadrant (yn > 0 xe < 0):  pi/2  <= angle <=  pi
    ###      Rule -1* + 5pi/2 maps to: 2pi               3pi/2
    ####################################################################
    angle = numpy.arctan2(yn, xe) # yes, in that order
    bearing = [theta2bearing(i) for i in angle]

    # If bearing in degrees isexpected on return:
    if wantdeg:
        bearing = numpy.array([math.degrees(i) for i in bearing], 'f')

    return bearing, length

def bear2LatLon(bearing, distance, oLon, oLat):
    """
    Calculate the longitude and latitude of a new point from an origin
    point given a distance and bearing
    """
    radius = 6367.0 # Earth radius (km)
    oLon = math.radians(oLon)
    oLat = math.radians(oLat)
    bear = math.radians(bearing)

    nLat = math.asin(numpy.sin(oLat)*numpy.cos(distance/radius) + \
            numpy.cos(oLat)*numpy.sin(distance/radius)*numpy.cos(bear))
    aa = numpy.sin(bear)*numpy.sin(distance/radius)*numpy.cos(oLat)
    bb = numpy.cos(distance/radius) - numpy.sin(oLat)*numpy.sin(nLat)

    nLon = oLon + numpy.arctan2(aa, bb)

    return math.degrees(nLon), math.degrees(nLat)

def latLon2XY(xr, yr, lat, lon, ieast=1, azimuth=0):
    """
    latLon2XY: calculate the cartesian distance between consecutive
               lat,lon points

    USE: [xn,ye] = ll2xy(xr,yr,latr,lonr,azimuth,ieast)
           ieast  =  1 if longitude increases toward the East
                     (normal case)
                  = -1 if longitude increases toward the West
           xr,yr  = scalars, normally 0
           lat    = array, reference point lat
           lon    = array, reference point lon
           azimuth= local coordinate system constructed with origin at
                    latr,lonr, X axis ('North') in direction of azimuth,
                    and Y axis such that X x Y = Z(down)
                    when going from (lat,lon) to (x,y)
                    scalar or array
    OUTPUT: [xn, ye] = array of northward and eastward distances between
                       consecutive points. Use xy2r to convert to a
                       distance between consecutive points.

    NOTES:  * IF INPUT LATR LONR AZIMUTH & XORLAT YORLON ARE ARRAYS,
              THEY MUST BE OF THE SAME DIMENSION!!!
              NO ERROR CHECKING DONE!!!

            * assumes geographical coordinates and azimuth in decimal
              degrees local Cartesian coordinates in km

    BOMB NOTES: will bomb at North and South poles
    Matlab code by AM & SP
    Python code by Craig Arthur 2006
    """

    #if len(lat) != len(lon):
    #   raise ArrayMismatch, "Input array sizes do not match"

    radius = 6367.0 # Earth radius (km)

    lat = numpy.radians(lat)
    lon = numpy.radians(lon)

    # Is azimuth fixed or variable?
    if numpy.size(azimuth) == 1:
        angle = numpy.radians(azimuth)*numpy.ones(lat.size - 1)
    else:
        angle = numpy.radians(azimuth)

    cosazi = numpy.cos(angle)
    sinazi = numpy.sin(angle)

    xntru = xr + radius*(numpy.diff(lat))
    yetru = yr + ieast*radius*(numpy.diff(lon))*numpy.cos(lat[1:])
    xn = xntru*cosazi + yetru*sinazi
    ye = -xntru*sinazi + yetru*cosazi

    return xn, ye

def distGC(lat, lon):
    """
    Distance based on the great circle navigation between pairs of points

    Input:
    lat - a pair of latitude values for the two points 
    lon - a pair of longitude values for the two points

    Output:
    distance (in kilometres) between the two points, based on
    great circle navigation

    Example:

    dist = distGC([-20, -40],[120,190])
    
    """
    radius = 6367.0 # Earth radius (km)

    lat = numpy.radians(lat)
    lon = numpy.radians(lon)

    angular_distance = math.acos(math.sin(lat[0])*math.sin(lat[1]) + \
                       math.cos(lat[0])*math.cos(lat[1])*math.cos(lon[0] - lon[1]))

    return radius*180.0*angular_distance



def gridLatLonDist(cLon, cLat, lonArray, latArray, units=None):
    """
    Generate a grid containing the spherical earth distance
    of the points defined by (lonarray, latarray) from the
    point defined by (clon, clat).
    (lonarray,latarray) and (clon,clat) are in degrees.
    Returns distances in km by default, other units specified by the
    'units' kwarg.

    Based on m_lldist.m by Rich Pawlowicz (rich@ocgy.ubc.ca)
    Modified by Craig Arthur 2006-11-13

    Input:
    cLon - longitude of the point to measure the distance from
    cLat - latitude of the point to measure the distance from
    lonArray - 1-d array of longitude values that will define the grid over
               which distances will be calculated
    latArray - 1-d array of latitude values that will define the grid over
               which distances will be calculated
    units - units of distance to be returned (default is kilometre)

    Output:
    dist - 2-d array containing the distance of the points defined in lonArray 
           and latArray from the point (cLon, cLat)

    Example:
    lonArray = numpy.arange(90.,100.,0.1)
    latArray = numpy.arange(-20.,-10.,0.1)
    dist = gridLatLonDist( 105., -15., lonArray, latArray,'km') 

    """

    # #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    # cLat_cos = 0.0
    # cLat_sin = 0.0
    # lat = empty(len(latArray), 'd')
    # lon = empty(len(lonArray), 'd')

    # dLon_sin = empty(len(lonArray), 'd')
    # dLat_sin = empty(len(latArray), 'd')
    # lat_cos = empty(len(latArray), 'd')

    # dist = empty([len(latArray), len(lonArray)], 'd')

    # code = """
        # #include <math.h>

        # double radius = 6367.0;
        # double toRads = 0.017453292519943295;

        # double cLon_ = cLon;
        # double cLat_ = cLat;

        # cLon_ = cLon_*toRads;
        # cLat_ = cLat_*toRads;

        # cLat_cos = cos(cLat_);

        # for (int i = 0; i < NlonArray[0]; ++i)
        # {
            # lon(i) = lonArray(i)*toRads;
            # double dLon = (lon(i) - cLon_)/2.0;
            # dLon_sin(i) = sin(dLon);
        # }

        # for (int i = 0; i < NlatArray[0]; ++i)
        # {
            # lat(i) = latArray(i)*toRads;
            # lat_cos(i) = cos(lat(i));

            # double dLat = (lat(i) - cLat_)/2.0;
            # dLat_sin(i) = sin(dLat);
        # }

        # for (int j = 0; j < NlatArray[0]; ++j)
        # {
            # for (int i = 0; i < NlonArray[0]; ++i)
            # {
                 # double a = pow(dLat_sin(j), 2) + \
                            # cLat_cos*lat_cos(j)*pow(dLon_sin(i), 2);
                 # double c = 2.0*atan2(sqrt(fabs(a)), sqrt(1 - a));

                 # dist(j, i) = radius*c;
            # }
        # }
    # """
    # err = weave.inline(code,
                       # ['cLon', 'cLat', 'lonArray', 'latArray', 'lat', 'lon',
                        # 'dLon_sin', 'dLat_sin', 'lat_cos', 'dist', 'cLat_cos'],
                       # type_converters=converters.blitz,
                       # compiler = 'gcc')
    # #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    radius = 6367.0

    lat = numpy.radians(latArray)
    lon = numpy.radians(lonArray)

    cLon = math.radians(cLon)
    cLat = math.radians(cLat)
    lon_, lat_ = numpy.meshgrid(lon, lat)

    dLon= lon_ - cLon
    dLat= lat_ - cLat

    a = numpy.square(numpy.sin(dLat/2.0)) + numpy.cos(cLat)*numpy.cos(lat_)*numpy.square(numpy.sin(dLon/2.0))
    c = 2.0*numpy.arctan2(numpy.sqrt(numpy.absolute(a)),numpy. sqrt(1 - a))
    dist = radius*c

    dist = metutils.convert(dist, "km", units)

    return dist

def gridLatLonBear(cLon, cLat, lonArray, latArray):
    """
    Generate a grid containing the bearing of the points defined by
    (lonArray,latArray) from the point defined by (cLon,cLat).
    (lonArray,latArray) and (cLon,cLat) are in degrees.
    Returns bearing in radians.

    Input:
    cLon - longitude of the point to measure the bearing from
    cLat - latitude of the point to measure the bearing from
    lonArray - 1-d array of longitude values that will define the grid over
               which bearings will be calculated
    latArray - 1-d array of latitude values that will define the grid over
               which bearingss will be calculated
    
    Output:
    bear - 2-d array containing the bearing (direction) of the points defined in 
           lonArray and latArray from the point (cLon, cLat)

    Example:
    lonArray = numpy.arange(90.,100.,0.1)
    latArray = numpy.arange(-20.,-10.,0.1)
    bear = gridLatLonBear( 105., -15., lonArray, latArray )
    """

    # #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    # lat = empty(len(latArray), 'd')
    # lon = empty(len(lonArray), 'd')

    # dLon_sin = empty(len(lonArray), 'd')
    # dLon_cos = empty(len(lonArray), 'd')
    # lat_sin = empty(len(latArray), 'd')
    # lat_cos = empty(len(latArray), 'd')

    # bearing = empty([len(latArray), len(lonArray)], 'd')

    # code = """
        # #include <math.h>

        # double toRads = 0.017453292519943295;

        # double cLon_ = cLon;
        # double cLat_ = cLat;

        # cLon_ = cLon_*toRads;
        # cLat_ = cLat_*toRads;

        # double cLat_cos = cos(cLat_);
        # double cLat_sin = sin(cLat_);

        # for (int i = 0; i < NlonArray[0]; ++i)
        # {
            # lon(i) = lonArray(i)*toRads;
            # double dLon = lon(i) - cLon_;
            # dLon_sin(i) = sin(dLon);
            # dLon_cos(i) = cos(dLon);
        # }

        # for (int i = 0; i < NlatArray[0]; ++i)
        # {
            # lat(i) = latArray(i)*toRads;
            # lat_sin(i) = sin(lat(i));
            # lat_cos(i) = cos(lat(i));
        # }

        # for (int j = 0; j < NlatArray[0]; ++j)
        # {
            # for (int i = 0; i < NlonArray[0]; ++i)
            # {
                # double alpha = dLon_sin(i)*lat_cos(j);
                # double beta = (cLat_cos*lat_sin(j)) - (cLat_sin*lat_cos(j)*dLon_cos(i));

                # bearing(j, i) = atan2(alpha, beta);
            # }
        # }
    # """
    # err = weave.inline(code,
                       # ['cLon', 'cLat', 'lonArray', 'latArray', 'lat', 'lon',
                        # 'dLon_sin', 'dLon_cos', 'lat_sin', 'lat_cos', 'bearing'],
                       # type_converters=converters.blitz,
                       # compiler = 'gcc')
    # #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    lat = numpy.radians(latArray)
    lon = numpy.radians(lonArray)

    cLon = math.radians(cLon)
    cLat = math.radians(cLat)
    lon_, lat_ = numpy.meshgrid(lon, lat)

    dLon= lon_ - cLon
    #dLat= lat_ - cLat

    alpha = numpy.sin(dLon)*numpy.cos(lat_)
    beta = numpy.cos(cLat)*numpy.sin(lat_) - \
           numpy.sin(cLat)*numpy.cos(lat_)*numpy.cos(dLon)

    bearing = numpy.arctan2(alpha,beta)

    return bearing

def bearing2theta(bearing):
    """
    Converts bearing in azimuth coordinate system
    into theta in cartesian coordinate system
    Assumes -2*pi <= bearing <= 2*pi
    """
    theta = numpy.pi/2. - bearing
    theta = numpy.mod(theta, 2.*numpy.pi)
    
    return theta

def theta2bearing(theta):
    """
    Converts a cartesian angle (in radians) to an azimuthal
    bearing (in radians).
    Assumes -2*pi <= theta <= 2*pi
    """
    bearing = 2.*numpy.pi - (theta - numpy.pi/2.)
    bearing = numpy.mod(bearing, 2.*numpy.pi)
    
    return bearing

def makeGrid(cLon, cLat, margin=2, resolution=0.01, minLon=None, maxLon=None,
             minLat=None, maxLat=None):
    """
    Generate a grid of the distance and angle of a grid of points
    surrounding a storm centre given the location of the storm.
    The grid margin and grid size can be set in configuration files.
    xMargin, yMargin and gridSize are in degrees
    """
    if type(cLon)==list or type(cLat)==list or type(cLon)==numpy.ndarray or type(cLat)==numpy.ndarray:
        raise TypeError, "Input values must be scalar values"

    gridSize = int(resolution*1000)

    if minLon:
        minLon_ = int(1000*(minLon)) - int(1000*margin)
    else:
        minLon_ = int(1000*(cLon)) - int(1000*margin)
    if maxLon:
        maxLon_ = int(1000*(maxLon)) + int(1000*margin) + 1
    else:
        maxLon_ = int(1000*(cLon)) + int(1000*margin) + 1
    if minLat:
        minLat_ = int(1000*(minLat)) - int(1000*margin)
    else:
        minLat_= int(1000*(cLat)) - int(1000*margin)
    if maxLat:
        maxLat_ = int(1000*(maxLat)) + int(1000*margin) + 1
    else:
        maxLat_= int(1000*(cLat)) + int(1000*margin) + 1

    xGrid = numpy.array(numpy.arange(minLon_, maxLon_, gridSize), dtype=int)
    yGrid = numpy.array(numpy.arange(minLat_, maxLat_, gridSize), dtype=int)

    R = gridLatLonDist(cLon, cLat, xGrid/1000., yGrid/1000.)
    numpy.putmask(R, R==0, 1e-30)
    theta = numpy.pi/2. - gridLatLonBear(cLon, cLat, xGrid/1000., yGrid/1000.)

    return R, theta

def makeGridDomain(cLon, cLat, minLon, maxLon, minLat, maxLat, margin=2, resolution=0.01):
    """
    Generate a grid of the distance and angle of a grid of points
    surrounding a storm centre given the location of the storm.
    The grid margin and grid size can be set in configuration files.
    xMargin, yMargin and gridSize are in degrees


    """
    if (type(cLon)==list or type(cLat)==list or 
        type(cLon)==numpy.ndarray or type(cLat)==numpy.ndarray):
        raise TypeError, "Input values must be scalar values"
    gridSize = int(resolution*1000)
    minLon_ = int(1000*(minLon)) - int(1000*margin)
    maxLon_ = int(1000*(maxLon)) + int(1000*margin) + 1
    minLat_ = int(1000*(minLat)) - int(1000*margin)
    maxLat_ = int(1000*(maxLat)) + int(1000*margin) + 1

    xGrid = numpy.array(numpy.arange(minLon_, maxLon_, gridSize), dtype=int)
    yGrid = numpy.array(numpy.arange(minLat_, maxLat_, gridSize), dtype=int)

    R = gridLatLonDist(cLon, cLat, xGrid/1000., yGrid/1000.)
    numpy.putmask(R, R==0, 1e-30)
    theta = numpy.pi/2. - gridLatLonBear(cLon, cLat, xGrid/1000., yGrid/1000.)
    return R, theta

def meshLatLon(cLon, cLat, margin=2, resolution=0.01):
    """
    Create a meshgrid of the lon/lat grid.
    """
    if (type(cLon)==list or type(cLat)==list or 
        type(cLon)==numpy.ndarray or type(cLat)==numpy.ndarray):
        raise TypeError, "Input values must be scalar values"
    gridSize = int(1000*resolution)

    minLon = int(1000*(cLon-margin))
    maxLon = int(1000*(cLon+margin))+gridSize
    minLat = int(1000*(cLat-margin))
    maxLat = int(1000*(cLat+margin))+gridSize

    xx = numpy.array(numpy.arange(minLon, maxLon, gridSize))
    yy = numpy.array(numpy.arange(minLat, maxLat, gridSize))

    xGrid,yGrid = numpy.meshgrid(xx, yy)
    return xGrid/1000., yGrid/1000.

def meshLatLonDomain(minLon, maxLon, minLat, maxLat, margin=2, resolution=0.01):
    """
    Create a meshgrid of the lon/lat grid.
    """
    gridSize = int(1000*resolution)

    minLon_ = int(1000*(minLon-margin))
    maxLon_ = int(1000*(maxLon+margin))+gridSize
    minLat_ = int(1000*(minLat-margin))
    maxLat_ = int(1000*(maxLat+margin))+gridSize

    xx = numpy.array(numpy.arange(minLon_, maxLon_, gridSize))
    yy = numpy.array(numpy.arange(minLat_, maxLat_, gridSize))

    xGrid,yGrid = numpy.meshgrid(xx, yy)
    return xGrid/1000., yGrid/1000.

def dist2GC(cLon1, cLat1, cLon2, cLat2, lonArray, latArray, units="km"):
    """
    Calculate the distance between an array of points and the great circle
    joining two (other) points.
    All input values are in degrees.
    By default returns distance in km, other units specified by the
    'units' kwarg.

    Based on a cross-track error formulation from:
    http://williams.best.vwh.net/avform.htm#XTE
    """

    # Calculate distance and bearing from first point to array of points:
    dist_ = gridLatLonDist(cLon1, cLat1, lonArray, latArray, units="rad")
    bear_ = gridLatLonBear(cLon1, cLat1, lonArray, latArray)

    #bearing of the cyclone:
    cyc_bear_ = latLon2Azi([cLon1, cLon2], [cLat1, cLat2])

    dist2GC_ = numpy.asin(sin(dist_)*numpy.sin(bear_-cyc_bear_))

    dist2GC = metutils.convert(dist2GC_, "rad", units)
    return dist2GC

def coriolis(lat):
    """Calculate the Coriolis factor
    Calculate the Coriolis factor (f) for a given latitude (degrees).
    If a list is passed, return a list, else return a single value.
    """
    omega = 2*numpy.pi/24./3600.
    f = 2*omega*numpy.sin(numpy.radians(lat))

    return f

def find_index( array, value ):
    """
    Find the index of 'array' with a value closest to 'value'
    Input:
    array - numpy array of data values
    value - a value to search the array for (or find the index of the nearest value to 'value'
    
    Output:
    idx - integer representing the index of 'array' that most closely matches 'value'
    
    Example:
    idx = find_index( numpy.arange(0.,100.,0.5), 15.25 )
    """
    if type( value )==numpy.ndarray or type( value )==list:
        raise ValueError, "Value cannot be an array"

    if ( value > array.max() ): 
        # Value is above the largest value in the array - return the last index:
        return len(array) - 1
    elif ( value < array.min() ):
        # Value is below minimum value in the array - return the first index:
        return 0
    else:
        # argmin gives us the index corresponding to the minimum value of the array.
        idx = ( abs( array - value ) ).argmin( )
        return idx

def find_nearest( array, value ):
    """
    Find the closest value in 'array' to 'value'
    
    Input:
    array - numpy array of data values
    value - a value to search the array for (or find the index of the nearest value to 'value'
    
    Output:
    array[idx], where idx is the index of array that corresponds to the value closest to 'value'

    Example:
    n = find_nearest( numpy.arange(0,100.,0.5), 15.25 )
    """
    if type( value )==numpy.ndarray or type( value )==list:
        raise ValueError, "Value cannot be an array"
    idx = find_index( array, value )

    try:
        v = array[idx]
    except IndexError:
        raise
    else:
        return v
