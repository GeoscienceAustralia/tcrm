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


#-----------------------------------------------------------------------
# Imports:
#-----------------------------------------------------------------------
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
from numpy import *
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
    return sqrt(x**2 + y**2)

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
    angle = arctan2(yn, xe) # yes, in that order
    bearing = [theta2bearing(i) for i in angle]

    # If bearing in degrees isexpected on return:
    if wantdeg:
        bearing = array([math.degrees(i) for i in bearing], 'f')

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

    nLat = math.asin(sin(oLat)*cos(distance/radius) + \
            cos(oLat)*sin(distance/radius)*cos(bear))
    aa = sin(bear)*sin(distance/radius)*cos(oLat)
    bb = cos(distance/radius) - sin(oLat)*sin(nLat)

    nLon = oLon + arctan2(aa, bb)

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

    lat = radians(lat)
    lon = radians(lon)

    # Is azimuth fixed or variable?
    if size(azimuth) == 1:
        angle = radians(azimuth)*ones(lat.size - 1)
    else:
        angle = radians(azimuth)

    cosazi = cos(angle)
    sinazi = sin(angle)

    xntru = xr + radius*(diff(lat))
    yetru = yr + ieast*radius*(diff(lon))*cos(lat[1:])
    xn = xntru*cosazi + yetru*sinazi
    ye = -xntru*sinazi + yetru*cosazi

    return xn, ye

def distGC(lat, lon):
    """Distance based on the great circle navigation
    """
    radius = 6367.0 # Earth radius (km)

    lat = radians(lat)
    lon = radians(lon)

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

    lat = radians(latArray)
    lon = radians(lonArray)

    cLon = math.radians(cLon)
    cLat = math.radians(cLat)
    lon_, lat_ = meshgrid(lon, lat)

    dLon= lon_ - cLon
    dLat= lat_ - cLat

    a = square(sin(dLat/2.0)) + cos(cLat)*cos(lat_)*square(sin(dLon/2.0))
    c = 2.0*arctan2(sqrt(absolute(a)), sqrt(1 - a))
    dist = radius*c

    dist = metutils.convert(dist, "km", units)

    return dist

def gridLatLonBear(cLon, cLat, lonArray, latArray):
    """
    Generate a grid containing the bearing of the points defined by
    (lonArray,latArray) from the point defined by (cLon,cLat).
    (lonArray,latArray) and (cLon,cLat) are in degrees.
    Returns bearing in radians.
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

    lat = radians(latArray)
    lon = radians(lonArray)

    cLon = math.radians(cLon)
    cLat = math.radians(cLat)
    lon_, lat_ = meshgrid(lon, lat)

    dLon= lon_ - cLon
    #dLat= lat_ - cLat

    alpha = sin(dLon)*cos(lat_)
    beta = cos(cLat)*sin(lat_) - sin(cLat)*cos(lat_)*cos(dLon)

    bearing = arctan2(alpha,beta)

    return bearing

def bearing2theta(bearing):
    """
    Converts bearing in azimuth coordinate system
    into theta in cartesian coordinate system
    Assumes -2*pi <= bearing <= 2*pi
    """
    theta = pi/2. - bearing
    if (type(bearing) is ndarray) or (type(bearing) is list):
        for i in xrange(len(bearing)):
            if theta[i] > 2.*pi:
                theta[i] -= 2.*pi
            else:
                pass

    else:
        if theta > 2.*pi:
            theta = theta - 2.*pi

    return theta

def theta2bearing(theta):
    """
    Converts a cartesian angle (in radians) to an azimuthal
    bearing (in radians).
    Assumes -2*pi <= theta <= 2*pi
    """
    bearing = 2.*pi - (theta - pi/2.)
    if (type(theta) is ndarray) or (type(theta) is list):
        for i in xrange(len(theta)):
            if bearing[i] >= 2.*pi:
                bearing[i] -= 2.*pi
            elif bearing[i] < 0:
                bearing[i] += 2.*pi
            else:
                pass
    else:
        if bearing >= 2.*pi:
            bearing -= 2.*pi
        elif bearing < 0:
            bearing += 2.*pi
        else:
            pass

    return bearing

def makeGrid(cLon, cLat, margin=2, resolution=0.01):
    """
    Generate a grid of the distance and angle of a grid of points
    surrounding a storm centre given the location of the storm.
    The grid margin and grid size can be set in configuration files.
    xMargin, yMargin and gridSize are in degrees
    """

    gridSize = int(resolution*100)
    minLon = int(100*(cLon)) - int(100*margin)
    maxLon = int(100*(cLon)) + int(100*margin) + 1
    minLat = int(100*(cLat)) - int(100*margin)
    maxLat = int(100*(cLat)) + int(100*margin) + 1

    xGrid = array(arange(minLon, maxLon, gridSize), dtype=int)
    yGrid = array(arange(minLat, maxLat, gridSize), dtype=int)

    R = gridLatLonDist(cLon, cLat, xGrid/100., yGrid/100.)
    putmask(R, R==0, 1e-30)
    theta = pi/2. - gridLatLonBear(cLon, cLat, xGrid/100., yGrid/100.)
    return R, theta

def makeGridDomain(cLon, cLat, minLon, maxLon, minLat, maxLat, margin=2, resolution=0.01):
    """
    Generate a grid of the distance and angle of a grid of points
    surrounding a storm centre given the location of the storm.
    The grid margin and grid size can be set in configuration files.
    xMargin, yMargin and gridSize are in degrees
    """

    gridSize = int(resolution*100)
    minLon_ = int(100*(minLon)) - int(100*margin)
    maxLon_ = int(100*(maxLon)) + int(100*margin) + 1
    minLat_ = int(100*(minLat)) - int(100*margin)
    maxLat_ = int(100*(maxLat)) + int(100*margin) + 1

    xGrid = array(arange(minLon_, maxLon_, gridSize), dtype=int)
    yGrid = array(arange(minLat_, maxLat_, gridSize), dtype=int)

    R = gridLatLonDist(cLon, cLat, xGrid/100., yGrid/100.)
    putmask(R, R==0, 1e-30)
    theta = pi/2. - gridLatLonBear(cLon, cLat, xGrid/100., yGrid/100.)
    return R, theta

def meshLatLon(cLon, cLat, margin=2, resolution=0.01):
    """
    Create a meshgrid of the lon/lat grid.
    """
    gridSize = int(100*resolution)

    minLon = int(100*(cLon-margin))
    maxLon = int(100*(cLon+margin))+gridSize
    minLat = int(100*(cLat-margin))
    maxLat = int(100*(cLat+margin))+gridSize

    xx = array(arange(minLon, maxLon, gridSize))
    yy = array(arange(minLat, maxLat, gridSize))

    xGrid,yGrid = meshgrid(xx, yy)
    return xGrid/100., yGrid/100.

def meshLatLonDomain(minLon, maxLon, minLat, maxLat, margin=2, resolution=0.01):
    """
    Create a meshgrid of the lon/lat grid.
    """
    gridSize = int(100*resolution)

    minLon_ = int(100*(minLon-margin))
    maxLon_ = int(100*(maxLon+margin))+gridSize
    minLat_ = int(100*(minLat-margin))
    maxLat_ = int(100*(maxLat+margin))+gridSize

    xx = array(arange(minLon_, maxLon_, gridSize))
    yy = array(arange(minLat_, maxLat_, gridSize))

    xGrid,yGrid = meshgrid(xx, yy)
    return xGrid/100., yGrid/100.

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

    dist2GC_ = asin(sin(dist_)*sin(bear_-cyc_bear_))

    dist2GC = metutils.convert(dist2GC_, "rad", units)
    return dist2GC

def coriolis(lat):
    """Calculate the Coriolis factor
    Calculate the Coriolis factor (f) for a given latitude (degrees).
    If a list is passed, return a list, else return a single value.
    """
    omega = 2*math.pi/24./3600.
    f = 2*omega*sin(radians(lat))

    return f
