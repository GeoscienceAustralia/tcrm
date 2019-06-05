"""
:mod:`maputils` -- mapping functions
====================================

.. module:: maputils
    :synopsis: Mapping functions.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

Contains mapping functions that operate on arrays. Supercedes some of
the functions lying around in some unusual places (like DataProcess).

"""

import logging

import numpy as np
import math
from . import metutils


# C weave code disabled for now.  The code speeds up the windfield interface module by ~6% but
# does not appear to work on some systems.
#from scipy import weave
#from scipy.weave import converters

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def xy2r(x, y):
    """
    Given x and y arrays, returns the distance between consecutive
    elements.

    :param x: x-coordinate of points.
    :param y: y-coordinate of points.
    :type x: :class:`numpy.ndarray`
    :type y: :class:`numpy.ndarray`

    :returns: Distance (in native units) between consecutive points.
    :rtype: :class:`numpy.ndarray`
    """

    #if len(x) != len(y):
    #   raise ArrayMismatch, "Input array sizes do not match"
    return np.sqrt(x**2 + y**2)

def latLon2Azi(lat, lon, ieast=1, azimuth=0, wantdeg=True):
    """
    Returns the bearing and distance (in km) between consecutive
    members of the array pair (lat,lon).

    :param lat: Latitudes of positions.
    :param lon: Longitudes of positions.
    :param int ieast: 1 for longitudes increasing towards East, -1 for
                      longitudes increasing towards West (default 1).
    :param float azimuth: Local coordinate system constructed with origin at
                          latr,lonr, X axis ('North') in direction of azimuth,
                          and Y axis such that X x Y = Z(down)
                          when going from (lat,lon) to (x,y) (default 0).
    :param boolean wantdeg: If ``True`` return bearings as degrees, not radians.

    :returns: azimuth (+ve clockwise from north) and distance (in km).

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
    angle = np.arctan2(yn, xe) # yes, in that order
    bearing = [theta2bearing(i) for i in angle]

    # If bearing in degrees isexpected on return:
    if wantdeg:
        bearing = np.array([math.degrees(i) for i in bearing], 'f')

    return bearing, length

def bear2LatLon(bearing, distance, oLon, oLat):
    """
    Calculate the longitude and latitude of a new point from an origin
    point given a distance and bearing.

    :param bearing: Direction to new position (degrees, +ve clockwise
                    from north).
    :param distance: Distance to new position (km).
    :param oLon: Initial longitude.
    :param oLat: Initial latitude.

    :returns: new longitude and latitude (in degrees)
    """
    radius = 6367.0 # Earth radius (km)
    oLon = math.radians(oLon)
    oLat = math.radians(oLat)
    bear = math.radians(bearing)

    nLat = math.asin(np.sin(oLat) * np.cos(distance / radius) + \
            np.cos(oLat) * np.sin(distance / radius) * np.cos(bear))
    aa = np.sin(bear) * np.sin(distance / radius) * np.cos(oLat)
    bb = np.cos(distance / radius) - np.sin(oLat) * np.sin(nLat)

    nLon = oLon + np.arctan2(aa, bb)

    return math.degrees(nLon), math.degrees(nLat)

def latLon2XY(xr, yr, lat, lon, ieast=1, azimuth=0):
    """
    Calculate the cartesian distance between consecutive lat,lon
    points. Will bomb at North and South Poles. Assumes geographical
    coordinates and azimuth in decimal degrees, local Cartesian
    coordinates in km.

    :param xr: Reference longitude, normally 0.
    :param yr: Reference latitude, normally 0.
    :param lat: Array of latitudes.
    :param lon: Array of longitudes.
    :param int ieast: 1 if longitude increases toward the East
                     (normal case), -1 if longitude increases
                     toward the West.
    :param int azimuth: local coordinate system constructed with
                        origin at latr,lonr, X axis ('North') in
                        direction of azimuth, and Y axis such that X x
                        Y = Z(down) when going from (lat,lon) to (x,y)
                        scalar or array.

    :returns: Array of northward and eastward distances between
              consecutive points. use :func:`xy2r` to convert to a
              distance between consecutive points.
    """

    #if len(lat) != len(lon):
    #   raise ArrayMismatch, "Input array sizes do not match"

    radius = 6367.0 # Earth radius (km)

    lat = np.radians(lat)
    lon = np.radians(lon)

    # Is azimuth fixed or variable?
    if np.size(azimuth) == 1:
        angle = np.radians(azimuth)*np.ones(lat.size - 1)
    else:
        angle = np.radians(azimuth)

    cosazi = np.cos(angle)
    sinazi = np.sin(angle)

    xntru = xr + radius * (np.diff(lat))
    yetru = yr + ieast * radius * (np.diff(lon)) * np.cos(lat[1:])
    xn = xntru * cosazi + yetru * sinazi
    ye = -xntru * sinazi + yetru * cosazi

    return xn, ye

def distGC(lat, lon):
    """
    Distance based on the great circle navigation between pairs of points.

    :param lat: A pair of latitude values for the two points.
    :param lon: A pair of longitude values for the two points.

    :returns: Distance (in kilometres) between the two points, based on
              great circle navigation.

    Example::

        >>> dist = distGC([-20, -40],[120,190])
        6914.42

    """
    radius = 6367.0 # Earth radius (km)

    lat = np.radians(lat)
    lon = np.radians(lon)

    angular_distance = math.acos(math.sin(lat[0]) * math.sin(lat[1]) + \
                       math.cos(lat[0]) * math.cos(lat[1]) * \
                       math.cos(lon[0] - lon[1]))

    return radius*angular_distance



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


    :param float cLon: Longitude of the point to measure the distance from.
    :param float cLat: Latitude of the point to measure the distance from.
    :param lonArray: 1-d array of longitude values that will define the
                     grid over which distances will be calculated.
    :param latArray: 1-d array of latitude values that will define the
                     grid over which distances will be calculated.
    :param str units: Units of distance to be returned (default is kilometre)

    :returns: 2-d array containing the distance of the points defined in
             ``lonArray`` and ``latArray`` from the point
             (``cLon``, ``cLat``).

    Example::

        >>> lonArray = np.arange(90.,100.,0.1)
        >>> latArray = np.arange(-20.,-10.,0.1)
        >>> dist = gridLatLonDist( 105., -15., lonArray, latArray, 'km')

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

    lat = np.radians(latArray)
    lon = np.radians(lonArray)

    cLon = math.radians(cLon)
    cLat = math.radians(cLat)
    lon_, lat_ = np.meshgrid(lon, lat)

    dLon = lon_ - cLon
    dLat = lat_ - cLat

    a = np.square(np.sin(dLat / 2.0)) + \
        np.cos(cLat) * np.cos(lat_) * np.square(np.sin(dLon / 2.0))
    c = 2.0 * np.arctan2(np.sqrt(np.absolute(a)), np.sqrt(1 - a))
    dist = radius * c

    dist = metutils.convert(dist, "km", units)

    return dist

def gridLatLonBear(cLon, cLat, lonArray, latArray):
    """
    Generate a grid containing the bearing of the points defined by
    (lonArray,latArray) from the point defined by (cLon,cLat).
    (lonArray,latArray) and (cLon,cLat) are in degrees.
    Returns bearing in radians.

    :param float cLon: Longitude of the point to measure the distance from.
    :param float cLat: Latitude of the point to measure the distance from.
    :param lonArray: 1-d array of longitude values that will define the
                     grid over which distances will be calculated.
    :param latArray: 1-d array of latitude values that will define the
                     grid over which distances will be calculated.
    :returns: 2-d array containing the bearing (direction) of the points
              defined in ``lonArray`` and ``latArray`` from the point
              (``cLon``, ``cLat``)

    Example::

        >>> from maputils import gridLatLonBear
        >>> import numpy as np
        >>> lonArray = np.arange(90.,100.,0.1)
        >>> latArray = np.arange(-20.,-10.,0.1)
        >>> gridLatLonBear( 105., -15., lonArray, latArray)
        array([[-1.94475949, -1.94659552, -1.94845671, ..., -2.36416927,
                -2.37344337, -2.38290081],
               [-1.93835542, -1.94015859, -1.94198663, ..., -2.35390045,
                -2.36317282, -2.37263233],
               [-1.93192776, -1.93369762, -1.93549204, ..., -2.34343069,
                -2.35269718, -2.36215458],
                ...,
               [-1.29066433, -1.28850464, -1.28632113, ..., -0.84374983,
                -0.83405688, -0.82416555],
               [-1.28446304, -1.28227062, -1.28005406, ..., -0.83332654,
                -0.82361918, -0.813717  ],
               [-1.27828819, -1.27606348, -1.27381433, ..., -0.82310335,
                -0.81338586, -0.80347714]])

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

    lat = np.radians(latArray)
    lon = np.radians(lonArray)

    cLon = math.radians(cLon)
    cLat = math.radians(cLat)
    lon_, lat_ = np.meshgrid(lon, lat)

    dLon = lon_ - cLon
    #dLat= lat_ - cLat

    alpha = np.sin(dLon) * np.cos(lat_)
    beta = np.cos(cLat) * np.sin(lat_) - \
           np.sin(cLat) * np.cos(lat_) * np.cos(dLon)

    bearing = np.arctan2(alpha, beta)

    return bearing

def bearing2theta(bearing):
    """
    Converts bearing in azimuth coordinate system into theta in
    cartesian coordinate system. Assumes -2*pi <= bearing <= 2*pi

    :param bearing: Bearing to convert (in radians) (+ve clockwise
                    from north).

    :returns: Angle in cartesian coordinate system (+ve anticlockwise
              from east).

    """
    theta = np.pi / 2. - bearing
    theta = np.mod(theta, 2.*np.pi)

    return theta

def theta2bearing(theta):
    """
    Converts a cartesian angle (in radians) to an azimuthal bearing
    (in radians). Assumes -2*pi <= theta <= 2*pi

    :param theta: Angle in cartesian coordinate system (+ve anticlockwise
                  from east).

    :returns: Bearing in azimuth coordinate system (+ve clockwise
                    from north).

    """
    bearing = 2. * np.pi - (theta - np.pi / 2.)
    bearing = np.mod(bearing, 2. * np.pi)

    return bearing

def makeGrid(cLon, cLat, margin=2, resolution=0.01, minLon=None, maxLon=None,
             minLat=None, maxLat=None):
    """
    Generate a grid of the distance and angle of a grid of points
    surrounding a storm centre given the location of the storm. The
    grid margin and grid size can be set in configuration
    files. xMargin, yMargin and gridSize are in degrees.

    :param float cLon: Reference longitude.
    :param float cLat: Reference latitude.
    :param float margin: Distance (in degrees) around the centre to fit the
                         grid.
    :param float resolution: Resolution of the grid (in degrees).
    :param float minLon: Minimum longitude of points to include in the grid.
    :param float maxLon: Maximum longitude of points to include in the grid.
    :param float minLat: Minimum latitude of points to include in the grid.
    :param float maxLat: Maximum latitude of points to include in the grid.

    :returns: 2 2-d arrays containing the distance (km) and bearing (azimuthal)
              of all points in a grid from the ``cLon``, ``cLat``.
    """
    if (type(cLon)==list or type(cLat)==list or
        type(cLon)==np.ndarray or type(cLat)==np.ndarray):
        raise TypeError("Input values must be scalar values")

    gridSize = int(resolution * 1000)

    if minLon:
        minLon_ = int(1000 * (minLon)) - int(1000 * margin)
    else:
        minLon_ = int(1000 * (cLon)) - int(1000 * margin)
    if maxLon:
        maxLon_ = int(1000 * (maxLon)) + int(1000 * margin) + 1
    else:
        maxLon_ = int(1000 * (cLon)) + int(1000 * margin) + 1
    if minLat:
        minLat_ = int(1000 * (minLat)) - int(1000 * margin)
    else:
        minLat_ = int(1000 * (cLat)) - int(1000 * margin)
    if maxLat:
        maxLat_ = int(1000 * (maxLat)) + int(1000 * margin) + 1
    else:
        maxLat_ = int(1000 * (cLat)) + int(1000 * margin) + 1

    xGrid = np.array(np.arange(minLon_, maxLon_, gridSize), dtype=int)
    yGrid = np.array(np.arange(minLat_, maxLat_, gridSize), dtype=int)

    R = gridLatLonDist(cLon, cLat, xGrid / 1000., yGrid / 1000.)
    np.putmask(R, R==0, 1e-30)
    theta = np.pi/2. - gridLatLonBear(cLon, cLat, xGrid / 1000., yGrid / 1000.)

    return R, theta

def makeGridDomain(cLon, cLat, minLon, maxLon, minLat, maxLat,
                   margin=2, resolution=0.01):
    """
    Generate a grid of the distance and angle of a grid of points
    across a complete model domain, given the location of the storm.

    :param float cLon: Reference longitude.
    :param float cLat: Reference latitude.
    :param float minLon: Minimum longitude of points to include in the grid.
    :param float maxLon: Maximum longitude of points to include in the grid.
    :param float minLat: Minimum latitude of points to include in the grid.
    :param float maxLat: Maximum latitude of points to include in the grid.
    :param float margin: Distance (in degrees) around the centre to fit the
                         grid.
    :param float resolution: Resolution of the grid (in degrees).

    :returns: 2 2-d arrays containing the distance (km) and bearing (azimuthal)
              of all points in a grid from the ``cLon``, ``cLat``, spanning the
              complete region.

    """
    if (type(cLon)==list or type(cLat)==list or
        type(cLon)==np.ndarray or type(cLat)==np.ndarray):
        raise TypeError("Input values must be scalar values")
    gridSize = int(resolution * 1000)
    minLon_ = int(1000 * (minLon)) - int(1000 * margin)
    maxLon_ = int(1000 * (maxLon)) + int(1000 * margin) + 1
    minLat_ = int(1000 * (minLat)) - int(1000 * margin)
    maxLat_ = int(1000 * (maxLat)) + int(1000 * margin) + 1

    xGrid = np.array(np.arange(minLon_, maxLon_, gridSize), dtype=int)
    yGrid = np.array(np.arange(minLat_, maxLat_, gridSize), dtype=int)

    R = gridLatLonDist(cLon, cLat, xGrid / 1000., yGrid / 1000.)
    np.putmask(R, R==0, 1e-30)
    theta = np.pi / 2. - gridLatLonBear(cLon, cLat,
                                        xGrid / 1000., yGrid / 1000.)
    return R, theta

def meshLatLon(cLon, cLat, margin=2, resolution=0.01):
    """
    Create a meshgrid of the longitudes and latitudes of a grid.

    :param float cLon: Longitude of centre of grid.
    :param float cLat: Latitude of centre of grid.
    :param float margin: Distance (in degrees) around the centre to
                         build the grid.
    :param float resolution: Resolution of the grid (degrees).

    :returns: Coordinate matrices for the longitude and latitude
              vectors, covering the region within ``margin``
              degrees of (``cLon``, ``cLat``).

    """
    if (type(cLon)==list or type(cLat)==list or
        type(cLon)==np.ndarray or type(cLat)==np.ndarray):
        raise TypeError("Input values must be scalar values")
    gridSize = int(1000 * resolution)

    minLon = int(1000 * (cLon - margin))
    maxLon = int(1000 * (cLon + margin)) + gridSize
    minLat = int(1000 * (cLat - margin))
    maxLat = int(1000 * (cLat + margin)) + gridSize

    xx = np.array(np.arange(minLon, maxLon, gridSize))
    yy = np.array(np.arange(minLat, maxLat, gridSize))

    xGrid, yGrid = np.meshgrid(xx, yy)
    return xGrid / 1000., yGrid / 1000.

def meshLatLonDomain(minLon, maxLon, minLat, maxLat,
                     margin=2, resolution=0.01):
    """
    Create a meshgrid of the lon/lat grid across th full model domain.

    :param float minLon: Minimum longitude of the domain.
    :param float maxLon: Maximum longitude of the domain.
    :param float minLat: Minimum latitude of the domain.
    :param float maxLat: Maximum latitude of the domain.
    :param float margin: Distance (in degrees) around the centre to
                         build the grid.
    :param float resolution: Resolution of the grid (degrees).

    :returns: Coordinate matrices for the longitude and latitude
              vectors, covering the full domain, plus an additional
              margin of ``margin`` degrees.
    """
    gridSize = int(1000 * resolution)

    minLon_ = int(1000 * (minLon - margin))
    maxLon_ = int(1000 * (maxLon + margin)) + gridSize
    minLat_ = int(1000 * (minLat - margin))
    maxLat_ = int(1000 * (maxLat + margin)) + gridSize

    xx = np.array(np.arange(minLon_, maxLon_, gridSize))
    yy = np.array(np.arange(minLat_, maxLat_, gridSize))

    xGrid, yGrid = np.meshgrid(xx, yy)
    return xGrid / 1000., yGrid / 1000.

def dist2GC(cLon1, cLat1, cLon2, cLat2, lonArray, latArray, units="km"):
    """
    Calculate the distance between an array of points and the great
    circle joining two (other) points. All input values are in
    degrees. By default returns distance in km, other units specified
    by the 'units' kwarg.

    Based on a cross-track error formulation from:
    http://williams.best.vwh.net/avform.htm#XTE

    :param float cLon1: Longitude of first point.
    :param float cLat1: Latitude of first point.
    :param float cLon2: Longitude of second point.
    :param float cLat2: Latitude of second point.
    :param lonArray: :class:`numpy.ndarray` of longitudes for which
                     the distance to the line joining the two points
                     will be calculated.
    :param latArray: :class:`numpy.ndarray` of latitudes for which the
                      distance to the line joining the two points will
                      be calculated.

    :returns: 2-d array of distances between the array points and the
              line joining two points.
    :rtype: :class:`numpy.ndarray`
    """

    # Calculate distance and bearing from first point to array of points:
    dist_ = gridLatLonDist(cLon1, cLat1, lonArray, latArray, units="rad")
    bear_ = gridLatLonBear(cLon1, cLat1, lonArray, latArray)

    #bearing of the cyclone:
    cyc_bear_ = latLon2Azi([cLon1, cLon2], [cLat1, cLat2])

    dist2GC_ = np.arcsin(np.sin(dist_) * np.sin(bear_ - cyc_bear_))

    distance = metutils.convert(dist2GC_, "rad", units)
    return distance

def coriolis(lat):
    """
    Calculate the Coriolis factor (f) for a given latitude (degrees).
    If a list is passed, return a list, else return a single value.

    :param lat: Latitude (degrees).
    :type  lat: Array-like.

    :returns: Coriolis factor.

    """
    omega = 2 * np.pi / 24. / 3600.
    f = 2 * omega * np.sin(np.radians(lat))

    return f

def find_index(array, value):
    """
    Find the index of 'array' with a value closest to 'value'

    :param array: array of data values.
    :type array: :class:`numpy.ndarray` or `list`.
    :param value: a value to search `array` for
                 (or find the index of the nearest value to `value`).

    :param int idx: index of `array` that most closely matches
                    `value`
    :raises: ValueError if `value` is a :class:`numpy.ndarray` or
             a list

    Example::

        >>> find_index(np.arange(0., 100., 0.5), 15.25)
        30

    """
    if type(value) == np.ndarray or type(value) == list:
        raise ValueError("Value cannot be an array")

    if (value > array.max()):
        # Value is above the largest value in the array - return the last index:
        return len(array) - 1
    elif (value < array.min()):
        # Value is below minimum value in the array - return the first index:
        return 0
    else:
        # argmin gives us the index corresponding to the minimum value of the array.
        idx = (abs(array - value)).argmin()
        return idx

def find_nearest(array, value):
    """
    Find the closest value in 'array' to 'value'

    :param array: array of data values.
    :type array: :class:`numpy.ndarray`
    :param value: a value to search the array for (or find the index of the
                  nearest value to 'value'
    :type value: int or float

    :returns: array[idx], where idx is the index of array that corresponds
              to the value closest to 'value'.

    :raises ValueError: If `value` is a :class:`numpy.ndarray` or
                        a list.
    :raises IndexError: If the `value` cannot be found in the `array`

    Example::

        >>> n = find_nearest( np.arange(0,100.,0.5), 15.25 )
        15.0

    """
    if type(value) == np.ndarray or type(value) == list:
        raise ValueError("Value cannot be an array")
    idx = find_index(array, value)

    try:
        v = array[idx]
    except IndexError:
        raise
    else:
        return v
