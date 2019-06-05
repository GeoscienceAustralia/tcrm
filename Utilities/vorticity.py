import numpy as np
import numpy.ma as ma
import sys, io

from Utilities.metutils import coriolis

def deriv(*args, **kwargs):
    """Calculate the derivative along a single dimension.


    Calling Sequence:
        Result = deriv([x_in,] y_in, missing=1e+20, algorithm='default')


    Positional Input Arguments:
    * x_in:  Abscissa values of y_in to take with respect to.  If 
      not set, the derivative of y_in is take with respect to unit 
      abscissa intervals.  Numeric array of same shape and size as 
      y_in.  Must be monotonic and with no duplicate values.
      Optional.  First positional argument out of two, if present.

    * y_in:  Ordinate values, to take the derivative with respect 
      to.  Numeric array vector of rank 1.  Required.  Second posi-
      tional argument, if x_in is present; only positional argument 
      if x_in is absent.


    Keyword Input Arguments:
    * missing:  If y_in and/or x_in has missing values, this is the 
      missing value value.  Scalar.  Default is 1e+20.

    * algorithm:  Name of the algorithm to use.  String scalar.
      Default is 'default'.  Possible values include:
      + 'default':  Default method (currently set to 'order1').
      + 'order1':  First-order finite-differencing (backward and
        forward differencing used at the endpoints, and centered
        differencing used everywhere else).  If abscissa intervals
        are irregular, differencing will be correspondingly asym-
        metric.


    Output Result:
    * Derivative of y_in with respect to x_in (or unit interval 
      abscissa, if x_in is not given).  Numeric array of same shape 
      and size as y_in.  If there are missing values, those elements 
      in the output are set to the value in |missing|.  For instance, 
      if y_in is only one element, a one-element vector is returned 
      as the derivative with the value of |missing|.  If there are 
      missing values in the output due to math errors and |missing| 
      is set to None, output will fill those missing values with the 
      MA default value of 1e+20.  


    References:
    * Press, W. H., et al. (1992):  Numerical Recipes in Fortran 
      77:  The Art of Scientific Computing.  New York, NY:  Cambridge
      University Press, pp. 180-184.

    * Wang, Y. (1999):  "Numerical Differentiation," Introduction to 
      MHD Numerical Simulation in Space, ESS265: Instrumentation, 
      Data Processing and Data Analysis in Space Physics (UCLA).
      URL:  http://www-ssc.igpp.ucla.edu/personnel/russell/ESS265/
      Ch10/ylwang/node21.html.


    Example with one argument, no missing values, using the default
    method:
    >>> from deriv import deriv
    >>> import Numeric as N
    >>> y = N.sin(N.arange(8))
    >>> dydx = deriv(y)
    >>> ['%.7g' % dydx[i] for i in range(4)]
    ['0.841471', '0.4546487', '-0.3501755', '-0.83305']
    >>> true = N.cos(N.arange(8))  #- Compare with exact solution
    >>> ['%.7g' % true[i] for i in range(4)]  
    ['1', '0.5403023', '-0.4161468', '-0.9899925']

    Example with two arguments with missing values, using first-
    order differencing:
    >>> x = N.arange(8)/(2.*N.pi)
    >>> y = N.sin(x)
    >>> y[3] = 1e20            #- Set an element to missing value
    >>> dydx = deriv(x, y, missing=1e20, algorithm='order1')
    >>> ['%.7g' % dydx[i] for i in range(5)]
    ['0.9957836', '0.9831985', '1e+20', '0.8844179', '1e+20']
    >>> true = N.cos(x)       #- Compare with exact solution
    >>> ['%.7g' % true[i] for i in range(5)]  
    ['1', '0.9873616', '0.9497657', '0.8881628', '0.8041098']
    """

    #- Establish y_in and x_in from *args:

    if len(args) == 1:
        y_in = args[0]
        x_in = np.arange(len(y_in), dtype=y_in.dtype)
    elif len(args) == 2:
        x_in = args[0]
        y_in = args[1]
    else:
        raise ValueError("deriv:  Bad inputs")


    #- Establish missing and algorithm from *kwargs:

    if ('missing' in kwargs) == 1:
        missing = kwargs['missing']
    else:
        missing = 1e+20

    if ('algorithm' in kwargs) == 1:
        algorithm = kwargs['algorithm']
    else:
        algorithm = 'default'


    #- Check positional and keyword inputs for possible errors:

    if (len(y_in.shape) != 1) or (len(x_in.shape) != 1):
        raise ValueError("deriv:  Inputs not a vector")
    if type(algorithm) != type(''):
        raise ValueError("deriv:  algorithm not str")


    #- Set algorithm_to_use variable, based on the algorithm keyword.
    #  The algorithm_to_use tells which algorithm below to actually
    #  use (so here is where we set what algorithm to use for default):

    if algorithm == 'default':
        algorithm_to_use = 'order1'
    else:
        algorithm_to_use = algorithm


    #- Change input to MA:  just set to input value unless there are
    #  missing values, in which case add mask:

    if missing == None:
        x = ma.masked_array(x_in)
        y = ma.masked_array(y_in)
    else:
        x = ma.masked_values(x_in, missing, copy=0)
        y = ma.masked_values(y_in, missing, copy=0)


    #- Calculate and return derivative:

    #  * Create working arrays that are consistent with a 3-point
    #    stencil in the interior and 2-point stencil on the ends:
    #    *im1 means the point before index i, *ip1 means the point 
    #    after index i, and the i index array is just plain x or 
    #    y; the endpadded arrays replicate the ends of x and y.
    #    I use an MA array filled approach instead of concatentation
    #    because the MA concatenation routine doesn't work right
    #    when the endpoint element is a missing value:

    x_endpadded = ma.zeros(x.size+2, dtype=x.dtype)
    x_endpadded[0]    = x[0]
    x_endpadded[1:-1] = x 
    x_endpadded[-1]   = x[-1]

    y_endpadded = ma.zeros(y.size+2, dtype=y.dtype)
    y_endpadded[0]    = y[0]
    y_endpadded[1:-1] = y
    y_endpadded[-1]   = y[-1]

    y_im1 = y_endpadded[:-2]
    y_ip1 = y_endpadded[2:]
    x_im1 = x_endpadded[:-2]
    x_ip1 = x_endpadded[2:]


    #  * Option 1:  First-order differencing (interior points use
    #    centered differencing, and end points use forward or back-
    #    ward differencing, as applicable):

    if algorithm_to_use == 'order1':
        dydx = (y_ip1 - y_im1) / (x_ip1 - x_im1) 


    #  * Option 2:  Bad algorithm specified:

    else:
        raise ValueError("deriv:  bad algorithm")


    #- Return derivative as Numeric array:

    return ma.filled( dydx, missing )

def has_close(data, value, rtol=1.e-5, atol=1.e-8):
    """Test if data has any values "equal" to value.

    Returns 1 if any of the elements of argument data has a value 
    "equal" to argument value; returns 0 otherwise.  If data or value
    is floating point, "equal" means where abs(data-value) <= atol + 
    rtol * abs(value).  This is essentially the same algorithm used 
    in the Numeric function allclose.  If data and value are integer,
    "equal" means strict equality.

    Positional Input Arguments:
    * data:   Data.  Scalar or Numeric array, Python list/tuple of
              any size and shape.  Floating or integer type.
    * value:  Test value.  Scalar or 1 element Numeric array, list,
              or tuple.  Floating or integer type.

    Keyword Input Arguments:
    * rtol:   "Relative" tolerance.  Default is 1.e-5.  Used in the
              comparison between data and value only if the two are 
              floating point.  Floating or integer type.
    * atol:   "Absolute" tolerance.  Default is 1.e-8.  Used in the
              comparison between data and value only if the two are 
              floating point.  Floating or integer type.

    Examples:
    >>> from has_close import has_close
    >>> data = [20., -32., -1., 2., 5., 29.]
    >>> has_close(data, -1.)
    1
    >>> has_close(data, 10.)
    0
    """


    #- Make sure data is Numeric type and value is a scalar within the
    #  function:

    dataN  = np.array(data)
    valueS = np.array(value)


    #- Safe compare if floating.  Strict compare if integer.  Any other
    #  type returns an error:

    if (dataN.dtype == float) or (valueS.dtype == type(1.)):
        closemask = np.less_equal(np.abs(dataN-valueS), 
                                  atol+rtol*np.abs(valueS))
    elif (dataN.dtype == np.int32) and (valueS.dtype == type(1)):
        closemask = np.where(dataN == valueS, 1, 0)
    else:
        print((dataN.dtype))
        print((type(valueS)))
        raise ValueError("has_close:  Inputs must be float or integer")


    #- Return true if any elements of data has value:

    if ma.maximum(closemask) == 1:
        return 1
    else:
        return 0

def can_use_sphere(longitude, latitude):
    """Test if can use sphere package.
   
    Calling Sequence:
        Result = can_use_sphere(longitude, latitude)

    Note that Result is a 3-element tuple, not a scalar.  See the
    description below for details.

    Test if the longitude-latitude domain and Python package configur-
    ation allows use of the NCAR SPHEREPACK 3.0 package.  Specifically 
    the function tests:

    (1) Can you import sphere?
    (2) Is the longitude vector evenly spaced?
    (3) Does the longitude vector entirely span the globe, but doesn't
        repeat?
    (4) Is the latitude vector gaussian?  If not, is it evenly spaced 
        and includes the poles?

    If the answer to (1) is no, then there's no use to go through the 
    other tests and the function return tuple element [0] is 0.  If (2) 
    or (3) is no, return tuple element [0] is 0.  If (4) is both not 
    gaussian and not evenly spaced with the poles, return tuple 
    element [0] is 0.  Otherwise, return tuple element [0] is 1.  Note 
    that (4) uses the latitude checker in the sphere package.


    Method Arguments:
    * longitude:  Vector of longitudes of the domain [deg].  Can be 
      in any order, and a Numeric array or regular list/tuple.

    * latitude:  Vector of latitudes of the domain [deg].  Can be in
      any order, and a Numeric array or regular list/tuple.


    Output Result:
    * A 3-element tuple:
      [0]:  1 if the longitude-latitude domain and Python package 
            configuration passes tests to allow use of the NCAR 
            SPHEREPACK 3.0 package.  0 if does not pass those tests.
      [1]:  String containing messages written to stdout by calls to 
            module sphere used for checking latitude (e.g. if it is
            gaussian, etc.).  If there are no messages or there were 
            no calls needed to sphere for checking latitude, value is 
            empty string.
      [2]:  Same as [1] except contains stderr messages.


    Examples:
    >>> from can_use_sphere import can_use_sphere
    >>> import numpy as np
    >>> lon = np.arange(36)*10
    >>> lat = [-90, -60, -30, 0, 30, 60, 90]
    >>> can_use_sphere(lon, lat)[0]
    1
    >>> lat = [-90, -60, -30, 0, 30, 60, 87]
    >>> can_use_sphere(lon, lat)[0]
    0
    >>> can_use_sphere(lon, lat)[1].splitlines()[1]
    'CANNOT PROCESS THE DATA - Latitude values are incorrect'
    """


    #- Test if the sphere package exists:

    try:  import sphere
    except ImportError:  return (0,'','')


    #- Convert input to Numeric and sort:
   
    lon = np.sort(np.array(longitude))
    lat = np.sort(np.array(latitude))


    #- Check if either vector is less than 2 elements.  If so,
    #  return false:

    if (len(lon) < 2) or (len(lat) < 2):  return (0,'','')


    #- Is the longitude vector evenly spaced?  If not, return false:

    diff_lon = lon[1:] - lon[0:-1]
    if not np.allclose(diff_lon, diff_lon[0]):  return (0,'','')


    #- Does the longitude vector exactly spans the globe but without
    #  repeating?  If not, return false:

    if not np.allclose( (lon[-1]+diff_lon[0]-lon[0]), 360.0 ):
       return (0,'','')


    #- Check latitude (e.g. whether it is gaussian, whether it includes
    #  the pole points if they are evenly spaced) and any other bad 
    #  conditions for sphere using the sphere package checker.  If 
    #  sphere can't execute correctly, return 0.  Also return standard 
    #  error and standard out:

    temp_out = sys.stdout
    temp_err = sys.stderr
    sys.stdout = io.StringIO()
    sys.stderr = io.StringIO()

    try:
        try:
            sph_obj = sphere.Sphere(lon, lat)
            sph_obj_exception = 0
        except:
            sph_obj_exception = 1
    finally:
        stdout = sys.stdout.getvalue()
        stderr = sys.stderr.getvalue()
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = temp_out
        sys.stderr = temp_err
        if sph_obj_exception == 1:  return (0, stdout, stderr)


    #- Return true if haven't returned yet:

    return (1, stdout, stderr)

def curl_2d( x, y, Fx, Fy, missing=1e+20 \
           , algorithm='default', R_sphere=6.37122e+6):
    """Curl of a vector F on a 2-D "rectangular" grid.

    The 2-D grid F is defined on is rectangular, meaning that while
    the grid spacings do not have to be even, each grid box is
    rectangular and grid boundaries are rectilinear and parallel to
    each other.  If the 2-D grid is on a sphere, we assume that lines
    of constant longitude and latitude form grids that are rectangular
    (though in reality this is not strictly true).

    Since F does not have a z-component, the curl of F is positive
    pointing as a normal out of the plane defined by the basis for 
    x and y.


    Positional Input Arguments:
    * x:  x-coordinate of each F.  1-D Numeric array with number of 
      elements equal to the number of columns in array F.  Typically 
      x is in km, or on a sphere the longitude in deg.  Floating or
      integer type.

    * y:  y-coordinate of each F.  1-D Numeric array with number of 
      elements equal to the number of rows in array F.  Typically 
      y is in km, or on a sphere the latitude in deg.  Floating or
      integer type.

    * Fx:  x-component of vector quantity F.  2-D Numeric array of 
      shape (len(y), len(x)).  Floating or integer type.

    * Fy:  y-component of vector quantity F.  2-D Numeric array of 
      same shape and size as Fx.  Floating or integer type.


    Keyword Input Arguments:
    * missing:  If Fx and/or Fy has missing values, this is the 
      missing value value.  Scalar.  Default is 1e+20.

    * algorithm:  Name of the algorithm to use.  String scalar.
      Default is 'default'.  Possible values include:

      + 'default':  Default method.  Algorithm 'order1_cartesian' 
        is used.

      + 'default_spherical':  If 'spherepack' can be used, chooses
        that method.  Otherwise 'order1_spherical' is used.  Either
        way this option assumes the domain is on a sphere of radius
        equal to keyword R_sphere.

      + 'order1_cartesian':  First-order finite-differencing (back-
        ward and forward differencing used at the endpoints, and 
        centered differencing used everywhere else).  If abscissa 
        intervals are irregular, differencing will be correspon-
        dingly asymmetric.  Grid is assumed Cartesian.  The curl
        returned is in units of F divided by the units of x and y
        (x and y must be in the same units).

      + 'order1_spherical':  First-order finite-differencing (back-
        ward and forward differencing used at the endpoints, and 
        centered differencing used everywhere else).  If abscissa 
        intervals are irregular, differencing will be correspon-
        dingly asymmetric.  Grid is assumed spherical and x and y
        are longitude and latitude (in deg).  The curl returned is 
        in units of F divided by units of R_sphere.  Note because
        of singularities in the algorithm near the poles, the
        curl for all points north of 88N or south of 88S latitude 
        are set to missing.

      + 'spherepack':  The NCAR package SPHEREPACK 3.0 is used for 
        calculating the curl.  The spatial domain is the global 
        sphere; x and y are longitude and latitude (in deg), 
        respectively; x must be evenly spaced and y must be gauss-
        ian or evenly spaced.  The domain must cover the entire 
        sphere and x and y must be longitude and latitude, respec-
        tively, in deg.  There can be no missing values.  The 
        algorithm also assumes that x and y are monotonically in-
        creasing from index 0.  The curl returned is in units of F 
        divided by m.  Thus, if F is velocity in m/s, the curl of 
        F is in units 1/s (and is the vorticity).  Note that in 
        this function's implementation, the final curl output using 
        algorithm 'spherepack' is adjusted for R_sphere not equal 
        to the default mean earth radius, if appropriate, so this 
        algorithm can be used on other planets.

        Note that SPHEREPACK operates on single precision floating
        point (Float32) values.  This function silently converts
        input to Float32 for the computations, as needed (the input 
        parameters are not altered, however).

    * R_sphere:  If the grid is the surface of a sphere, this is 
      the radius of the sphere.  This keyword is only used when 
      the algorithm is 'order1_spherical' or 'spherepack' (either 
      chosen explicitly by the algorithm keyword or automatically 
      when 'default_spherical' is used).  Default value is the 
      "mean" radius of the Earth, in meters.  "Mean" is in quotes
      because this value changes depending on what you consider the
      mean; the default in this function is the same Earth radius
      used in module sphere.  Note that if the level is a distance 
      z above the "mean" radius of the earth, R_sphere should be 
      set to the "mean" radius plus z.  
      
      Value can be a scalar or a Numeric or MA array of the same 
      size and shape as Fx.  If keyword is an array, the value of
      R_sphere for each element corresponds to the same location
      as in the Fx and Fy arrays; there should also be no missing 
      values in R_sphere (function doesn't check for this, however).


    Output:
    * Curl of F.  Units depend on algorithm used to calculate curl
      (see above discussion of keyword algorithm).  Numeric array of 
      same shape as Fx and Fy inputs.  
      
      If there are missing values, those elements in the output that 
      used missing values in the calculation of the curl are set to 
      the value in keyword missing.  If there are missing values in 
      the output due to math errors and keyword missing is set to 
      None, output will fill those missing values with the MA default 
      value of 1e+20.


    References:
    * Glickman, T. S. (Ed.) (2000):  "Curl," Glossary of Meteorology,
      Boston, MA:  American Meteorological Society, ISBN 1-878220-
      34-9.
    * Haltiner, G. J, and R. T. Williams (1980):  Numerical Prediction
      and Dynamic Meteorology, New York:  John Wiley & Sons, pp. 8-9.
    * Holton, J. R. (1992):  An Introduction to Dynamic Meteorology,
      San Diego, CA:  Academic Press, p. 482.
    * Kauffman, B. G., and W. G. Large (2002):  The CCSM Coupler,
      Version 5.0, User's Guide, Source Code Reference, and Scientific
      Description.  Boulder, CO:  NCAR.  Value for earth radius is 
      given at URL:
      http://www.ccsm.ucar.edu/models/ccsm2.0/cpl5/users_guide/node10.html
    * Overview of the CDAT Interface to the NCAR SPHEREPACK 3.0.  
      SPHEREPACK is by J. C. Adams and P. N. Swarztrauber.  More
      information on SPHEREPACK can be found at:  
         http://www.scd.ucar.edu/css/software/spherepack/
      while information on the CDAT implementation used here is at:
         http://esg.llnl.gov/cdat/getting_started/cdat.htm


    Example of a synthetic dataset of winds.  The curl is the vorti-
    city.  Images of the wind field and the vorticity are online at
    http://www.johnny-lin.com/py_pkgs/gemath/doc/test_curl_2d_add.html.
    Velocity at y greater than 60 and less than -60 are set to 0, to
    prevent spectral algorithms from giving errors at the poles.  Tiny 
    values less than 1e-10 are also set to 0:

    (a) Import statements and create data:

    >>> import Numeric as N
    >>> from curl_2d import curl_2d
    >>> nx = 181
    >>> ny = 91
    >>> x = N.arange(nx) * 2.0 - 180.0
    >>> y = N.arange(ny) * 2.0 - 90.0
    >>> y_2d = N.reshape( N.repeat(y, nx), (ny, nx) )
    >>> x_2d = N.reshape( N.repeat(N.reshape(x,(1,nx)), ny), (ny, nx) )
    >>> u = N.sin(x_2d*N.pi/180.)*N.sin((y_2d-90.)*N.pi/30.)
    >>> v = N.sin((x_2d-30.)*N.pi/180.)*N.sin((y_2d-90.)*N.pi/30.)
    >>> N.putmask( u, N.where(N.absolute(u) < 1e-10, 1, 0), 0. )
    >>> N.putmask( v, N.where(N.absolute(v) < 1e-10, 1, 0), 0. )
    >>> N.putmask( u, N.where(N.absolute(y_2d) > 60., 1, 0), 0. )
    >>> N.putmask( v, N.where(N.absolute(y_2d) > 60., 1, 0), 0. )

    (b) Use order 1 Cartesian algorithm:  If x and y are in meters,
        and u and v are in m/s, the curl is in 1/s:

    >>> curl = curl_2d(x, y, u, v, algorithm='order1_cartesian')
    >>> ['%.7g' % curl[45,i] for i in range(88,92)]
    ['-0.007251593', '-0.003628007', '0', '0.003628007']
    >>> ['%.7g' % curl[0,i]  for i in range(8,12)]
    ['0', '0', '0', '0']

    (c) Use spectral method from SPHEREPACK (need to first remove the
        repeating dateline point from the domain).  Here we take x and
        y in deg, while u and v are in m/s.  If so the curl is in 1/s.
        These results as much less than for example (b) above since
        the domain is much larger:

    >>> x = x[0:-1]
    >>> u = u[:,0:-1]
    >>> v = v[:,0:-1]
    >>> curl = curl_2d(x, y, u, v, algorithm='spherepack')
    >>> ['%.7g' % curl[45,i] for i in range(88,92)]
    ['-6.436402e-08', '-3.220152e-08', '1.33852e-13', '3.220165e-08']
    >>> ['%.7g' % curl[0,i]  for i in range(8,12)]
    ['-2.299736e-20', '-4.806512e-21', '-9.283145e-22', '-1.368445e-20']

    (d) Use order 1 spherical algorithm:  x and y are in deg and u and
        v are in m/s.  The curl is in 1/s.  Note that these results are
        nearly identical to those from SPHEREPACK, except near the poles
        the values are set to missing:

    >>> curl = curl_2d(x, y, u, v, algorithm='order1_spherical')
    >>> ['%.7g' % curl[45,i] for i in range(88,92)]
    ['-6.517317e-08', '-3.260645e-08', '0', '3.260645e-08']
    >>> ['%.7g' % curl[0,i]  for i in range(8,12)]
    ['1e+20', '1e+20', '1e+20', '1e+20']

    (e) Use "default" algorithm, which for this domain chooses the 
        order 1 cartesian algorithm:

    >>> curl = curl_2d(x, y, u, v, algorithm='default')
    >>> ['%.7g' % curl[45,i] for i in range(88,92)]
    ['-0.007251593', '-0.003628007', '0', '0.003628007']

    (f) Use "default_spherical" algorithm with missing data (equal to
        1e+20), which because of the missing data will choose the order
        1 spherical algorithm.  We also do the calculation on Mars, and
        illustrate passing in the radius as a constant as well as an
        array of the same size and shape as u and v:

    >>> u[30:46,90:114] = 1e+20
    >>> curl = curl_2d( x, y, u, v, algorithm='default_spherical' \
                      , R_sphere=3390e+3)
    >>> ['%.7g' % curl[45,i] for i in range(88,92)]
    ['-1.224875e-07', '-6.128107e-08', '1e+20', '1e+20']

    >>> R_mars = N.zeros(u.shape, typecode=N.Float) + 3390e+3
    >>> curl = curl_2d( x, y, u, v, algorithm='default_spherical' \
                      , R_sphere=R_mars )
    >>> ['%.7g' % curl[45,i] for i in range(88,92)]
    ['-1.224875e-07', '-6.128107e-08', '1e+20', '1e+20']
    """
    #import MA
    #import Numeric as N
    #from has_close import has_close
    #from can_use_sphere import can_use_sphere



#------- Nested Function:  Curl By Order 1 Cartesian Finite Diff. ------

    def _order1_cartesian_curl():
        """
        Calculate curl using Cartesian first-order differencing.

        An ma array is returned.

        Algorithm:  First-order differencing (interior points use 
        centered differencing, and end points use forward or backward 
        differencing, as applicable) in Cartesian coordinates (see 
        Glickman [2000], p. 194).
        """
        dFy_dx_N = np.zeros( (len(y), len(x)), dtype=float )
        dFx_dy_N = np.zeros( (len(y), len(x)), dtype=float )

        for iy in range(len(y)):
            dFy_dx_N[iy,:] = deriv( x, np.ravel(Fy[iy,:]) \
                                  , missing=missing, algorithm='order1')

        for ix in range(len(x)):
            dFx_dy_N[:,ix] = deriv( y, np.ravel(Fx[:,ix]) \
                                  , missing=missing, algorithm='order1')

        dFy_dx = ma.masked_values(dFy_dx_N, missing, copy=0)
        dFx_dy = ma.masked_values(dFx_dy_N, missing, copy=0)

        return dFy_dx - dFx_dy




#------- Nested Function:  Curl By Order 1 Spherical Finite Diff. ------

    def _order1_spherical_curl():
        """
        Calculate curl using spherical first-order differencing.

        An ma array is returned.

        Algorithm:  First-order differencing (interior points use 
        centered differencing, and end points use forward or backward 
        differencing, as applicable) in spherical coordinates (see 
        Holton 1992).

        Note because of singularities in the algorithm near the poles, 
        the return value for all points north of 88N or south of 88S 
        latitude are set to mask true.

        Key to some variables:
        * xrad is x in radians (x is assumed in deg)
        * yrad is y in radians (y is assumed in deg)
        * ddFy is d(Fy) / d(xrad) as a masked array.  The version with
          a "_N" suffix means the Numeric version.
        * ddFx is d(Fx * cos(yrad)) / d(yrad) as a masked array.  The 
          version with a "_N" suffix means the Numeric version.
        """

        xrad = x * np.pi/180.0
        yrad = y * np.pi/180.0


        #- Derivative preliminaries:

        ddFy_N = np.zeros((len(yrad), len(xrad)), dtype=float)
        ddFx_N = np.zeros((len(yrad), len(xrad)), dtype=float)

        for iy in range(len(yrad)):
            ddFy_N[iy, :] = deriv(xrad, np.ravel(Fy[iy, :]),
                                  missing=missing, algorithm='order1')

        for ix in range(len(xrad)):
            tmp_MA = ma.masked_values(np.ravel(Fx[:, ix]), missing, copy=0)
            tmp_N  = ma.filled(tmp_MA * np.cos(yrad), missing)
            ddFx_N[:, ix] = deriv(yrad, tmp_N,
                                  missing=missing,
                                  algorithm='order1')

        ddFy = ma.masked_values(ddFy_N, missing, copy=0)
        ddFx = ma.masked_values(ddFx_N, missing, copy=0)


        #- Calculate the curl:

        tmpr = (ddFy - ddFx) / \
               (R_sphere * np.reshape(np.repeat(np.cos(yrad), len(xrad)),
                                      (len(yrad), len(xrad))))


        #- Make points "near" poles be missing and return:

        np_mask_lat =  88.0
        sp_mask_lat = -88.0
        yarr = np.reshape(np.repeat(y,len(x)), (len(y), len(x)))

        np_mask = ma.make_mask(np.where(yarr > np_mask_lat, 1, 0))
        sp_mask = ma.make_mask(np.where(yarr < sp_mask_lat, 1, 0))

        tmpr_mask = tmpr.mask
        tmpr_mask = ma.mask_or(tmpr_mask, np_mask)
        tmpr_mask = ma.mask_or(tmpr_mask, sp_mask)

        return ma.masked_array(tmpr, mask=tmpr_mask)




#--------- Nested Function:  Curl By SPHEREPACK Spectral Method --------

    def _spherepack_curl():
        """
        Calculate curl using SPHEREPACK spectral methods.
        
        If the value of R_sphere is not the same as the earth radius 
        used in package sphere, adjust curl to reflect R_sphere.  A 
        Numeric array is returned.

        SPHEREPACK outputs a number of diagnostic messages to
        stdout (and perhaps to stderr).  The one that's irritating
        is the one saying data will be converted to Float32.  Thus,
        we silently ensure that calculations are used only on Float32 
        data, to surpress the warnings.  I did this because redirect-
        ing stdout/stderr turned out to be too difficult; just
        altering the sys attributes didn't work.
        """
        import sphere

        sph_obj = sphere.Sphere( x.astype(ma.Float32) \
                               , y.astype(ma.Float32) )
        if np.allclose(np.array(R_sphere), sphere.radius):
            curl = sph_obj.vrt( Fx.astype(ma.Float32) \
                              , Fy.astype(ma.Float32) )
        else:
            curl = sph_obj.vrt( Fx.astype(ma.Float32) \
                              , Fy.astype(ma.Float32) ) \
                 * sphere.radius / R_sphere

        return curl




#-------------- Overall Function:  Settings, Body, Return --------------

    #- Choose algorithm to compute curl:
    
    if algorithm == 'default':
        _calculate_curl = _order1_cartesian_curl

    elif algorithm == 'default_spherical':
        if (can_use_sphere(x,y)[0] == 1) and \
           (not has_close(Fx, missing)) and \
           (not has_close(Fy, missing)):
            _calculate_curl = _spherepack_curl
        else:
            _calculate_curl = _order1_spherical_curl

    elif algorithm == 'order1_cartesian':
        _calculate_curl = _order1_cartesian_curl

    elif algorithm == 'order1_spherical':
        _calculate_curl = _order1_spherical_curl

    elif algorithm == 'spherepack':
        if has_close(Fx, missing) or has_close(Fy, missing):
            raise ValueError("curl_2d:  has missing values")
        else:
            _calculate_curl = _spherepack_curl

    else:
        raise ValueError("curl_2d:  bad algorithm")


    #- Calculate curl and return from function:

    return ma.filled( _calculate_curl(), missing )

def relative(u, v, x, y):
    return curl_2d(x, y, u, v, 
                   algorithm='order1_spherical', 
                   R_sphere=6378137.)

def absolute(u, v, x, y):
    rel = curl_2d(x, y, u, v, 
                  algorithm='order1_spherical', 
                  R_sphere=6378137.)

    xx, yy = np.meshgrid(x, y)
    f = coriolis(yy)
    zeta = rel + f
    return zeta
