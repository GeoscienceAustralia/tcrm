"""
:py:mod:`loadData` - load TC track data from formatted files
============================================================

.. module:: loadData
   :synopsis: Load formatted csv file containing tropical
              cyclone track information

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

Load a formatted csv file that contains tropical cyclone track
information and return a collection of :class:`Track` objects.

The format of the files is defined in a named section of the input
configuration file. See the :ref:`Source Formats <sourceformats>` section
of the User Guide for details on specifying the format.

TODO:
Modify the source to be a dict containing the kwargs to
:mod:`numpy.genfromtxt` to make reading additional formats possible
(e.g. include converters for some data fields, as in B-deck format
track files.

"""

import sys
import logging
import numpy as np
from . import metutils
from . import maputils
from . import nctools
from . import interp3d

from datetime import datetime, timedelta
from .columns import colReadCSV
from Utilities.config import ConfigParser, cnfGetIniValue
from Utilities.track import Track, trackFields, trackTypes

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

"""

TRACKFILE_COLS = ('Indicator', 'CycloneNumber', 'Year', 'Month',
                  'Day', 'Hour', 'Minute', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'WindSpeed', 'rMax', 'EnvPressure')

TRACKFILE_FMTS = ('i', 'i', 'i', 'i',
                  'i', 'i', 'i', 'f',
                  'f', 'f', 'f', 'f', 'f',
                  'f', 'f', 'f')

TRACKFILE_OUTFMT = ('%i,%i,%i,%i,'
                    '%i,%i,%i,%5.1f,'
                    '%8.3f,%8.3f,%6.2f,%6.2f,%7.2f,'
                    '%6.2f,%6.2f,%7.2f')

class Track(object):

    ""
    A single tropical cyclone track.

    The object exposes the track data through the object attributes.
    For example, If `data` contains the tropical cyclone track data
    (`numpy.array`) loaded with the :meth:`readTrackData` function,
    then the central pressure column can be printed out with the
    code::

        t = Track(data)
        print(t.CentralPressure)


    :type  data: numpy.ndarray
    :param data: the tropical cyclone track data.
    ""

    def __init__(self, data):
        self.data = data
        self.trackId = None
        self.trackfile = None
        self.trackMinPressure = None
        self.trackMaxWind = None

    def __getattr__(self, key):
        ""
        Get the `key` from the `data` object.

        :type  key: str
        :param key: the key to lookup in the `data` object.
        ""
        if key.startswith('__') and key.endswith('__'):
            return super(Track, self).__getattr__(key)
        return self.data[key]

#FORMATS:
bdeck = {
    "delimiter": ",",
    "names" : ("basin", "num", "date", "lat", "lon", "vmax", "pressure", "poci", "rmax"),
    "dtype" : ("|U2", "i", "object", "f8", "f8", "f8", "f8", "f8", "f8"),
    "usecols" : (0, 1, 2, 6, 7, 8, 9, 17, 19),
    "converters" : {
                0: lambda s: s.strip(),
                1: lambda s: s.strip(),
                2: lambda s: datetime.strptime(s.strip(), "%Y%m%d%H"),
                6: lambda s: float(s.strip("NSEW")) / 10.,
                7: lambda s: float(s.strip("NSEW")) / 10.,
                8: lambda s: float(s.strip()),
                9: lambda s: float(s.strip()),
                11: lambda s: metutils.convert(float(s.strip()), "nm", "km")
            },
    "autostrip" : True
    }

ibtracs = {
    "delimiter" : ",",
    "names" : ("tcserialno", "season", "num", "date", "lat", "lon", "pressure"),
    "dtype" : ("|U13", "i", "i", "object", "f8", "f8", "f8"),
    "usecols" : (0, 1, 2, 6, 8, 9, 11),
    "converters" : {
                0: lambda s: s.strip(),
                1: lambda s: int(s.strip()),
                2: lambda s: int(s.strip()),
                6: lambda s: datetime.strptime(s.strip(), "%Y-%m-%d %H:%M:%S"),
                8: lambda s: float(s.strip()),
                9: lambda s: float(s.strip()),
                11: lambda s: float(s.strip())
             },
    "skip_header" : 3,
    "autostrip" : True,
    }

tcrm = {
    "delimiter" : ",",
    "names" : ('CycloneNumber', 'Datetime', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'EnvPressure', 'rMax'),
    "dtype" : ('i', 'object', 'f', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'),
    "converters" : {
        0: lambda s: int(float(s.strip() or 0)),
        1: lambda s: datetime.strptime(s.strip(), "%Y-%m-%d %H:%M:%S"),
        5: lambda s: metutils.convert(float(s.strip() or 0), 'kph', 'mps'),
        6: lambda s: maputils.bearing2theta(float(s.strip() or 0) * np.pi / 180.),
        7: lambda s: metutils.convert(float(s.strip() or 0), 'hPa', 'Pa'),
        8: lambda s: metutils.convert(float(s.strip() or 0), 'hPa', 'Pa'),
        },
    "autostrip" : True
    }

"""

def getSpeedBearing(index, lon, lat, deltatime, ieast=1,
                    missingValue=sys.maxsize):
    """
    Calculate the speed and bearing of a TC.

    :param index: Array of 0/1 indicating start of new TC (1)
    :type index: :class:`numpy.ndarray`
    :param lon: Longitudes of TC positions.
    :type  lon: :class:`numpy.ndarray`
    :param lat: Latitudes of TC positions.
    :type  lat: :class:`numpy.ndarray`
    :param deltatime: Time difference (hours) between
                      consecutive TC observations.
    :type  deltatime: :class:`numpy.ndarray`
    :param int ieast: Indicate which direction has positive
                      longitude. 1 = positive longitude eastwards
                      -1 = positive longiture westwards.

    :param missingValue: Replace questionable values with `missingValue`.
    :type missingValue: int or float, default = `sys.maxsize`


    :returns: speed and bearing : :class:`numpy.ndarray`

    Example::

        >>> speed, bearing = getSpeedBearing(index, lon, lat, deltatime)

    """

    bear_, dist_ = maputils.latLon2Azi(lat, lon, ieast, azimuth=0)
    assert bear_.size == index.size - 1
    assert dist_.size == index.size - 1
    bearing = np.zeros(index.size, 'f')
    bearing[1:] = bear_
    np.putmask(bearing, index, missingValue)

    dist = np.zeros(index.size, 'f')
    dist[1:] = dist_
    speed = dist / deltatime
    # Delete speeds less than 0, greated than 200,
    # or where indicator == 1.
    np.putmask(speed, (speed < 0), missingValue)
    np.putmask(speed, (speed > 200), missingValue)
    np.putmask(speed, index, missingValue)
    np.putmask(speed, np.isnan(speed), missingValue)

    return speed, bearing


def maxWindSpeed(index, deltatime, lon, lat, pressure, penv,
                 gustfactor=0.9524):
    """
    Calculate the 10-minute-mean maximum wind speed from the central
    pressure deficit, using the method described in Holland et al. (2010).

    :param indicator: Array (values of 1 or 0) indicating the beginning of
                      a new TC in the input dataset.
    :param deltatime: Time difference (in hours) between each point in the
                      record.
    :param lon: Longitudes of TC postions.
    :param lat: Latitudes of TC positions.
    :param pressure: Central pressure estimate of TCs (hPa).
    :param penv: Environmental pressure estimates for each TC postion (hPa).
    :param float gf: Gust factor - default value represents converting from a
                     1-minute sustained wind speed to a 10-minute mean wind
                     speed. Based on Harper et al. 2010, WMO-TD1555.
    :type indicator: :class:`numpy.ndarray`
    :type deltatime: :class:`numpy.ndarray`
    :type lon: :class:`numpy.ndarray`
    :type lat: :class:`numpy.ndarray`
    :type pressure: :class:`numpy.ndarray`
    :type penv: :class:`numpy.ndarray`

    :returns: :class:`numpy.ndarray` of estimated wind speed based on
              central pressure deficit.

    Example::

      >>> v = maxWindSpeed(indicator, dt, lon, lat, pressure, penv)

    """

    # Speed and bearing:
    speed, bearing = getSpeedBearing(index, lon, lat, deltatime)
    speed = metutils.convert(speed, 'kmh', 'mps')
    np.putmask(speed, speed > 10e+3, 0)

    # Pressure deficit:
    deltap = penv - pressure

    # Pressure rate of change
    dpt = np.zeros(index.size, 'f')
    dpt[1:] = np.diff(pressure)
    dpdt = dpt / deltatime
    np.putmask(dpdt, index, 0)
    np.putmask(dpdt, np.isnan(dpdt) |
               np.isinf(dpdt) |
               (np.abs(dpdt) > 5.), 0)

    # Estimated pressure at the radius of maximum wind:
    prmw = pressure + deltap / 3.7

    # Calculate thermodynamic variables at RMW:
    tsurf = 28.0 - 3 * (np.abs(lat) - 10.) / 20.
    qmix = 0.9 * (3.802 / prmw) * np.exp(17.67 * tsurf / (243.5 + tsurf))
    tvs = (tsurf + 273.15) * (1. + 0.81 * qmix)
    rho = prmw * 100. / (tvs * 287.04)

    chi = 0.6 * (1.0 - deltap / 215.)
    beta = -0.000044 * np.power(deltap, 2.) + \
        0.01 * deltap + 0.03 * dpdt - 0.014 * np.abs(lat) + \
        0.15 * np.power(speed, chi) + 1.

    # Holland's P-W relation derives a 1-minute mean wind speed, so we often
    # need to convert to some other averaging period. I use the recommendations
    # of Harper et al. (2010) WMO TD-1555:
    # Common values are( Assuming "At-sea" conditions):
    # 10-min mean: 0.95 (default)
    # 3-second gust: 1.11

    v = gustfactor * np.sqrt(deltap * 100 * beta / (rho * np.exp(1.)))
    np.putmask(v, (np.isnan(v) |
                   np.isinf(v) |
                   (pressure >= 10e+7) |
                   (pressure <= 0) |
                   (speed >= 10e+7)), 0)

    return v


def getInitialPositions(data):
    """

    Get the array indices corresponding to the initial position of TCs in
    the input dataset. This is done through examining the data for a number
    of specific fields to see when they change, or if the data has a field
    that indicates such an instance.

    :param dict data: :class:`dict` of arrays that contains the data loaded
                      from the input file

    :returns: :class:`numpy.ndarray` of indices that can be used to slice
              the observations and return those corresponding to an
              initial TC position.

    .. note::
       Using only the 'num' field will result in different results
       than using 'num' and 'season'. From experience, the 'num' field
       in best-track datasets refers to the sequential number of the
       storm for that season - i.e. it starts at 1 for each season and
       increments for each new storm. The use of 'num' only should be
       reserved for those situations where the dataset is known to
       have unique numbers for each storm (e.g in simulated data).

    """
    try:
        indicator = np.array(data['index'], 'i')
        LOG.info("Using index contained in file to "
                 "determine initial TC positions")
        return indicator
    except ValueError:
        LOG.error("'index' field cannot be converted to integer")
    except KeyError:
        pass

    try:
        tcSerialNo = data['tcserialno']
        LOG.info("Using TC serial number to determine initial "
                 "TC positions")
        indicator = np.ones(len(tcSerialNo), 'i')
        for i in range(1, len(tcSerialNo)):
            if tcSerialNo[i] == tcSerialNo[i - 1]:
                indicator[i] = 0
        return indicator
    except (ValueError, KeyError):
        pass

    try:
        num = np.array(data['num'], 'i')
        season = np.array(data['season'], 'i')
        LOG.info("Using season and TC number to determine initial "
                 "TC positions")
        indicator = np.ones(num.size, 'i')
        for i in range(1, len(num)):
            if (season[i] == season[i - 1]) and (num[i] == num[i - 1]):
                indicator[i] = 0
        return indicator
    except KeyError:
        pass
    except ValueError:
        LOG.error("'num' field cannot be converted to an integer")


    try:
        num = np.array(data['num'], 'i')
        LOG.info("Using TC number to determine initial TC positions "
                 "(no season information)")
        indicator = np.ones(num.size, 'i')
        ind_ = np.diff(num)
        ind_[np.where(ind_ != 0)] = 1
        indicator[1:] = ind_
        return indicator
    except KeyError:
        pass
    except ValueError:
        LOG.error("'num' field cannot be converted to an integer")

    raise KeyError(('Insufficient input file columns have been specified'
                    'Check the input file has enough fields to determine'
                    'TC starting positions'))


def date2ymdh(dates, datefmt='%Y-%m-%d %H:%M:%S'):
    """
    Convert date strings to arrays for date components.

    :param dates: Array of str objects that describe a date.
    :type  dates: :class:`numpy.ndarray`

    :param str datefmt: Format string of the date values in `dates`,
                        default='%Y-%m-%d %H:%M:%S'

    :returns: :class:`numpy.ndarray` arrays containing the year,
              month, day, hour and minute of the dates contained
              in the input array `dates`.

    """

    import re
    pattern = re.compile("%y")
    if pattern.search(datefmt):
        raise ValueError("Cannot use 2-digit year formats in date format")

    year = np.empty(len(dates), 'i')
    month = np.empty(len(dates), 'i')
    day = np.empty(len(dates), 'i')
    hour = np.empty(len(dates), 'i')
    minute = np.empty(len(dates), 'i')
    datetimes = np.empty(len(dates), datetime)

    for i in range(len(dates)):
        try:
            d = datetime.strptime(str(dates[i]), datefmt)
        except ValueError as e:
            LOG.exception("Error in date information for record {0}".format(i))
            LOG.exception(repr(e))
            raise
        else:
            year[i] = d.year
            month[i] = d.month
            day[i] = d.day
            hour[i] = d.hour
            minute[i] = d.minute
            datetimes[i] = d

    return year, month, day, hour, minute, datetimes


def parseDates(data, indicator, datefmt='%Y-%m-%d %H:%M:%S'):
    """
    Parse the date/time information to extract year, month, day, hour and
    minute details for the input dataset.

    :param dict data: :class:`dict` of arrays that contains the data loaded
                      from the input file.
    :param indicator: Array (values of 1 or 0) indicating the beginning of
                      a new TC in the input dataset.
    :param str datefmt: Format string of the date values in `dates`,
                        default='%Y-%m-%d %H:%M:%S'

    :type indicator: :class:`numpy.ndarray`

    :returns: :class:`numpy.ndarray`s of year, month, day, hour, minute
              and :class:`datetime.datetime` objects.

    """
    try:
        year, month, day, hour, minute, datetimes = date2ymdh(data['date'], datefmt)
    except (KeyError):
        # Sort out date/time information:
        month = np.array(data['month'], 'i')
        day = np.array(data['day'], 'i')
        hour = np.array(data['hour'], 'i')
        try:
            year = np.array(data['year'], 'i')
        except (ValueError, KeyError):
            # Create dummy variable year - applicable for datasets
            # such as WindRiskTech which contain no year information.
            year = np.zeros(month.size, 'i')
            for i in range(len(year)):
                if indicator[i] > 0:
                    fill_year = 2000
                if month[i] == 1:
                    fill_year = 2001
                year[i] = fill_year

        try:
            minute = np.array(data['minute'], 'i')
        except (ValueError, KeyError):
            if hour.max() >= 100:
                minute = np.mod(hour, 100)
                hour = hour // 100
            else:
                LOG.warning("Missing minute data from input data" + \
                               "- setting minutes to 00 for all times")
                minute = np.zeros((hour.size), 'i')

        datetimes = np.array([datetime(y, m, d, h, mn) for y, m, d, h, mn
                              in zip(year, month, day, hour, minute)])

    return year, month, day, hour, minute, datetimes

def parseAge(data, indicator):
    """
    Parse the TC age information to get a proxy date record. Assumes every TC
    starts at 2000-01-01 00:00 and calculates year, month, day, hour and
    minute values based on the age field.

    Assumes the age is given in hours since the initial timestep.

    :param dict data: :class:`dict` of arrays that contains the data loaded
                      from the input file.
    :param indicator: Array (values of 1 or 0) indicating the beginning of
                      a new TC in the input dataset.
    :type indicator: :class:`numpy.ndarray`

    :returns: :class:`numpy.ndarray`s of year, month, day, hour, minute
              and :class:`datetime.datetime` objects.
    """

    start_time = datetime(2000, 1, 1, 0, 0)
    delta = np.array([timedelta(i) for i in data['age']/24.])
    times_ = start_time + delta
    year = np.ones(indicator.size, 'i')
    month = np.ones(indicator.size, 'i')
    day = np.ones(indicator.size, 'i')
    hour = np.ones(indicator.size, 'i')
    minute = np.zeros(indicator.size, 'i')

    for i, dt in enumerate(times_):
        year[i] = dt.year
        month[i] = dt.month
        day[i] = dt.day
        hour[i] = dt.hour
        minute[i] = dt.minute

    return year, month, day, hour, minute, times_


def getTimeDelta(year, month, day, hour, minute):
    """
    Calculate the time difference between consecutive observations.

    :param year: :class:`numpy.ndarray` of the year of all observations.
    :param month: As for year, but for the month of observation.
    :param day: As for year, but for the day of observation.
    :param hour: As for year, but for the hour of the observation.
    :param minutes: As for year, but for the hour of the observation.
    :param seconds: As for year, but for the hour of the observation.

    :returns: :class:`numpy.ndarray` of time difference between
              observations in hours.

    """
    dates = [datetime(*x) for x in zip(year, month, day, hour, minute)]
    diffs = [d1 - d0 for d0, d1 in zip(dates[:-1], dates[1:])]
    return np.array([0.0] + [round(d.days * 24. + d.seconds / 3600.)
                             for d in diffs], 'f')

def getTimeElapsed(indicator, year, month, day, hour, minute):
    """
    Calculate the age in hours of each event.


    :param indicator: Array (values of 1 or 0) indicating the beginning of
                      a new TC in the input dataset.
    :type indicator: :class:`numpy.ndarray`
    :param year: :class:`numpy.ndarray` of the year of all observations.
    :param month: As for year, but for the month of observation.
    :param day: As for year, but for the day of observation.
    :param hour: As for year, but for the hour of the observation.
    :param minutes: As for year, but for the hour of the observation.

    :returns: :class:`numpy.ndarray` of time since the initial observation
              for each TC (in hours).
    """
    dates = [datetime(*x) for x in zip(year, month, day, hour, minute)]

    timeElapsed = []
    for i, d in zip(indicator, dates):
        if i == 1:
            start = d
            timeElapsed.append(0.0)
        else:
            delta = d - start
            timeElapsed.append(delta.total_seconds()/3600.)

    return np.array(timeElapsed)


def getTime(year, month, day, hour, minute):
    """
    Calculate the number of days since 0001-01-01 00:00:00 UTC + 1

    :param year: Year values.
    :param month: Month values.
    :param day: Day values.
    :param hour: Hour values.
    :param minute: Minute values.

    :type year: :class:`numpy.ndarray` or list
    :type month: :class:`numpy.ndarray` or list
    :type day: :class:`numpy.ndarray` or list
    :type hour: :class:`numpy.ndarray` or list
    :type minute: :class:`numpy.ndarray` or list

    :return: :class:`numpy.ndarray` of days since 0001-01-01 00:00:00 UTC + 1

    """
    dates = [datetime(*x) for x in zip(year, month, day, hour, minute)]
    return np.array([d.toordinal() + d.hour / 24. for d in dates], 'f')


def julianDays(year, month, day, hour, minute):
    """
    Calculate the julian day (day of year) based on the known date/time
    information.

    :param year: :class:`numpy.ndarray` of the year of all observations.
    :param month: As for year, but for the month of observation.
    :param day: As for year, but for the day of observation.
    :param hour: As for year, but for the hour of the observation.
    :param minute: As for year, but for the hour of the observation.

    :returns: :class:`numpy.ndarray` of julian day values for each
              observation.

    """
    LOG.debug("Calculating julian day (day of year) values")

    if np.any(year < 0):
        raise ValueError("Error in input year information - check input file")
    if np.any(month >= 13):
        raise ValueError("Error in input month information - check input file")
    if np.any(day > 31):
        raise ValueError("Error in input day information - check input file")
    if np.any(hour > 24):
        raise ValueError("Error in input hour information - check input file")
    if np.any(minute > 60):
        raise ValueError(
            "Error in input minute information - check input file")

    # set all years prior to 1900 to 1904 - strftime() requires year >=1900;
    # and in the Gregorian calendar, 1900 is not a leap year (and there are
    # many years prior to 1900 that are!).
    second = np.zeros((hour.size), 'i')
    jyear = np.copy(year)
    jyear[np.where(jyear < 1900)] = 1904
    day = [datetime(jyear[i], month[i], day[i], hour[i], minute[i],
                    second[i]) for i in range(year.size)]

    jdays = np.array([int(day[i].strftime("%j")) for
                      i in range(year.size)])
    return jdays


def ltmPressure(jdays, time, lon, lat, ncfile, ncvar='slp'):
    """
    Extract pressure value from a daily long-term mean SLP dataset at the
    given day of year and lon,lat position
    To use this function (and hence some form of daily LTM SLP data) requires
    knowledge of the day of year.

    :param jdays: Julian day (day of year) values.
    :param time: Time of day for each observation (fraction of a day).
    :param lon: Longitude of TC position.
    :param lat: Latitude of TC position.
    :param str ncfile: Path to netCDF file containing daily long-term mean
                       sea level pressure data.
    :param str ncvar: Name of the netcdf variable that holds the SLP data

    :type  jdays: :class:`numpy.ndarray`
    :type  time: :class:`numpy.ndarray`
    :type  lon: :class:`numpy.ndarray`
    :type  lat: :class:`numpy.ndarray`

    :returns: :class:`numpy.ndarray` of long-term mean sea level pressure
              values at the day of year and positions given.
    """
    jtime = jdays + np.modf(time)[0]
    coords = np.array([jtime, lat, lon])

    LOG.debug("Sampling data from MSLP data in {0}".format(ncfile))
    ncobj = nctools.ncLoadFile(ncfile)
    if ncvar not in ncobj.variables:
        raise KeyError(f"{ncfile} does not contain a variable called '{ncvar}'")
    
    slpunits = getattr(ncobj.variables[ncvar], 'units')

    data = nctools.ncGetData(ncobj, ncvar)
    # Get the MSLP by interpolating to the location of the TC:
    penv = interp3d.interp3d(data, coords, scale=[365., 180., 360.],
                             offset=[0., -90., 0.])
    penv = metutils.convert(penv, slpunits, 'hPa')
    del data
    ncobj.close()
    del ncobj

    return penv

def getPoci(penv, pcentre, lat, jdays, eps,
            coeffs=[2324.1564738613392, -0.6539853183796136,
                    -1.3984456535888878, 0.00074072928008818927,
                    0.0044469231429346088, -1.4337623534206905],
            missingValue=sys.maxsize):
    """
    Calculate a modified pressure for the outermost closed isobar, based
    on a model of daily long-term mean SLP values, central pressure,
    latitude and day of year.

    :param penv: environmental pressure estimate (from long term mean pressure
                 dataset, hPa).
    :param pcentre: Central pressure of storm (hPa).
    :param lat: Latitude of storm (degrees).
    :param jdays: Julian day (day of year).
    :param eps: random variate. Retained as a constant for a single storm.
    :param list coeffs: Coefficients of the model. Defaults based on 
                        Southern Hemisphere data 
                        (IBTrACS v03r06, 1981-2014). 

    :returns: Revised estimate for the pressure of outermost closed isobar.
    """

    if len(coeffs) < 6:
        LOG.warn("Insufficient coefficients for poci calculation")
        LOG.warn("Using default values")
        coeffs=[2324.1564738613392, -0.6539853183796136,
                -1.3984456535888878, 0.00074072928008818927,
                0.0044469231429346088, -1.4337623534206905]

    if isinstance(penv, (np.ndarray, list)) and \
      isinstance(pcentre, (np.ndarray, list)) and \
      isinstance(lat, (np.ndarray, list)) and \
      isinstance(jdays, (np.ndarray, list)):
        assert len(penv) == len(pcentre)
        assert len(penv) == len(lat)
        assert len(penv) == len(jdays)
      
    poci_model = coeffs[0] + coeffs[1]*penv + coeffs[2]*pcentre \
      + coeffs[3]*pcentre*pcentre + coeffs[4]*lat*lat + \
        coeffs[5]*np.cos(np.pi*2*jdays/365) + eps

    if isinstance(poci_model, (np.ndarray, list)):
        nvidx = np.where(pcentre == missingValue)
        poci_model[nvidx] = np.nan

        nvidx = np.where(penv < pcentre)
        poci_model[nvidx] = np.nan

    elif penv < pcentre:
        poci_model = np.nan
    elif pcentre == missingValue:
        poci_model = np.nan
    
    return poci_model
    

def filterPressure(pressure, inputPressureUnits='hPa',
                   missingValue=sys.maxsize):
    """
    Filter pressure values to remove any non-physical values.

    :param pressure: input pressure values to check.
    :param str inputPressureUnits: The units of the pressure values.
                     Can be one of ``hPa``, ``Pa``, ``kPa``,
                     ``Pascals`` or ``mmHg``.
    :param missingValue: replace all null values in the input data
                         with this value.
    :type pressure: :class:`numpy.ndarray`
    :type missingValue: int or float (default ``sys.maxsize``)

    :returns: :class:`numpy.ndarray` with only valid pressure values.

    """

    novalue_index = np.where(pressure == missingValue)
    pressure = metutils.convert(pressure, inputPressureUnits, "hPa")
    pressure[novalue_index] = missingValue

    # Convert any non-physical central pressure values to maximum integer
    # This is required because IBTrACS has a mix of missing value codes
    # (i.e. -999, 0, 9999) in the same global dataset.
    pressure = np.where((pressure < 600) | (pressure > 1100),
                        missingValue, pressure)
    return pressure

def getMinPressure(track, missingValue=sys.maxsize):
    """
    Determine the minimum pressure of a :class:`Track` instance

    :param track: A :class:`Track` instance
    :param missingValue: Replace missing values with this value
                         (default ``sys.maxsize``).

    :returns: :class:`Track.trackMinPressure` attribute updated

    """

    p = track.CentralPressure
    if np.all(p == missingValue):
        track.trackMinPressure = missingValue
    else:
        track.trackMinPressure = p[p != missingValue].min()

def getMaxWind(track, missingValue=sys.maxsize):
    """
    Determine the maximum wind speed of a :class:`Track` instance

    :param track: A :class:`Track` instance
    :param missingValue: replace all null values in the input data
                         with this value.
    :type missingValue: int or float (default ``sys.maxsize``)

    :returns: :class:`Track.trackMaxWind` attribute updated with calculated
              wind speed updated.

    """

    w = track.WindSpeed
    if np.all(w == missingValue):
        track.trackMaxWind = missingValue
    else:
        track.trackMaxWind = w[w != missingValue].max()


def loadTrackFile(configFile, trackFile, source, missingValue=0,
                  calculateWindSpeed=True):
    """
    Load TC track data from the given input file, from a specified source.
    The configFile is a configuration file that contains a section called
    'source' that describes the data.
    This returns a collection of :class:`Track` objects that contains
    the details of the TC tracks in the input file.

    :param str configFile: Configuration file with a section ``source``.
    :param str trackFile: Path to a csv-formatted file containing TC data.
    :pararm str source: Name of the source format of the TC data. There
                        *must* be a section in ``configFile`` matching
                        this string, containing the details of the format
                        of the data.
    :param missingValue: Replace all null values in the input data with
                         this value (default=0).
    :param boolean calculateWindSpeed: Calculate maximum wind speed using
                                       a pressure-wind relation described
                                       in :func:`maxWindSpeed`

    :returns: A collection of :class:`Track` objects.
              If any of the variables are not present in the input
              dataset, they are (where possible) calculated
              (date/time/windspeed), sampled from default datasets
              (e.g. environmental pressure) or set to the missing value.

    Example::

      >>> tracks = loadTrackFile('tcrm.ini', 'IBTRaCS.csv', 'IBTrACS' )

    """

    LOG.info("Loading %s" % trackFile)
    inputData = colReadCSV(configFile, trackFile, source) #,
                          #nullValue=missingValue)

    config = ConfigParser()
    config.read(configFile)

    inputSpeedUnits = config.get(source, 'SpeedUnits')
    inputPressureUnits = config.get(source, 'PressureUnits')
    inputLengthUnits = config.get(source, 'LengthUnits')
    inputDateFormat = config.get(source, 'DateFormat')

    if config.getboolean('DataProcess', 'FilterSeasons'):
        startSeason = config.getint('DataProcess', 'StartSeason')
        idx = np.where(inputData['season'] >= startSeason)[0]
        inputData = inputData[idx]

    # Determine the initial TC positions...
    indicator = getInitialPositions(inputData)


    # Sort date/time information
    if 'age' in inputData.dtype.names:
        year, month, day, hour, minute, datetimes = parseAge(inputData,
                                                             indicator)
        timeElapsed = inputData['age']
    else:
        year, month, day, hour, minute, datetimes = parseDates(inputData,
                                                               indicator,
                                                               inputDateFormat)
        timeElapsed = getTimeElapsed(indicator, year, month, day, hour, minute)

    # Time between observations:
    dt = getTimeDelta(year, month, day, hour, minute)

    # Calculate julian days
    jdays = julianDays(year, month, day, hour, minute)

    lat = np.array(inputData['lat'], 'd')
    lon = np.mod(np.array(inputData['lon'], 'd'), 360)
    delta_lon = np.diff(lon)
    delta_lat = np.diff(lat)

    # Split into separate tracks if large jump occurs (delta_lon > 10 degrees
    # or delta_lat > 5 degrees)
    # This avoids two tracks being accidentally combined when seasons and track
    # numbers match but basins are different as occurs in the IBTrACS dataset.
    # This problem can also be prevented if the 'tcserialno' column is
    # specified.
    indicator[np.where(delta_lon > 10)[0] + 1] = 1
    indicator[np.where(delta_lat > 5)[0] + 1] = 1

    pressure = filterPressure(np.array(inputData['pressure'], 'd'),
                              inputPressureUnits, missingValue)
    try:
        windspeed = np.array(inputData['vmax'], 'd')
        novalue_index = np.where(windspeed == sys.maxsize)
        windspeed = metutils.convert(windspeed, inputSpeedUnits, "mps")
        windspeed[novalue_index] = missingValue
    except (ValueError, KeyError):
        LOG.debug("No max wind speed data - all values will be zero")
        windspeed = np.zeros(indicator.size, 'f')
    assert lat.size == indicator.size
    assert lon.size == indicator.size
    assert pressure.size == indicator.size

    try:
        rmax = np.array(inputData['rmax'])
        novalue_index = np.where(rmax == missingValue)
        rmax = metutils.convert(rmax, inputLengthUnits, "km")
        rmax[novalue_index] = missingValue

    except (ValueError, KeyError):
        LOG.debug("No radius to max wind data - all values will be zero")
        rmax = np.zeros(indicator.size, 'f')

    if 'penv' in inputData.dtype.names:
        penv = np.array(inputData['penv'], 'd')
    else:
        LOG.debug("No ambient MSLP data in this input file")
        LOG.debug("Sampling data from MSLP data defined in "
                  "configuration file")
        # Warning: using sampled data will likely lead to some odd behaviour
        # near the boundary of the MSLP grid boundaries - higher resolution
        # MSLP data will decrease this unusual behaviour.

        try:
            ncfile = cnfGetIniValue(configFile, 'Input', 'MSLPFile')
        except:
            LOG.exception("No input MSLP file specified in configuration")
            raise

        try:
            ncvar = cnfGetIniValue(configFile, 'Input', 'MSLPVariableName')
        except:
            LOG.debug("Using default variable name of 'slp' for sea level pressure data")
            ncvar = 'slp'

        time = getTime(year, month, day, hour, minute)
        penv = ltmPressure(jdays, time, lon, lat, ncfile, ncvar)

    if 'poci' in inputData.dtype.names:
        poci = np.array(inputData['poci'], 'd')
    else:
        LOG.debug("Determining poci")
        eps = np.random.normal(0, scale=2.5717)
        poci = getPoci(penv, pressure, lat, jdays, eps)


    speed, bearing = getSpeedBearing(indicator, lon, lat, dt,
                                     missingValue=missingValue)

    if calculateWindSpeed:
        windspeed = maxWindSpeed(indicator, dt, lon, lat, pressure, poci)

    TCID = np.cumsum(indicator)

    data = np.empty(len(indicator),
                    dtype={
                        'names': trackFields,
                        'formats': trackTypes
                    })
    for key, value in zip(trackFields,
                          [indicator, TCID, year, month,
                           day, hour, minute, timeElapsed,
                           datetimes, lon, lat, speed, bearing,
                           pressure, windspeed, rmax, poci]):
        data[key] = value

    tracks = []
    n = np.max(TCID)
    for i in range(1, n + 1):
        track = Track(data[TCID == i])
        track.trackId = (i, n)
        track.trackfile = trackFile
        getMinPressure(track, missingValue)
        getMaxWind(track, missingValue)
        tracks.append(track)

    return tracks
