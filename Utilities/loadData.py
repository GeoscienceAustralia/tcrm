#!/usr/bin/env python
"""
 Title: loadData.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2011-10-27 2:21:PM
 Description: Load a csv-format file of TC track data and return a series
              of arrays containing the data. Expected fields which have no
              data are replaced with numpy.zero arrays matching the length
              of the other fields.

 Version :$Rev$

 $Id$
"""

#import os
import sys
#import pdb
import logging

import numpy
from datetime import datetime
import metutils
import maputils
import nctools
import interp3d
from matplotlib.dates import date2num, num2date
from columns import colReadCSV
#from config import cnfGetIniValue
#from grid import SampleGrid
from config import ConfigParser

__version__ = "$Id$"

logger = logging.getLogger(__name__)

def getSpeedBearing(index, lon, lat, deltatime, ieast=1):
    """
    Calculate the speed and bearing of a TC based on the position and time
    difference between observations.
    """
    bear_, dist_ = maputils.latLon2Azi(lat, lon, ieast, azimuth=0)
    assert bear_.size == index.size - 1
    assert dist_.size == index.size - 1
    bearing = numpy.zeros(index.size, 'f')
    bearing[1:] = bear_
    numpy.putmask(bearing, index, sys.maxint)

    dist = numpy.zeros(index.size, 'f')
    dist[1:] = dist_
    speed = dist/deltatime
    # Delete speeds less than 0, greated than 200,
    # or where indicator == 1.
    numpy.putmask(speed, (speed < 0) | (speed > 200) | index, sys.maxint)
    numpy.putmask(speed, numpy.isnan(speed), sys.maxint)

    return speed, bearing

def maxWindSpeed(index, deltatime, lon, lat, pressure, penv,
                 gustfactor=0.9524):
    """
    Calculate the 10-minute-mean maximum wind speed from the central
    pressure deficit, using the method described in Holland et al. (2010)

    Input: indicator - array (values of 1 or 0) indicating the beginning of
               a new TC in the input dataset;
           deltatime - time difference (in hours) between each point in the
               record;
           lon - array of longitudes of TC positions
           lat - array of latitides of TC positions
           pressure - central pressure estimates of TCs
           penv - environmental pressure values
           gf - gust factor - default value represents converting from a
               1-minute sustained wind speed to a 10-minute mean wind speed.
               Based on Harper et al. 2010, WMO-TD1555

    Output: array of estimated wind speed based on estimated central pressure.

    Example: v = maxWindSpeed(indicator, dt, lon, lat, pressure, penv)

    """

    # Speed and bearing:
    speed, bearing = getSpeedBearing(index, lon, lat, deltatime)
    speed = metutils.convert(speed, 'kmh', 'mps')

    # Pressure deficit:
    deltap = penv - pressure

    # Pressure rate of change
    dpt = numpy.zeros(index.size,'f')
    dpt[1:] = numpy.diff(pressure)
    dpdt = dpt/deltatime
    numpy.putmask(dpdt, index, 0)
    numpy.putmask(dpdt, numpy.isnan(dpdt) | numpy.isinf(dpdt) , 0)


    # Estimated pressure at the radius of maximum wind:
    prmw = pressure + deltap/3.7

    # Calculate thermodynamic variables at RMW:
    tsurf = 28.0 - 3*(numpy.abs(lat) - 10.)/20.
    qmix = 0.9*(3.802/prmw)*numpy.exp(17.67*tsurf/(243.5 + tsurf))
    tvs = (tsurf + 273.15)*(1. + 0.81*qmix)
    rho = prmw*100./(tvs*287.04)

    chi = 0.6*(1.0 - deltap/215.)
    beta = -0.000044*numpy.power(deltap, 2.) + \
             0.01*deltap + 0.03*dpdt - 0.014*numpy.abs(lat) + \
             0.15*numpy.power(speed, chi)+ 1.

    # Holland's P-W relation derives a 1-minute mean wind speed, so we often
    # need to convert to some other averaging period. I use the recommendations
    # of Harper et al. (2010) WMO TD-1555:
    # Common values are( Assuming "At-sea" conditions):
    # 10-min mean: 0.95 (default)
    # 3-second gust: 1.11

    v = gustfactor*numpy.sqrt(deltap*100*beta/(rho*numpy.exp(1.)))

    numpy.putmask(v, (numpy.isnan(v) |
                      numpy.isinf(v) |
                      (pressure>=sys.maxint) ), 0)

    return v

def getInitialPositions(data):
    """
    Get the array indices corresponding to the initial position of TCs in
    the input dataset. This is done through examining the data for a number
    of specific fields to see when they change, or if the data has a field
    that indicates such an instance.
    Input:
        data - dict of arrays that contains the data loaded from the input
                file

    Output:
        numpy.ndarray of indexes that can be used to slice the observations
        and return those corresponding to an initial TC position

    NOTE: using only the 'num' field will result in different results than
    using 'num' and 'season'. From experience, the 'num' field in best-track
    datasets refers to the sequential number of the storm for that season -
    i.e. it starts at 1 for each season and increments for each new storm.

    The use of 'num' only should be reserved for those situations where the
    dataset is known to have unique numbers for each storm (e.g in simulated
    data).

    """

    if data.has_key('index'):
        logger.debug("Using index contained in file to determine \
                            initial TC positions")
        indicator = numpy.array(data['index'], 'i')
    else:
        if data.has_key('tcserialno'):
            logger.debug("Using TC serial number to determine \
                                initial TC positions")
            tcSerialNo = data['tcserialno']
            indicator = numpy.ones(len(tcSerialNo), 'i')
            for i in range(1, len(tcSerialNo)):
                if tcSerialNo[i] == tcSerialNo[i-1]:
                    indicator[i] = 0
        elif data.has_key('season') and data.has_key('num'):
            logger.debug("Using season and TC number to determine \
                                initial TC positions")
            num = numpy.array(data['num'], 'i')
            season = numpy.array(data['season'], 'i')
            indicator = numpy.ones(num.size, 'i')
            for i in range(1, len(num)):
                if (season[i] == season[i-1]) and (num[i] == num[i-1]):
                    indicator[i] = 0

        elif data.has_key('num'):
            logger.debug("Using TC number to determine initial \
                                TC positions (no season information)")
            num = numpy.array(data['num'], 'i')
            indicator = numpy.ones(num.size, 'i')
            ind_ = numpy.diff(num)
            ind_[numpy.where(ind_ != 0)] = 1
            indicator[1:] = ind_
        else:
            raise KeyError ("Insufficient input file columns have \
                                    been specified to run TCRM.")

    return indicator

def date2ymdh(dates, datefmt='%Y-%m-%d %H:%M:%S'):

    import re
    pattern = re.compile("%y")
    if pattern.search(datefmt):
        raise ValueError("Cannot use 2-digit year formats in date format")

    year = numpy.empty(len(dates), 'i')
    month = numpy.empty(len(dates), 'i')
    day = numpy.empty(len(dates), 'i')
    hour = numpy.empty(len(dates), 'i')
    minute = numpy.empty(len(dates), 'i')

    for i in xrange(len(dates)):
        try:
            d = datetime.strptime(dates[i], datefmt)
        except ValueError:
            raise ValueError("Error in date information for record %d"%i)
        else:
            year[i] = d.year
            month[i] = d.month
            day[i] = d.day
            hour[i] = d.hour
            minute[i] = d.minute

    return year, month, day, hour, minute

def parseDates(data, indicator, datefmt='%Y-%m-%d %H:%M:%S'):
    """
    Parse the date/time information to extract year, month, day, hour and
    minute details for the input dataset
    """
    if data.has_key('date'):
        year, month, day, hour, minute = date2ymdh(data['date'], datefmt)

    else:
        # Sort out date/time information:
        month = numpy.array(data['month'], 'i')
        day = numpy.array(data['day'], 'i')
        hour = numpy.array(data['hour'], 'i')
        try:
            year = numpy.array(data['year'], 'i')
        except KeyError:
            # Create dummy variable year - applicable for datasets
            # such as WindRiskTech which contain no year information.
            year = numpy.zeros(month.size, 'i')
            for i in range(len(year)):
                if indicator[i] > 0:
                    fill_year = 2000
                if month[i] == 1:
                    fill_year = 2001
                year[i] = fill_year

        if data.has_key('minute'):
            minute = numpy.array(data['minute'], 'i')
        elif hour.max() >= 100:
            minute = numpy.mod(hour, 100)
            hour = hour/100
        else:
            logger.warning("Missing minute data from input data - \
                                    setting minutes to 00 for all times")
            minute = numpy.zeros((hour.size), 'i')


    return year, month, day, hour, minute

def parseAge(data, indicator):
    """
    Parse the TC age information to get a proxy date record. Assumes every TC
    starts at 2000-01-01 00:00 and calculates year, month, day, hour and
    minute values based on the age field.
    """

    start_time = date2num(datetime(2000, 1, 1, 0, 0))
    times_ = start_time + data['age']/24.
    d = num2date(times_)
    year = 2000*numpy.ones(indicator.size, 'i')
    month = numpy.ones(indicator.size, 'i')
    day = numpy.ones(indicator.size, 'i')
    hour = 12*numpy.ones(indicator.size, 'i')
    minute = numpy.zeros(indicator.size, 'i')

    for i in xrange(len(d)):
        year[i] = d[i].year
        month[i] = d[i].month
        day[i] = d[i].day
        hour[i] = d[i].hour
        minute[i] = d[i].minute

    return year, month, day, hour, minute

def getTimeDelta(year, month, day, hour, minute):
    """
    Calculate the time difference between consecutive observations

    Input:
        year - numpy.ndarray of the year of all observations
        month - as for year, but for the month of observation
        day - as for year, but for the day of observation
        hour - as for year, but for the hour of the observation
        minutes - as for year, but for the hour of the observation
        seconds - as for year, but for the hour of the observation

    Output:
        dt - numpy.ndarray of time difference between observations in hours
    """
    if numpy.any(year < 0):
        raise ValueError("Error in input year information - check input file")
    if numpy.any(month >= 13):
        raise ValueError("Error in input month information - check input file")
    if numpy.any(day > 31):
        raise ValueError("Error in input day information - check input file")
    if numpy.any(hour > 24):
        raise ValueError("Error in input hour information - check input file")
    if numpy.any(minute > 60):
        raise ValueError("Error in input minute information - check input file")

    logger.debug("Calculating time difference between observations")
    second = numpy.zeros((hour.size), 'i')
    try:
        day_ = [datetime(year[i], month[i], day[i], hour[i],
                                  minute[i], second[i])
                for i in xrange(year.size)]
    except ValueError:
        raise ValueError("Error in date information - check your input file")


    try:
        time_ = date2num(day_)
    except ValueError:
        raise ValueError("Error in day values - check your input file")

    dt = numpy.zeros(year.size, 'f')
    dt[1:] = 24.0*numpy.diff(time_)

    return dt

def getTime(year, month, day, hour, minute):
    """
    Calculate the number of days since 0001-01-01 00:00:00 UTC + 1

    """

    logger.debug("Calculating time in hours since epoch")
    second = numpy.zeros((hour.size), 'i')

    if numpy.any(year < 0):
        raise ValueError("Error in input year information - check input file")
    if numpy.any(month >= 13):
        raise ValueError("Error in input month information - check input file")
    if numpy.any(day > 31):
        raise ValueError("Error in input day information - check input file")
    if numpy.any(hour > 24):
        raise ValueError("Error in input hour information - check input file")
    if numpy.any(minute > 60):
        raise ValueError("Error in input minute information - check input file")

    try:
        day_ = [datetime(year[i], month[i], day[i], hour[i],
                                  minute[i], second[i])
                for i in xrange(year.size)]
    except ValueError:
        raise ValueError("Error in date information - check input file")

    try:
        time = date2num(day_)
    except ValueError:
        raise ValueError("Error in day values - check input file")

    else:
        return numpy.array(time,'f')


def julianDays(year, month, day, hour, minute):
    """
    Calculate the julian day (day of year) based on the known date/time
    information
    """
    logger.debug("Calculating julian day (day of year) values")

    if numpy.any(year < 0):
        raise ValueError("Error in input year information - check input file")
    if numpy.any(month >= 13):
        raise ValueError("Error in input month information - check input file")
    if numpy.any(day > 31):
        raise ValueError("Error in input day information - check input file")
    if numpy.any(hour > 24):
        raise ValueError("Error in input hour information - check input file")
    if numpy.any(minute > 60):
        raise ValueError("Error in input minute information - check input file")

    # set all years prior to 1900 to 1904 - strftime() requires year >=1900;
    # and in the Gregorian calendar, 1900 is not a leap year (and there are
    # many years prior to 1900 that are!).
    second = numpy.zeros((hour.size), 'i')
    jyear = numpy.copy(year)
    jyear[numpy.where(jyear<1900)]=1904
    day = [datetime(jyear[i], month[i], day[i], hour[i], minute[i],
                              second[i]) for i in xrange(year.size)]

    jdays = numpy.array([int(day[i].strftime("%j")) for
                            i in xrange(year.size)])
    return jdays

def ltmPressure(jdays, time, lon, lat, ncfile):
    """
    Extract pressure value from a daily long-term mean SLP dataset at the
    given day of year and lon,lat position
    To use this function (and hence some form of daily LTM SLP data) requires
    knowledge of the day of year.
    """
    jtime = jdays + numpy.modf(time)[0]
    coords = numpy.array([jtime, lat, lon])

    logger.debug("Sampling data from MSLP data in {0}".format(ncfile))
    ncobj = nctools.ncLoadFile(ncfile)
    data = nctools.ncGetData(ncobj, 'mslp')
    # Get the MSLP by interpolating to the location of the TC:
    penv = interp3d.interp3d(data, coords)
    del data
    ncobj.close()
    del ncobj

    return penv

def filterPressure(pressure, inputPressureUnits='hPa'):
    """
    Filter pressure values to remove any non-physical values
    """

    novalue_index = numpy.where(pressure==sys.maxint)
    pressure = metutils.convert(pressure, inputPressureUnits, "hPa")
    pressure[novalue_index] = sys.maxint

    # Convert any non-physical central pressure values to maximum integer
    # This is required because IBTrACS has a mix of missing value codes
    # (i.e. -999, 0, 9999) in the same global dataset.
    pressure = numpy.where((pressure < 600) | (pressure > 1100),
                           sys.maxint, pressure)
    return pressure

def loadTrackFile(configFile, trackFile, source, missingValue=0,
                  calculateWindSpeed=True):
    """
    Load TC track data from the given input file, from a specified source.
    The configFile is a configuration file that contains a section called
    'source' that describes the data.
    This returns a series of arrays containing the data. See the return line for the
    common names of the data returned.

    Input:
    configFile: configuration file with a section 'source'
    trackFile:  path to a csv-formatted file containing TC data
    source:     string describing the source format of the TC data. There *must* be
                a section in 'configFile' matching this string, containing the
                details of the format of the data
    missingValue: replace all null values in the input data with this value (default=0)
    calculateWindSpeed: Boolean (default True), calculate maximum wind speed using
                a pressure-wind relation described in maxWindSpeed()

    Output:
    A series of arrays with the required variables for use in TCRM:
        indicator, year, month, day, hour, minute, lon, lat, pressure, speed,
        bearing, windspeed, rmax, penv
    If any of these variables are not present in the input dataset, they
    are (where possible) calculated (date/time/windspeed), sampled from default
    datasets (e.g. environmental pressure) or set to the missing value.

    Example:
    indicator,year,month,day,hour,minute,lon,lat,pressure,speed,bearing,\
    windspeed,rmax,penv = loadTrackFile('tcrm.ini', 'IBTRaCS.csv', 'IBTrACS' )

    """
    logger.debug("Loading %s"%trackFile)
    inputData = colReadCSV(configFile, trackFile, source,
                           nullValue=missingValue)

    config = ConfigParser()
    config.read(configFile)

    inputSpeedUnits = config.get(source, 'SpeedUnits')
    inputPressureUnits = config.get(source, 'PressureUnits')
    inputLengthUnits = config.get(source, 'LengthUnits')
    inputDateFormat = config.get(source, 'DateFormat')

    # Determine the initial TC positions...
    indicator = getInitialPositions(inputData)

    # Sort date/time information
    if 'age' in inputData:
        year, month, day, hour, minute = parseAge(inputData, indicator)
    else:
        year, month, day, hour, minute = parseDates(inputData, indicator,
                                                    inputDateFormat)

    # Time between observations:
    dt = getTimeDelta(year, month, day, hour, minute)

    # Calculate julian days
    jdays = julianDays(year, month, day, hour, minute)


    lat = numpy.array(inputData['lat'], 'd')
    lon = numpy.mod(numpy.array(inputData['lon'], 'd'), 360)
    delta_lon = numpy.diff(lon)
    delta_lat = numpy.diff(lat)

    # Split into separate tracks if large jump occurs (delta_lon > 10 degrees
    # or delta_lat > 5 degrees)
    # This avoids two tracks being accidentally combined when seasons and track
    # numbers match but basins are different as occurs in the IBTrACS dataset.
    # This problem can also be prevented if the 'tcserialno' column is
    # specified.
    indicator[numpy.where(delta_lon > 10)[0] + 1] = 1
    indicator[numpy.where(delta_lat > 5)[0] + 1] = 1

    pressure = filterPressure(numpy.array(inputData['pressure'], 'd'),
                              inputPressureUnits)
    try:
        windspeed = numpy.array(inputData['vmax'], 'd')
        novalue_index = numpy.where(windspeed==sys.maxint)
        windspeed = metutils.convert(windspeed, inputSpeedUnits, "mps")
        windspeed[novalue_index] = sys.maxint
    except KeyError:
        logger.debug("No max wind speed data - all values will be zero")
        windspeed = numpy.zeros(indicator.size, 'f')
    assert lat.size == indicator.size
    assert lon.size == indicator.size
    assert pressure.size == indicator.size
    #assert vmax.size == indicator.size

    try:
        rmax = numpy.array(inputData['rmax'])
        novalue_index = numpy.where(rmax==sys.maxint)
        rmax = metutils.convert(rmax, inputLengthUnits, "km")
        rmax[novalue_index] = sys.maxint

    except KeyError:
        logger.debug("No radius to max wind data - all values will be zero")
        rmax = numpy.zeros( indicator.size, 'f' )

    if 'penv' in inputData:
        penv = numpy.array( inputData['penv'], 'd' )
    else:
        logger.info("No ambient MSLP data in this input file")
        logger.info("Sampling data from MSLP data defined in configuration file")
        # Warning: using sampled data will likely lead to some odd behaviour near
        # the boundary of the MSLP grid boundaries - higher resolution MSLP
        # data will decrease this unusual behaviour.

        try:
            ncfile = cnfGetIniValue( configFile, 'Input', 'MSLPFile')
        except:
            logger.exception("No input MSLP file specified in configuration")
            raise
        time = getTime(year, month, day, hour, minute)
        penv = ltmPressure(jdays, time, lon, lat, ncfile)

    speed, bearing = getSpeedBearing(indicator, lon, lat, dt)

    if calculateWindSpeed:
        windspeed = maxWindSpeed(indicator, dt, lon, lat, pressure, penv)

    return indicator, year, month, day, hour, minute, lon, lat, pressure, \
            speed, bearing, windspeed, rmax, penv
