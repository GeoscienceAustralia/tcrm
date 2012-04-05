#!/usr/bin/env python
"""
 Title: loadData.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2011-10-27 2:21:PM
 Description: Load a csv-format file of TC track data and return a series
              of arrays containing the data. Expected fields which have no
              data are replaced with numpy.zero arrays matching the length
              of the other fields.

 Version :$Rev: 686 $

 $Id: loadData.py 686 2012-03-29 04:24:59Z carthur $
"""

import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

import numpy
from datetime import datetime
import metutils
import maputils
from matplotlib.dates import date2num, num2date
from columns import colReadCSV
from config import cnfGetIniValue
from grid import SampleGrid

logger = logging.getLogger()

def loadTrackFile(configFile, trackFile, source, missingValue=0):
    """
    Load TC track data from the given input file, from a specified source.
    The configFile is a configuration file that contains a section called
    'source' that describes the data.
    This returns a series of arrays containing the data. See the return line
    """
    logger.debug("Loading %s"%trackFile)
    inputData = colReadCSV(configFile, trackFile, source, nullValue=missingValue )
    inputSpeedUnits = cnfGetIniValue(configFile, source, 'SpeedUnits', 'mps')
    inputPressureUnits = cnfGetIniValue(configFile, source, 'PressureUnits', 'hPa')
    inputLengthUnits = cnfGetIniValue(configFile, source, 'LengthUnits', 'km')

    # Determine the initial TC positions...
    if 'index' in inputData:
        logger.debug("Using index contained in file to determine initial TC positions")
        indicator = numpy.array(inputData['index'], 'i')
    else:
        if 'tcserialno' in inputData:
            tcSerialNo = inputData['tcserialno']
            indicator = numpy.ones(len(tcSerialNo), 'i')
            for i in xrange(1, len(tcSerialNo)):
                if tcSerialNo[i] == tcSerialNo[i-1]:
                    indicator[i] = 0
        elif ('season' in inputData) and ('num' in inputData):
            logger.debug("Using season and TC number to determine initial TC positions")
            num = numpy.array(inputData['num'], 'i')
            season = numpy.array(inputData['season'], 'i')
            indicator = numpy.ones(num.size, 'i')
            for i in xrange(1, len(num)):
                if season[i] == season[i-1] and num[i] == num[i-1]:
                    indicator[i] = 0
        elif 'num' in inputData:
            logger.debug("Using TC number to determine initial TC positions (no season information)")
            num = numpy.array(inputData['num'],'i')
            indicator = numpy.ones(num.size,'i')
            ind_ = numpy.diff(num)
            ind_[numpy.where(ind_ > 0)] = 1
            indicator[1:] = ind_
        else:
            logger.critical("Insufficient input file columns have been specified to run TCRM.")
            sys.exit(2)
    # Sort date/time information
    if 'age' in inputData:
        dt_ = numpy.diff(inputData['age'])
        dt = numpy.empty(indicator.size, 'f')
        dt[1:] = dt_
        start_time = date2num(datetime(2000,1,1,0,0))
        times_ = start_time + inputData['age']/24.
        d = num2date(times_)
        year = 2000*numpy.ones(indicator.size, 'i')
        month = numpy.ones(indicator.size,'i')
        day = numpy.ones(indicator.size,'i')
        hour = 12*numpy.ones(indicator.size,'i')
        minute = numpy.zeros(indicator.size,'i')
        for i in xrange(len(d)):
            year[i] = d[i].year
            month[i] = d[i].month
            day[i] = d[i].day
            hour[i] = d[i].hour
            minute[i] = d[i].minute
    else:
        if 'date' in inputData:
            year = numpy.empty(len(indicator), 'i')
            month = numpy.empty(len(indicator), 'i')
            day = numpy.empty(len(indicator), 'i')
            hour = numpy.empty(len(indicator), 'i')
            minute = numpy.empty(len(indicator), 'i')
            datefmt = cnfGetIniValue(configFile, source, 'DateFormat', '%Y-%m-%d %H:%M:%S')
            for i in xrange(len(inputData['date'])):
                try:
                    d = datetime.strptime(inputData['date'][i], datefmt)
                except ValueError:
                    logger.critical("Error in date information for record %d"%i)
                    logger.critical(sys.exc_info()[1])
                    logger.critical("Check the date/time records in your input file")

                year[i] = d.year
                month[i] = d.month
                day[i] = d.day
                hour[i] = d.hour
                minute[i] = d.minute
        else:
            # Sort out date/time information:
            month = numpy.array(inputData['month'], 'i')
            day = numpy.array(inputData['day'], 'i')
            hour = numpy.array(inputData['hour'], 'i')
            try:
                year = numpy.array(inputData['year'], 'i')
            except:
                # Create dummy variable year - applicable for datasets
                # such as WindRiskTech which contain no year information.
                year = numpy.zeros(indicator.size, 'i')
                for i in range(len(year)):
                    if xindicator[i] > 0:
                        fill_year = 2000
                    if month[i] == 1:
                        fill_year = 2001
                    year[i] = fill_year

            try:
                minute = numpy.array(inputData['minute'], 'i')
                assert minute.size == indicator.size
            except KeyError:
                # Create dummy variable minute:
                logger.debug("Missing minute data from input data - setting minutes to 00 for all times")
                minute = numpy.zeros((hour.size), 'i')

        # Create the dummy variable second for use in function datenum
        second = numpy.zeros((hour.size), 'i')

        # Time between observations:
        try:
            day_ = [datetime(year[i], month[i], day[i], hour[i], minute[i], second[i]) for i in xrange(year.size)]
        except ValueError:
            logger.critical("Error in date information")
            logger.critical(sys.exc_info()[1])
            logger.critical("Check the date/time records in your input file")
            sys.exit(2)
        try:
            time_ = date2num(day_)
        except ValueError:
            logger.critical("Error in day values")
            logger.critical(sys.exc_info()[1])
            logger.critical("Check the date/time records in your input file")
            sys.exit(2)

        dt_ = 24.0*numpy.diff(time_)
        dt = numpy.empty(indicator.size, 'f')
        dt[1:] = dt_
        # Calculate julian days - set all years prior to 1900 to 1904
        # strftime() requires year >=1900; and
        # in the Gregorian calendar, 1900 is not a leap year (and there are many years prior to 1900
        # that are!).
        jyear = numpy.copy(year)
        jyear[numpy.where(jyear<1900)]=1904
        jday_ = [datetime(jyear[i], month[i], day[i], hour[i], minute[i], second[i]) for i in xrange(jyear.size)]
        jdays = numpy.array([int(jday_[i].strftime("%j")) for i in xrange(jyear.size)])

    lat = numpy.array(inputData['lat'], 'd')
    lon = numpy.mod(numpy.array(inputData['lon'], 'd'), 360)
    delta_lon = numpy.diff(lon)
    delta_lat = numpy.diff(lat)
    # Split into separate tracks if large jump occurs (delta_lon > 10 degrees or delta_lat > 5 degrees)
    # This avoids two tracks being accidentally combined when seasons and track numbers match but
    # basins are different as occurs in the IBTrACS dataset.  This problem can also be prevented if
    # the 'tcserialno' column is specified.
    indicator[numpy.where(delta_lon > 10)[0] + 1] = 1
    indicator[numpy.where(delta_lat > 5)[0] + 1] = 1

    # Save information required for frequency auto-calculation
    if 'season' in inputData:
        origin_seasonOrYear = numpy.array(inputData['season'], 'i').compress(indicator)
        header = 'Season'
    else:
        origin_seasonOrYear = year.compress(indicator)
        header = 'Year'
    origin_lon = lon.compress(indicator)
    origin_lat = lat.compress(indicator)

    pressure = numpy.array(inputData['pressure'], 'd')
    novalue_index = numpy.where(pressure==sys.maxint)
    pressure = metutils.convert(pressure, inputPressureUnits, "hPa")
    pressure[novalue_index] = sys.maxint

    # Convert any non-physical central pressure values to maximum integer
    # This is required because IBTrACS has a mix of missing value codes (i.e. -999, 0, 9999)
    # in the same global dataset.
    pressure = numpy.where((pressure < 600) | (pressure > 1100), sys.maxint, pressure)
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
        # Warning: using sampled data will likely lead to some odd behaviour near the boundary of the
        # MSLP grid boundaries - higher resolution MSLP data will decrease this unusual behaviour.
        penv = numpy.zeros( len( lon ) )
        mslp = SampleGrid( cnfGetIniValue( configFile, 'Input', 'MSLPGrid' ) )
        for i in xrange( len( lon ) ):
            penv[i] = mslp.sampleGrid(lon[i], lat[i])

    # ieast : parameter used in latLon2Azi --> should be a config
    # setting describing the input data.
    ieast = 1

    # Determine the index of initial cyclone observations, excluding
    # those cyclones that have only one observation. This is used
    # for calculating initial bearing and speed
    indicator2 = numpy.where(indicator > 0, 1, 0)   # ensure indicator is only ones and zeros
    initIndex = numpy.concatenate([numpy.where(numpy.diff(indicator2) == -1, 1, 0), [0]])

    # Calculate the bearing and distance (km) of every two
    # consecutive records using ll2azi
    bear_, dist_ = maputils.latLon2Azi(lat, lon, ieast, azimuth=0)
    assert bear_.size == indicator.size - 1
    assert dist_.size == indicator.size - 1
    bearing = numpy.zeros(indicator.size, 'f')
    bearing[1:] = bear_
    numpy.putmask(bearing, indicator, sys.maxint)

    dist = numpy.zeros(indicator.size, 'f')
    dist[1:] = dist_
    speed = dist/dt
    # Delete speeds less than 0, greated than 200,
    # or where indicator == 1.
    numpy.putmask(speed, (speed < 0) | (speed > 200) | indicator, sys.maxint)
    numpy.putmask(speed, numpy.isnan(speed), sys.maxint)


    return indicator,year,month,day,hour,minute,lon,lat,pressure,speed,bearing,windspeed,rmax,penv
