#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

 Title: interpTrack.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2007-10-25 9:51:AM
 Description: Interpolates an historic TC track to a shorter time
 interval.  As an example, the csv data for TC Tracy is included in at
 the bottom of this file.
 
 Version :$Rev: 512 $

 ModifiedBy: Nicholas Summons
 ModifiedDate: 2011-02-08
 Modification: Replaced lat/lon interpolation with different spline 
               algorithm that exhibits less overshoot.
               Pressure and Rmax interpolation was replaced with linear
               interpolation.

 $Id: interpTrack.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, sys, pdb, logging, getopt

from scipy import interpolate
import numpy
import datetime
import pylab

from files import flConfigFile, flModuleName, flSaveFile, flStartLog
from config import cnfGetIniValue, logger
from columns import colReadCSV
from grid import SampleGrid


import maputils
import metutils

__version__ = '$Id: interpTrack.py 512 2011-10-31 07:20:38Z nsummons $'

def usage():
    print "interpTrack.py:"
    print "Interpolate the observed points of a tropical cyclone temporally"
    print "for use in modelling a scenario event in TCRM"
    print "Example usage:"
    print "interpTrack.py -c <config file> -l <log file> -v"
    print ""
    print "Options:"
    print "-h, --help:   prints this help message"
    print "-c, --config: configuration path (default value is interpTrack.ini)"
    print "-l --logfile: path to log file to record actions (default value is interpTrack.log)"
    print "-v --verbose: True|False - print all logging messages to the screen"
    print ""
    print "Created by Craig Arthur, 2007-10-25 9:51:AM"
    print __version__

def main(argv):
    "Main part of the program"
    gConfigFile = flConfigFile()
    logFIle = flConfigFile(".log")
    verbose = False
    logger = logging.getLogger()

    try:
        opts, args = getopt.getopt(argv,"hc:l:v",["help","config=","logfile=","verbose",])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt,arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(2)
        elif opt in ("-c", "--config"):
            gConfigFile = arg
        elif opt in ("-l","--logfile"):
            logFile = arg
        elif opt in ("-v","--verbose"):
            verbose = True
    #pdb.set_trace()

    flStartLog(cnfGetIniValue(gConfigFile, 'Logging', 'LogFile', flConfigFile('.log')),
               cnfGetIniValue(gConfigFile, 'Logging', 'LogLevel', 'INFO'),
               cnfGetIniValue(gConfigFile, 'Logging', 'Verbose', False))


    inputFile = cnfGetIniValue(gConfigFile, 'Input', 'File')
    logger.info("Processing %s"%inputFile)
    source = cnfGetIniValue(gConfigFile, 'Input', 'Source')
    inputData = colReadCSV(gConfigFile,inputFile, source)
    #pdb.set_trace()
    logger.info("Extracting date/time information")
    if inputData.has_key('index'):
        indicator = numpy.array(inputData['index'], 'i')
    # Sort date/time information
    if inputData.has_key('date'):
        logger.debug("File contains date information")
        year = numpy.empty(len(indicator), 'i')
        month = numpy.empty(len(indicator), 'i')
        day = numpy.empty(len(indicator), 'i')
        hour = numpy.empty(len(indicator), 'i')
        minute = numpy.empty(len(indicator), 'i')
        datefmt = cnfGetIniValue(gConfigFile, source, 'DateFormat', '%Y-%m-%d %H:%M:%S')
        for i in range(len(inputData['date'])):
            d = datetime.datetime.strptime(inputData['date'][i], datefmt)
            year[i] = d.year
            month[i] = d.month
            day[i] = d.day
            hour[i] = d.hour
            minute[i] = d.minute
    else:
        # Sort out date/time information:
        logger.debug("File contains year/month/day/hour information")
        month = numpy.array(inputData['month'], 'i')
        day = numpy.array(inputData['day'], 'i')
        hour = numpy.array(inputData['hour'], 'i')
        try:
            year = numpy.array(inputData['year'], 'i')
        except:
            # Create dummy variable year - applicable for datasets
            # such as WindRiskTech which contain no year information.
            year = 2000*numpy.ones(indicator.size, 'i')
            ii = numpy.where(month < 7)
            year[ii] = 2001
        try:
            minute = numpy.array(inputData['minute'], 'i')
            assert minute.size == hour.size
        except KeyError:
            # Create dummy variable minute:
            logger.warning("Missing minute data from input data - setting minutes to 00 for all times")
            minute = numpy.zeros((hour.size), 'i')

    # Create the dummy variable second for use in function datenum
    second = numpy.zeros((hour.size), 'i')

    # Time between observations:
    day_ = [datetime.datetime(year[i], month[i], day[i], hour[i],minute[i], second[i]) for i in xrange(year.size)]
    time_ = pylab.date2num(day_)
    dt_ = 24.0*numpy.diff(time_)
    dt = numpy.empty(hour.size, 'f')
    dt[1:] = dt_

    # At this stage, convert all times to a time after initial observation:
    timestep = 24.0*(time_ - time_[0])
    #pdb.set_trace()

    lat = numpy.array(inputData['lat'], 'd')
    lon = numpy.array(inputData['lon'], 'd')
    pressure = numpy.array(inputData['pressure'], 'd')
    rmax = numpy.array(inputData['rmax'],'d')
    bear_, dist_ = maputils.latLon2Azi(lat, lon, 1, azimuth=0)
    bear = numpy.empty(hour.size, 'f')
    bear[1:] = bear_
    bear[0] = bear_[1]
    dist = numpy.empty(hour.size, 'f')
    dist[1:] = dist_
    speed = dist/dt
    #pdb.set_trace()
    if inputData.has_key('penv'):
        penv = numpy.array(inputData['penv'], 'd')
    else:
        logger.info("No ambient MSLP data in this input file")
        logger.info("Sampling data from MSLP data defined in configuration file")
        # Warning: using sampled data will likely lead to some odd behaviour near the boundary of the
        # MSLP grid boundaries - higher resolution MSLP data will decrease this unusual behaviour.
        penv = numpy.zeros(len(lon))
        mslp = SampleGrid(cnfGetIniValue(gConfigFile, 'Input', 'MSLPGrid'))
        for i in range(len(lon)):
            penv[i] = mslp.sampleGrid(lon[i], lat[i])

    delta = cnfGetIniValue(gConfigFile,'Output','Delta',0.1)
    newtime = numpy.arange(timestep[0], timestep[-1], delta)
    nid = numpy.ones(newtime.size)

    logger.info("Interpolating data...")
    nLon = interpolate.splev(newtime, interpolate.splrep(timestep, lon, s=0), der=0)
    nLat = interpolate.splev(newtime, interpolate.splrep(timestep, lat, s=0), der=0)    
    #nvFm = interpolate.spline(timestep, speed, newtime,kind='smoothest')
    #nthetaFm = interpolate.spline(timestep, bear, newtime,kind='smoothest')
    #nthetaFm = numpy.mod(nthetaFm,360)
    npCentre = interpolate.interp1d(timestep, pressure, kind='linear')(newtime)
    npEnv = interpolate.interp1d(timestep, penv, kind='linear')(newtime)
    nrMax = interpolate.interp1d(timestep, rmax, kind='linear')(newtime)
    
    bear_, dist_ = maputils.latLon2Azi(nLat, nLon, 1, azimuth=0)
    nthetaFm = numpy.empty(newtime.size, 'f')
    nthetaFm[1:] = bear_
    dist = numpy.empty(newtime.size, 'f')
    dist[1:] = dist_
    nvFm = dist/delta
    #pdb.set_trace()

    newTrack = numpy.transpose(numpy.concatenate(([nid], [newtime], [nLon], [nLat], [nthetaFm],
                                      [nvFm], [npCentre], [npEnv], [nrMax]), axis=0))
    header=''
    outputFile = cnfGetIniValue(gConfigFile,'Output','File')
    logger.info("Saving interpolated data to %s"%(outputFile))
    flSaveFile(outputFile,newTrack,header,',', fmt='%f')
    logger.info("Completed %s"%(sys.argv[0]))


#####################################################
"""
def interpTrack(cfgFile):
    ""interpTrack(cfgFile):
    Interpolates TC track information to the time resolution given in
    fractions of an hour.  cfgFile: configuration file containing
    details of the track file to interpolate
    ""
    cfg = myutils.loadConfig(open(cfgFile))
    delta = cfg['Timestep.delta']
    trackData = myutils.textread(cfgFile)
    iD = array(trackData['cycloneId'], int)
    Age = array(trackData['cycloneAge'], float)
    Lon = array(trackData['cLon'], float)
    Lat = array(trackData['cLat'], float)
    vFm = array(trackData['vFm'], float)
    thetaFm = array(trackData['thetaFm'])
    pCentre = array(trackData['pCentre'], float)
    pEnv = array(trackData['pEnv'], float)
    rMax = array(trackData['rMax'], float)
    year = array(trackData['Yr'], int)
    month = array(trackData['Mon'], int)
    day = array(trackData['Day'], int)
    hour = array(trackData['Hr'], int)
    min = array(trackData['Min'], int)

    datestamp = [datetime.datetime(year[i], month[i], day[i], hour[i], min[i])
                 for i in xrange(year.size)]
    time = 24.0*pylab.date2num(datestamp)      # Time in hours since 0001-01-01 00:00:00UTC
    newtime = arange(time[0], time[-1], delta)

    nid = interp(newtime, time, iD)
    nage = interp(newtime, time, Age)
    nLon = interp(newtime, time, Lon)
    nLat = interp(newtime, time, Lat)
    nvFm = interp(newtime, time, vFm)
    nthetaFm = interp(newtime, time, thetaFm)
    npCentre = interp(newtime, time, pCentre)
    npEnv = interp(newtime, time, pEnv)
    nrMax = interp(newtime, time, rMax)

    newTrack = transpose(concatenate(([nid], [nage], [nLon], [nLat], [nthetaFm],
                                      [nvFm], [npCentre], [npEnv], [nrMax]), axis=0))
    header=''
    myutils.save(cfg['File.Output'], newTrack,header,',', fmt='%f')

"""

# Call the main program:
if __name__ == "__main__":
    main(sys.argv[1:])

"""
Test data for interpTrack:
ID, Age, Lon, Lat, Bearing,Speed, central pressure, Environmental Pressure, rMax, Year, Month, Day, Hour, Minute, Seconds
1,0,132.1,-9.2,360,0,99000,100800,15,1974,12,21,11,30,0
1,6.5,131.5,-9.6,236,12.2,99000,100800,15,1974,12,21,18,0,0
1,12,131.2,-9.9,224.6,8.5,99000,100800,12,1974,12,21,23,30,0
1,18.5,130.9,-10.5,206.2,11.4,98500,100800,12,1974,12,22,6,0,0
1,24,130.6,-10.8,224.5,8.5,98500,100800,12,1974,12,22,11,30,0
1,28.5,130.5,-11.1,198.1,7.8,98500,100800,12,1974,12,22,16,0,0
1,36,130.3,-11.2,243,3.3,98000,100800,12,1974,12,22,23,30,0
1,48,129.9,-11.6,224.5,5.2,97500,100800,12,1974,12,23,11,30,0
1,60,129.9,-11.9,180,2.8,97000,100800,12,1974,12,23,23,30,0
1,64,130.1,-12,117,6.1,96500,100800,10,1974,12,24,3,30,0
1,69.5,130.3,-12.2,135.6,5.7,96000,100800,10,1974,12,24,9,0,0
1,72,130.5,-12.2,90,8.7,95500,100800,10,1974,12,24,11,30,0
1,75.5,130.6,-12.3,135.6,4.4,95000,100800,10,1974,12,24,15,0,0
1,78,130.8,-12.4,117.1,9.8,95000,100800,8,1974,12,24,17,30,0
1,78.5,130.8,-12.4,117.1,0,95000,100800,8,1974,12,24,18,0,0
1,84,131.3,-12.5,101.5,10.1,95000,100800,8,1974,12,24,23,30,0
1,87,131.6,-12.6,108.8,11.5,97000,100800,8,1974,12,25,2,30,0

"""
