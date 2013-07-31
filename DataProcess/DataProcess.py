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


Title: DataProcess.py
Author: Craig Arthur, craig.arthur@ga.gov.au

CreationDate: 2006-10-23
Description: Processes the database of historical TC's into
suitably formatted text files. Data is written to plain text files for
ease of access (this may be upgraded in future versions to netCDF format
files). Currently extracts fields containing the value of cyclone
parameters, but no information on the change of parameters.
This is a revision of the original DataProcess class written by
Geoff Xu (2006).

Version: $Rev: 832 $

ModifiedBy: C. Arthur
ModifiedDate: 2006-11-07
Modification: Eliminated convoluted processes to calculate bearing and speed.
              Also removed self.Nan references (replaced with None and a
              check on values before writing them to file).

ModifiedBy: N. Habili
ModifiedDate: 2006-11-13
Modification: Style changes. Upgraded to numpy. Created submethods.

Version: 528
ModifiedBy: C. Arthur
ModifiedDate: 2006-11-15
Modification: Added several submethods to calculate changes in cyclone parameters.

Version: 534
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-10-02 10:13:AM
Modification: Removed need for individual output files to be listed in the config file. Only the path
              to the output directory is needed. The path given is checked for the correct
              path separator at the end and is appended if necessary.

Version: 58
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-10-03 4:01:PM
Modification: Corrected bug which was resulting in spurious rates of bearing change.

Version: 75
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 10/04/08 11:42:AM
Modification: Changed logging method

Version: 108
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 17/09/08 11:34:AM
Modification: Added rMax to values that are extracted. The process will
              first try to extract the values from the inputData dict. If
              a KeyError is returned it's assumed the data is not available,
              issues a warning and moves on.

Version: 161
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-01-20
Modification: Corrected bug in determining initIndex that was returning the
              first (roughly) 1000 indices, rather than indices corresponding to
              the required initial TC locations.

Version: 210
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-04-23 4:03:PM
Modification: Determines a land/sea flag for all TC positions. This is used
              to select obs when generating statistical properties for a
              grid cell.

Version: 277
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-09-11
Modification: Made vMax an optional input parameter (i.e. it's non-essential). The
              updated BoM dataset caused problems because it contains no wind speed
              information, only estimated central pressures.

Version: 292
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-05-15
Modification: Generate a histogram of the frequency of events

Version: $Rev: 832 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-06-24 5:45:PM
Modification: Altered the calculation of dt to make allowance for TCRM data
              to be processsed (for validation purposes).

SeeAlso: (related programs)
Constraints:

$Id: DataProcess.py 832 2012-03-28 07:23:32Z nsummons $
"""

import os, sys, logging
import numpy as np
from os.path import join as pjoin

import Utilities.maputils as maputils
import Utilities.metutils as metutils
import Utilities.stats as stats
import Utilities.loadData as loadData

from Utilities.grid import SampleGrid
from Utilities.files import flModuleName, flSaveFile, flStartLog
from Utilities.columns import colReadCSV
from Utilities import pathLocator
from CalcTrackDomain import CalcTrackDomain
from Utilities.config import ConfigParser

# Switch off minor warning messages
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


class DataProcess:
    """
    Processes the database of historical TCs into suitably formatted
    text files.  Data is written to plain text files for ease of
    access (this may be upgraded in future versions to netCDF files).
    Currently extracts fields containing the value of cyclone
    parameters, but no information on the change of parameters.

    Parameters
    -------
    files: Dictionary containing input and output file names.

    Members
    -------
    :type  configFile: string
    :param configFile: Configuration file containing simulation settings

    :type  progressbar: :class:`progressbar`
    :param progressbar: :attr:`progressbar` object to print progress to 
                        STDOUT

    Methods
    -------
    processdata()
        process initial raw data into ASCII data that is readable by the
        program

    Internal Methods
    -------
    _lonLat(lon, lat, indicator)
        Extract longitudes and latitudes
    _bearing(bear, indicator)
        Extract bearings
    _speed(dist, dt, indicator)
        Extract speeds
    _pressure(pressure, indicator)
        Extract pressures
    _pressureRate(pressure, dt, indicator)
        Extract rate of presssure change
    _speedRate(dist, dt, indicator)
        Extract the acceleration (rate of change of speed)
    _bearingChange(bear, dt, indicator)
        Extract the rate of change of bearing
    _windSpeed(vmax, indicator)
        Extract the maximum sustained wind speed
    """

    def __init__(self, configFile, progressbar=None):
        """
        Initialize the data include tool instance, Nan value and all
        full path names of the files in which data will be stored.
        """
        #CalcTD = CalcTrackDomain(configFile)
        #self.domain = CalcTD.calc()
        self.configFile = configFile
        self.progressbar = progressbar
        self.logger = logging.getLogger(__name__)
        self.logger.info("Initialising DataProcess")

        config = ConfigParser()
        config.read(configFile)

        self.outputPath = config.get('Output', 'Path')
        self.processPath = pjoin(self.outputPath, 'process')

        # Determine TCRM input directory
        tcrm_dir = pathLocator.getRootDirectory()
        self.tcrm_input_dir = pjoin(tcrm_dir, 'input')

        landmask = config.get('Input', 'LandMask')

        self.landmask = SampleGrid(landmask)

        fmt = config.get('Output', 'Format')

        self.ncflag = False
        if fmt.startswith("nc"):
            self.logger.debug("Output format is netcdf")
            self.ncflag = True
            self.data = {}
            #dimensions = {records}
            #variables = {init_index(records),
            #             genesis_index(records),
            #             non_init_index(records),
            #             lon(records), lat(records),
            #             year(records), month(records),
            #             day(records), hour(records),
            #             minute(records), julianday(records),
            #             bearing(records), speed(records),
            #             pressure(records), lsflag(records), }
            #global_attributes = dict(description=
            #                         source_file=,
            #                         source_file_moddate,
            #                         landmask_file=,
            #                         version=,)
        elif fmt.startswith("txt"):
            self.logger.debug("Output format is text")
            self.origin_year = pjoin(self.processPath, 'origin_year')



    def __doc__(self):
        return 'Processes the database of historical TCs into \
            suitably formatted text files. Data is written to \
            plain text files for ease of access (this may be  \
            upgraded in future versions to netCDF format files).\
            Currently extracts fields containing the value of \
            cyclone parameters, but no information on the \
            change of parameters.'

    def extractTracks(self, index, lon, lat):
        """
        Extract tracks that only *start* within the pre-defined domain.
        The function returns the indices of the tracks that begin within
        the domain.
        """
        outIndex = []
        flag = 0
        for i in xrange(len(index)):
            if index[i] == 1:
                # A new track:
                if ( stats.between(lon[i], self.domain['xMin'], self.domain['xMax']) &
                     stats.between(lat[i], self.domain['yMin'], self.domain['yMax']) ):
                    # We have a track starting within the spatial domain:
                    outIndex.append(i)
                    flag = 1
                else:
                    flag = 0
            else:
                if flag == 1:
                    outIndex.append(i)
                else:
                    pass

        return outIndex

    def filterData(self, data, index):
        return data[:,index]

    def processData(self, restrictToWindfieldDomain=False):
        """
        Process raw data into ASCII files that can be read by the main
        components of the system

        """
        config = ConfigParser()
        config.read(self.configFile)

        self.logger.info("Running %s" % flModuleName())
        inputFile = config.get('DataProcess', 'InputFile')

        # If input file has no path information, default to tcrm input folder
        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.tcrm_input_dir, inputFile)

        self.logger.info("Processing %s" % inputFile)
        
        self.source = config.get('DataProcess', 'Source')

        inputData = colReadCSV(self.configFile, inputFile, self.source)

        inputSpeedUnits = config.get(self.source, 'SpeedUnits')
        inputPressureUnits = config.get(self.source, 'PressureUnits')
        inputLengthUnits = config.get(self.source, 'LengthUnits')
        startSeason = config.getint('DataProcess', 'StartSeason')

        indicator = loadData.getInitialPositions(inputData)
        lat = np.array(inputData['lat'], 'd')
        lon = np.mod(np.array(inputData['lon'], 'd'), 360)

        if restrictToWindfieldDomain:
            # Filter the input arrays to only retain the tracks that
            # pass through the windfield domain.
            CD = CalcTrackDomain(self.configFile)
            self.domain = CD.calcDomainFromTracks(indicator, lon, lat)
            domainIndex = self.extractTracks(indicator, lon, lat)
            inputData = self.filterData(inputData, domainIndex)
            indicator = indicator[domainIndex]
            lon = lon[domainIndex]
            lat = lat[domainIndex]

        if self.progressbar is not None:
            self.progressbar.update(0.125)

        # Sort date/time information
        try:
            dt = np.empty(indicator.size, 'f')
            dt[1:] = np.diff(inputData['age'])
        except (ValueError, KeyError):

            try:
                # Find indicies that satisfy minimum season filter
                idx = np.where(inputData['season'] >= startSeason)[0]
                # Filter records:
                inputData = self.filterData(inputData, idx)
                indicator = indicator[idx]
                lon = lon[idx]
                lat = lat[idx]
            except (ValueError, KeyError):
                pass

            year, month, day, hour, minute = loadData.parseDates(inputData,
                                                             indicator)

            # Time between observations:
            dt = loadData.getTimeDelta(year, month, day, hour, minute)

            # Calculate julian days:
            jdays = loadData.julianDays(year, month, day, hour, minute)

        delta_lon = np.diff(lon)
        delta_lat = np.diff(lat)

        # Split into separate tracks if large jump occurs (delta_lon >
        # 15 degrees or delta_lat > 5 degrees) This avoids two tracks
        # being accidentally combined when seasons and track numbers
        # match but basins are different as occurs in the IBTrACS
        # dataset.  This problem can also be prevented if the
        # 'tcserialno' column is specified.
        indicator[np.where(delta_lon > 15)[0] + 1] = 1
        indicator[np.where(delta_lat > 5)[0] + 1] = 1

        # Save information required for frequency auto-calculation
        try:
            origin_seasonOrYear = np.array(inputData['season'], 'i').compress(indicator)
            header = 'Season'
        except (ValueError, KeyError):
            origin_seasonOrYear = year.compress(indicator)
            header = 'Year'

        flSaveFile(self.origin_year, np.transpose(origin_seasonOrYear),
                   header, ',', fmt='%d')

        pressure = np.array(inputData['pressure'], 'd')
        novalue_index = np.where(pressure==sys.maxint)
        pressure = metutils.convert(pressure, inputPressureUnits, "hPa")
        pressure[novalue_index] = sys.maxint

        # Convert any non-physical central pressure values to maximum integer
        # This is required because IBTrACS has a mix of missing value codes
        # (i.e. -999, 0, 9999) in the same global dataset.
        pressure = np.where((pressure < 600) | (pressure > 1100),
                               sys.maxint, pressure)

        if self.progressbar is not None:
            self.progressbar.update(0.25)

        try:
            vmax = np.array(inputData['vmax'], 'd')
        except (ValueError, KeyError):
            self.logger.warning("No max wind speed data")
            vmax = np.empty(indicator.size, 'f')
        else:
            novalue_index = np.where(vmax==sys.maxint)
            vmax = metutils.convert(vmax, inputSpeedUnits, "mps")
            vmax[novalue_index] = sys.maxint

        assert lat.size == indicator.size
        assert lon.size == indicator.size
        assert pressure.size == indicator.size
        #assert vmax.size == indicator.size

        try:
            rmax = np.array(inputData['rmax'])
            novalue_index = np.where(rmax==sys.maxint)
            rmax = metutils.convert(rmax, inputLengthUnits, "km")
            rmax[novalue_index] = sys.maxint

            self._rmax(rmax, indicator)
            self._rmaxRate(rmax, dt, indicator)
        except (ValueError, KeyError):
            self.logger.warning("No rmax data available")

        if self.ncflag:
            self.data['index'] = indicator

        # ieast : parameter used in latLon2Azi 
        # FIXME: should be a config setting describing the input data.
        ieast = 1

        # Determine the index of initial cyclone observations, excluding
        # those cyclones that have only one observation. This is used
        # for calculating initial bearing and speed
        indicator2 = np.where(indicator > 0, 1, 0)
        initIndex = np.concatenate([np.where(np.diff(indicator2) == \
                                            -1, 1, 0), [0]])

        # Calculate the bearing and distance (km) of every two
        # consecutive records using ll2azi
        bear_, dist_ = maputils.latLon2Azi(lat, lon, ieast, azimuth=0)
        assert bear_.size == indicator.size - 1
        assert dist_.size == indicator.size - 1
        bear = np.empty(indicator.size, 'f')
        bear[1:] = bear_
        dist = np.empty(indicator.size, 'f')
        dist[1:] = dist_

        self._lonLat(lon, lat, indicator, initIndex)
        self._bearing(bear, indicator, initIndex)
        self._bearingRate(bear, dt, indicator)
        if self.progressbar is not None:
            self.progressbar.update(0.375)
        self._speed(dist, dt, indicator, initIndex)
        self._speedRate(dist, dt, indicator)
        self._pressure(pressure, indicator)
        self._pressureRate(pressure, dt, indicator)
        self._windSpeed(vmax)

        try:
            self._frequency(year, indicator)
            self._juliandays(jdays, indicator, year)
        except (ValueError, KeyError):
            pass

        self.logger.info("Completed %s" % flModuleName())
        if self.progressbar is not None:
            self.progressbar.update(0.5)

    def _lonLat(self, lon, lat, indicator, initIndex):
        """
        Extract longitudes and latitudes for all obs, initial obs, TC
        origins and determine a land/sea flag indicating if the TC
        position is over land or sea.

        Input: lon - array of TC longitudes
               lat - array of TC latitudes
               indicator - array of ones/zeros representing initial TC
                           observations, including TCs with a single
                           observation
               initIndex - array of ones/zeros representing initial TC
                           observations, excluding TCs with a single
                           observation

        Output: None - data is written to file

        """

        self.logger.info('Extracting longitudes and latitudes')
        lsflag = np.zeros(len(lon))
        for i, [x, y] in enumerate(zip(lon, lat)):
            if self.landmask.sampleGrid(x, y) > 0:
                lsflag[i] = 1

        lonOne = lon.compress(indicator)
        latOne = lat.compress(indicator)
        lsflagOne = lsflag.compress(indicator)
        lonInit = lon.compress(initIndex)
        latInit = lat.compress(initIndex)
        lsflagInit = lsflag.compress(initIndex)

        origin_lon_lat = pjoin(self.processPath, 'origin_lon_lat')
        init_lon_lat = pjoin(self.processPath, 'init_lon_lat')
        all_lon_lat = pjoin(self.processPath, 'all_lon_lat')

        # Output the lon & lat of cyclone origins
        self.logger.debug('Outputting data into %s' % init_lon_lat)
        self.logger.debug('Outputting data into %s' % origin_lon_lat)
        self.logger.debug('Outputting data into %s' % all_lon_lat)

        header = 'Longitude, Latitude, LSFlag'
        if self.ncflag:
            self.data['longitude'] = lon
            self.data['latitude'] = lat
            self.data['lsflag'] = lsflag
        else:
            flSaveFile(origin_lon_lat,
                       np.transpose([lonOne, latOne, lsflagOne]),
                       header, ',', fmt='%6.2f')
            flSaveFile(init_lon_lat,
                       np.transpose([lonInit, latInit, lsflagInit]),
                       header, ',', fmt='%6.2f')
            flSaveFile(all_lon_lat,
                       np.transpose([lon, lat, lsflag]),
                       header, ',', fmt='%6.2f')

            # Output all cyclone positions:
            cyclone_tracks = pjoin(self.processPath, 'cyclone_tracks')
            self.logger.debug('Outputting data into %s' % cyclone_tracks)
            header = 'Cyclone Origin,Longitude,Latitude, LSflag'
            flSaveFile(cyclone_tracks,
                       np.transpose([indicator, lon, lat, lsflag]),
                       header, ',', fmt='%6.2f')

    def _bearing(self, bear, indicator, initIndex):
        """
        Extract bearings for all obs, initial obs and TC origins
        Input: bear - array of bearing of TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
               initIndex - array of ones/zeros representing initial TC
                           observations (excluding TCs with a single
                           observation)
        Output: None - data is written to file
        """

        self.logger.info('Extracting bearings')

        #extract all bearings
        np.putmask(bear, indicator, sys.maxint)

        # extract initial bearings
        initBearingIndex = np.flatnonzero(initIndex[:-1]) + 1
        initBearing = bear.take(initBearingIndex)

        # extract non-initial bearings
        indicator_ = indicator.copy()
        indicator_.put(initBearingIndex, 1)
        bearingNoInit = bear.compress(indicator_ == 0)

        if self.ncflag:
            self.data['bearing'] = bear
            self.data['init_bearing'] = initBearing
            self.data['bearing_no_init'] = bearingNoInit
        else:
            all_bearing = pjoin(self.processPath, 'all_bearing')
            self.logger.debug('Outputting data into %s' % all_bearing)
            header = 'all cyclone bearing in degrees'
            flSaveFile(all_bearing, bear, header, fmt='%6.2f')

            init_bearing = pjoin(self.processPath, 'init_bearing')
            self.logger.debug('Outputting data into %s' % init_bearing)
            header = 'initial cyclone bearing in degrees'
            flSaveFile(init_bearing, initBearing, header, fmt='%6.2f')
            bearing_no_init = pjoin(self.processPath, 'bearing_no_init')
            self.logger.debug('Outputting data into %s' % bearing_no_init)
            header = 'cyclone bearings without initial ones in degrees'
            flSaveFile(bearing_no_init, bearingNoInit, header, fmt='%6.2f')

    def _speed(self, dist, dt, indicator, initIndex):
        """
        Extract speeds for all obs, initial obs and TC origins
        Input: dist - array of distances between consecutive TC
                      observations
               dt - array of times between consecutive TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
               initIndex - array of ones/zeros representing initial TC
                           observations (excluding TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting speeds')
        speed = dist/dt
        # Delete speeds less than 0, greated than 200,
        # or where indicator == 1.
        np.putmask(speed, (speed < 0) | (speed > 200) | indicator,
                      sys.maxint)
        np.putmask(speed, np.isnan(speed), sys.maxint)

        initSpeedIndex = np.flatnonzero(initIndex[:-1]) + 1
        initSpeed = speed.take(initSpeedIndex)
        indicator_ = indicator.copy()
        indicator_.put(initSpeedIndex, 1)
        speedNoInit = speed.compress(indicator_ == 0)

        if self.ncflag:
            self.data['speed'] = speed
            self.data['init_speed'] = initSpeed
            self.data['speed_no_init'] = speedNoInit
        else:
            init_speed = pjoin(self.processPath, 'init_speed')
            all_speed = pjoin(self.processPath, 'all_speed')
            speed_no_init = pjoin(self.processPath, 'speed_no_init')
            # Extract all speeds
            self.logger.debug('Outputting data into %s'%all_speed)
            header = 'all cyclone speed in km/hour'
            flSaveFile(all_speed, speed, header, fmt='%6.2f')

            # Extract initial speeds
            self.logger.debug('Outputting data into %s'%init_speed)
            header = 'initial cyclone speed in km/hour'
            flSaveFile(init_speed, initSpeed, header, fmt='%f')

            # Extract speeds, excluding initial speeds
            self.logger.debug('Outputting data into %s'%speed_no_init)
            header = 'cyclone speed without initial ones in km/hour'
            flSaveFile(speed_no_init, speedNoInit, header, fmt='%6.2f')

    def _pressure(self, pressure, indicator):
        """Extract pressure for all obs, initial obs and TC origins
        Input: pressure - array of central pressure observations for TC
                          observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
               initIndex - array of ones/zeros representing initial TC
                           observations (excluding TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting pressures')
        initPressure = pressure.compress(indicator)
        pressureNoInit = pressure.compress(indicator == 0)
        pressureNoInit = pressureNoInit.compress(pressureNoInit < sys.maxint)

        if self.ncflag:
            self.data['pressure'] = pressure
            self.data['init_pressure'] = initPressure
            self.data['pressure_no_init'] = pressureNoInit
        else:
            init_pressure = pjoin(self.processPath, 'init_pressure')
            all_pressure = pjoin(self.processPath, 'all_pressure')
            pressure_no_init = pjoin(self.processPath, 'pressure_no_init')
            # Extract all pressure
            self.logger.debug('Outputting data into %s'%all_pressure)
            header = 'all cyclone pressure in hPa'
            flSaveFile(all_pressure, pressure, header, fmt='%7.2f')

            # Extract initial pressures
            self.logger.debug('Outputting data into %s'%init_pressure)
            header = 'initial cyclone pressure in hPa'
            flSaveFile(init_pressure, initPressure, header, fmt='%7.2f')

            # Extract pressures, excluding initial times
            self.logger.debug('Outputting data into %s'%pressure_no_init)
            header = 'cyclone pressure without initial ones in hPa'
            flSaveFile(pressure_no_init, pressureNoInit, header, fmt='%7.2f')

    def _pressureRate(self, pressure, dt, indicator):
        """Extract the rate of pressure change from the pressure values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in pressure from the previous observation
        is undefined. Entries corresponding to records with no pressure
        observation are also set to maxint.
        Input: pressure - array of central pressure observations for TC
               observations
               dt - array of times between consecutive TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting the rate of pressure change')

        #Change in pressure:
        pressureChange_ = np.diff(pressure)
        pressureChange = np.empty(indicator.size, 'f')
        pressureChange[1:] = pressureChange_

        # Rate of pressure change:
        pressureRate = pressureChange/dt

        # Mask rates corresponding to initial times, times when
        # the pressure is known to be missing, and when the
        # pressure rate is greater than 10 hPa/hour (a sanity check).
        # The highest rate of intensification on record is
        # Typhoon Forrest (Sept 1983) 100 mb in 24 hrs.
        
        np.putmask(pressureRate, indicator, sys.maxint)
        np.putmask(pressureRate, pressure >= sys.maxint, sys.maxint)
        np.putmask(pressureRate, np.isnan(pressureRate), sys.maxint)
        np.putmask(pressureRate, np.abs(pressureRate)>10, sys.maxint)

        if self.ncflag:
            self.data['pressureRate'] = pressureRate
        else:
            pressure_rate = pjoin(self.processPath, 'pressure_rate')
            self.logger.debug('Outputting data into %s'%pressure_rate)
            header = 'All pressure change rates (hPa/hr)'
            flSaveFile(pressure_rate, pressureRate, header, fmt='%6.2f')

    def _bearingRate(self, bear, dt, indicator):
        """Extract the rate of bearing change for each cyclone:
        Entries corresponding to initial position reports and the
        second observation are set to maxint. The first entry is set
        to maxint as there is no bearing associated with it and the
        second entry is therefore non-sensical.
        Input: bear - array of bearings between consecutive TC
                      observations
               dt - array of times between consecutive TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting the rate of bearing change')

        bearingChange_ = np.diff(bear)
        ii = np.where((bearingChange_ > 180.))
        jj = np.where((bearingChange_ < -180.))
        bearingChange_[ii] -= 360.
        bearingChange_[jj] += 360.
        bearingChange = np.empty(indicator.size, 'd')
        bearingChange[1:] = bearingChange_

        bearingRate = bearingChange/dt

        np.putmask(bearingRate, indicator, sys.maxint)
        np.putmask(bearingRate[1:], indicator[:-1], sys.maxint)
        np.putmask(bearingRate, (bearingRate >= sys.maxint) |
                                    (bearingRate <= -sys.maxint),
                                      sys.maxint)

        np.putmask(bearingRate, np.isnan(bearingRate), sys.maxint)

        if self.ncflag:
            self.data['bearingRate'] = bearingRate
        else:
            bearing_rate = pjoin(self.processPath, 'bearing_rate')
            self.logger.debug('Outputting data into %s'%bearing_rate)
            header = 'All bearing change rates (degrees/hr)'
            flSaveFile(bearing_rate, bearingRate, header, fmt='%6.2f')

    def _speedRate(self, dist, dt, indicator):
        """Extract the rate of speed change for each cyclone:
        Note this results in some odd values for the accelerations,
        propagated from odd position reports. Entries corresponding to
        initial position reports and the second observation are set to
        maxint. The first entry is set to maxint as there is no speed
        associated with it and the second is therefore non-sensical.
        Input: dist - array of distances between consecutive TC
                      observations
               dt - array of times between consecutive TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting the rate of speed change for each cyclone')

        speed = dist/dt
        speedChange_ = np.diff(speed)
        speedChange = np.empty(indicator.size, 'd')
        speedChange[1:] = speedChange_

        indicator_ = indicator.copy()
        np.putmask(indicator_, (speed < 0) | (speed > 200), 1)

        speedRate = speedChange/dt

        np.putmask(speedRate, indicator_, sys.maxint)
        np.putmask(speedRate[1:], indicator_[:-1], sys.maxint)
        np.putmask(speedRate, (speedRate >= sys.maxint) |
                                 (speedRate <= -sys.maxint), sys.maxint)

        np.putmask(speedRate, np.isnan(speedRate), sys.maxint)

        if self.ncflag:
            self.data['speedRate'] = speedRate
        else:
            speed_rate = pjoin(self.processPath, 'speed_rate')
            self.logger.debug('Outputting data into %s'%speed_rate)
            header = 'All speed change rates (km/hr/hr)'
            flSaveFile(speed_rate, speedRate, header, fmt='%6.2f')

    def _windSpeed(self, windSpeed):
        """Extract maximum sustained wind speeds
        Input: windSpeed - array of windspeeds for TC observations

        Output: None - data is written to file
        """
        self.logger.info('Extracting maximum sustained wind speeds')
        np.putmask(windSpeed, windSpeed > 200., sys.maxint)
        if self.ncflag:
            self.data['windspeed'] = windSpeed
        else:
            wind_speed = pjoin(self.processPath, 'wind_speed')
            self.logger.debug('Outputting data into %s'%wind_speed)
            header = 'Maximum wind speed (m/s)'
            flSaveFile(wind_speed, windSpeed, header, fmt='%6.2f')

    def _rmax(self, rmax, indicator):
        """Extract radii to maximum wind:
        Input: rmax - array of radii to maximum winds for TC
                      observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info("Extracting radii to maximum winds")
        initrmax = rmax.compress(indicator)
        rmaxNoInit = rmax.compress(indicator == 0)
        rmaxNoInit = rmaxNoInit.compress(rmaxNoInit < sys.maxint)
        if self.ncflag:
            self.data['rmax'] = rmax
            self.data['init_rmax'] = initrmax
            self.data['rmax_no_init'] = rmaxNoInit
        else:
            init_rmax = pjoin(self.processPath, 'init_rmax')
            all_rmax = pjoin(self.processPath, 'all_rmax')
            rmax_no_init = pjoin(self.processPath, 'rmax_no_init')
            #extract all rmax
            self.logger.debug('Outputting data into %s'%all_rmax)
            header = 'rMax (km)'
            flSaveFile(all_rmax, rmax, header, fmt='%6.2f')

            #extract initial rmax
            self.logger.debug('Outputting data into %s'%init_rmax)
            header = 'initial rmax (km)'
            flSaveFile(init_rmax, initrmax, header, fmt='%6.2f')

            #extract rmax no init
            self.logger.debug('Outputting data into %s'%rmax_no_init)
            header = 'rmax excluding initial ones (km)'
            flSaveFile(rmax_no_init, rmaxNoInit, header, fmt='%6.2f')

    def _rmaxRate(self, rmax, dt, indicator):
        """Extract the rate of size change from the rmax values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in rmax from the previous observation is
        undefined. Entries corresponding to records with no rmax
        observation are also set to maxint.
        Input: rmax - array of radii to maximum winds for TC
                      observations
               dt - array of times between consecutive TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting the rate of size change')

        #Change in rmax:
        rmaxChange_ = np.diff(rmax)
        rmaxChange = np.empty(indicator.size, 'f')
        rmaxChange[1:] = rmaxChange_

        # Rate of rmax change:
        rmaxRate = rmaxChange/dt

        # Mask rates corresponding to initial times and times when
        # the rmax is known to be missing.
        self.logger.debug('Outputting data into %s'%self.rmax_rate)
        np.putmask(rmaxRate, indicator, sys.maxint)
        np.putmask(rmaxRate, rmax >= sys.maxint, sys.maxint)
        np.putmask(rmaxRate, (rmaxRate >= sys.maxint) |
                                (rmaxRate <= -sys.maxint), sys.maxint)
        np.putmask(rmaxRate, np.isnan(rmaxRate), sys.maxint)

        if self.ncflag:
            self.data['rmaxRate'] = rmaxRate
        else:
            rmax_rate = pjoin(self.processPath, 'rmax_rate') 
            header = 'All rmax change rates (km/hr)'
            flSaveFile(rmax_rate, rmaxRate, header, fmt='%6.2f')

    def _frequency(self, years, indicator):
        # Generate a histogram of the annual frequency of events from the input data
        self.logger.info('Extracting annual frequency of events')
        minYr = years.min()
        maxYr = years.max()
        genesisYears = years.compress(indicator)
        if minYr == maxYr:
            self.logger.info("First and last year of input data are the same")
            self.logger.info("Cannot generate histogram of frequency")
        else:
            frequency = pjoin(self.processPath, 'frequency')
            bins = np.arange(minYr,maxYr+2,1)
            n, b = np.histogram(genesisYears, bins)
            header = 'Year,count'
            flSaveFile(frequency, np.transpose( [bins[:-1], n] ),
                       header, fmt='%6.2f' )

            self.logger.info( "Mean annual frequency: %5.1f"%np.mean( n ) )
            self.logger.info( "Standard deviation: %5.1f"%np.std( n ) )

    def _juliandays(self, jdays, indicator, years ):
        """
        Generate a distribution of the formation day of
        year from observations
        """

        self.logger.info( "Calculating annual distribution of observations" )

        # Do a bodgy job of addressing 29th of February (there surely
        # must be a recommended way of accounting for leap years)

        for i in range( len( jdays ) ):
            if ( years[i]%4 == 0 ) and ( jdays[i] >= 60 ):
                jdays[i] -= 1

        bins = np.arange( 1, 367 )
        n,b = np.histogram( jdays.compress( indicator ), bins )
        header = 'Day,count'
        jday_genesis = pjoin(self.processPath, 'jday_genesis')
        jday_observations = pjoin(self.processPath, 'jday_obs')
        jday = pjoin( self.processPath, 'jdays' )
        # Distribution of genesis days (histogram):
        flSaveFile( jday_genesis, np.transpose( [bins[:-1], n] ),
                    header, fmt='%d', delimiter=',' )
        n,b = np.histogram( jdays, bins)
        # Distribution of all days (histogram):
        flSaveFile( jday_observations, np.transpose( [bins[:-1],n] ),
                    header, fmt='%d', delimiter=',' )
        # All days:
        flSaveFile( jday, np.transpose( jdays.compress( indicator ) ),
                    header='Day', fmt='%d' )

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename does not exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file does not exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg

    config = ConfigParser()
    config.read(configFile)

    logFile = config.read('Logging', 'LogFile')
    logLevel = config.read('Logging', 'LogLevel')
    verbose = config.read('Logging', 'Verbose')

    flStartLog(logFile, logLevel, verbose)

    dp = DataProcess(configFile)
    dp.processData()
    logging.shutdown()
