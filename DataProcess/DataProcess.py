"""
:mod:`DataProcess` -- basic processing of input database
========================================================

.. module:: DataProcess
    :synopsis: Processes the database of historical TCs into
               suitably formatted text files.

Processes the database of historical TCs into simply formatted
text files. Extracts locations, including genesis location, and
calculates speed, bearing, genesis day (as the day of the year)
annual frequency (of the complete database), pressure rates of
change, speed rates of change and bearing rates of change.


"""

import os
import sys
import logging
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
from .CalcTrackDomain import CalcTrackDomain
from Utilities.config import ConfigParser

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


class DataProcess(object):
    """
    Processes the database of historical TCs into suitably formatted
    text files.  Data is written to plain text files for ease of
    access (this may be upgraded in future versions to netCDF files).
    Currently extracts fields containing the value of cyclone
    parameters, but no information on the change of parameters.

    :type  configFile: str
    :param configFile: Configuration file containing simulation settings

    :type  progressbar: :class:`progressbar`
    :param progressbar: :attr:`progressbar` object to print progress to
                        STDOUT

    Internal Methods:
    _lonLat(lon, lat, indicator) Extract longitudes and latitudes
    _bearing(bear, indicator) Extract bearings
    _speed(dist, dt, indicator) Extract speeds
    _pressure(pressure, indicator) Extract pressures
    _pressureRate(pressure, dt, indicator) Extract rate of presssure change
    _speedRate(dist, dt, indicator) Extract the acceleration (rate of change of speed)
    _bearingChange(bear, dt, indicator) Extract the rate of change of bearing
    _windSpeed(vmax, indicator) Extract the maximum sustained wind speed

    """

    def __init__(self, configFile, progressbar=None):

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
            # variables = {init_index(records),
            #             genesis_index(records),
            #             non_init_index(records),
            #             lon(records), lat(records),
            #             year(records), month(records),
            #             day(records), hour(records),
            #             minute(records), julianday(records),
            #             bearing(records), speed(records),
            #             pressure(records), lsflag(records), }
            # global_attributes = dict(description=
            #                         source_file=,
            #                         source_file_moddate,
            #                         landmask_file=,
            #                         version=,)
        elif fmt.startswith("txt"):
            self.logger.debug("Output format is text")
            self.origin_year = pjoin(self.processPath, 'origin_year')

    def extractTracks(self, index, lon, lat):
        """
        Extract tracks that only *start* within the pre-defined domain.
        The function returns the indices of the tracks that begin within
        the domain.

        :param index: indicator of initial positions of TCs (1 = new TC,
                      0 = continuation of previous TC).
        :param lon: longitudes of TC positions
        :param lat: latitudes of TC positions

        :type index: `numpy.ndarray`
        :type lon: `numpy.ndarray`
        :type lat: `numpy.ndarray`

        :returns: indices corresponding to all points of all tracks that
                  form within the model domain
        :rtype: `numpy.ndarray`

        """
        outIndex = []
        flag = 0
        for i in range(len(index)):
            if index[i] == 1:
                # A new track:
                if (stats.between(lon[i],
                                  self.domain['xMin'],
                                  self.domain['xMax']) &
                        stats.between(lat[i],
                                      self.domain['yMin'],
                                      self.domain['yMax'])):
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

    def processData(self, restrictToWindfieldDomain=False):
        """
        Process raw data into ASCII files that can be read by the main
        components of the system

        :param bool restrictToWindfieldDomain: if True, only process data
            within the wind field domain, otherwise, process data from
            across the track generation domain.

        """
        config = ConfigParser()
        config.read(self.configFile)

        self.logger.info("Running {0}".format(flModuleName()))

        if config.has_option('DataProcess', 'InputFile'):
            inputFile = config.get('DataProcess', 'InputFile')
            self.logger.info(f"Input file from DataProcess: {inputFile}")

        if config.has_option('DataProcess', 'Source'):
            source = config.get('DataProcess', 'Source')
            self.logger.info(f"Loading {source} dataset")
            fn = config.get(source, 'Filename')
            path = config.get(source, 'Path')
            inputFile = pjoin(path, fn)
            self.logger.info(f"Input file set to {inputFile}")

        # If input file has no path information, default to tcrm input folder
        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.tcrm_input_dir, inputFile)

        self.logger.info(f"Processing {inputFile}")

        self.source = config.get('DataProcess', 'Source')

        inputData = colReadCSV(self.configFile, inputFile, self.source)

        inputSpeedUnits = config.get(self.source, 'SpeedUnits')
        inputPressureUnits = config.get(self.source, 'PressureUnits')
        inputLengthUnits = config.get(self.source, 'LengthUnits')
        startSeason = config.getint('DataProcess', 'StartSeason')

        indicator = loadData.getInitialPositions(inputData)
        lat = np.array(inputData['lat'], 'float32')
        lon = np.mod(np.array(inputData['lon'], 'float32'), 360.)

        if restrictToWindfieldDomain:
            # Filter the input arrays to only retain the tracks that
            # pass through the windfield domain.
            CD = CalcTrackDomain(self.configFile)
            self.domain = CD.calcDomainFromTracks(indicator, lon, lat)
            domainIndex = self.extractTracks(indicator, lon, lat)
            inputData = inputData[domainIndex]
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
                self.logger.info(("Filtering input data by season: "
                                  f"season >= {startSeason}"))
                # Find indicies that satisfy minimum season filter
                idx = np.where(inputData['season'] >= startSeason)[0]
                # Filter records:
                inputData = inputData[idx]
                indicator = indicator[idx]
                lon = lon[idx]
                lat = lat[idx]
            except (ValueError, KeyError):
                pass

            year, month, day, hour, minute, datetimes \
                = loadData.parseDates(inputData, indicator)

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
            origin_seasonOrYear = np.array(
                inputData['season'], 'i').compress(indicator)
            header = 'Season'
        except (ValueError, KeyError):
            origin_seasonOrYear = year.compress(indicator)
            header = 'Year'

        flSaveFile(self.origin_year, np.transpose(origin_seasonOrYear),
                   header, ',', fmt='%d')

        pressure = np.array(inputData['pressure'], 'float32')
        novalue_index = np.where(pressure == sys.maxsize)
        pressure = metutils.convert(pressure, inputPressureUnits, "hPa")
        pressure[novalue_index] = sys.maxsize

        # Convert any non-physical central pressure values to maximum integer
        # This is required because IBTrACS has a mix of missing value codes
        # (i.e. -999, 0, 9999) in the same global dataset.
        pressure = np.where((pressure < 600) | (pressure > 1100),
                            sys.maxsize, pressure)

        if self.progressbar is not None:
            self.progressbar.update(0.25)

        try:
            vmax = np.array(inputData['vmax'], 'float32')
        except (ValueError, KeyError):
            self.logger.warning("No max wind speed data")
            vmax = np.empty(indicator.size, 'f')
        else:
            novalue_index = np.where(vmax == sys.maxsize)
            vmax = metutils.convert(vmax, inputSpeedUnits, "mps")
            vmax[novalue_index] = sys.maxsize

        assert lat.size == indicator.size
        assert lon.size == indicator.size
        assert pressure.size == indicator.size
        #assert vmax.size == indicator.size

        try:
            rmax = np.array(inputData['rmax'])
            novalue_index = np.where(rmax == sys.maxsize)
            rmax = metutils.convert(rmax, inputLengthUnits, "km")
            rmax[novalue_index] = sys.maxsize

            self._rmax(rmax, indicator)
            self._rmaxRate(rmax, dt, indicator)
        except (ValueError, KeyError):
            self.logger.warning("No rmax data available")
            rmax = np.empty(indicator.size, 'f')

	try:
            r34 = np.array(inputData['r34'], 'd')
	    r34 = np.where((r34 < 0) | (r34 > 500), sys.maxint, r34)
            novalue_index = np.where(r34 == sys.maxint)
            r34           = metutils.convert(r34, inputLengthUnits, "km")
 	    r34[novalue_index] = sys.maxint
            self._r34    (r34,     indicator)
            self._r34Rate(r34, dt, indicator)
        except (ValueError, KeyError):
            self.logger.warning("No r34 data available")
            r34 = np.empty(indicator.size, 'f')

        try:
            r64 = np.array(inputData['r64'], 'd')
	    r64 = np.where((r64 < 0) | (r64 > 300), sys.maxint, r64)
            novalue_index = np.where(r64 == sys.maxint)
            r64           = metutils.convert(r64, inputLengthUnits, "km")
	    r64[novalue_index] = sys.maxint
            self._r64    (r64,     indicator)
            self._r64Rate(r64, dt, indicator)
        except (ValueError, KeyError):
            self.logger.warning("No r64 data available")
            r64 = np.empty(indicator.size, 'f')

        if self.ncflag:
            self.data['index'] = indicator

        # ieast : parameter used in latLon2Azi
        # FIXME: should be a config setting describing the input data.
        ieast = 1

        # Determine the index of initial cyclone observations, excluding
        # those cyclones that have only one observation. This is used
        # for calculating initial bearing and speed
        indicator2 = np.where(indicator > 0, 1, 0)
        initIndex = np.concatenate([np.where(np.diff(indicator2) ==
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
        self._windSpeed(vmax, indicator)
        self._windSpeedRate(vmax, dt, indicator)

        try:
            self._frequency(year, indicator)
            self._juliandays(jdays, indicator, year)
        except (ValueError, KeyError):
            pass

        self.logger.info("Completed {0}".format(flModuleName()))
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
        self.logger.debug('Outputting data into {0}' .format(init_lon_lat))
        self.logger.debug('Outputting data into {0}'.format(origin_lon_lat))
        self.logger.debug('Outputting data into {0}'.format(all_lon_lat))

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
            self.logger.debug('Outputting data into {0}'.format(cyclone_tracks))
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

        # extract all bearings
        np.putmask(bear, indicator, sys.maxsize)

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
            self.logger.debug('Outputting data into {0}'.format(all_bearing))
            header = 'all cyclone bearing in degrees'
            flSaveFile(all_bearing, bear, header, fmt='%6.2f')

            init_bearing = pjoin(self.processPath, 'init_bearing')
            self.logger.debug('Outputting data into {0}'.format(init_bearing))
            header = 'initial cyclone bearing in degrees'
            flSaveFile(init_bearing, initBearing, header, fmt='%6.2f')
            bearing_no_init = pjoin(self.processPath, 'bearing_no_init')
            self.logger.debug('Outputting data into {0}'.format(bearing_no_init))
            header = 'cyclone bearings without initial ones in degrees'
            flSaveFile(bearing_no_init, bearingNoInit, header, fmt='%6.2f')

    def _speed(self, dist, dt, indicator, initIndex):
        """
        Extract speeds for all obs, initial obs and TC origins
        Input: dist - array of distances between consecutive TC
                      observations (km)
               dt - array of times between consecutive TC observations (hours)
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
               initIndex - array of ones/zeros representing initial TC
                           observations (excluding TCs with a single
                           observation)
        Output: None - data is written to file
        """
        self.logger.info('Extracting speeds')
        speed = dist / dt
        # Delete speeds less than 0, greated than 200,
        # or where indicator == 1.
        np.putmask(speed, (speed < 0) | (speed > 200) | indicator,
                   sys.maxsize)
        np.putmask(speed, np.isnan(speed), sys.maxsize)

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
            self.logger.debug('Outputting data into {0}'.format(all_speed))
            header = 'all cyclone speed in km/h'
            flSaveFile(all_speed, speed, header, fmt='%6.2f')

            # Extract initial speeds
            self.logger.debug('Outputting data into {0}'.format(init_speed))
            header = 'initial cyclone speed in km/h'
            flSaveFile(init_speed, initSpeed, header, fmt='%f')

            # Extract speeds, excluding initial speeds
            self.logger.debug('Outputting data into {0}'.format(speed_no_init))
            header = 'cyclone speed without initial ones in km/h'
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
        pressureNoInit = pressureNoInit.compress(pressureNoInit < sys.maxsize)

        if self.ncflag:
            self.data['pressure'] = pressure
            self.data['init_pressure'] = initPressure
            self.data['pressure_no_init'] = pressureNoInit
        else:
            init_pressure = pjoin(self.processPath, 'init_pressure')
            all_pressure = pjoin(self.processPath, 'all_pressure')
            pressure_no_init = pjoin(self.processPath, 'pressure_no_init')
            # Extract all pressure
            self.logger.debug('Outputting data into {0}'.format(all_pressure))
            header = 'all cyclone pressure in hPa'
            flSaveFile(all_pressure, pressure, header, fmt='%7.2f')

            # Extract initial pressures
            self.logger.debug('Outputting data into {0}'.format(init_pressure))
            header = 'initial cyclone pressure in hPa'
            flSaveFile(init_pressure, initPressure, header, fmt='%7.2f')

            # Extract pressures, excluding initial times
            self.logger.debug('Outputting data into {0}'.format(pressure_no_init))
            header = 'cyclone pressure without initial ones in hPa'
            flSaveFile(pressure_no_init, pressureNoInit, header, fmt='%7.2f')

    def _pressureRate(self, pressure, dt, indicator):
        """Extract the rate of pressure change from the pressure values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in pressure from the previous observation
        is undefined. Entries corresponding to records with no pressure
        observation are also set to maxint.

        :param pressure: :class:`numpy.ndarray` of central pressure 
                         observations for TCs
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation

        Output: None - data is written to file
        """
        self.logger.info('Extracting the rate of pressure change')

        # Change in pressure:
        pressureChange_ = np.diff(pressure)
        pressureChange = np.empty(indicator.size, 'f')
        pressureChange[1:] = pressureChange_

        # Rate of pressure change:
        pressureRate = pressureChange / dt

        # Mask rates corresponding to initial times, times when
        # the pressure is known to be missing, and when the
        # pressure rate is greater than 10 hPa/hour (a sanity check).
        # The highest rate of intensification on record is
        # Typhoon Forrest (Sept 1983) 100 mb in 24 hrs.

        np.putmask(pressureRate, indicator, sys.maxsize)
        np.putmask(pressureRate, pressure >= sys.maxsize, sys.maxsize)
        np.putmask(pressureRate, np.isnan(pressureRate), sys.maxsize)
        np.putmask(pressureRate, np.abs(pressureRate) > 10, sys.maxsize)

        if self.ncflag:
            self.data['pressureRate'] = pressureRate
        else:
            pressure_rate = pjoin(self.processPath, 'pressure_rate')
            self.logger.debug('Outputting data into {0}'.format(pressure_rate))
            header = 'All pressure change rates (hPa/hr)'
            flSaveFile(pressure_rate, pressureRate, header, fmt='%6.2f')

    def _bearingRate(self, bear, dt, indicator):
        """Extract the rate of bearing change for each cyclone:
        Entries corresponding to initial position reports and the
        second observation are set to maxint. The first entry is set
        to maxint as there is no bearing associated with it and the
        second entry is therefore non-sensical.

        :param bear :class:`numpy.ndarray` of bearing observations
                     for TCs
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation

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

        bearingRate = bearingChange / dt

        np.putmask(bearingRate, indicator, sys.maxsize)
        np.putmask(bearingRate[1:], indicator[:-1], sys.maxsize)
        np.putmask(bearingRate, (bearingRate >= sys.maxsize) |
                   (bearingRate <= -sys.maxsize),
                   sys.maxsize)

        np.putmask(bearingRate, np.isnan(bearingRate), sys.maxsize)

        if self.ncflag:
            self.data['bearingRate'] = bearingRate
        else:
            bearing_rate = pjoin(self.processPath, 'bearing_rate')
            self.logger.debug('Outputting data into {0}'.format(bearing_rate))
            header = 'All bearing change rates (degrees/hr)'
            flSaveFile(bearing_rate, bearingRate, header, fmt='%6.2f')

    def _speedRate(self, dist, dt, indicator):
        """Extract the rate of speed change for each cyclone:
        Note this results in some odd values for the accelerations,
        propagated from odd position reports. Entries corresponding to
        initial position reports and the second observation are set to
        maxint. The first entry is set to maxint as there is no speed
        associated with it and the second is therefore non-sensical.

        :param dist: :class:`numpy.ndarray` of distances between consecutive 
                     TC observations
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation

        Output: None - data is written to file
        """
        self.logger.info(
            'Extracting the rate of speed change for each cyclone')

        speed = dist / dt
        speedChange_ = np.diff(speed)
        speedChange = np.empty(indicator.size, 'd')
        speedChange[1:] = speedChange_

        indicator_ = indicator.copy()
        np.putmask(indicator_, (speed < 0) | (speed > 200), 1)

        speedRate = speedChange / dt

        np.putmask(speedRate, indicator_, sys.maxsize)
        np.putmask(speedRate[1:], indicator_[:-1], sys.maxsize)
        np.putmask(speedRate, (speedRate >= sys.maxsize) |
                   (speedRate <= -sys.maxsize), sys.maxsize)

        np.putmask(speedRate, np.isnan(speedRate), sys.maxsize)

        if self.ncflag:
            self.data['speedRate'] = speedRate
        else:
            speed_rate = pjoin(self.processPath, 'speed_rate')
            self.logger.debug('Outputting data into {0}'.format(speed_rate))
            header = 'All speed change rates (km/hr/hr)'
            flSaveFile(speed_rate, speedRate, header, fmt='%6.2f')

     def _windSpeed(self, vmax, indicator):
        """
        Extract maximum sustained wind speeds. Data is written to file.

        :param vmax: :class:`numpy.ndarray` of wind speeds from TC obs.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation
        
        """
        self.logger.info('Extracting maximum sustained wind speeds')
        np.putmask(windSpeed, windSpeed > 200., sys.maxsize)
        if self.ncflag:
            self.data['vmax'] = vmax
	    self.data['init_vmax'] = initVmax
            self.data['vmax_no_init'] = vmaxNoInit
        else:
	    init_vmax = pjoin(self.processPath, 'init_vmax')
            all_vmax = pjoin(self.processPath, 'all_vmax')
            vmax_no_init = pjoin(self.processPath, 'vmax_no_init')
	
	    # Extract all vmax
            self.logger.debug('Outputting data into %s' % all_vmax)
            header = 'all cyclone vmax  '
            flSaveFile(all_vmax, vmax, header, fmt='%7.2f')

            # Extract initial vmax
            self.logger.debug('Outputting data into %s' % init_vmax)
            header = 'initial cyclone vmax'
            flSaveFile(init_vmax, initVmax, header, fmt='%7.2f')

            # Extract vmax, excluding initial times
            self.logger.debug('Outputting data into %s' % vmax_no_init)
            header = 'cyclone vmax without initial ones '
            flSaveFile(vmax_no_init, vmaxNoInit, header, fmt='%7.2f')

    def _windSpeedRate(self, vmax, dt, indicator):
        """Extract the rate of vmax change from the vmax values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in vmax from the previous observation
        is undefined. Entries corresponding to records with no vmax
        observation are also set to maxint.

        :param vmax: :class:`numpy.ndarray` of central vmax observations
                     for TCs
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation
        """
        self.logger.info('Extracting the rate of vmax change')

        # Change in vmax:
        vmaxChange_ = np.diff(vmax)
        vmaxChange = np.empty(indicator.size, 'f')
        vmaxChange[1:] = vmaxChange_

        # Rate of vmax change:
        vmaxRate = vmaxChange / dt

        # Mask rates corresponding to initial times, times when
        # the vmax is known to be missing, and when the
        # vmax rate is greater than 10 m/s (a sanity check).
        # The highest rate of intensification on record is
        # Typhoon Forrest (Sept 1983) 100 mb in 24 hrs.

        np.putmask(vmaxRate, indicator, sys.maxint)
        np.putmask(vmaxRate, vmax >= sys.maxint, sys.maxint)
        np.putmask(vmaxRate, (vmaxRate >= sys.maxint) | \
                   (vmaxRate <= -sys.maxint), sys.maxint)
        np.putmask(vmaxRate, np.isnan(vmaxRate), sys.maxint)
        np.putmask(vmaxRate, np.abs(vmaxRate) > 80, sys.maxint)

        if self.ncflag:
            self.data['vmaxRate'] = vmaxRate
        else:
            vmax_rate = pjoin(self.processPath, 'vmax_rate')
            self.logger.debug('Outputting data into {0}'.format(vmax_rate))
            header = 'All vmax change rates '
            flSaveFile(vmax_rate, vmaxRate, header, fmt='%6.2f')

    def _rmax(self, rmax, indicator):
        """
        Extract radii to maximum wind.

        :param rmax: :class:`numpy.ndarray` of radii to maximum winds 
                     for TC observatons
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with a
                          single observation.

        Output: None - data is written to file
        """
        self.logger.info("Extracting radii to maximum winds")
        initrmax = rmax.compress(indicator)
        rmaxNoInit = rmax.compress(indicator == 0)
        rmaxNoInit = rmaxNoInit.compress(rmaxNoInit < sys.maxsize)
        if self.ncflag:
            self.data['rmax'] = rmax
            self.data['init_rmax'] = initrmax
            self.data['rmax_no_init'] = rmaxNoInit
        else:
            init_rmax = pjoin(self.processPath, 'init_rmax')
            all_rmax = pjoin(self.processPath, 'all_rmax')
            rmax_no_init = pjoin(self.processPath, 'rmax_no_init')
            # extract all rmax
            self.logger.debug('Outputting data into {0}'.format(all_rmax))
            header = 'rMax (km)'
            flSaveFile(all_rmax, rmax, header, fmt='%6.2f')

            # extract initial rmax
            self.logger.debug('Outputting data into {0}'.format(init_rmax))
            header = 'initial rmax (km)'
            flSaveFile(init_rmax, initrmax, header, fmt='%6.2f')

            # extract rmax no init
            self.logger.debug('Outputting data into {0}'.format(rmax_no_init))
            header = 'rmax excluding initial ones (km)'
            flSaveFile(rmax_no_init, rmaxNoInit, header, fmt='%6.2f')

    def _rmaxRate(self, rmax, dt, indicator):
        """Extract the rate of size change from the rmax values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in rmax from the previous observation is
        undefined. Entries corresponding to records with no rmax
        observation are also set to maxint.

        :param rmax: :class:`numpy.ndarray` of radii to maximum winds 
                     for TC observatons.
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with a
                          single observation.

        Output: None - data is written to file

        """
        self.logger.info('Extracting the rate of size change')

        # Change in rmax:
        rmaxChange_ = np.diff(rmax)
        rmaxChange = np.empty(indicator.size, 'f')
        rmaxChange[1:] = rmaxChange_

        # Rate of rmax change:
        rmaxRate = rmaxChange / dt

        # Mask rates corresponding to initial times and times when
        # the rmax is known to be missing.
        self.logger.debug('Outputting data into {0}'
                          .format(pjoin(self.processPath, 'rmax_rate')))
        np.putmask(rmaxRate, indicator, sys.maxsize)
        np.putmask(rmaxRate, rmax >= sys.maxsize, sys.maxsize)
        np.putmask(rmaxRate, (rmaxRate >= sys.maxsize) |
                   (rmaxRate <= -sys.maxsize), sys.maxsize)
        np.putmask(rmaxRate, np.isnan(rmaxRate), sys.maxsize)

        if self.ncflag:
            self.data['rmaxRate'] = rmaxRate
        else:
            rmax_rate = pjoin(self.processPath, 'rmax_rate')
            header = 'All rmax change rates (km/hr)'
            flSaveFile(rmax_rate, rmaxRate, header, fmt='%6.2f')

    def _r34(self, r34, indicator):
        """
        Extract radii of 34kt wind:
        :param r34: :class:`numpy.ndarray` of radii to 34-kt winds 
                      for TC observatons.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with a
                          single observation.

        Output: None - data is written to file
        """
        self.logger.info("Extracting radii to 34kt winds")
        initR34   = r34.compress(indicator)
        r34NoInit = r34.compress(indicator == 0)
        r34NoInit = r34NoInit.compress(r34NoInit < sys.maxint)
        if self.ncflag:
            self.data['r34']         = r34
            self.data['init_r34']    = initR34
            self.data['r34_no_init'] = r34NoInit
        else:
            init_r34    = pjoin(self.processPath, 'init_r34')
            all_r34     = pjoin(self.processPath, 'all_r34')
            r34_no_init = pjoin(self.processPath, 'r34_no_init')

            # extract all r34
            self.logger.debug('Outputting data into %s' % all_r34)
            header = 'r34 (km)'
            flSaveFile(all_r34, r34, header, fmt='%6.2f')

            # extract initial r34
            self.logger.debug('Outputting data into %s' % init_r34)
            header = 'initial r34 (km)'
            flSaveFile(init_r34, initR34, header, fmt='%6.2f')

            # extract r34 no init
            self.logger.debug('Outputting data into %s' % r34_no_init)
            header = 'r34 excluding initial ones (km)'
            flSaveFile(r34_no_init, r34NoInit, header, fmt='%6.2f')

    def _r34Rate(self, r34, dt, indicator):
        """Extract the rate of size change from the r34 values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in r34 from the previous observation is
        undefined. Entries corresponding to records with no r34
        observation are also set to maxint.
        :param vmax: :class:`numpy.ndarray` of central vmax observations
                     for TCs
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation
        Input: r34 - array of r34
               dt - array of times between consecutive TC observations
               indicator - array of ones/zeros representing initial TC
                           observations (including TCs with a single
                           observation)
        Output: None - data is written to file
        
        """
        self.logger.info('Extracting the rate of r34 size change')

        # Change in r34:
        r34Change_ = np.diff(r34)
        r34Change = np.empty(indicator.size, 'f')
        r34Change[1:] = r34Change_

        # Rate of r34 change:
        r34Rate = r34Change / dt

        # Mask rates corresponding to initial times and times when
        # the r34 is known to be missing.
        self.logger.debug('Outputting data into %s' %
                          pjoin(self.processPath, 'r34_rate'))
        np.putmask(r34Rate, indicator, sys.maxint)
        np.putmask(r34Rate, r34 >= sys.maxint, sys.maxint)
        np.putmask(r34Rate, (r34Rate >= sys.maxint) |
                            (r34Rate <= -sys.maxint), sys.maxint)
        np.putmask(r34Rate, np.isnan(r34Rate), sys.maxint)
	np.putmask(r34Rate, np.abs(r34Rate) > 150, sys.maxint)

        if self.ncflag:
            self.data['r34Rate'] = r34Rate
        else:
            r34_rate = pjoin(self.processPath, 'r34_rate')
            header = 'All r34 change rates (km/hr)'
            flSaveFile(r34_rate, r34Rate, header, fmt='%6.2f')

    def _r64(self, r64, indicator):
        """Extract radii of 64kt wind:

        :param r64: :class:`numpy.ndarray` of radius to 64-kt winds 
                     from TC observations.
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation

        Output: None - data is written to file
        """
        self.logger.info("Extracting radii to 64kt winds")
        initR64   = r64.compress(indicator)
        r64NoInit = r64.compress(indicator == 0)
        r64NoInit = r64NoInit.compress(r64NoInit < sys.maxint)
        if self.ncflag:
            self.data['r64'] = r64
            self.data['init_r64'] = initR64
            self.data['r64_no_init'] = r64NoInit
        else:
            init_r64 = pjoin(self.processPath, 'init_r64')
            all_r64 = pjoin(self.processPath, 'all_r64')
            r64_no_init = pjoin(self.processPath, 'r64_no_init')

            # extract all r64
            self.logger.debug('Outputting data into {0}'.format(all_r64))
            header = 'r64 (km)'
            flSaveFile(all_r64, r64, header, fmt='%6.2f')

            # extract initial r64
            self.logger.debug('Outputting data into {0}'.format(init_r64))
            header = 'initial r64 (km)'
            flSaveFile(init_r64, initR64, header, fmt='%6.2f')

            # extract r64 no init
            self.logger.debug('Outputting data into {0}'.format(r64_no_init))
            header = 'r64 excluding initial ones (km)'
            flSaveFile(r64_no_init, r64NoInit, header, fmt='%6.2f')

    def _r64Rate(self, r64, dt, indicator):
        """Extract the rate of size change from the r64 values.

        Entries corresponding to initial cyclone reports are set to
        maxint, as the change in r64 from the previous observation is
        undefined. Entries corresponding to records with no r64
        observation are also set to maxint.

        :param r64: :class:`numpy.ndarray` of radius to 64-kt winds 
                     from TC observations.
        :param dt: :class:`numpy.ndarray` of times between consecutive
                   TC observations.
        :param indicator: :class:`numpy.ndarray` of ones/zeros representing
                          initial TC observations, including TCs with
                          a single observation

        Output: None - data is written to file
        
        """
        self.logger.info('Extracting the rate of r64 size change')

        # Change in r64:
        r64Change_    = np.diff(r64)
        r64Change     = np.empty(indicator.size, 'f')
        r64Change[1:] = r64Change_

        # Rate of r64 change:
        r64Rate = r64Change / dt
        
        # Mask rates corresponding to initial times and times when
        # the r64 is known to be missing.
        self.logger.debug('Outputting data into %s' %
                          pjoin(self.processPath, 'r64_rate'))
        np.putmask(r64Rate, indicator, sys.maxint)
        np.putmask(r64Rate, r64 >= sys.maxint,        sys.maxint)
        np.putmask(r64Rate, (r64Rate >= sys.maxint) |
                            (r64Rate <= -sys.maxint), sys.maxint)
        np.putmask(r64Rate, np.isnan(r64Rate),        sys.maxint)
        np.putmask(r64Rate, np.abs(r64Rate) > 150,    sys.maxint)

        if self.ncflag:
            self.data['r64Rate'] = r64Rate
        else:
            r64_rate = pjoin(self.processPath, 'r64_rate')
            header = 'All r64 change rates (km/hr)'
            flSaveFile(r64_rate, r64Rate, header, fmt='%6.2f')



    def _frequency(self, years, indicator):
        """
        Generate a histogram of the annual frequency of events from the input
        data
        """
        self.logger.info('Extracting annual frequency of events')
        minYr = years.min()
        maxYr = years.max()
        genesisYears = years.compress(indicator)
        if minYr == maxYr:
            self.logger.info("First and last year of input data are the same")
            self.logger.info("Cannot generate histogram of frequency")
        else:
            frequency = pjoin(self.processPath, 'frequency')
            bins = np.arange(minYr, maxYr + 2, 1)
            n, b = np.histogram(genesisYears, bins)
            header = 'Year,count'
            flSaveFile(frequency, np.transpose([bins[:-1], n]),
                       header, fmt='%6.2f')

            self.logger.info("Mean annual frequency: {0}".format(np.mean(n)))
            self.logger.info("Standard deviation: {0}".format(np.std(n)))

    def _juliandays(self, jdays, indicator, years):
        """
        Generate a distribution of the formation day of
        year from observations
        """

        self.logger.info("Calculating annual distribution of observations")

        # Do a bodgy job of addressing 29th of February (there surely
        # must be a recommended way of accounting for leap years)

        for i in range(len(jdays)):
            if (years[i] % 4 == 0) and (jdays[i] >= 60):
                jdays[i] -= 1

        bins = np.arange(1, 367)
        n, b = np.histogram(jdays.compress(indicator), bins)
        header = 'Day,count'
        jday_genesis = pjoin(self.processPath, 'jday_genesis')
        jday_observations = pjoin(self.processPath, 'jday_obs')
        jday = pjoin(self.processPath, 'jdays')
        # Distribution of genesis days (histogram):
        flSaveFile(jday_genesis, np.transpose([bins[:-1], n]),
                   header, fmt='%d', delimiter=',')
        n, b = np.histogram(jdays, bins)
        # Distribution of all days (histogram):
        flSaveFile(jday_observations, np.transpose([bins[:-1], n]),
                   header, fmt='%d', delimiter=',')
        # All days:
        flSaveFile(jday, np.transpose(jdays.compress(indicator)),
                   header='Day', fmt='%d')

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename does not exist =>
        # raise error
        if not os.path.exists(configFile):
            error_msg = ("No configuration file specified, please type: "
                         "python tcrm.py -c {config filename}.ini")
            raise IOError(error_msg)
    # If config file does not exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile + "' not found"
        raise IOError(error_msg)

    config = ConfigParser()
    config.read(configFile)

    logFile = config.read('Logging', 'LogFile')
    logLevel = config.read('Logging', 'LogLevel')
    verbose = config.read('Logging', 'Verbose')

    flStartLog(logFile, logLevel, verbose)

    dp = DataProcess(configFile)
    dp.processData()
    logging.shutdown()
