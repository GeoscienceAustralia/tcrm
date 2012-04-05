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


Title: TrackGenerator.py

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-11-28
Description: Generate synthetic cyclone tracks

SeeAlso: (related programs)
Constraints:
Version: $Rev: 811 $

ModifiedDate: 2007-02-19
ModifiedBy: N. Habili
Modifications: Modifications to __init__, generatePath, singleTrack,
getFileSampleParameters, and eliminateStep in order to make the class
operable.
Conformance with style guide

Version: 520
ModifiedDate: 2007-05-14
ModifiedBy: N. Habili
Modifications: Made singleTrack, getFileSampleParameters, and
eliminateStep private members.

Version: 522
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-07-11
Modification: Output directory is set separately in config file.
              Sample origins using SamplingOrigin() method.
              Origin can also be set manually.

Version: 551
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-09-10 3:40:PM
Modification: Can optionally set initial self.pressure, self.speed and
self.bearing for the cyclone

Version: 640
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-12-07 3:24:PM
Modification: The generation of synthetic tracks now uses a method based
              on Hall, T.M. & Jewson, S. (2007): Statistical modelling
              of North Atlantic tropical cyclone tracks.  Tellus 59A,
              486-498. The changes in parameters are a combination of
              the previous timestep value and a random variation.
              Provides improved distributions of the parameters.

Version: 46
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-04-01 11:29:AM
Modification: Corrected error in calculation of decay rate following
              landfall.

Version: Rev: 131
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-07-02 9:42:AM
Modification: Replaced references to 'files' dictionary with
              configuration file. Allows this to be integrated as a set
              of functions into the rest of TCRM, rather than running as
              a stand-alone feature.

Version: 152
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-04-02 10:26:14
Modification: Optionally writes tracks to shapefiles, set by suffix on
              trackFile.
              i.e. if the suffix is .shp, then a shapefile (and
              corresponding dbf file) is created, otherwise a csv file
              is created.
Version: 217
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-05-01 11:52:AM
Modification: Uses threading to save data in separate threads

Version: 276
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-10-08 12:05:AM
Modification: Changed the calculation of the size parameter to ensure no
              negative sizes are generated

Version: 312
ModifiedBy: Nicholas Summons, nicholas.summons@ga.gov.au
ModifiedDate: 2010-06-04 05:17:PM
Modification: Added option to generate MSLP seasonal averages during run time

Version: 370
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-09-03
Modification: Separated the generation of the autocorrelations for each parameter
              (speed, bearing, central pressure change and size change) into separate
              functions. Because of the loops testing resulting values for each parameter
              were calling the same generic function, all parameters would be updated.
              This would be equivalent to stepping some variables multiple time steps
              rather than one time step.

Version: $Rev: 811 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2011-06-07 10:53:AM
Modification: Save coefficients to netCDF file to permit further analysis.
              Saves mu, sig, alpha and min for the bearing, speed, pressure and
              pressure rate of change variables, over water and over land.

$Id: TrackGenerator.py 811 2012-02-24 05:10:11Z carthur $
"""
import os, sys, pdb, logging

from scipy import array, empty, flatnonzero, transpose, concatenate
import math
from numpy import mod, random, arange

import Utilities.stats as stats
import StatInterface.SamplingOrigin as SamplingOrigin
import StatInterface.SamplingParameters as SamplingParameters
import StatInterface.generateStats as generateStats
import trackLandfall
import Utilities.nctools as nctools

import Utilities.Cmap as Cmap
import Utilities.Cstats as Cstats
from Utilities.AsyncRun import AsyncRun
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile, flSaveFile
from Utilities.grid import SampleGrid
from MSLP.mslp_seasonal_clim import MSLPGrid
from Utilities.shptools import shpSaveTrackFile

__version__ = '$Id: TrackGenerator.py 811 2012-02-24 05:10:11Z carthur $'

class TrackGenerator:
    """
    Description: Generates a random cyclone track based on the
    distributions of self.speed, self.bearing, and central self.pressure
    given a cyclone origin.

    Parameters:
    dt : float timestep (in hours) for cyclone simulation
    tsteps : integer maximum number of timesteps to simulate.
    Other information is determined from the configuration file

    Members:

    Methods:

    Internal Methods:

    """

    def __init__(self, configFile, dt=1, tsteps=360, autoCalc_gridLimit=None, progressbar=None):
        """Initialise required fields"""
        self.configFile = configFile
        self.logger = logging.getLogger()
        self.logger.info("Initialising TrackGenerator")

        self.outputPath = cnfGetIniValue(self.configFile, 'Output', 'Path')
        self.processPath = os.path.join(self.outputPath, 'process')

        gridLimitStr = cnfGetIniValue(self.configFile, 'TrackGenerator', 'gridLimit', '')
        if gridLimitStr is not '':
            try:
                self.gridLimit = eval(gridLimitStr)
            except SyntaxError:
                self.logger.exception('Error! gridLimit is not a dictionary' )
        else:
            self.gridLimit = autoCalc_gridLimit
            self.logger.info('No gridLimit specified - using automatic selection: ' + str(self.gridLimit))

        self.setInnerDomainFreq = cnfGetIniValue(configFile, 'TrackGenerator', 'SetInnerDomainFreq', False)

        if self.setInnerDomainFreq:
            # Load gridLimit for inner domain
            self.logger.info('Setting TC frequency on inner domain')
            gridLimitStr = cnfGetIniValue(self.configFile, 'Region', 'gridLimit', '')
            if gridLimitStr is not '':
                try:
                    self.glInner = eval(gridLimitStr)
                except SyntaxError:
                    self.logger.exception('Error! gridLimit is not a dictionary')
            else:
                self.logger.error('Error! No gridLimit for inner domain found in config file.')

        self.gridSpace = eval(cnfGetIniValue(self.configFile, 'TrackGenerator',
                                             'gridSpace'))
        self.gridInc = eval(cnfGetIniValue(self.configFile, 'TrackGenerator',
                                           'gridInc'))
        self.missingValue = cnfGetIniValue(self.configFile, 'StatInterface',
                                           'MissingValue', sys.maxint)



        self.dt = dt
        self.tsteps = tsteps
        self.timeOverflow = dt * tsteps

        self.trackLandfall = trackLandfall.LandfallDecay(self.configFile, self.dt)
        mslp_configsetting = cnfGetIniValue(self.configFile, 'Input', 'MSLPGrid')

        # Check whether MSLPGrid config setting is a list of month numbers for run time generation or
        # a pre-generated file.
        try:
            mnth_sel = set(mslp_configsetting.strip('[]{}() ').replace(',', ' ').split(' '))
            mnth_sel.discard('')
            if mnth_sel.issubset([str(k) for k in range(1,13)]):
                MSLP_type = 1
                mnth_sel_int = [int(k) for k in mnth_sel]
            else:
                MSLP_type = 2
        except:
            MSLP_type = 2

        if MSLP_type == 1:
            # Calculate MSLP seasonal average during run time
            self.logger.info("Generating MSLP seasonal average")
            self.mslp = MSLPGrid(mnth_sel_int)
        elif MSLP_type == 2:
            # Load MSLP seasonal average from file
            self.logger.info("Loading MSLP seasonal average from file")
            self.mslp = SampleGrid(mslp_configsetting)
        else:
            self.logger.info("Critical Error")

        self.Origin = SamplingOrigin.SamplingOrigin(os.path.join(self.processPath, 'originPDF.nc'), None, None)

        self.lon = empty(self.tsteps, 'd')
        self.lat = empty(self.tsteps, 'd')
        self.index = empty(self.tsteps, 'd')
        self.bearing = empty(self.tsteps, 'd')
        self.speed = empty(self.tsteps, 'd')
        self.dist = empty(self.tsteps, 'd')
        self.pressure = empty(self.tsteps, 'd')
        self.penv = empty(self.tsteps,'d')
        self.rmax = empty(self.tsteps, 'd')
        self.elapsedTime = empty(self.tsteps, 'i')

        try:
            #obtain cyclone parameter data for the cell number
            self.allCDFInitBearing = flLoadFile(os.path.join(self.processPath, 'all_cell_cdf_init_bearing'), '%', ',')
            self.allCDFInitSpeed = flLoadFile(os.path.join(self.processPath, 'all_cell_cdf_init_speed'), '%', ',')
            self.allCDFInitPressure = flLoadFile(os.path.join(self.processPath, 'all_cell_cdf_init_pressure'), '%', ',')

            """
            self.allCDFBearingRate = flLoadFile(self.processPath+'all_cell_cdf_bearing_rate','%',',')
            self.allCDFSpeedRate = flLoadFile(self.processPath+'all_cell_cdf_speed_rate','%',',')
            self.allCDFPressureRate = flLoadFile(self.processPath+'all_cell_cdf_pressure_rate','%',',')
            """
        except IOError:
            self.logger.critical('Error! CDF distribution files does not exist!')
            self.logger.critical('Run AllDistribution option in main to generate those files.')
            raise

        try:
            self.allCDFInitSize = flLoadFile(os.path.join(self.processPath, 'all_cell_cdf_init_rmax'), '%', ',')
        except IOError:
            self.allCDFInitSize = None
            distMean = cnfGetIniValue(configFile, 'RMW', 'mean', -1)
            distSigma = cnfGetIniValue(configFile, 'RMW', 'sigma', -1)
            if (distMean == -1) or (distSigma == -1):
                distMean = 57.0
                distSigma = 0.6
                # Ref: McConochie, J.D., T.A. Hardy and L.B. Mason, 2004: Modelling tropical cyclone over-water wind and pressure fields.
                #      Ocean Engineering, 31, 1757-1782
                self.logger.warning('No RMW distribution data or lognormal parameters given. Using published distributions from McConochie et al. 2004.')

            self.cdfSize = transpose(array(stats.rMaxDist(distMean, distSigma, maxrad=120.0)))
            pass
        else:
            #self.sStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath, 'all_rmax'),
            #                                          os.path.join(self.processPath, 'all_lon_lat'),
            #                                          self.gridLimit, self.gridSpace,
            #                                          self.gridInc, minSample=100,
            #                                          missingValue=self.missingValue)

            self.dsStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath,'rmax_rate'),
                                                      os.path.join(self.processPath,'all_lon_lat'),
                                                      self.gridLimit, self.gridSpace,
                                                      self.gridInc, minSample=100,
                                                      missingValue=self.missingValue)

        self.vStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath, 'all_speed'),
                                                  os.path.join(self.processPath, 'all_lon_lat'),
                                                  self.gridLimit, self.gridSpace,
                                                  self.gridInc, minSample=100,
                                                  missingValue=self.missingValue, 
                                                  progressbar=progressbar, prgStartValue=0, prgEndValue=0.075)

        self.pStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath, 'all_pressure'),
                                                  os.path.join(self.processPath, 'all_lon_lat'),
                                                  self.gridLimit, self.gridSpace,
                                                  self.gridInc, minSample=100,
                                                  missingValue=self.missingValue, 
                                                  progressbar=progressbar, prgStartValue=0.075, prgEndValue=0.15)

        self.bStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath, 'all_bearing'),
                                                  os.path.join(self.processPath, 'all_lon_lat'),
                                                  self.gridLimit, self.gridSpace, self.gridInc,
                                                  minSample=100, angular=True,
                                                  missingValue=self.missingValue,
                                                  progressbar=progressbar, prgStartValue=0.15, prgEndValue=0.225)

        # These instances of GenerateStats are for the rates of change.
        #self.dvStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath,'speed_rate'),
        #                                           os.path.join(self.processPath,'all_lon_lat'),
        #                                           self.gridLimit, self.gridSpace, self.gridInc,
        #                                           minSample=100,missingValue=self.missingValue)

        self.dpStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath,'pressure_rate'),
                                                   os.path.join(self.processPath,'all_lon_lat'),
                                                   self.gridLimit, self.gridSpace, self.gridInc,
                                                   minSample=100, missingValue=self.missingValue, 
                                                   progressbar=progressbar, prgStartValue=0.225, prgEndValue=0.30)

        #self.dbStats = generateStats.GenerateStats(self.configFile, os.path.join(self.processPath,'bearing_rate'),
        #                                           os.path.join(self.processPath,'all_lon_lat'),
        #                                           self.gridLimit, self.gridSpace, self.gridInc,
        #                                           minSample=100,angular=True, missingValue=self.missingValue)
        self._saveData()

    def __doc__(self):
        """Documentation on the function of the class"""
        return 'Generate synthetic cyclone track and other parameters'

    def generatePath(self, nEvents=1, outputFile=None, initLon=None,
                     initLat=None, initSpeed=None, initBearing=None,
                     initPressure=None):
        """generatePath(self, [outputFile, [oLon, oLat]])
        Cyclone movement path in a particular cyclone origin based on
        random walk.  By specifying [initLon, initLat], the origin of
        all self.numCyclones events will be the defined point
        """
        self.logger.debug('Generating %d tropical cyclone tracks'%nEvents)
        
        # Initialise array to prevent error when nEvents=0
        results = []
        
        for j in xrange(1, nEvents + 1):
            if (j % (nEvents/100.)) == 0. and j != 0.:
                self.logger.debug('Processed %i cyclones'%j)

            index = []
            #while age[-1] < 12:
            while len(index) < 2:
                if initLon is None:
                    oLon, oLat = self.Origin.generateOneSample()
                    while (oLon <= self.gridLimit['xMin'] or
                           oLon >= self.gridLimit['xMax'] or
                           oLat <= self.gridLimit['yMin'] or
                           oLat >= self.gridLimit['yMax']):
                        oLon, oLat = self.Origin.generateOneSample()
                else:
                    oLon = initLon
                    oLat = initLat
                self.currentCell = None
                self.logger.debug("Cyclone origin: %6.2f, %6.2f"%(oLon,oLat))
                index, age, lon, lat, speed, bearing, pressure, \
                penv, rmax = self._singleTrack(j, oLon, oLat, initSpeed,
                                               initBearing, initPressure)
                if age[-1] < 12:
                    index = []

                if self.setInnerDomainFreq:
                    insideInnerDomain = [lon[k] > self.glInner['xMin'] and \
                                         lon[k] < self.glInner['xMax'] and \
                                         lat[k] > self.glInner['yMin'] and \
                                         lat[k] < self.glInner['yMax'] \
                                         for k in range(len(lon))]
                    if not any(insideInnerDomain):
                        # Reject track if no points pass within inner domain
                        index = []


            if j == 1:
                results = transpose(array([index, age, lon, lat, speed,
                                           bearing, pressure, penv, rmax]))
            else:
                results = concatenate((results,
                                       transpose(array([index, age, lon, lat,
                                                        speed, bearing,
                                                        pressure, penv,
                                                        rmax]))))

            self.logger.debug("Completed track generation for simulation %d"%j)
        if outputFile is None:
            return results
        elif outputFile.endswith("shp"):
            self.logger.debug('Outputting data into %s'%outputFile)
            fields = {}
            fields['Index'] = {'Type':1, 'Length':5, 'Precision':0, 'Data':results[:, 0]}
            fields['Time'] = {'Type':2, 'Length':7, 'Precision':1, 'Data':results[:, 1]}
            fields['Longitude'] = {'Type':2, 'Length':7, 'Precision':2, 'Data':results[:, 2]}
            fields['Latitude'] = {'Type':2, 'Length':7, 'Precision':2, 'Data':results[:, 3]}
            fields['Speed'] = {'Type':2, 'Length':6, 'Precision':1, 'Data':results[:, 4]}
            fields['Bearing'] = {'Type':2, 'Length':6, 'Precision':1, 'Data':results[:, 5]}
            fields['Pressure'] = {'Type':2, 'Length':6, 'Precision':1, 'Data':results[:, 6]}
            fields['pEnv'] = {'Type':2, 'Length':6, 'Precision':1, 'Data':results[:, 7]}
            fields['rMax'] = {'Type':2, 'Length':5, 'Precision':1, 'Data':results[:, 8]}

            args = {"filename":outputFile, "lon":results[:,2], "lat":results[:,3], "fields":fields}
            thr = AsyncRun(shpSaveTrackFile, args)
            try:
                thr.start()
            except:
                raise
            #shpSaveTrackFile(outputFile, results[:,2], results[:,3], fields)

        else:
            self.logger.debug('Outputting data into %s'%outputFile)
            header = 'CycloneNumber,TimeElapsed(hr),Longitude(degree),Latitude(degree),Speed(km/hr),Bearing(degrees),CentralPressure(hPa),EnvPressure(hPa),rMax(km)'
            args = {"filename":outputFile, "data":results, "header":header,
                    "delimiter":',', "fmt":'%7.2f'}
            fl = AsyncRun(flSaveFile, args)
            fl.start()


    def _singleTrack(self, index, oLon, oLat, initSpeed=None,
                     initBearing=None, initPressure=None):
        """Generate one cyclone path using random walk method"""
        self.logger.debug('Generating a single TC track')
        try:
            cellNum = Cstats.getCellNum(oLon, oLat, self.gridLimit,
                                        self.gridSpace)
        except ValueError:
            self.logger.critical("TC origin is outside of grid: %6.2fE, %6.2fS"%(oLon, oLat))
            raise
        self.logger.debug("Cyclone starts in cell %d"%cellNum)
        # We need to get this first:
        self.penv[0] = self.mslp.sampleGrid(oLon, oLat)

        if initBearing is None:
            ind = self.allCDFInitBearing[:, 0] == cellNum
            cdfInitBearing = self.allCDFInitBearing[ind, 1:3]
            initBearing = SamplingParameters.SamplingParameters(cdfInitBearing).generateOneSample()

        if initSpeed is None:
            ind = self.allCDFInitSpeed[:, 0] == cellNum
            cdfInitSpeed = self.allCDFInitSpeed[ind, 1:3]
            initSpeed = SamplingParameters.SamplingParameters(cdfInitSpeed).generateOneSample()

        if initPressure is None:
            ind = self.allCDFInitPressure[:, 0] == cellNum
            cdfInitPressure = self.allCDFInitPressure[ind, 1:3]
            initPressure = self.penv[0]
            # A little gotcha - the initial self.pressure sampled was greater than the environmental self.pressure.
            # While the TC would be immediately terminated, the initial timestep would still be included
            # in the list of events.
            while initPressure >= self.penv[0]:
                initPressure = SamplingParameters.SamplingParameters(cdfInitPressure).generateOneSample()

        if self.allCDFInitSize is None:
            # Sample from an artificial distribution of cyclone sizes:
            initRmax = SamplingParameters.SamplingParameters(self.cdfSize[:, [0, 2]]).generateOneSample()
        else:
            ind = self.allCDFInitSize[:, 0] == cellNum
            cdfSize = self.allCDFInitSize[ind, 1:3]
            initRmax = SamplingParameters.SamplingParameters(cdfSize).generateOneSample()

        self.lon[0] = oLon
        self.lat[0] = oLat
        self.index[0] = index
        self.elapsedTime[0] = 0
        self.bearing[0] = initBearing
        self.speed[0] = initSpeed
        self.pressure[0] = initPressure
        self.dist[0] = self.dt * initSpeed
        self.rmax[0] = initRmax

        self.vChi = 0.0
        self.bChi = 0.0
        self.pChi = 0.0
        self.sChi = 0.0
        self.dpChi = 0.0
        self.dsChi = 0.0

        self.theta = initBearing
        self.dp = 0.0
        self.v = initSpeed
        self.ds = 0.0
        self.offshorePressure = initPressure

        for i in xrange(1, self.tsteps):
            self.index[i] = index
            self.lon[i], self.lat[i] = Cmap.bear2LatLon(self.bearing[i-1], self.dist[i-1], self.lon[i-1], self.lat[i-1])
            self.penv[i] = self.mslp.sampleGrid(self.lon[i], self.lat[i])

            if (self.lon[i] < self.gridLimit['xMin'] or
                self.lon[i] >= self.gridLimit['xMax'] or
                self.lat[i] <= self.gridLimit['yMin'] or
                self.lat[i] > self.gridLimit['yMax']):
                return self.index[:i], self.elapsedTime[:i], self.lon[:i], \
                       self.lat[:i], self.speed[:i], self.bearing[:i], \
                       self.pressure[:i], self.penv[:i], self.rmax[:i]

            # Cell number:
            cellNum = Cstats.getCellNum(self.lon[i], self.lat[i], self.gridLimit,self.gridSpace)
            # Is the cyclone over land?
            on_land = self.trackLandfall.onLand(self.lon[i], self.lat[i])

            #self._getCoeffs(self.lon[i], self.lat[i], i)

            self._getPressureCoeffs(cellNum, i, on_land)
            self._getBearingCoeffs(cellNum, i, on_land)
            self._getSpeedCoeffs(cellNum, i, on_land)
            if self.allCDFInitSize is not None:
                self._getSizeCoeffs(cellNum, i, on_land)

            self.bearing[i] = self.theta
            self.speed[i] = self.v
            self.pressure[i] = self.pressure[i-1] + self.dp*self.dt

            # Check for negative self.speed:
            while self.speed[i] < 0.0:
                self.logger.debug("Speed is less than zero - recalculating")
                self._getSpeedCoeffs(cellNum, i, on_land)
                self.speed[i] = self.v

            # Central pressure:
            if on_land:
                self.pressure[i] = self.trackLandfall.pChange(self.offshorePressure, self.penv[i])
                self.logger.debug( "Central pressure: %7.2f"%self.pressure[i])
            else:
                self.pressure[i] = self.pressure[i-1] + self.dp*self.dt
                if (self.pressure[i] < (self.pStats.coeffs.min[cellNum]-4.*self.pStats.coeffs.sig[cellNum])):
                    # If the central pressure of the synthetic storm is more than 4 std deviations lower than the
                    # minimum observed central pressure, automatically start raising the central pressure.
                    self.logger.debug("Pressure is extremely low - recalculating")
                    #pdb.set_trace()
                    #self._getPressureCoeffs(cellNum, i, on_land)
                    self.pressure[i] = self.pressure[i-1] + abs(self.dp)*self.dt
                self.offshorePressure = self.pressure[i]

            # Radius to maximum winds:
            if self.allCDFInitSize is not None:
                self.rmax[i] = self.rmax[i-1] + self.ds*self.dt
                while self.rmax[i] <= 1.0:
                    self.logger.debug("Size is less than one - recalculating")
                    self._getSizeCoeffs(cellNum, i, on_land)
                    self.rmax[i] = self.rmax[i-1] + self.ds*self.dt

            else:
                self.rmax[i] = self.rmax[i-1]

            self.dist[i] = self.dt * self.speed[i]
            self.elapsedTime[i] = self.elapsedTime[i-1] + self.dt
            if self._eliminateStep(self.pressure[i], self.penv[i],
                                   self.elapsedTime[i], self.lon[0],
                                   self.lat[0], self.lon[i], self.lat[i]):
                return self.index[:i], self.elapsedTime[:i], self.lon[:i], \
                       self.lat[:i], self.speed[:i], self.bearing[:i], \
                       self.pressure[:i], self.penv[:i], self.rmax[:i]

        return self.index, self.elapsedTime, self.lon, self.lat, self.speed, \
               self.bearing, self.pressure, self.penv, self.rmax

    def _getCoeffs(self, lon, lat, i):
        cellNum = Cstats.getCellNum(lon, lat, self.gridLimit, self.gridSpace)
        if self.trackLandfall.onLand(self.lon[i], self.lat[i]):
            self.bChi = self.bStats.coeffs.lalpha[cellNum]*self.bChi + \
                        random.normal()*self.bStats.coeffs.lphi[cellNum]
            self.dpChi = self.dpStats.coeffs.lalpha[cellNum]*self.dpChi + \
                        random.normal()*self.dpStats.coeffs.lphi[cellNum]
            self.vChi = self.vStats.coeffs.lalpha[cellNum]*self.vChi + \
                        random.normal()*self.vStats.coeffs.lphi[cellNum]
            if i == 1:
                self.theta += math.degrees(self.bStats.coeffs.lsig[cellNum]*self.bChi)
                self.v += abs(self.vStats.coeffs.lsig[cellNum]*self.vChi)
                self.dp += self.dpStats.coeffs.lsig[cellNum]*self.dpChi

            else:
                self.v = abs(self.vStats.coeffs.lmu[cellNum] + \
                         self.vStats.coeffs.lsig[cellNum]*self.vChi)
                self.theta = math.degrees(self.bStats.coeffs.lmu[cellNum] + \
                             self.bStats.coeffs.lsig[cellNum]*self.bChi)
                self.dp = self.dpStats.coeffs.lmu[cellNum] + \
                         self.dpStats.coeffs.lsig[cellNum]*self.dpChi
        else:
            # TC is over water
            self.bChi = self.bStats.coeffs.alpha[cellNum]*self.bChi + \
                        random.normal()*self.bStats.coeffs.phi[cellNum]
            self.vChi = self.vStats.coeffs.alpha[cellNum]*self.vChi + \
                        random.normal()*self.vStats.coeffs.phi[cellNum]
            self.dpChi = self.dpStats.coeffs.alpha[cellNum]*self.dpChi + \
                        random.normal()*self.dpStats.coeffs.phi[cellNum]
            if i == 1:
                self.theta += math.degrees(self.bStats.coeffs.sig[cellNum]*self.bChi)
                self.v += abs(self.vStats.coeffs.sig[cellNum]*self.vChi)
                self.dp += self.dpStats.coeffs.sig[cellNum]*self.dpChi
            else:
                self.v = abs(self.vStats.coeffs.mu[cellNum] + \
                         self.vStats.coeffs.sig[cellNum]*self.vChi)
                self.theta = math.degrees(self.bStats.coeffs.mu[cellNum] + \
                             self.bStats.coeffs.sig[cellNum]*self.bChi)
                self.dp = self.dpStats.coeffs.mu[cellNum] + \
                         self.dpStats.coeffs.sig[cellNum]*self.dpChi

        self.theta = mod(self.theta, 360)


        # Only change the value of the rMax coefficients if we have data:
        if self.allCDFInitSize is not None:
            if self.trackLandfall.onLand(self.lon[i], self.lat[i]):
                self.dsChi = self.dsStats.coeffs.lalpha[cellNum]*self.dsChi + \
                            random.normal()*self.dsStats.coeffs.lphi[cellNum]
                if i == 1:
                    self.ds += self.dsStats.coeffs.lsig[cellNum]*self.dsChi
                else:
                    self.ds = self.dsStats.coeffs.lmu[cellNum] + \
                             self.dsStats.coeffs.lsig[cellNum]*self.dsChi
            else:
                self.dsChi = self.dsStats.coeffs.alpha[cellNum]*self.dsChi + \
                            random.normal()*self.dsStats.coeffs.phi[cellNum]
                if i == 1:
                    self.ds += self.dsStats.coeffs.sig[cellNum]*self.dsChi
                else:
                    self.ds = self.dsStats.coeffs.mu[cellNum] + \
                             self.dsStats.coeffs.sig[cellNum]*self.dsChi


    def _getPressureCoeffs(self, cellNum, i, on_land):
        """
        Retrieve the coefficients for pressure in cell specified by (lon, lat)
        """
        if on_land:
            self.dpChi = self.dpStats.coeffs.lalpha[cellNum]*self.dpChi + \
                         random.normal()*self.dpStats.coeffs.lphi[cellNum]
            if i == 1:
                self.dp += self.dpStats.coeffs.lsig[cellNum]*self.dpChi
            else:
                self.dp = self.dpStats.coeffs.lmu[cellNum] + \
                         self.dpStats.coeffs.lsig[cellNum]*self.dpChi
        else:
            # TC is over water
            self.dpChi = self.dpStats.coeffs.alpha[cellNum]*self.dpChi + \
                        random.normal()*self.dpStats.coeffs.phi[cellNum]
            if i == 1:
                self.dp += self.dpStats.coeffs.sig[cellNum]*self.dpChi
            else:
                self.dp = self.dpStats.coeffs.mu[cellNum] + \
                         self.dpStats.coeffs.sig[cellNum]*self.dpChi

    def _getBearingCoeffs(self, cellNum, i, on_land):
        """
        Retrieve the coefficients for bearing in cell specified by (lon, lat)
        """
        if on_land:
            # TC is over land:
            self.bChi = self.bStats.coeffs.lalpha[cellNum]*self.bChi + \
                        random.normal()*self.bStats.coeffs.lphi[cellNum]
            if i == 1:
                self.theta += math.degrees(self.bStats.coeffs.lsig[cellNum]*self.bChi)
            else:
                self.theta = math.degrees(self.bStats.coeffs.lmu[cellNum] + \
                             self.bStats.coeffs.lsig[cellNum]*self.bChi)

        else:
            # TC is over water:
            self.bChi = self.bStats.coeffs.alpha[cellNum]*self.bChi + \
                        random.normal()*self.bStats.coeffs.phi[cellNum]
            if i == 1:
                self.theta += math.degrees(self.bStats.coeffs.sig[cellNum]*self.bChi)
            else:
                self.theta = math.degrees(self.bStats.coeffs.mu[cellNum] + \
                             self.bStats.coeffs.sig[cellNum]*self.bChi)
        self.theta = mod(self.theta,360.)

    def _getSpeedCoeffs(self, cellNum, i, on_land):
        """
        Retrieve coefficients for speed in cell specified by (lon, lat).
        """
        if on_land:
            # TC is over land:
            self.vChi = self.vStats.coeffs.lalpha[cellNum]*self.vChi + \
                        random.normal()*self.vStats.coeffs.lphi[cellNum]
            if i == 1:
                self.v += abs(self.vStats.coeffs.lsig[cellNum]*self.vChi)
            else:
                self.v = abs(self.vStats.coeffs.lmu[cellNum] + \
                         self.vStats.coeffs.lsig[cellNum]*self.vChi)
        else:
            # TC is over water:
            self.vChi = self.vStats.coeffs.alpha[cellNum]*self.vChi + \
                        random.normal()*self.vStats.coeffs.phi[cellNum]
            if i == 1:
                self.v += abs(self.vStats.coeffs.sig[cellNum]*self.vChi)
            else:
                self.v = abs(self.vStats.coeffs.mu[cellNum] + \
                         self.vStats.coeffs.sig[cellNum]*self.vChi)

    def _getSizeCoeffs(self, cellNum, i, on_land):
        """
        Retrieve coefficients for size in cell specified by (lon, lat).
        """
        if on_land:
            self.dsChi = self.dsStats.coeffs.lalpha[cellNum]*self.dsChi + \
                            random.normal()*self.dsStats.coeffs.lphi[cellNum]
            if i == 1:
                self.ds += self.dsStats.coeffs.lsig[cellNum]*self.dsChi
            else:
                self.ds = self.dsStats.coeffs.lmu[cellNum] + \
                             self.dsStats.coeffs.lsig[cellNum]*self.dsChi
        else:
            self.dsChi = self.dsStats.coeffs.alpha[cellNum]*self.dsChi + \
                            random.normal()*self.dsStats.coeffs.phi[cellNum]
            if i == 1:
                self.ds += self.dsStats.coeffs.sig[cellNum]*self.dsChi
            else:
                self.ds = self.dsStats.coeffs.mu[cellNum] + \
                             self.dsStats.coeffs.sig[cellNum]*self.dsChi


    def _getFileSampleParameters(self, lon0, lat0):
        """
        Obtains a random sample from cdf distribution files for cyclone
        parameters in each random walk step
        """
        # find the cell number given the cyclone location
        cellNum = Cstats.getCellNum(lon0, lat0, self.gridLimit, self.gridSpace)

        if self.currentCell != cellNum:  #the cyclone has stepped into a different grid
            # extract the cdf data for a particular cell from the files
            ind = self.allCDFBearingRate[:, 0] == cellNum
            self.cdfBearingRate = self.allCDFBearingRate[ind, 1:3]

            ind = self.allCDFSpeedRate[:, 0] == cellNum
            self.cdfSpeedRate = self.allCDFSpeedRate[ind, 1:3]

            ind = self.allCDFPressureRate[:, 0] == cellNum
            self.cdfPressureRate = self.allCDFPressureRate[ind, 1:3]

            self.currentCell = cellNum

        bearingRate = SamplingParameters.SamplingParameters(self.cdfBearingRate).generateOneSample()
        speedRate = SamplingParameters.SamplingParameters(self.cdfSpeedRate).generateOneSample()
        pressureRate = SamplingParameters.SamplingParameters(self.cdfPressureRate).generateOneSample()

        return bearingRate, speedRate, pressureRate

    def _eliminateStep(self, pressure, penv, elapsedTime, lon0, lat0, nextlon,
                       nextlat):
        """Cyclone life over if certain conditions are met assumptions
        that have been placed to end a cyclone life
        1. central self.pressure risen over ambient self.pressure value
           specified by user
        2. cyclone has existed more than time_overflow value specified
           by user
        3. cyclone eye has moved more than distance_overflow value from
           origin specified by user
        """
        if penv - pressure < 5.0:
            return True

        if pressure > 1000.0:
            return True

        return False

    def _saveData(self):
        """
        Save coefficients to a netcdf file to permit further analysis.

        """
        lon = arange(self.gridLimit['xMin'],self.gridLimit['xMax'],self.gridSpace['x'])
        lat = arange(self.gridLimit['yMax'],self.gridLimit['yMin'],-1*self.gridSpace['y'])

        nx = len(lon)
        ny = len(lat)

        dimensions = {0:{'name':'lat','values':lat,'dtype':'f',
                         'atts':{'long_name':'Latitude','units':'degrees_north'} },
                      1:{'name':'lon','values':lon,'dtype':'f',
                         'atts':{'long_name':'Longitude','units':'degrees_east'} } }

        variables = {0:{'name':'vmu','dims':('lat','lon'),
                        'values':self.vStats.coeffs.mu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean forward speed',
                                'units':'m/s'} },
                     1:{'name':'valpha','dims':('lat','lon'),
                        'values':self.vStats.coeffs.alpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of forward speed',
                                'units':''} },
                     2:{'name':'vsig','dims':('lat','lon'),
                        'values':self.vStats.coeffs.sig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation forward speed',
                                'units':'m/s'} },
                     3:{'name':'vmin','dims':('lat','lon'),
                        'values':self.vStats.coeffs.min.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum forward speed',
                                'units':'m/s'} },
                     4:{'name':'vlmu','dims':('lat','lon'),
                        'values':self.vStats.coeffs.lmu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean forward speed (over land)',
                                'units':'m/s'} },
                     5:{'name':'vlalpha','dims':('lat','lon'),
                        'values':self.vStats.coeffs.lalpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of forward speed (over land)',
                                'units':''} },
                     6:{'name':'vlsig','dims':('lat','lon'),
                        'values':self.vStats.coeffs.lsig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of forward speed (over land)',
                                'units':'m/s'} },
                     7:{'name':'vlmin','dims':('lat','lon'),
                        'values':self.vStats.coeffs.lmin.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum forward speed (over land)',
                                'units':'m/s'} },

                     8:{'name':'bmu','dims':('lat','lon'),
                        'values':self.bStats.coeffs.mu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean bearing',
                                'units':'degrees'} },
                     9:{'name':'balpha','dims':('lat','lon'),
                        'values':self.bStats.coeffs.alpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of bearing',
                                'units':''} },
                     10:{'name':'bsig','dims':('lat','lon'),
                        'values':self.bStats.coeffs.sig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of bearing',
                                'units':'degrees'} },
                     11:{'name':'bmin','dims':('lat','lon'),
                        'values':self.bStats.coeffs.min.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum bearing',
                                'units':'degrees'} },
                     12:{'name':'blmu','dims':('lat','lon'),
                        'values':self.bStats.coeffs.lmu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean bearing(over land)',
                                'units':'degrees'} },
                     13:{'name':'blalpha','dims':('lat','lon'),
                        'values':self.bStats.coeffs.lalpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of bearing (over land)',
                                'units':''} },
                     14:{'name':'blsig','dims':('lat','lon'),
                        'values':self.bStats.coeffs.lsig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of bearing (over land)',
                                'units':'degrees'} },
                     15:{'name':'blmin','dims':('lat','lon'),
                        'values':self.bStats.coeffs.lmin.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum bearing (over land)',
                                'units':'degrees'} },

                     16:{'name':'pmu','dims':('lat','lon'),
                        'values':self.pStats.coeffs.mu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean central pressure',
                                'units':'hPa'} },
                     17:{'name':'palpha','dims':('lat','lon'),
                        'values':self.pStats.coeffs.alpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of central pressure',
                                'units':''} },
                     18:{'name':'psig','dims':('lat','lon'),
                        'values':self.pStats.coeffs.sig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of central pressure',
                                'units':'hPa'} },
                     19:{'name':'pmin','dims':('lat','lon'),
                        'values':self.pStats.coeffs.min.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum central pressure',
                                'units':'hPa'} },
                     20:{'name':'plmu','dims':('lat','lon'),
                        'values':self.pStats.coeffs.lmu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean central pressure (over land)',
                                'units':'hPa'} },
                     21:{'name':'plalpha','dims':('lat','lon'),
                        'values':self.pStats.coeffs.lalpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of central pressure (over land)',
                                'units':''} },
                     22:{'name':'plsig','dims':('lat','lon'),
                        'values':self.pStats.coeffs.lsig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of central pressure (over land)',
                                'units':'hPa'} },
                     23:{'name':'plmin','dims':('lat','lon'),
                        'values':self.pStats.coeffs.lmin.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum central pressure (over land)',
                                'units':'hPa'} },

                     24:{'name':'dpmu','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.mu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean rate of pressure change',
                                'units':'hPa/h'} },
                     25:{'name':'dpalpha','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.alpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of rate of pressure change',
                                'units':''} },
                     26:{'name':'dpsig','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.sig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of rate of pressure change',
                                'units':'hPa/h'} },
                     27:{'name':'dpmin','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.min.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum rate of pressure change',
                                'units':'hPa/h'} },
                     28:{'name':'dplmu','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.lmu.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Mean rate of pressure change (over land)',
                                'units':'hPa/h'} },
                     29:{'name':'dplalpha','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.lalpha.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Lag-1 autocorrelation of rate of pressure change (over land)',
                                'units':''} },
                     30:{'name':'dplsig','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.lsig.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Standard deviation of rate of pressure change (over land)',
                                'units':'hPa/h'} },
                     31:{'name':'dplmin','dims':('lat','lon'),
                        'values':self.dpStats.coeffs.lmin.reshape((ny,nx)),
                        'dtype':'d',
                        'atts':{'long_name':'Minimum rate of pressure change (over land)',
                                'units':'hPa/h'} } }

        outputFile = os.path.join(self.processPath,'coefficients.nc')
        nctools._ncSaveGrid(outputFile, dimensions, variables,
                            nodata=self.missingValue,datatitle=None,dtype='f',
                            writedata=True, keepfileopen=False)

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg

    flStartLog(cnfGetIniValue(configFile, 'Logging', 'LogFile', __file__.rstrip('.py') + '.log'),
               cnfGetIniValue(configFile, 'Logging', 'LogLevel', 'DEBUG'),
               cnfGetIniValue(configFile, 'Logging', 'Verbose', True))

    outputPath = cnfGetIniValue(configFile, 'Output', 'Path', os.getcwd())
    trackFile = os.path.join(outputPath, 'sample_tracks.200.shp')
    numCyclones = 50
    dt = 1
    tsteps = 360

    initLon = 155.
    initLat = -25.
    initBearing = 225.
    initSpeed = 5.
    initPressure = 950.
    tracks = TrackGenerator(configFile)
    tracks.generatePath(numCyclones, trackFile)  #, initLon,initLat,initSpeed,initBearing,initPressure)
"""

    gridLimit = eval(cnfGetIniValue(configFile,'Parameters','gridLimit'))
    gridSpace = eval(cnfGetIniValue(configFile,'Parameters','gridSpace'))
    gridInc = eval(cnfGetIniValue(configFile,'Parameters','gridInc'))
    gridLimit = eval(cfgMain['Parameters.gridLimit'])
    gridSpace = eval(cfgMain['Parameters.gridSpace'])
    gridInc = eval(cfgMain['Parameters.gridInc'])
""""""
    import profile
    import pstats
    #profile.run("tracks.generatePath(cfgFiles['Output.sample_tracks'])", 'fooprof')

    profile.run("tracks.generatePath(cfgMain['Output.path'] + 'sample_tracks_200.txt')", 'fooprof')
    p = pstats.Stats('fooprof')
    p.sort_stats('cumulative').print_stats(10)
    p.sort_stats('cumulative').strip_dirs().print_callees()

    print 'total time = ', time.time() - start
"""
