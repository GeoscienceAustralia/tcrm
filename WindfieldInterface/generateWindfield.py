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


Title: generateWindfield.py - Generate the wind field associated with an
       instance of a TC.

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-12-01
Description:

SeeAlso: (related programs)
Constraints:

Version: 314
ModifiedBy: N. Habili, nariman.habili@ga.gov.au
ModifiedDate: 2007-01-12
Modification: Speed data written to a super-array.
          Lat/Longs are now integer in order to avoid rounding errors.

Version: 608
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2007-10-18
Modification: Converted to class structure.

Version: 285
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: $LastChangedDate: 2012-02-21 18:52:50 +1100 (Tue, 21 Feb 2012) $
Modification: Super-array limits set by passing dictionary of
              {'xmin','xmax','ymin','ymax'}
              Modified logging method

Version: $Rev: 810 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-10-14 3:32:PM
Modification: If set to save all times, the data (gust wind speed, u and
              v components and sea level pressure) are written to a netcdf
              file at each timestep when the cyclone is within the
              specified region.

$Id: generateWindfield.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

import numpy
import windProfile
import windVorticity
import PressureInterface.pressureProfile as pressureProfile
import windField as wF
import Utilities.maputils as maputils
import Utilities.metutils as metutils
import Utilities.nctools as nctools

import math
import time

from Utilities.grid import grdSave
from Utilities.config import cnfGetIniValue
from Utilities.columns import colReadCSV
from Utilities.timeseries import timeseries
from PlotInterface.plotWindfield import plotWindfield


class generateWindfield:

    def __init__(self, config_file):

        self.config_file = config_file

        self.logger = logging.getLogger()
        self.logger.info("Initialising GenerateWindfield")

        self.profileType = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                          'profileType', 'holland')
        self.windFieldType = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                            'windFieldType', 'kepert')

        self.logger.debug('Using %s radial profile and %s wind field'
                           %(self.profileType, self.windFieldType))

        self.beta = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                   'beta', 1.3)
        self.beta1 = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                    'beta1', self.beta)
        self.beta2 = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                    'beta2', self.beta)
        self.thetaMax = math.radians(cnfGetIniValue(self.config_file,
                                                    'WindfieldInterface',
                                                    'thetaMax', 70.0))
        self.title = 'Cyclone wind field model'
        self.info = 'Profile: %s \nAsymmetry: %s' \
                    %(self.profileType, self.windFieldType)
        self.margin = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                     'Margin', 2)
        self.resolution = cnfGetIniValue(self.config_file, 'WindfieldInterface',
                                         'Resolution', 0.01)

        gridLimitStr = cnfGetIniValue(self.config_file, 'WindfieldInterface', 
                                      'gridLimit', '')

        if gridLimitStr is '':
            gridLimitStr = cnfGetIniValue(self.config_file, 'Region', 'gridLimit')

        try:
            self.gL = eval(gridLimitStr)
        except SyntaxError:
            self.logger.exception('Error! gridLimit is not a dictionary' )

        self.outputPath = os.path.join(cnfGetIniValue(self.config_file, 'Output', 'Path'),
                                       'windfield')
                                       
        if cnfGetIniValue(self.config_file, 'Timeseries', 'Extract', False):
            tsPath = os.path.join(cnfGetIniValue(self.config_file, 'Output', 'Path'), 
                                  'process', 'timeseries')
                                  
            self.logger.info('Timeseries data will be collected')
            self.ts = timeseries(config_file, outputPath=tsPath)

    def calcwind(self, data, saveAllTimes=False):
        """
        Calculate the windfield around a track given the parameters of
        position, motion, intensity and size.
        """

        [index, age, cLon, cLat, vFm, thetaFm, pCentre, pEnv, rMax] = \
            [data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5],
             data[:,6], data[:,7], data[:,8]]

        # Correct format of TC position to eliminate rounding errors
        # when matching the parent grid and each timestep grid:
        cLon_ = 100*cLon
        cLat_ = 100*cLat
        cLon_ = cLon_.astype(int)
        cLat_ = cLat_.astype(int)

        # Determine bounds of parent grid:
        #minLat = self.gL['yMin']
        #maxLat = self.gL['yMax']
        #minLon = self.gL['xMin']
        #maxLon = self.gL['xMax']

        minLat_ = (self.gL['yMin']-self.margin)*100
        minLon_ = (self.gL['xMin']-self.margin)*100
        maxLat_ = (self.gL['yMax']+self.margin)*100
        maxLon_ = (self.gL['xMax']+self.margin)*100

        cGridLat = numpy.arange(minLat_, maxLat_+self.resolution,
                                self.resolution*100, dtype=int)
        cGridLon = numpy.arange(minLon_, maxLon_+self.resolution,
                                self.resolution*100, dtype=int)

        [cGridX, cGridY] = numpy.meshgrid(cGridLon, cGridLat)

        # We determine the times for which the cyclone is within the domain
        # so that we minimise the size of the resulting output file.
        # This is due to a short-coming of Scientific.IO.NetCDF, which does
        # not permit files greater than 2GB.
        ii=[]
        for i in range(len(index)):
            if ((cLon[i] >= self.gL['xMax']) or
                (cLon[i] <= self.gL['xMin']) or
                (cLat[i] >= self.gL['yMax']) or
                (cLat[i] <= self.gL['yMin'])):
                continue
            else:
                ii.append(i)
        if saveAllTimes:
            self.logger.info("Saving output at all timesteps")
            self.logger.info("Output file will contain fields at all times when cyclone is within domain")
            missingValue = -9999.
            dimensions = {0:{'name':'time', 'values':age[ii], 'dtype':'f',
                                'atts':{'long_name':'Time', 'units':'hours'} },
                          1:{'name':'lat','values':cGridLat/100., 'dtype':'f',
                                'atts':{'long_name':'Latitude',
                                        'units':'degrees_north'} },
                          2:{'name':'lon','values':cGridLon/100., 'dtype':'f',
                                'atts':{'long_name':'Longitude',
                                        'units':'degrees_east'} } }

            variables = {0:{'name':'vmax', 'dims':('time', 'lat', 'lon'),
                            'values':numpy.array(missingValue), 'dtype':'d',
                            'atts':{'long_name':'3-second gust wind speed',
                                    'units':'m/s'} },
                        1:{'name':'ua', 'dims':('time', 'lat', 'lon'),
                            'values':numpy.array(missingValue), 'dtype':'d',
                            'atts':{'long_name':'Eastward wind',
                                    'units':'m/s'} },
                         2:{'name':'va', 'dims':('time', 'lat', 'lon'),
                            'values':numpy.array(missingValue), 'dtype':'d',
                            'atts':{'long_name':'Northward wind',
                                    'units':'m/s'} },
                        3:{'name':'psl', 'dims':('time', 'lat', 'lon'),
                            'values':numpy.array(missingValue), 'dtype':'d',
                            'atts':{'long_name':'Air pressure at sea level',
                                    'units':'hPa'} } }

            # Create the netcdf file, and return the netcdf object so we can
            # write data to the file:
            outputFile = os.path.join(self.outputPath,'tcrm.nc')
            self.logger.info("Creating output file: %s"%outputFile)
            ncobj = nctools.ncSaveGrid(outputFile, dimensions, variables,
                                       nodata=missingValue, writedata=False,
                                       dtype='d', keepfileopen=True)

            # Create the variable objects:
            gust_varobj = ncobj.variables['vmax']
            ua_varobj = ncobj.variables['ua']
            va_varobj = ncobj.variables['va']
            psl_varobj = ncobj.variables['psl']


        UU = numpy.empty(numpy.shape(cGridX), float)
        VV = numpy.empty(numpy.shape(cGridY), float)
        bearing = numpy.empty(numpy.shape(cGridX), float)
        speed = numpy.zeros(numpy.shape(cGridX), float)
        gust = numpy.zeros(numpy.shape(cGridX), float)
        

        # The 'if' statement is required to prevent an error when there are no storms for a given simulation
        if len(pEnv) > 0:
            pressure = pEnv[0]*numpy.ones(numpy.shape(cGridX), float)
        else:
            pressure = numpy.NaN

        self.logger.debug('Shape of speed array is: %4d by %4d ' % numpy.shape(speed))
        self.logger.debug('Minimum longitude: %7.3f' % (cGridLon.min()/100.))
        self.logger.debug('Minimum latitude: %7.3f' % (cGridLat.min()/100.))

        times = numpy.zeros(len(index))
        ii = 0
        for i in range(len(index)):
            times[i] = time.time()
            if ((cLon[i] >= self.gL['xMax']) or
                (cLon[i] <= self.gL['xMin']) or
                (cLat[i] >= self.gL['yMax']) or
                (cLat[i] <= self.gL['yMin'])):
                self.logger.debug("Cyclone centre is outside grid. Skipping timestep")
                if hasattr(self,'ts'):
                    self.ts.extract(None, None, None, pEnv[i], cGridLon,
                                    cGridLat, age[i])

                continue
            else:
                self.logger.debug('Calculating wind field at timestep: %s' %i )
                self.logger.debug('Forward speed: %03f | direction %03f | Pcentre %06.1f | Penv %06.1f | rMax %04f' %(vFm[i], thetaFm[i], pCentre[i], pEnv[i], rMax[i]))
                R, lam = maputils.makeGrid(cLon[i], cLat[i], self.margin,
                                           self.resolution)
                xgrid, ygrid = maputils.meshLatLon(cLon[i], cLat[i],
                                                   self.margin,
                                                   self.resolution)
                f = metutils.coriolis(cLat[i])

                profile = windProfile.WindProfile(R, pEnv[i], pCentre[i],
                                                  rMax[i], cLat[i], cLon[i],
                                                  self.beta, beta1=self.beta1,
                                                  beta2=self.beta2)
                vorticity = windVorticity.WindVorticity(R, pEnv[i], pCentre[i],
                                                        rMax[i], cLat[i],
                                                        cLon[i], self.beta,
                                                        beta1=self.beta1,
                                                        beta2=self.beta2)
                prsprofile = pressureProfile.PrsProfile(R, pEnv[i], pCentre[i], 
                                                        rMax[i], cLat[i], cLon[i],
                                                        self.beta, beta1=self.beta1,
                                                        beta2=self.beta2)
                try:
                    wind_ = getattr(profile, self.profileType)
                except AttributeError:
                    self.logger.exception('%s not implemented in windProfile'
                                           %self.windFieldType)

                try:
                    vort_ = getattr(vorticity, self.profileType)
                except AttributeError:
                    self.logger.exception('%s not implemented in windVorticity'
                                           %self.windFieldType)

                try:
                    pressure_ = getattr(prsprofile, self.profileType)
                except AttributeError:
                    self.logger.exception( '%s not implemented in pressureProfile' %self.windFieldType )

                V_ = wind_()
                Z_ = vort_(abs(V_).max())
                P_ = pressure_()
                self.logger.debug('Maximum speed: %05.1f'%(abs(V_).max()))

                field_ = wF.WindField(R, lam, rMax[i], f, V_, Z_, vFm[i],
                                      thetaFm[i], self.thetaMax)
                try:
                    field = getattr(field_, self.windFieldType)
                except AttributeError:
                    self.logger.exception('%s not implemented in windField'
                                           %self.windFieldType)

                Ux, Vy = field()
                # Add gust factor (very elementary at this stage):
                Ux *= 1.38
                Vy *= 1.38

                sp = numpy.sqrt(Ux**2 + Vy**2)
                bear = ((numpy.arctan2(-Ux, -Vy))*180./numpy.pi)

                jmin = int(int(cLat_[i]-minLat_-int(100*self.margin)) \
                           /int(self.resolution*100))
                jmax = int(int(cLat_[i]-minLat_+int(100*self.margin)) \
                           /int(self.resolution*100))+1
                imin = int(int(cLon_[i]-minLon_-int(100*self.margin)) \
                           /int(self.resolution*100))
                imax = int(int(cLon_[i]-minLon_+int(100*self.margin)) \
                           /int(self.resolution*100))+1

                if saveAllTimes:
                    self.logger.debug("Saving data to file...")
                    gust_varobj[ii, jmin:jmax, imin:imax] = sp
                    ua_varobj[ii, jmin:jmax, imin:imax] = Ux
                    va_varobj[ii, jmin:jmax, imin:imax] = Vy
                    psl_varobj[ii, jmin:jmax, imin:imax] = P_


                speed_ = (speed[jmin:jmax, imin:imax]).copy()

                mask = sp > speed_

                numpy.putmask(speed_, mask, sp)
                speed[jmin:jmax, imin:imax] = speed_

                bear_ = (bearing[jmin:jmax, imin:imax]).copy()
                numpy.putmask(bear_, mask, bear)
                bearing[jmin:jmax, imin:imax] = bear_

                UU_ = (UU[jmin:jmax, imin:imax]).copy()
                numpy.putmask(UU_, mask, Ux)
                UU[jmin:jmax, imin:imax] = UU_

                VV_ = (VV[jmin:jmax, imin:imax]).copy()
                numpy.putmask(VV_, mask, Vy)
                VV[jmin:jmax, imin:imax] = VV_

                prs_ = (pressure[jmin:jmax, imin:imax]).copy()
                pmask = P_ < prs_
                numpy.putmask(prs_, pmask, P_)
                pressure[jmin:jmax, imin:imax] = prs_

                self.logger.debug("Timestep took %.3f seconds"
                                   %(time.time() - times[i]))

                if hasattr(self, 'ts'):
                    self.ts.extract(sp, Ux, Vy, P_, xgrid[0, :], ygrid[:, 0], age[i])


                ii += 1

        if saveAllTimes:
            # Set the actual_range attribute and close the netcdf file:
            setattr(gust_varobj, 'actual_range', 
                    [gust_varobj[:].min(), gust_varobj[:].max()])
            setattr(ua_varobj, 'actual_range', 
                    [ua_varobj[:].min(), ua_varobj[:].max()])
            setattr(va_varobj, 'actual_range', 
                    [va_varobj[:].min(), va_varobj[:].max()])
            setattr(psl_varobj, 'actual_range', 
                    [psl_varobj[:].min(), psl_varobj[:].max()])
            ncobj.close()

        # The 'if' statement is required to prevent an error when there are no storms for a given simulation
        if len(times) > 0:
            self.logger.debug('Windfield calculation took %.2f seconds'%(time.time()-times[0]))

        if hasattr(self, 'ts'):
            self.ts.shutdown()

        return speed, UU, VV, pressure, cGridLon/100., cGridLat/100.

if __name__ == "__main__":
    try:
        config_file = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        config_file = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(config_file):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError, error_msg
    # If config file doesn't exist => raise error
    if not os.path.exists(config_file):
        error_msg = "Configuration file '" + config_file +"' not found"
        raise IOError, error_msg

    logging.basicConfig(level=getattr(logging, cnfGetIniValue(config_file, 
                                                              'Logging', 
                                                              'LogLevel', 
                                                              'INFO')),
                        format='%(asctime)s %(name)-15s: %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=cnfGetIniValue(config_file, 'Logging', 'LogFile',
                                                __file__.rstrip('.py') + '.log'),
                        filemode='w')
                        
    simulationName = cnfGetIniValue(config_file, 'Simulation', 'Name')
    outputPath = cnfGetIniValue(config_file, 'Output', 'Path')
    GW = generateWindfield(config_file)

    source = cnfGetIniValue(config_file, 'WindfieldInterface', 'Source')
    input_file = cnfGetIniValue(config_file, 'WindfieldInterface', 'TrackFile')
    track_data = colReadCSV(config_file, input_file, source)

    index = numpy.array(track_data['index'], int)
    age = numpy.array(track_data['age'], float)
    lon = numpy.array(track_data['lon'], float)
    lat = numpy.array(track_data['lat'], float)
    speed = numpy.array(track_data['speed'], float)
    penv = numpy.array(track_data['penv'], float)
    pressure = numpy.array(track_data['pressure'], float)
    rmax = numpy.array(track_data['rmax'], float)
    try:
        speedUnits = cnfGetIniValue(config_file, source, 'SpeedUnits', 'mps')
    except NameError:
        pass
    else:
        speed = metutils.convert(speed, speedUnits, 'mps')

    try:
        prUnits = cnfGetIniValue(config_file, source, 'PressureUnits', 'hPa')
    except NameError:
        pass
    else:
        pressure = metutils.convert(pressure, prUnits, 'Pa')
        penv = metutils.convert(penv, prUnits, 'Pa')

    # At this time, assume the storm motion direction is given in
    # degrees true.  Thus the storm motion direction must be converted
    # to radians.
    bearing = numpy.array(track_data['bearing'])*numpy.pi/180.
    bearing = maputils.bearing2theta(bearing)


    data = numpy.transpose(numpy.array([index, age, lon, lat, speed, bearing,
                            pressure, penv, rmax]))
    speed, UU, VV, lon, lat = GW.calcwind(data, False, False)
    [gridX, gridY] = numpy.meshgrid(lon, lat)
    plotWindfield(gridX, gridY, numpy.flipud(speed), title="Windfield",
                  fileName=os.path.join(outputPath, 'windfield.png'))
    grdSave(os.path.join(outputPath, 'gust.txt'), speed, lon, lat,
            0.01, fmt='%4.1f')
    grdSave(os.path.join(outputPath, 'uu.txt'), UU, lon, lat,
            0.01, fmt='%4.1f')
    grdSave(os.path.join(outputPath, 'vv.txt'), VV, lon, lat,
            0.01, fmt='%4.1f')
