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

Title: autoPlotHazard.py
Author: Nicholas Summons
Email: nicholas.summons@ga.gov.au
CreationDate: 2011-08-05
Description: Automatically plots hazard maps for each return period and 
             hazard curves for each WMO station within domain.
             Adapted from compareGrids.py code and plotHazardCurves.py developed by Craig Arthur.
"""

import os, sys, pdb, logging
import numpy
import pylab
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma
from matplotlib import pyplot

from Utilities.files import flConfigFile, flStartLog
from Utilities.config import cnfGetIniValue
from Utilities.colours import colourMap
import Utilities.nctools as nctools
from Utilities.smooth import smooth
from Utilities import pathLocator
from Utilities import metutils
from Utilities.progressbar import ProgressBar
import csv

class AutoPlotHazard:

    def __init__(self, configFile):
        self.logger = logging.getLogger()
        self.outputPath = cnfGetIniValue(configFile,'Output','Path')
        self.inputFile = os.path.join(self.outputPath, 'hazard', 'hazard.nc')
        self.plotPath = os.path.join(self.outputPath, 'plots', 'hazard')
        self.margin = cnfGetIniValue(configFile, 'WindfieldInterface', 'Margin', 2.0)
        self.plotSpeedUnits = cnfGetIniValue(configFile, 'HazardInterface', 'PlotSpeedUnits', 'mps')
        if self.plotSpeedUnits == 'mps':
            self.speedunitlabel = 'metres / second'
        elif self.plotSpeedUnits == 'mph':
            self.speedunitlabel = 'miles / hour'
        elif self.plotSpeedUnits == 'kph':
            self.speedunitlabel = 'kilometres / hour'
        elif self.plotSpeedUnits == 'kts':
            self.speedunitlabel = 'knots'
        else:
            self.speedunitlabel = self.plotSpeedUnits        
        self.pbar = ProgressBar('(6/6) Plotting results:      ')

    def plot(self):
            lon, lat, years, inputData = self._loadFile(self.inputFile, 'wspd')
            
            # Plot return period wind speed maps
            for record in range(inputData.shape[0]):
                self._plotHazardMap(inputData[record,:,:], lon, lat, int(years[record]), self.plotPath)
                self.pbar.update((record + 1) / float(inputData.shape[0]), 0.0, 0.9)
            
            # Plot return period curves
            self._plotHazardCurves(self.inputFile, self.plotPath, self.margin)
            self.pbar.update(1.0)

    def _loadFile(self, inputFile, varname):
        try:
            ncobj = nctools.ncLoadFile(inputFile)
            lon = nctools.ncGetDims(ncobj,'lon')
            lat = nctools.ncGetDims(ncobj,'lat')
            years = nctools.ncGetDims(ncobj,'years')
            data = nctools.ncGetData(ncobj,varname)
            mv = getattr(ncobj.variables[varname],'_FillValue')
            ncobj.close()
        except:
            self.logger.critical("Cannot load input file: %s"%inputFile)
            try:
                ncobj.close()
            except:
                pass
            raise

        # Create a masked array:
        mask = (data==mv)
        mdata = ma.array(data,mask=mask)
        return lon,lat,years,mdata

    def _plotHazardMap(self, inputData, lon, lat, year, plotPath):
        llLon = lon[0]
        llLat = lat[0]
        urLon = lon[-1]
        urLat = lat[-1]
        res = 'i'
        if (urLon - llLon) > 20:
            dl = 10.
        else:
            dl = 5.
        ilon = numpy.where(((lon>llLon) & (lon<urLon)))[0]
        ilat = numpy.where(((lat>llLat) & (lat<urLat)))[0]
        [x,y] = numpy.meshgrid(lon[ilon],lat[ilat])
        if self.plotSpeedUnits == 'mps':
            levels = pylab.arange(30, 101., 5.)
        elif self.plotSpeedUnits == 'kts':
            levels = pylab.arange(60, 201., 10.)
        elif self.plotSpeedUnits == 'kph':
            levels = pylab.arange(80, 361., 20.)
        elif self.plotSpeedUnits == 'mph':
            levels = pylab.arange(80, 221., 10.)
        eps = 10e-3
        dmask = (inputData<eps)
        inputData.mask[:] = dmask
        inputData.data[:] = metutils.convert(inputData.data, 'mps', self.plotSpeedUnits)
        inputData.data[:] = smooth(inputData.data, 40)
        meridians=numpy.arange(dl*numpy.floor(llLon/dl),dl*numpy.ceil(urLon/dl),dl)
        parallels=numpy.arange(dl*numpy.floor(llLat/dl),dl*numpy.ceil(urLat/dl),dl)
        pylab.clf()
        m = Basemap(projection='cyl',
                    resolution=res,
                    llcrnrlon=llLon,
                    urcrnrlon=urLon,
                    llcrnrlat=llLat,
                    urcrnrlat=urLat)
        m.contourf(x,y,inputData[ilat][:,ilon],levels,extend='both',cmap=colourMap('pccspvar','stretched'))
        cb = pylab.colorbar(shrink=0.5,orientation='horizontal',extend='both',ticks=levels[::2],pad=0.05)
        cb.set_label('Maximum gust wind speed (' + self.speedunitlabel + ')',fontsize=10)
        if cb.orientation=='horizontal':
            for t in cb.ax.get_xticklabels():
                t.set_fontsize(8)
        else:
            for t in cb.ax.get_yticklabels():
                t.set_fontsize(8)
        pylab.title(str(year) + '-Year Return Period Cyclonic Wind Hazard')
        coastlinewidth = 0.5
        m.drawcoastlines(linewidth=coastlinewidth)
        m.drawparallels(parallels,labels=[1,0,0,1],fontsize=9,linewidth=0.2)
        m.drawmeridians(meridians,labels=[1,0,0,1],fontsize=9,linewidth=0.2)
        pylab.grid(True)
        imageFilename = str(year) + 'yrRP_hazard_map' + '.png'
        pylab.savefig(os.path.join(plotPath, imageFilename))

    def _plotHazardCurves(self, inputFile, plotPath, margin):
        """
        Plot the hazard values stored in hazardFile, at the stns
        stored in stnFile.
        """

        # Open data file
        try:
            ncobj = nctools.ncLoadFile(inputFile)
            lon = nctools.ncGetDims(ncobj, 'lon')
            lat = nctools.ncGetDims(ncobj, 'lat')
            years = nctools.ncGetDims(ncobj, 'years')
        except:
            self.logger.critical("Cannot load input file: %s"%inputFile)
            raise

        # Load data
        wspd = nctools.ncGetData(ncobj, 'wspd')
        try:
            w05  = nctools.ncGetData(ncobj, 'wspd05')
            w95 = nctools.ncGetData(ncobj, 'wspd95')
            ciBounds = True
        except:
            ciBounds = False
        ncobj.close()

        # Crop region since windfields are not generated when storms
        # are within the margin distance from domain edges
        minLon = min(lon) + (margin * 2.0)
        maxLon = max(lon) - (margin * 2.0)
        minLat = min(lat) + (margin * 2.0)
        maxLat = max(lat) - (margin * 2.0)
        stnDict = self._getStations()

        # Use the same maximum value for all stations to simplify station intercomparisons
        defaultMax = numpy.ceil(metutils.convert(100.0, 'mps', self.plotSpeedUnits)/10.0)*10.0

        for stn in stnDict:
            stnlon = stnDict[stn][0]
            stnlat = stnDict[stn][1]

            if (stnlon>minLon) and (stnlon<maxLon) and (stnlat>minLat) and (stnlat<maxLat):
                self.logger.info("Plotting return period curve for %s"%stn)
                i = numpy.where(stnlon<=lon)[0][0]
                j = numpy.where(stnlat<=lat)[0][0]
                stnWspd = metutils.convert(wspd[:,j,i], 'mps', self.plotSpeedUnits)
                maxWspd = stnWspd.max()
                if ciBounds:
                    stnWspd05 = metutils.convert(w05[:,j,i], 'mps', self.plotSpeedUnits)
                    stnWspd95  = metutils.convert(w95[:,j,i], 'mps', self.plotSpeedUnits)

                pyplot.clf()
                if stnWspd[0] > 0:
                    pyplot.semilogx(years, stnWspd, 'b-', linewidth=2, subsx=years, label='Return period wind speed')
                    if stnWspd05[0] > 0 and ciBounds:
                        pyplot.semilogx(years, stnWspd05, 'k--', linewidth=1, subsx=years, label='5th percentile')
                    if stnWspd95[0] > 0 and ciBounds:
                        pyplot.semilogx(years,stnWspd95, 'k--', linewidth=1, subsx=years, label='95th percentile')
                        maxWspd = numpy.max(maxWspd, stnWspd95.max())
                else:
                    continue
                pyplot.xlabel('Return period (years)')
                pyplot.ylabel('Wind speed (' + self.speedunitlabel + ')')
                years2 = numpy.array([25.0 * (2 ** k) for k in range(7)])
                pyplot.xticks(years2,years2.astype(int))
                
                pyplot.xlim(years.min(),years.max())

                # Override default maximum if exceeded by hazard values
                pyplot.ylim(0.0, max(numpy.ceil(maxWspd/10.0)*10.0, defaultMax))
                pyplot.title("Return period wind speeds at " + stnDict[stn][2] + ", " + stnDict[stn][3] + "\n(%5.1f,%5.1f)"%(stnlon,stnlat))
                pyplot.grid(True)
                pyplot.savefig(os.path.join(plotPath, 'RP_curve_Station%s.%s'%(stn,"png")))

    def _getStations(self):
        # Use station file located in TCRM input directory
        tcrm_dir = pathLocator.getRootDirectory()
        tcrm_input_dir = os.path.join(tcrm_dir, 'input')
        stnFile = os.path.join(tcrm_input_dir, 'station_list.txt')
        
        csvReader = csv.reader(open(stnFile, 'rb'), delimiter=';')
        stnDict = {}  
        for row in csvReader:
            stnlat_str = row[7]
            stnlat_str = stnlat_str.replace('-', ':')
            if stnlat_str[-1] == 'N':
                stnlat_str = stnlat_str[0:-1]
            elif stnlat_str[-1] == 'S':
                stnlat_str = '-' + stnlat_str[0:-1]
            else:
                continue

            stnlon_str = row[8]
            stnlon_str = stnlon_str.replace('-', ':')
            if stnlon_str[-1] == 'E':
                stnlon_str = stnlon_str[0:-1]
            elif stnlon_str[-1] == 'W':
                stnlon_str = '-' + stnlon_str[0:-1]
            else:
                continue

            degminsec = numpy.double(stnlat_str.split(':'))
            stnlat = 0
            for k in range(len(degminsec)):
                stnlat = stnlat + abs(degminsec[k] / 60.0**k)
            stnlat = stnlat * numpy.sign(degminsec[0])

            degminsec = numpy.double(stnlon_str.split(':'))
            stnlon = 0
            for k in range(len(degminsec)):
                stnlon = stnlon + abs(degminsec[k] / 60.0**k)
            stnlon = stnlon * numpy.sign(degminsec[0])
            stnlon = numpy.mod(stnlon, 360.0)
            if int(row[0] + row[1]) > 1:
                stnDict[int(row[0] + row[1])] = [stnlon, stnlat, row[3], row[5]]
        return stnDict