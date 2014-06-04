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

Title: autoPlotHazard.py
Author: Nicholas Summons
Email: nicholas.summons@ga.gov.au
CreationDate: 2011-08-05
Description: Automatically plots hazard maps for each return period and 
             hazard return curves for each locality within the domain.
             Adapted from compareGrids.py code and plotHazardCurves.py 
             developed by Craig Arthur.
"""

import logging
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot

try:
    from mpl_toolkits.basemap import Basemap
    NO_BASEMAP = False
except ImportError:
    NO_BASEMAP = True
    logging.warn('Basemap package not installed. Disabling some plots')

from os.path import join as pjoin

from Utilities.config import ConfigParser

from Utilities.maputils import find_index
import Utilities.nctools as nctools
from Utilities.smooth import smooth
from Utilities import pathLocator
from Utilities import metutils
from Utilities import colours
#from Utilities.progressbar import ProgressBar

from PlotInterface.maps import saveHazardMap

import sqlite3
import unicodedata

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

class PlotUnits(object):
    
    def __init__(self, units):
        labels = {
            'mps': 'metres/second',
            'mph': 'miles/hour',           
            'kts': 'knots',
            'kph': 'kilometres/hour',
            'kmh': 'kilomtres/hour'
        }
        
        levels = {
            'mps': np.arange(30, 101., 5.),
            'mph': np.arange(80, 221., 10.),           
            'kts': np.arange(60, 201., 10.),
            'kph': np.arange(80, 361., 20.),
            'kmh': np.arange(80, 361., 20.)
        }
        
        self.units = units
        self.label = labels[units]
        self.levels = levels[units]

class AutoPlotHazard(object):

    def __init__(self, configFile, progressbar=None):
        
        config = ConfigParser()
        config.read(configFile)
        
        outputPath = config.get('Output','Path')
        
        try:
            self.localityID = config.get('Region', 'LocalityID')
        except Exception:
            self.localityID = -999999
            
        self.inputFile = pjoin(outputPath, 'hazard', 'hazard.nc')
        self.plotPath = pjoin(outputPath, 'plots', 'hazard')
        self.plotUnits = PlotUnits(config.get('Hazard', 'PlotSpeedUnits'))
        
        self.progressbar = progressbar
        
        

    def plotMap(self):
        """Plot return period wind speed maps"""
        
        lon, lat, years, inputData = self.loadFile(self.inputFile, 'wspd')
        [xgrid, ygrid] = np.meshgrid(lon, lat)
        inputData = metutils.convert(inputData, 'mps', self.plotUnits.units)

        map_kwargs = dict(llcrnrlon=xdata.min(),
                          llcrnrlat=ydata.min(),
                          urcrnrlon=xdata.max(),
                          urcrnrlat=ydata.max(),
                          projection='merc',
                          resolution='i')

        for i, year in enumerate(years):
            title = '%d-Year Return Period Cyclonic Wind Hazard' % (year)
            imageFilename = '%d_yrRP_hazard_map.png' % (year)
            filename = pjoin(self.plotPath, imageFilename)
            cbarlab = "Wind speed (%s)"%self.plotUnits.units
            saveHazardMap(inputData[i, :, :], xgrid, ygrid, title, cbarlab, 
                            map_kwargs, filename)

            #self.plotHazardMap(inputData[i, :, :], lon, lat, 
            #                   int(year), self.plotPath)

            self.progressbar.update((i + 1) / float(len(years)), 0.0, 0.9)
            
    def plotCurves(self):
        """Plot hazard curves for speified locations"""                    

        tcrm_dir = pathLocator.getRootDirectory()
        localitiesDataFile = pjoin(tcrm_dir, 'input', 'localities.dat')
        self.sqlcon = sqlite3.connect(localitiesDataFile)
        self.sqlcur = self.sqlcon.cursor()
        
        self.plotHazardCurves(self.inputFile, self.plotPath)
        self.progressbar.update(1.0)

    def loadFile(self, inputFile, varname):
        """
        Load a variable from a netcdf file and return data as a masked array
        
        :param str inputFile: path to a netcdf file containing hazard data.
        :param str varname: name of the netcdf variable to plot.
        
        :returns: lon, lat, years and data (as a masked array)
        """
        
        try:
            ncobj = nctools.ncLoadFile(inputFile)
            lon = nctools.ncGetDims(ncobj, 'lon')
            lat = nctools.ncGetDims(ncobj, 'lat')
            years = nctools.ncGetDims(ncobj, 'years')
            data = nctools.ncGetData(ncobj, varname)
            mv = getattr(ncobj.variables[varname], '_FillValue')
            ncobj.close()
        except:
            self.logger.critical("Cannot load input file: %s"%inputFile)
            try:
                ncobj.close()
            except (IOError, KeyError, RuntimeError):
                pass
            raise

        # Create a masked array:
        mask = (data==mv)
        mdata = ma.array(data, mask=mask)
        return lon, lat, years, mdata

    def plotHazardMap(self, inputData, lon, lat, year, plotPath):
        """
        Plot a hazard map
        
        :param inputData: 2D array of data to plot on a map.
        :type inputData: :class:`numpy.ndarray`
        
        """
        if NO_BASEMAP:
            return

        llLon = lon[0]
        llLat = lat[0]
        urLon = lon[-1]
        urLat = lat[-1]
        res = 'i'
        if (urLon - llLon) > 20:
            dl = 10.
        else:
            dl = 5.
            
        ilon = np.where(((lon > llLon) & (lon < urLon)))[0]
        ilat = np.where(((lat > llLat) & (lat < urLat)))[0]
        [x, y] = np.meshgrid(lon[ilon], lat[ilat])
        
        # FIXME - data is already in a masked array - why mask here?
        eps = 10e-3
        dmask = (inputData < eps)
        inputData.mask[:] = dmask
        inputData.data[:] = metutils.convert(inputData.data, 'mps', 
                                             self.plotUnits.units)
                                             
        inputData.data[:] = smooth(inputData.data, 40)
        
        meridians=np.arange(dl * np.floor(llLon / dl), 
                            dl * np.ceil(urLon / dl), dl)
        parallels=np.arange(dl * np.floor(llLat / dl), 
                            dl * np.ceil(urLat / dl), dl)
        
        pyplot.clf()
        m = Basemap(projection='cyl',
                    resolution=res,
                    llcrnrlon=llLon,
                    urcrnrlon=urLon,
                    llcrnrlat=llLat,
                    urcrnrlat=urLat)
        m.contourf(x, y, inputData[ilat][:, ilon], self.plotUnits.levels, 
                   extend='both', cmap=pyplot.get_cmap('pccspvar'))
                   
        cb = pyplot.colorbar(shrink=0.5, orientation='horizontal',
                            extend='both', ticks=self.plotUnits.levels[::2], 
                            pad=0.05)
                            
        cb.set_label('Maximum gust wind speed (' + self.plotUnits.label + ')',
                     fontsize=10)
                     
        if cb.orientation=='horizontal':
            for t in cb.ax.get_xticklabels():
                t.set_fontsize(8)
        else:
            for t in cb.ax.get_yticklabels():
                t.set_fontsize(8)
        pyplot.title(str(year) + '-Year Return Period Cyclonic Wind Hazard')

        m.drawcoastlines(linewidth=0.5)
        m.drawparallels(parallels, labels=[1, 0, 0, 1], 
                        fontsize=9, linewidth=0.2)
        m.drawmeridians(meridians, labels=[1, 0, 0, 1],
                        fontsize=9, linewidth=0.2)
                        
        pyplot.grid(True)
        imageFilename = str(year) + 'yrRP_hazard_map' + '.png'
        pyplot.savefig(pjoin(plotPath, imageFilename))

    def plotHazardCurves(self, inputFile, plotPath):
        """
        Plot the hazard values stored in hazardFile, at the stns
        stored in stnFile.
        """

        log.info(("Plotting return period curves for locations within the "
                  "model domain"))
        # Open data file
        try:
            ncobj = nctools.ncLoadFile(inputFile)
            lon = nctools.ncGetDims(ncobj, 'lon')
            lat = nctools.ncGetDims(ncobj, 'lat')
            years = nctools.ncGetDims(ncobj, 'years')
        except (IOError, RuntimeError, KeyError):
            log.critical("Cannot load input file: %s"%inputFile)
            raise

        # Load data
        wspd = nctools.ncGetData(ncobj, 'wspd')
        try:
            wLower  = nctools.ncGetData(ncobj, 'wspdlower')
            wUpper = nctools.ncGetData(ncobj, 'wspdupper')
            ciBounds = True
        except KeyError:
            ciBounds = False
        ncobj.close()

        minLon = min(lon)
        maxLon = max(lon) 
        minLat = min(lat) 
        maxLat = max(lat) 

        # If locality is not found in domain, revert to plotting return 
        # curves for all localities in domain:
        self.sqlcur.execute(('select placename from localities where lon > ? '
                             'and lon < ? and lat > ? and lat < ? '
                             'and placeID = ?'), 
                             (minLon, maxLon, 
                              minLat, maxLat, 
                              str(self.localityID)))
                              
        if len([z[0] for z in self.sqlcur.fetchall()]) == 0:
            self.localityID = -99999

        if self.localityID == -99999:
            self.sqlcur.execute(('select placename, parentcountry, lat, lon '
                                 'from localities where lon > ? and lon < ? '
                                 'and lat > ? and lat < ?'), 
                                 (minLon, maxLon, minLat, maxLat))
        else:
            self.sqlcur.execute(('select placename, parentcountry, lat, lon '
                                 'from localities where placeID = ?'), 
                                 (str(self.localityID),))

        placeNames, parentCountries, placeLats, placeLons = \
            zip(*self.sqlcur.fetchall())
        placeNames = list(placeNames)
        parentCountries = list(parentCountries)
        placeLats = list(placeLats)
        placeLons = list(placeLons)

        # Use the same maximum value for all localities to simplify 
        # intercomparisons:
        defaultMax = np.ceil(metutils.convert(100.0, 'mps', 
                                              self.plotUnits.units)/10.0)*10.0

        for name, plat, plon, country in zip(placeNames, placeLats, 
                                             placeLons, parentCountries):

            log.debug("Plotting return period curve for %s"%name)
            
            i = find_index(lon, plon)
            j = find_index(lat, plat)
            
            placeWspd = metutils.convert(wspd[:, j, i], 'mps', 
                                         self.plotUnits.units)
            maxWspd = placeWspd.max()
            if ciBounds:
                placeWspdLower = metutils.convert(wLower[:,j,i], 'mps', 
                                                  self.plotUnits.units)
                placeWspdUpper  = metutils.convert(wUpper[:,j,i], 'mps', 
                                                   self.plotUnits.units)

            pyplot.clf()
            if placeWspd[0] > 0:
                pyplot.semilogx(years, placeWspd, 'r-',  
                                linewidth=2, subsx=years, 
                                label='Return period wind speed')
                if ciBounds:
                    if (placeWspdLower[0] > 0) and (placeWspdUpper[0] > 0):
                        pyplot.fill_between(years, placeWspdUpper, 
                                            placeWspdLower, facecolor='0.75', 
                                            edgecolor='0.99',
                                            alpha=0.7)
                        maxWspd = np.max([maxWspd, placeWspdUpper.max()])
                        
            else:
                continue
            pyplot.xlabel('Return period (years)')
            pyplot.ylabel('Wind speed (' + self.plotUnits.label + ')')
            #years2 = numpy.array([25.0 * (2 ** k) for k in range(7)])
            years2 = np.array([10., 25., 50., 100., 250., 500., 1000., 2500.])
            pyplot.xticks(years2, years2.astype(int))
            
            pyplot.xlim(years.min(), years.max())

            # Override default maximum if exceeded by hazard values
            pyplot.ylim(0.0, np.max([np.ceil(maxWspd/10.0)*10.0, defaultMax]))
            pyplot.title("Return period wind speeds at " + name + ", " \
                            + country + "\n(%5.1f,%5.1f)"%(plon, plat))
            pyplot.grid(True)
            
            # Convert unicode to ascii and remove spaces
            name = unicodedata.normalize('NFKD', name).encode('ascii', 'ignore')
            name.replace(' ', '')
            pyplot.savefig(pjoin(plotPath, 'RP_curve_%s.%s'%(name,"png")))
