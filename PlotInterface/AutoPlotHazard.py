"""
:mod:`AutoPlotHazard` -- Plot hazard maps and curves
====================================================

.. module:: AutoPlotHazard
    :synopsis: Plot hazard maps and curves for locations within the
               simulation domain.


.. moduleauthor:: Nicholas Summons <nicholas.summons@ga.gov.au>

Automatically plots hazard maps for each return period and hazard
return curves for each locality within the domain. Adapted from
compareGrids.py code and plotHazardCurves.py developed by Craig
Arthur.

"""

import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('shapely.geos').setLevel(logging.WARNING)

import database
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')

from os.path import join as pjoin

from Utilities.config import ConfigParser

from Utilities.maputils import find_index
from Utilities.smooth import smooth
import Utilities.nctools as nctools
from Utilities import pathLocator
from Utilities import metutils

from PlotInterface.maps import saveHazardMap
from PlotInterface.curves import saveHazardCurve

import sqlite3

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

class PlotUnits(object):

    def __init__(self, units):
        labels = {
            'mps': 'm/s',
            'mph': 'mi/h',
            'kts': 'kts',
            'kph': 'km/h',
            'kmh': 'km/h'
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

        outputPath = config.get('Output', 'Path')

        try:
            self.localityID = config.get('Region', 'LocalityID')
        except Exception:
            self.localityID = -999999

        self.inputFile = pjoin(outputPath, 'hazard', 'hazard.nc')
        self.plotPath = pjoin(outputPath, 'plots', 'hazard')
        self.plotUnits = PlotUnits(config.get('Hazard', 'PlotSpeedUnits'))
        self.ciBounds = config.getboolean('Hazard', 'CalculateCI')
        self.fit = config.get('Hazard', 'ExtremeValueDistribution')
        self.numsimulations = config.getint("TrackGenerator", "NumSimulations")
        self.smooth = config.getboolean("Hazard", "SmoothPlots")

        self.progressbar = progressbar

        self.db = database.HazardDatabase(configFile)

    def plotMap(self):
        """Plot return period wind speed maps"""

        lon, lat, years, inputData = self.loadFile(self.inputFile, 'wspd')
        if lon.min() > 180.:
            lon = lon - 360.

        [xgrid, ygrid] = np.meshgrid(lon, lat)
        dmask = inputData.mask
        inputData = metutils.convert(inputData, 'mps', self.plotUnits.units)
        inputData = ma.array(inputData, mask=dmask)
        map_kwargs = dict(llcrnrlon=xgrid.min(),
                          llcrnrlat=ygrid.min(),
                          urcrnrlon=xgrid.max(),
                          urcrnrlat=ygrid.max(),
                          projection='merc',
                          resolution='i')

        for i, year in enumerate(years):
            log.debug("Plotting %d-year return period hazard map", year)
            title = '%d-Year ARI Cyclonic Wind Hazard' % (year)
            imageFilename = '%d_yrRP_hazard_map.png' % (year)
            filename = pjoin(self.plotPath, imageFilename)
            cbarlab = "Wind speed (%s)"%self.plotUnits.units
            levels = self.plotUnits.levels
            if self.smooth:
                dx = np.mean(np.diff(xgrid))
                data = smooth(inputData[i, :, :], int(1/dx))
            else:
                data = inputData[i, :, :]

            saveHazardMap(data, xgrid, ygrid, title, levels,
                          cbarlab, map_kwargs, filename)

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
            years = nctools.ncGetDims(ncobj, 'ari')
            data = nctools.ncGetData(ncobj, varname)
            mv = getattr(ncobj.variables[varname], '_FillValue')
            ncobj.close()
            del ncobj
        except:
            log.critical("Cannot load input file: %s"%inputFile)
            try:
                ncobj.close()
            except (IOError, KeyError, RuntimeError):
                pass
            raise

        # Create a masked array:
        mask = (data == mv)
        mdata = ma.array(data, mask=mask)
        return lon, lat, years, mdata

    def getLocations(self):
        """
        Extract locations from the localities database

        :returns: Names, ID, latitude and longitude of all locations
                  within the model domain, plus the locations rec array
        """
        locations = self.db.getLocations()
        placeNames = list(locations['locName'])
        placeLats = list(locations['locLat'])
        placeLons = list(locations['locLon'])
        placeID = list(locations['locId'])

        return placeNames, placeID, placeLats, placeLons, locations

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
            years = nctools.ncGetDims(ncobj, 'ari')
        except (IOError, RuntimeError, KeyError):
            log.critical("Cannot load input file: %s"%inputFile)
            raise

        placeNames, placeID, placeLats, placeLons, locations = self.getLocations()

        for name, plat, plon, pID in zip(placeNames, placeLats, placeLons, placeID):
            pID = int(pID)

            log.debug("Plotting return period curve for %s"%name)
            i = find_index(lon, plon)
            j = find_index(lat, plat)

            xlabel = 'Average recurrence interval (years)'
            ylabel = 'Wind speed (%s)'%self.plotUnits.label
            title = "Return period wind speeds at " + name + ", \n(%5.1f,%5.1f)"%(plon, plat)

            name.replace(' ', '')
            log.debug("Working on {0}".format(name))
            filename = pjoin(plotPath, 'ARI_curve_%s.%s'%(pID, "png"))
            log.debug("Saving hazard curve for %s to %s"%(name, filename))
            wspd = ncobj.variables['wspd'][:, j, i]

            recs = database.queries.locationRecords(self.db, pID)
            data = np.zeros(int(self.numsimulations * 365.25))
            if len(recs) > 0:
                data[-len(recs):] = recs['wspd']

            allevents = np.sort(data)
            log.debug("allevents length = {0}".format(len(allevents)))

            placeWspd = metutils.convert(wspd, 'mps',
                                         self.plotUnits.units)
            if np.all(placeWspd.mask):
                log.debug("All values for {0} are null".format(name))
                continue

            if self.ciBounds:
                wspdLower = ncobj.variables['wspdlower'][:, j, i]
                wspdUpper = ncobj.variables['wspdupper'][:, j, i]
                placeWspdLower = metutils.convert(wspdLower, 'mps',
                                                  self.plotUnits.units)
                placeWspdUpper = metutils.convert(wspdUpper, 'mps',
                                                  self.plotUnits.units)
            else:
                placeWspdUpper = np.zeros(len(placeWspd))
                placeWspdLower = np.zeros(len(placeWspd))

            saveHazardCurve(years, allevents, placeWspd, placeWspdUpper, placeWspdLower,
                            xlabel, ylabel, title, filename, self.fit)

        ncobj.close()

