"""
:mod:`PressureDistributions` -- calculate central pressure distributions
========================================================================

.. module:: PressureDistribution
    :synopsis: calculate central pressure distributions based on historic and
               synthetic event sets

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import sys
import pdb
import logging

import numpy as np
import numpy.ma as ma

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile

from matplotlib import pyplot
from mpl_toolkits.basemap import Basemap


from Utilities.config import ConfigParser
from Utilities.metutils import convert
from Utilities.maputils import bearing2theta
from Utilities.loadData import loadTrackFile
from Utilities.track import Track
from Utilities.nctools import ncSaveGrid

# Importing :mod:`colours` makes a number of additional colour maps available:
from Utilities import colours

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

TRACKFILE_COLS = ('CycloneNumber', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'EnvPressure', 'rMax')

TRACKFILE_UNIT = ('', 'hr', 'degree', 'degree', 'kph', 'degrees',
                  'hPa', 'hPa', 'km')

TRACKFILE_FMTS = ('i', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f')

TRACKFILE_CNVT = {
    0: lambda s: int(float(s.strip() or 0)),
    4: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[4], 'mps'),
    5: lambda s: bearing2theta(float(s.strip() or 0) * np.pi / 180.),
    6: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[6], 'hPa'),
    7: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[7], 'hPa'),
}

def readTrackData(trackfile):
    """
    Read a track .csv file into a numpy.ndarray.

    The track format and converters are specified with the global variables

        TRACKFILE_COLS -- The column names
        TRACKFILE_FMTS -- The entry formats
        TRACKFILE_CNVT -- The column converters

    :param str trackfile: the track data filename.
    """
    try:
        return np.loadtxt(trackfile,
                          comments='%',
                          delimiter=',',
                          dtype={
                          'names': TRACKFILE_COLS,
                          'formats': TRACKFILE_FMTS},
                          converters=TRACKFILE_CNVT)
    except ValueError:
        # return an empty array with the appropriate `dtype` field names
        return np.empty(0, dtype={
                        'names': TRACKFILE_COLS,
                        'formats': TRACKFILE_FMTS})

def readMultipleTrackData(trackfile):
    """
    Reads all the track datas from a .csv file into a list of numpy.ndarrays.
    The tracks are seperated based in their cyclone id. This function calls
    `readTrackData` to read the data from the file.

    :type  trackfile: str
    :param trackfile: the track data filename.
    """

    log.debug("Reading multiple track data from {0}".format(trackfile))

    datas = []
    data = readTrackData(trackfile)
    if len(data) > 0:
        cycloneId = data['CycloneNumber']
        for i in range(1, np.max(cycloneId) + 1):
            if len(data[cycloneId == i]) > 0:
                datas.append(data[cycloneId == i])
            else:
                pass
    else:
        datas.append(data)
    return datas

def loadTracks(trackfile):
    """
    Read tracks from a track .csv file and return a list of :class:`Track`
    objects.

    This calls the function `readMultipleTrackData` to parse the track .csv
    file.

    :type  trackfile: str
    :param trackfile: the track data filename.
    """
    tracks = []
    datas = readMultipleTrackData(trackfile)
    n = len(datas)
    for i, data in enumerate(datas):
        track = Track(data)
        track.trackfile = trackfile
        track.trackId = (i, n)
        tracks.append(track)
    return tracks

class gridCell(object):
    def __init__(self, xmin, ymin, xmax, ymax, number, index):
        self.xmin = xmin
        self.ymin = ymin
        self.xmax = xmax
        self.ymax = ymax
        self.cell_number = number
        self.index = index

def plotDensity(x, y, data, llLon=None, llLat=None, urLon=None, urLat=None,
                res='i', dl=5., datarange=(-1.,1.), cmap='gist_heat_r', title=None,
                xlab='Longitude', ylab='Latitude', clabel=None, maskland=False,
                maskocean=False):
    """
    Plot a grid of values using pyplot.pcolormesh() to create a pseudocolour
    plot of a 2-dimensional array, with a basemap included.

    It updates the figure object, but does not have scope to save or display it.
    That way, the function can be used to generate a number of subplots, with
    the calling program handling the saving/display requirements.

    :param x: :class:`numpy.ndarray` of longitude points
    :param y: :class:`numpy.ndarray` of latitude points
    :param data: 2-d :class:`numpy.ndarray` of data values.

    """

    if (len(x.shape) < 2) and (len(y.shape) < 2):
        [xx, yy] = np.meshgrid(x, y)
    else:
        xx = x
        yy = y
    if llLon:
        llcrnrlon = llLon
    else:
        llcrnrlon = x.min()

    if llLat:
        llcrnrlat = llLat
    else:
        llcrnrlat = y.min()

    if urLon:
        urcrnrlon = urLon
    else:
        urcrnrlon = x.max()

    if urLat:
        urcrnrlat = urLat
    else:
        urcrnrlat = y.max()

    meridians = np.arange(dl*np.floor(llcrnrlon / dl),
                          dl*np.ceil(urcrnrlon / dl), dl)
    parallels = np.arange(dl*np.floor(llcrnrlat / dl),
                          dl*np.ceil(urcrnrlat / dl), dl)

    m = Basemap(projection='cyl',
                resolution=res,
                llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat)

    # Set the colour map:
    cmap = pyplot.get_cmap(cmap)

    if maskocean:
        try:
            from mpl_toolkits.basemap import maskoceans
        except ImportError:
            log.debug("Maskoceans module unavailable, skipping this command")
        else:
            datam = maskoceans(xx, yy, data, inlands=False)
            m.pcolormesh(xx, yy, datam, edgecolors='None',
                         vmin=datarange[0], vmax=datarange[1],
                         cmap=cmap)
    else:
        m.pcolormesh(xx, yy, data, edgecolors='None',
                     vmin=datarange[0], vmax=datarange[1],
                     cmap=cmap)

    m.drawcoastlines(linewidth=0.5)
    if maskland:
        m.fillcontinents(color='white')

    m.drawparallels(parallels, labels=[1, 0, 0, 0],
                    fontsize=7.5, linewidth=0.2)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1],
                    fontsize=7.5, linewidth=0.2)
    if ylab:
        pyplot.ylabel(ylab, fontsize=7.5, labelpad=25)
    if xlab:
        pyplot.xlabel(xlab, fontsize=7.5, labelpad=20)
    if title:
        pyplot.title(title)

    pyplot.grid(True)
    pyplot.tick_params(direction='out', length=4, right='off', top='off')

    cb = pyplot.colorbar(aspect=30, orientation='vertical',
                         extend='both', pad=0.1)
    cb.ax.tick_params(direction='in')

    if cb.orientation == 'horizontal':
        for t in cb.ax.get_xticklabels():
            t.set_fontsize(8)

    if clabel:
        cb.set_label(clabel)

    return

class PressureDistribution(object):
    def __init__(self, configFile):
        """
        Calculate central pressure distributions on a grid

        :param str configFile: path to a TCRM configuration file.
        """

        config = ConfigParser()
        config.read(configFile)
        self.configFile = configFile

        # Define the grid:
        gridLimit = config.geteval('Region', 'gridLimit')
        gridSpace = config.geteval('Region', 'GridSpace')

        self.lon_range = np.arange(gridLimit['xMin'],
                                   gridLimit['xMax'] + 0.1,
                                   gridSpace['x'])
        self.lat_range = np.arange(gridLimit['yMin'],
                                   gridLimit['yMax'] + 0.1,
                                   gridSpace['y'])

        outputPath = config.get('Output', 'Path')
        self.trackPath = pjoin(outputPath, 'tracks')
        self.plotPath = pjoin(outputPath, 'plots', 'stats')
        self.dataPath = pjoin(outputPath, 'process')

        self.synNumYears = config.getint('TrackGenerator',
                                         'yearspersimulation')
        cellnumber = 0
        self.gridCells = []
        for k in xrange(len(self.lon_range) - 1):
            for l in xrange(len(self.lat_range) - 1):
                ymin = self.lat_range[l]
                ymax = self.lat_range[l] + gridSpace['y']
                xmin = self.lon_range[k]
                xmax = self.lon_range[k] + gridSpace['x']
                self.gridCells.append(gridCell(xmin, ymin, xmax, ymax,
                                               cellnumber, (k, l)))
                cellnumber += 1

    def calculate(self, tracks):

        dataMean = ma.zeros((len(self.lon_range) - 1,
                             len(self.lat_range) - 1))
        dataMin = ma.zeros((len(self.lon_range) - 1,
                            len(self.lat_range) - 1))
        dataMax = ma.zeros((len(self.lon_range) - 1,
                            len(self.lat_range) - 1))
        dataMed = ma.zeros((len(self.lon_range) - 1,
                            len(self.lat_range) - 1))

        for cell in self.gridCells:
            vcell = np.array([])
            for t in tracks:
                ii = np.where(((t.Latitude >= cell.ymin) &
                               (t.Latitude < cell.ymax)) &
                              ((t.Longitude >= cell.xmin) &
                               (t.Longitude < cell.xmax)))[0]
                if len(ii) > 0:
                    vv = t.CentralPressure[ii].compress(t.CentralPressure[ii] < sys.maxint)
                    vcell = np.append(vcell, vv.compress(vv > 0.0))

            if len(vcell > 0):
                dataMean[cell.index[0], cell.index[1]] = np.mean(vcell)
                dataMin[cell.index[0], cell.index[1]] = np.min(vcell)
                dataMax[cell.index[0], cell.index[1]] = np.max(vcell)
                dataMed[cell.index[0], cell.index[1]] = np.median(vcell)

        dataMean = ma.masked_equal(dataMean, 0)
        dataMin = ma.masked_equal(dataMin, 0)
        dataMax = ma.masked_equal(dataMax, 0)
        dataMed = ma.masked_equal(dataMed, 0)
        return dataMean, dataMin, dataMax, dataMed

    def calcMinPressure(self, tracks):
        minCP = np.zeros(len(tracks))

        for i, t in enumerate(tracks):
            minCP[i] = t.CentralPressure.min()

        bins = np.arange(850., 1020., 5.)
        h, n = np.histogram(minCP, bins, normed=True)
        return h

    def historic(self):
        """Load historic data and calculate histogram"""
        log.info("Processing historical pressure distributions")
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        source = config.get('DataProcess', 'Source')

        path, base = os.path.split(inputFile)
        try:
            tracks = loadTrackFile(self.configFile, inputFile, source)
        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historical track file: {0}".format(inputFile))
            raise
        else:
            self.histMean, self.histMin, \
                self.histMax, self.histMed = self.calculate(tracks)

            self.histMinCP = self.calcMinPressure(tracks)

    def synthetic(self):
        """Load synthetic data and calculate histogram"""
        log.info("Processing synthetic pressure distributions")

        filelist = os.listdir(self.trackPath)
        trackfiles = [pjoin(self.trackPath, f) for f in filelist
                      if f.startswith('tracks')]
        synMean = -9999. * ma.ones((len(trackfiles),
                                    len(self.lon_range) - 1,
                                   len(self.lat_range) - 1))
        synMin = -9999. * ma.ones((len(trackfiles),
                                   len(self.lon_range) - 1,
                                   len(self.lat_range) - 1))
        synMax = -9999. * ma.ones((len(trackfiles),
                                   len(self.lon_range) - 1,
                                   len(self.lat_range) - 1))
        synMed = -9999. * ma.ones((len(trackfiles),
                                   len(self.lon_range) - 1,
                                   len(self.lat_range) - 1))

        bins = np.arange(850., 1020., 5.)
        synMinCP = np.empty((len(trackfiles),
                             len(bins) - 1))

        for n, trackfile in enumerate(trackfiles):
            tracks = loadTracks(trackfile)
            synMean[n, :, :],  synMin[n, :, :], \
               synMax[n, :, :], synMed[n, :, :] = self.calculate(tracks)
            synMinCP[n, :] = self.calcMinPressure(tracks)

        synMean = ma.masked_values(synMean, -9999.)
        synMin = ma.masked_values(synMin, -9999.)
        synMed = ma.masked_values(synMed, -9999.)
        synMax = ma.masked_values(synMax, -9999.)

        self.synMean = ma.mean(synMean, axis=0)
        self.synMed = ma.mean(synMed, axis=0)
        self.synMin = ma.mean(synMin, axis=0)
        self.synMax = ma.mean(synMax, axis=0)

        self.synMeanUpper = percentile(ma.compressed(synMean), per=95, axis=0)
        self.synMeanLower = percentile(ma.compressed(synMean), per=5, axis=0)
        self.synMinUpper = percentile(ma.compressed(synMin), per=95, axis=0)
        self.synMinLower = percentile(ma.compressed(synMin), per=5, axis=0)

        self.synMinCPDist = np.mean(synMinCP, axis=0)
        self.synMinCPLower = percentile(synMinCP, per=5, axis=0)
        self.synMinCPUpper = percentile(synMinCP, per=95, axis=0)

    def plotPressureMean(self):
        """
        Plot a map of observed and synthetic mean pressure values

        """

        datarange = (950, 1000)
        pyplot.figure()
        ax1 = pyplot.subplot(211)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1],
                    np.transpose(self.histMean), datarange=datarange,
                    clabel="Mean central pressure (hPa)", cmap='gist_heat')
        ax1.text(self.lon_range[1], self.lat_range[1],"Historic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))

        ax2 = pyplot.subplot(212)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1],
                    np.transpose(self.synMean), datarange=datarange,
                    clabel="Mean central pressure (hPa)", cmap='gist_heat')
        ax2.text(self.lon_range[1], self.lat_range[1],"Mean synthetic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))
        pyplot.savefig(pjoin(self.plotPath, 'meanPressure.png'))

    def plotPressureMin(self):
        """
        Plot a map of observed and synthetic minimum central pressure values.

        """

        datarange = (900, 1000)
        pyplot.figure()
        ax1 = pyplot.subplot(211)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1],
                    np.transpose(self.histMin), datarange=datarange,
                    clabel="Minimum central pressure (hPa)", cmap='gist_heat')
        ax1.text(self.lon_range[1], self.lat_range[1],"Historic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))

        ax2 = pyplot.subplot(212)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1],
                    np.transpose(self.synMin), datarange=datarange,
                    clabel="Minimum central pressure (hPa)", cmap='gist_heat')
        ax2.text(self.lon_range[1], self.lat_range[1],"Mean synthetic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))
        pyplot.savefig(pjoin(self.plotPath, 'minPressure.png'))

    def plotPressureMinDiff(self):
        """
        Plot a map of the difference between observed and synthetic minimum
        pressure values.

        """

        datarange = (-50, 50)
        pyplot.figure()
        ax1 = pyplot.subplot(111)
        data = self.histMin - self.synMin
        plotDensity(self.lon_range[:-1], self.lat_range[:-1],
                    np.transpose(data), datarange=datarange,
                    clabel="Minimum central pressure difference (hPa)", cmap='rwbcmap')
        ax1.text(self.lon_range[1], self.lat_range[1],"Historic - mean synthetic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))

        pyplot.savefig(pjoin(self.plotPath, 'minPressureDiff.png'))

    def plotPressureMeanDiff(self):
        """
        Plot a map of the difference between observed and synthetic mean
        pressure values.

        """

        datarange = (-50, 50)
        pyplot.figure()
        ax1 = pyplot.subplot(111)
        data = self.histMean - self.synMean
        plotDensity(self.lon_range[:-1], self.lat_range[:-1],
                    np.transpose(data), datarange=datarange,
                    clabel="Mean central pressure difference (hPa)", cmap='rwbcmap')
        ax1.text(self.lon_range[1], self.lat_range[1],"Historic - mean synthetic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))

        pyplot.savefig(pjoin(self.plotPath, 'meanPressureDiff.png'))

    def plotMinPressureDistribution(self):
        """
        Plot a pdf of observed minimum central pressure values, and the mean
        of the synthetic event sets (plus 90th percentile values).

        """

        pyplot.figure()
        ax1 = pyplot.subplot(111)
        bins = np.arange(850., 1020., 5.)
        ax1.plot(bins[:-1], self.histMinCP, color='r', lw=2)
        ax1.plot(bins[:-1], self.synMinCPDist, color='k', lw=2)
        ax1.fill_between(bins[:-1], self.synMinCPUpper,
                         self.synMinCPLower, facecolor='0.75',
                         edgecolor='0.99', alpha=0.7)

        pyplot.xlim(850., 1020.)
        pyplot.xticks()
        pyplot.xlabel('Minimum pressure')
        pyplot.yticks()
        pyplot.ylabel('Probability')
        pyplot.grid(True)

        outputFile = pjoin(self.plotPath, 'minPressureDist.png')
        pyplot.savefig(outputFile)

    def save(self):
        dataFile = pjoin(self.dataPath, 'pressureDistribution.nc')

        # Simple sanity check (should also include the synthetic data):
        if not hasattr(self, 'hist'):
            log.critical("No historical data available!")
            log.critical("Check that data has been processed before trying to save data")
            return

        log.info('Saving track density data to {0}'.format(dataFile))
        dimensions = {
            0: {
                'name': 'lat',
                'values': self.lat_range[:-1],
                'dtype': 'f',
                'atts': {
                    'long_name': 'Latitude',
                    'units': 'degrees_north',
                    'axis': 'Y'
                }
            },
            1: {
                'name': 'lon',
                'values': self.lon_range[:-1],
                'dtype': 'f',
                'atts': {
                    'long_name': 'Longitude',
                    'units':'degrees_east',
                    'axis': 'X'
                }
            }
        }

        # Define variables:
        variables = {
            0: {
                'name': 'hist_mean',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.histMean),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Historical mean central pressure',
                    'units':'hPa'
                }
            },
            1: {
                'name': 'syn_mean',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synMean),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean synthetic mean central pressure',
                    'units':'hPa'
                }
            },
            2: {
                'name': 'hist_min',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.histMin),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Historical minimum central pressure',
                    'units':'hPa'
                }
            },
            3: {
                'name': 'syn_min',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synMin),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean synthetic minimum central pressure',
                    'units': 'hPa'
                }
            }
        }

        ncSaveGrid(dataFile, dimensions, variables)

    def run(self):
        """Run the pressure distribution evaluation"""
        self.historic()
        self.synthetic()

        self.plotPressureMean()
        self.plotPressureMin()

        self.plotPressureMeanDiff()
        self.plotPressureMinDiff()

        self.plotMinPressureDistribution()

        self.save()
