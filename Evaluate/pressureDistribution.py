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
import logging

import pdb

import numpy as np
import numpy.ma as ma

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile
from datetime import datetime

import matplotlib
matplotlib.use('Agg', warn=False)

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from Utilities.config import ConfigParser
from Utilities.metutils import convert
from Utilities.maputils import bearing2theta
from Utilities.loadData import loadTrackFile
from Utilities.track import Track
from Utilities import pathLocator
from Utilities.nctools import ncSaveGrid
from Utilities.parallel import attemptParallel, disableOnWorkers

# Importing :mod:`colours` makes a number of additional colour maps available:
from Utilities import colours

from PlotInterface.maps import ArrayMapFigure, saveFigure
from PlotInterface.curves import saveDistributionCurve
from PlotInterface.figures import QuantileFigure


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

TRACKFILE_COLS = ('CycloneNumber', 'Datetime', 'TimeElapsed', 'Longitude',
                  'Latitude', 'Speed', 'Bearing', 'CentralPressure',
                  'EnvPressure', 'rMax')

TRACKFILE_UNIT = ('', '%Y-%m-%d %H:%M:%S', 'hr', 'degree', 'degree', 'kph', 'degrees',
                  'hPa', 'hPa', 'km')

TRACKFILE_FMTS = ('i', datetime, 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f')

TRACKFILE_CNVT = {
    0: lambda s: int(float(s.strip() or 0)),
    1: lambda s: datetime.strptime(s.strip(), TRACKFILE_UNIT[1]),
    5: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[5], 'mps'),
    6: lambda s: bearing2theta(float(s.strip() or 0) * np.pi / 180.),
    7: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[7], 'hPa'),
    8: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[8], 'hPa'),
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

class PressureDistribution(object):
    def __init__(self, configFile):
        """
        Calculate central pressure distributions on a grid

        :param str configFile: path to a TCRM configuration file.
        """

        config = ConfigParser()
        config.read(configFile)
        self.configFile = configFile

        # Determine TCRM input directory
        tcrm_dir = pathLocator.getRootDirectory()
        self.inputPath = pjoin(tcrm_dir, 'input')

        # Define the grid:
        gridLimit = config.geteval('Region', 'gridLimit')
        gridSpace = config.geteval('Region', 'GridSpace')

        self.lon_range = np.arange(gridLimit['xMin'],
                                   gridLimit['xMax'] + 0.1,
                                   gridSpace['x'])
        self.lat_range = np.arange(gridLimit['yMin'],
                                   gridLimit['yMax'] + 0.1,
                                   gridSpace['y'])

        self.gridLimit = gridLimit
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
        log.debug("Processing %d tracks" % (len(tracks)))
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
            if t.inRegion(self.gridLimit):
                minCP[i] = t.CentralPressure.min()

        bins = np.arange(850., 1020., 5.)
        h, n = np.histogram(minCP, bins, normed=True)
        return h, minCP


    def calculateMeans(self, synMean, synMin, synMed, synMax, synMinCP):
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
        
        r = list(np.random.uniform(high=synMean.shape[0], size=3).astype(int))
        self.synRandomMinima = synMean[r, :, :]

    @disableOnWorkers
    def historic(self):
        """Load historic data and calculate histogram"""
        log.info("Processing historical pressure distributions")
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        source = config.get('DataProcess', 'Source')
        
        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.inputPath, inputFile)
        
        try:
            tracks = loadTrackFile(self.configFile, inputFile, source)
        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historical track file: {0}".format(inputFile))
            raise
        else:
            self.histMean, self.histMin, \
                self.histMax, self.histMed = self.calculate(tracks)

            self.histMinCPDist, self.histMinCP = self.calcMinPressure(tracks)

    def synthetic(self):
        """Load synthetic data and calculate histogram"""
        log.info("Processing synthetic pressure distributions")

        work_tag = 0
        result_tag = 1
        filelist = os.listdir(self.trackPath)
        trackfiles = sorted([pjoin(self.trackPath, f) for f in filelist
                             if f.startswith('tracks')])
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
        synMinCPDist = np.empty((len(trackfiles), len(bins) - 1))
        self.synMinCP = np.array([])
        if (pp.rank() == 0) and (pp.size() > 1):

            w = 0
            n = 0
            for d in range(1, pp.size()):
                pp.send(trackfiles[w], destination=d, tag=work_tag)
                log.debug("Processing track file %d of %d" % (w + 1, len(trackfiles)))
                w += 1

            terminated = 0
            while (terminated < pp.size() - 1):
                results, status = pp.receive(pp.any_source, tag=result_tag,
                                             return_status=True)
                

                sMean, sMin, sMax, sMed, sMinCPDist, sMinCP = results
                synMean[n, :, :] = sMean
                synMin[n, :, :] = sMin
                synMax[n, :, :] = sMax
                synMed[n, :, :] = sMed
                synMinCPDist[n, :] = sMinCPDist
                self.synMinCP = np.append(self.synMinCP, sMinCP)
                n += 1

                d = status.source

                if w < len(trackfiles):
                    pp.send(trackfiles[w], destination=d, tag=work_tag)
                    log.debug("Processing track file %d of %d" % (w + 1, len(trackfiles)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    terminated += 1

            self.calculateMeans(synMean, synMin, synMed, synMax, synMinCPDist)

        elif (pp.size() > 1) and (pp.rank() != 0):
            while(True):
                trackfile = pp.receive(source=0, tag=work_tag)
                if trackfile is None:
                    break
                
                log.debug("Processing %s" % (trackfile))
                tracks = loadTracks(trackfile)
                sMean, sMin, sMax, sMed = self.calculate(tracks)
                sMinCPDist, sMinCP = self.calcMinPressure(tracks)
                results = (sMean, sMin, sMax, sMed, sMinCPDist, sMinCP)
                pp.send(results, destination=0, tag=result_tag)
                
        elif pp.size() == 1 and pp.rank() == 0:
            # Assumed no Pypar - helps avoid the need to extend DummyPypar()
            for n, trackfile in enumerate(sorted(trackfiles)):
                tracks = loadTracks(trackfile)
                synMean[n, :, :], synMin[n, :, :], \
                    synMax[n, :, :], synMed[n, :, :] = self.calculate(tracks)
                synMinCPDist[n, :], sMinCP = self.calcMinPressure(tracks)
                self.synMinCP = np.append(self.synMinCP, sMinCP)
            self.calculateMeans(synMean, synMin, synMed, synMax, synMinCPDist)

    @disableOnWorkers
    def plotPressureMean(self):
        """
        Plot a map of observed and synthetic mean pressure values

        """
        
        
        datarange = (950, 1000)
        figure = ArrayMapFigure() 

        map_kwargs = dict(llcrnrlon=self.lon_range[:-1].min(),
                          llcrnrlat=self.lat_range[:-1].min(),
                          urcrnrlon=self.lon_range[:-1].max(),
                          urcrnrlat=self.lat_range[:-1].max(),
                          projection='merc',
                          resolution='i')

        cbarlab = "Mean central pressure (hPa)"
        xgrid, ygrid = np.meshgrid(self.lon_range[:-1], self.lat_range[:-1])
        figure.add(np.transpose(self.histMean), xgrid, ygrid, "Historic", datarange, 
                   cbarlab, map_kwargs)
        figure.add(np.transpose(self.synMean), xgrid, ygrid, "Synthetic", datarange, 
                   cbarlab, map_kwargs)

        figure.plot()
        outputFile = pjoin(self.plotPath, 'meanPressure.png')
        saveFigure(figure, outputFile)

    @disableOnWorkers
    def plotPressureMin(self):
        """
        Plot a map of observed and synthetic minimum central pressure values.

        """

        datarange = (900, 1000)
        figure = ArrayMapFigure()

        map_kwargs = dict(llcrnrlon=self.lon_range[:-1].min(),
                          llcrnrlat=self.lat_range[:-1].min(),
                          urcrnrlon=self.lon_range[:-1].max(),
                          urcrnrlat=self.lat_range[:-1].max(),
                          projection='merc',
                          resolution='i')

        cbarlab = "Minimum central pressure (hPa)"
        xgrid, ygrid = np.meshgrid(self.lon_range[:-1], self.lat_range[:-1])
        figure.add(np.transpose(self.histMin), xgrid, ygrid, "Historic", datarange, 
                   cbarlab, map_kwargs)
        figure.add(np.transpose(self.synMin), xgrid, ygrid, "Synthetic", datarange, 
                   cbarlab, map_kwargs)

        figure.plot()
        outputFile = pjoin(self.plotPath, 'minPressure.png')
        saveFigure(figure, outputFile)

    @disableOnWorkers
    def plotPressureMinDiff(self):
        """
        Plot a map of the difference between observed and synthetic minimum
        pressure values.

        """

        datarange = (-50, 50)
        figure = ArrayMapFigure()

        map_kwargs = dict(llcrnrlon=self.lon_range[:-1].min(),
                          llcrnrlat=self.lat_range[:-1].min(),
                          urcrnrlon=self.lon_range[:-1].max(),
                          urcrnrlat=self.lat_range[:-1].max(),
                          projection='merc',
                          resolution='i')

        cbarlab = "Minimum central pressure difference (hPa)"
        xgrid, ygrid = np.meshgrid(self.lon_range[:-1], self.lat_range[:-1])
        data = self.histMin - self.synMin
        figure.add(np.transpose(data), xgrid, ygrid, "", datarange, 
                   cbarlab, map_kwargs)
        figure.plot()
        outputFile = pjoin(self.plotPath, 'minPressureDiff.png')
        saveFigure(figure, outputFile)


    @disableOnWorkers
    def plotPressureMeanDiff(self):
        """
        Plot a map of the difference between observed and synthetic mean
        pressure values.

        """

        datarange = (-25, 25)
        figure = ArrayMapFigure()

        map_kwargs = dict(llcrnrlon=self.lon_range[:-1].min(),
                          llcrnrlat=self.lat_range[:-1].min(),
                          urcrnrlon=self.lon_range[:-1].max(),
                          urcrnrlat=self.lat_range[:-1].max(),
                          projection='merc',
                          resolution='i')

        cbarlab = "Mean central pressure difference (hPa)"
        data = self.histMean - self.synMean
        xgrid, ygrid = np.meshgrid(self.lon_range[:-1], self.lat_range[:-1])
        data = self.histMin - self.synMin
        figure.add(np.transpose(data), xgrid, ygrid, "", datarange, 
                   cbarlab, map_kwargs)
        figure.plot()
        outputFile = pjoin(self.plotPath, 'meanPressureDiff.png')
        saveFigure(figure, outputFile)

    @disableOnWorkers
    def plotMinPressureDistribution(self):
        """
        Plot a pdf of observed minimum central pressure values, and the mean
        of the synthetic event sets (plus 90th percentile values).

        """
        x = np.arange(850., 1020., 5.)[:-1]
        y1 = self.histMinCPDist
        y2 = self.synMinCPDist
        y2min= self.synMinCPLower
        y2max = self.synMinCPUpper
        outputFile = pjoin(self.plotPath, 'minPressureDist.png')
        saveDistributionCurve(x, y1, y2, y2max, y2min, "Minimum pressure (hPa)", 
                              "Probability", "Minimum pressure distribution", 
                              outputFile)

    @disableOnWorkers
    def plotMinPressureQuantiles(self):
        x = self.histMinCP
        y = self.synMinCP
        lims = (850, 1000)
        fig = QuantileFigure()
        fig.add(x.compress(x>0), y.compress(y>0), lims, "Observed pressure (hPa)", "Simulated pressure (hPa)",
                "Q-Q plot of minimum central pressure")
        fig.plot()
        outputFile = pjoin(self.plotPath, 'minPressureQuantiles.png')
        saveFigure(fig, outputFile)

    @disableOnWorkers
    def save(self):
        dataFile = pjoin(self.dataPath, 'pressureDistribution.nc')

        # Simple sanity check (should also include the synthetic data):
        if not hasattr(self, 'histMin'):
            log.critical("No historical data available!")
            log.critical("Check that data has been processed before trying to save data")
            return

        log.info('Saving pressure distribution data to {0}'.format(dataFile))
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
        histCPFile = pjoin(self.dataPath, 'histCP.csv')
        synCPFile = pjoin(self.dataPath, 'synCP.csv')
        np.savetxt(histCPFile, self.histMinCP)
        np.savetxt(synCPFile, self.synMinCP)

    def run(self):
        """Run the pressure distribution evaluation"""
        global pp
        pp = attemptParallel()

        self.historic()

        pp.barrier()

        self.synthetic()

        pp.barrier()

        self.plotPressureMean()
        self.plotPressureMin()

        self.plotPressureMeanDiff()
        self.plotPressureMinDiff()

        self.plotMinPressureDistribution()
        self.plotMinPressureQuantiles()
        self.save()
