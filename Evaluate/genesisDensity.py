"""
:mod:`GenesisDensity` -- calculate TC genesis density
=====================================================

.. module:: GenesisDensity
   :synopsis: Calculate density of TC genesis points over a grid (TCs/year)

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

TODO: Use `statsmodels.nonparametric.kde.KDEMultivariate` for calculating
      genesis density for synthetic event sets.

"""

import os
import logging

import numpy as np
import numpy.ma as ma

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile
from statsmodels.nonparametric.kernel_density import KDEMultivariate

from Utilities.config import ConfigParser
from Utilities.track import ncReadTrackData
from Utilities.nctools import ncSaveGrid
from Utilities.loadData import loadTrackFile
from Utilities.parallel import attemptParallel, disableOnWorkers
from Utilities import pathLocator

from PlotInterface.maps import FilledContourMapFigure, saveFigure, levels

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

def loadTracks(trackfile):
    """
    Read tracks from a track .nc file and return a list of :class:`Track`
    objects.

    This calls the function `ncReadTrackData` to parse the track .nc
    file.

    :type  trackfile: str
    :param trackfile: the track data filename.
    """
    tracks = ncReadTrackData(trackfile)
    return tracks

def loadTracksFromFiles(trackfiles):
    """
    Provides an iterator over all tracks, which have been loaded from
    a list of track files.

    :param list trackfiles: A list of track files to load.

    :returns: an iterator that yields :class:`Track` objects.
    """
    for f in trackfiles:
        tracks = loadTracks(f)
        for track in tracks:
            yield track

def loadTracksFromPath(path):
    """
    Load a collection of track files from a given path.

    :param str path: Path to load track files from.

    :returns: iterator of :class:`Track` objects.
    """
    files = os.listdir(path)
    trackfiles = [pjoin(path, f) for f in files if f.startswith('tracks')]
    msg = 'Processing %d track files in %s' % (len(trackfiles), path)
    log.info(msg)
    return loadTracksFromFiles(sorted(trackfiles))

class GenesisDensity(object):
    def __init__(self, configFile):
        """
        Calculate density of TC genesis positions on a grid

        :param str configFile: path to a TCRM configuration file.
        """

        config = ConfigParser()
        config.read(configFile)
        self.configFile = configFile

        # Define the grid:
        self.gridLimit = config.geteval('Region', 'gridLimit')
        self.gridSpace = config.geteval('Region', 'GridSpace')

        self.lon_range = np.arange(self.gridLimit['xMin'],
                                   self.gridLimit['xMax'] + 0.1,
                                   1.)
        self.lat_range = np.arange(self.gridLimit['yMin'],
                                   self.gridLimit['yMax'] + 0.1,
                                   1.)

        self.X, self.Y = np.meshgrid(self.lon_range, self.lat_range)

        outputPath = config.get('Output', 'Path')
        self.trackPath = pjoin(outputPath, 'tracks')
        self.plotPath = pjoin(outputPath, 'plots', 'stats')
        self.dataPath = pjoin(outputPath, 'process')

        # Determine TCRM input directory
        tcrm_dir = pathLocator.getRootDirectory()
        self.inputPath = pjoin(tcrm_dir, 'input')

        self.synNumYears = config.getint('TrackGenerator',
                                         'yearspersimulation')


    def calculate(self, tracks):
        """
        Calculate a histogram of TC genesis occurrences given a set
        of tracks

        :param tracks: Collection of :class:`Track` objects.
        """

        lon = np.array([])
        lat = np.array([])

        for t in tracks:
            lon = np.append(lon, t.Longitude[0])
            lat = np.append(lat, t.Latitude[0])
        histogram, x, y = np.histogram2d(lon, lat,
                                         [self.lon_range,
                                          self.lat_range],
                                         normed=False)
        return histogram

    def calculatePDF(self, tracks):
        """
        Calculate a 2-d probability density surface using kernel density
        estimation.

        :param tracks: Collection of :class:`Track` objects.
        """

        if len(tracks) == 0:
            # No tracks:
            return np.zeros(self.X.shape)

        lon = np.array([])
        lat = np.array([])
        for t in tracks:
            lon = np.append(lon, t.Longitude)
            lat = np.append(lat, t.Latitude)

        xy = np.vstack([self.X.ravel(), self.Y.ravel()])
        data = np.array([[lon], [lat]])

        kde = KDEMultivariate(data, bw='cv_ml', var_type='cc')
        pdf = kde.pdf(data_predict=xy)

        return pdf.reshape(self.X.shape)

    def _calculate(self, tracks):
        """
        Calculate a histogram of TC genesis counts given a set of tracks.

        :param tracks: Collection of :class:`Track` objects.
        """
        log.debug("Calculating PDF for set of {0:d} tracks".format(len(tracks)))

        hist = ma.zeros((len(self.lon_range) - 1,
                         len(self.lat_range) - 1))

        xy = np.vstack([self.X.ravel(), self.Y.ravel()])

        x = []
        y = []

        for track in tracks:
            if len(track.Longitude) == 0:
                pass
            elif len(track.Longitude) == 1:
                x.append(track.Longitude)
                y.append(track.Latitude)
            else:
                x.append(track.Longitude[0])
                y.append(track.Latitude[0])

        xx = np.array(x)
        yy = np.array(y)
        ii = np.where((xx >= self.gridLimit['xMin']) &
                      (xx <= self.gridLimit['xMax']) &
                      (yy >= self.gridLimit['yMin']) &
                      (yy <= self.gridLimit['yMax']))

        values = np.vstack([xx[ii], yy[ii]])
        kernel = KDEMultivariate(values, bw='cv_ml', var_type='cc')
        pdf = kernel.pdf(data_predict=xy)
        Z = np.reshape(pdf, self.X.shape)
        return Z.T

    def calculateMeans(self):
        """
        Calculate mean, median and percentiles of the :attr:`self.synHist`
        attribute.

        """

        self.synHist = ma.masked_values(self.synHist, -9999.)
        self.synHistMean = ma.mean(self.synHist, axis=0)
        self.medSynHist = ma.median(self.synHist, axis=0)

        self.synHistUpper = percentile(self.synHist, per=95, axis=0)
        self.synHistLower = percentile(self.synHist, per=5, axis=0)

    @disableOnWorkers
    def historic(self):
        """Load historic data and calculate histogram"""
        log.info("Processing historic track records")
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.inputPath, inputFile)

        source = config.get('DataProcess', 'Source')

        try:
            tracks = loadTrackFile(self.configFile, inputFile, source)

        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historical track file: {0}".\
                         format(inputFile))
            raise
        else:
            startYr = 9999
            endYr = 0
            for t in tracks:
                startYr = min(startYr, min(t.Year))
                endYr = max(endYr, max(t.Year))
            numYears = endYr - startYr
            log.info("Range of years: %d - %d" % (startYr, endYr))
            try:
                self.hist = self._calculate(tracks)
            #self.hist = self._calculate(tracks) / numYears
            except (ValueError):
                log.critical("KDE error: The number of observations must be larger than the number of variables")
                raise

    def synthetic(self):
        """Load synthetic data and calculate histogram"""

        filelist = os.listdir(self.trackPath)
        trackfiles = [pjoin(self.trackPath, f) for f in filelist
                      if f.startswith('tracks')]
        self.synHist = -9999. * np.ones((len(trackfiles),
                                         len(self.lon_range),
                                         len(self.lat_range)))

        work_tag = 0
        result_tag = 1

        if (comm.rank == 0) and (comm.size > 1):
            w = 0
            n = 0
            for d in range(1, comm.size):
                comm.Send(trackfiles[w], dest=d, tag=work_tag)
                log.debug("Processing track file {0:d} of {1:d}".\
                          format(w, len(trackfiles)))
                w += 1

            terminated = 0
            while terminated < comm.size - 1:
                results, status = comm.Recv(MPI.ANY_SOURCE, tag=result_tag,
                                             status=True)
                self.synHist[n, :, :] = results
                n += 1

                d = status.source
                if w < len(trackfiles):
                    comm.Send(trackfiles[w], dest=d, tag=work_tag)
                    log.debug("Processing track file {0:d} of {1:d}".\
                              format(w, len(trackfiles)))
                    w += 1
                else:
                    comm.Send(None, dest=d, tag=work_tag)
                    terminated += 1

            self.calculateMeans()

        elif (comm.size > 1) and (comm.rank != 0):
            while True:
                trackfile = comm.Recv(source=0, tag=work_tag)
                if trackfile is None:
                    break

                log.debug("Processing {0}".format(trackfile))
                tracks = loadTracks(trackfile)
                results = self._calculate(tracks) #/ self.synNumYears
                comm.Send(results, dest=0, tag=result_tag)

        elif (comm.size == 1) and (comm.rank == 0):
            for n, trackfile in enumerate(trackfiles):
                log.debug("Processing track file {0:d} of {1:d}".\
                          format(n + 1, len(trackfiles)))
                tracks = loadTracks(trackfile)
                self.synHist[n, :, :] = self._calculate(tracks) #/ self.synNumYears

            self.calculateMeans()

    @disableOnWorkers
    def save(self):
        dataFile = pjoin(self.dataPath, 'genesis_density.nc')

        # Simple sanity check (should also include the synthetic data):
        if not hasattr(self, 'hist'):
            log.critical("No historical data available!")
            log.critical(("Check that data has been processed "
                          "before trying to save data"))
            return

        log.info('Saving genesis density data to {0}'.format(dataFile))
        dimensions = {
            0: {
                'name': 'lat',
                'values': self.lat_range,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Latitude',
                    'units': 'degrees_north',
                    'axis': 'Y'
                }
            },
            1: {
                'name': 'lon',
                'values': self.lon_range,
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
                'name': 'hist_density',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.hist),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Historical genesis density',
                    'units':'TCs per 1-degree grid per year'
                }
            },
            1: {
                'name': 'syn_density',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistMean),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Genesis density - synthetic events',
                    'units':'TCs per 1-degree grid per year'
                }
            },
            2: {
                'name': 'syn_density_upper',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistUpper),
                'dtype': 'f',
                'atts': {
                    'long_name': ('Genesis density - upper percentile '
                                  '- synthetic events'),
                    'units':'TCs per 1-degree grid per year',
                    'percentile': '95'
                }
            },
            3: {
                'name': 'syn_density_lower',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistLower),
                'dtype': 'f',
                'atts': {
                    'long_name': ('Genesis density - lower percentile '
                                  '- synthetic events'),
                    'units': 'TCs per 1-degree grid per year',
                    'percentile': '5'
                }
            }
        }

        ncSaveGrid(dataFile, dimensions, variables)


    @disableOnWorkers
    def plotGenesisDensity(self):
        """Plot genesis density information"""

        datarange = (0, self.hist.max())
        figure = FilledContourMapFigure()
        lvls, exponent = levels(self.hist.max())
        map_kwargs = dict(llcrnrlon=self.lon_range.min(),
                          llcrnrlat=self.lat_range.min(),
                          urcrnrlon=self.lon_range.max(),
                          urcrnrlat=self.lat_range.max(),
                          projection='merc',
                          resolution='i')
        cbarlab = "TCs/yr"
        xgrid, ygrid = np.meshgrid(self.lon_range, self.lat_range)
        figure.add(self.hist.T*(10.**-exponent), self.X, self.Y,
                   "Historic", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure.add(self.synHistMean.T*(10.**-exponent), xgrid, ygrid,
                   "Synthetic", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure.plot()
        outputFile = pjoin(self.plotPath, 'genesis_density.png')
        saveFigure(figure, outputFile)

        figure2 = FilledContourMapFigure()
        figure2.add(self.hist.T*(10.**-exponent), xgrid, ygrid,
                    "Historic", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure2.add(self.synHist[0, :, :].T*(10.**-exponent), xgrid, ygrid,
                    "Synthetic", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure2.add(self.synHist[1, :, :].T*(10.**-exponent), xgrid, ygrid,
                    "Synthetic", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure2.add(self.synHist[2, :, :].T*(10.**-exponent), xgrid, ygrid,
                    "Synthetic", lvls*(10.**-exponent), cbarlab, map_kwargs)

        figure2.plot()
        outputFile = pjoin(self.plotPath, 'genesis_density_samples.png')
        saveFigure(figure2, outputFile)

    @disableOnWorkers
    def plotGenesisDensityPercentiles(self):
        """
        Plot upper and lower percentiles of genesis density derived from
        synthetic event sets

        """

        datarange = (0, self.hist.max())
        figure = FilledContourMapFigure()
        lvls, exponent = levels(self.hist.max())

        map_kwargs = dict(llcrnrlon=self.lon_range.min(),
                          llcrnrlat=self.lat_range.min(),
                          urcrnrlon=self.lon_range.max(),
                          urcrnrlat=self.lat_range.max(),
                          projection='merc',
                          resolution='i')
        cbarlab = "TCs/yr"
        xgrid, ygrid = np.meshgrid(self.lon_range, self.lat_range)

        figure.add(self.synHistUpper.T*(10.**-exponent), xgrid, ygrid,
                   "Upper percentile", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure.add(self.synHistLower.T*(10.**-exponent), xgrid, ygrid,
                   "Lower percentile", lvls*(10.**-exponent), cbarlab, map_kwargs)
        figure.plot()
        outputFile = pjoin(self.plotPath, 'genesis_density_percentiles.png')
        saveFigure(figure, outputFile)



    def run(self):
        """Run the track density evaluation"""
        global MPI, comm
        MPI = attemptParallel()
        comm = MPI.COMM_WORLD
        self.historic()

        comm.barrier()

        self.synthetic()

        comm.barrier()

        self.plotGenesisDensity()
        self.plotGenesisDensityPercentiles()

        self.save()
