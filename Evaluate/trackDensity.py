"""
:mod:`TrackDensity` -- calculate track density
==============================================

.. module:: TrackDensity
   :synopsis: Calculate density of TC tracks over a grid (TCs/degree/year)

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import logging

import numpy as np
import numpy.ma as ma

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile

from Utilities.loadData import loadTrackFile
from Utilities.config import ConfigParser
from Utilities.track import ncReadTrackData
from Utilities.nctools import ncSaveGrid
from Utilities.parallel import attemptParallel, disableOnWorkers
from Utilities import pathLocator

from PlotInterface.maps import ArrayMapFigure, saveFigure

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

class TrackDensity(object):
    def __init__(self, configFile):
        """
        Calculate density of TC positions on a grid

        :param str configFile: path to a TCRM configuration file.
        """

        config = ConfigParser()
        config.read(configFile)
        self.configFile = configFile

        # Define the grid:
        self.gridLimit = config.geteval('Region', 'gridLimit')
        gridSpace = config.geteval('Region', 'GridSpace')

        self.lon_range = np.arange(self.gridLimit['xMin'],
                                   self.gridLimit['xMax'] + 0.1,
                                   gridSpace['x'])
        self.lat_range = np.arange(self.gridLimit['yMin'],
                                   self.gridLimit['yMax'] + 0.1,
                                   gridSpace['y'])
        self.X, self.Y = np.meshgrid(self.lon_range, self.lat_range)

        self.map_kwargs = dict(llcrnrlon=self.lon_range.min(),
                               llcrnrlat=self.lat_range.min(),
                               urcrnrlon=self.lon_range.max(),
                               urcrnrlat=self.lat_range.max(),
                               projection='merc',
                               resolution='i')
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
        Calculate a histogram of TC occurrences given a set
        of tracks

        :param tracks: Collection of :class:`Track` objects.
        """

        lon = []
        lat = []

        for t in tracks:
            #if t.inRegion(self.gridLimit):
            lon = np.append(lon, t.Longitude)
            lat = np.append(lat, t.Latitude)

        histogram, x, y = np.histogram2d(lon, lat,
                                         [self.lon_range,
                                          self.lat_range],
                                          normed=False)
        return histogram

    def calculateMeans(self):
        """
        Calculate the mean, median and percentiles of the synthetic values
        """
        self.synHist = ma.masked_values(self.synHist, -9999.)
        self.synHistMean = ma.mean(self.synHist, axis=0)
        self.medSynHist = ma.median(self.synHist, axis=0)

        self.synHistUpper = percentile(self.synHist, per=95, axis=0)
        self.synHistLower = percentile(self.synHist, per=5, axis=0)

    @disableOnWorkers
    def historic(self):
        """
        Load historic data and calculate histogram.
        Note that the input historical data is filtered by year
        when it's loaded in `interpolateTracks.parseTracks()`.

        The timestep to interpolate to is set to match that of the
        synthetic event set (normally set to 1 hour).
        """
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

            self.hist = self.calculate(tracks) / numYears

    def synthetic(self):
        """Load synthetic data and calculate histogram"""

        filelist = os.listdir(self.trackPath)
        trackfiles = [pjoin(self.trackPath, f) for f in filelist
                      if f.startswith('tracks')]
        self.synHist = -9999. * np.ones((len(trackfiles),
                                         len(self.lon_range)-1,
                                         len(self.lat_range)-1))

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
            while (terminated < comm.size - 1):
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
            while(True):
                trackfile = comm.Recv(source=0, tag=work_tag)
                if trackfile is None:
                    break

                log.debug("Processing %s" % (trackfile))
                tracks = loadTracks(trackfile)
                results = self.calculate(tracks) / self.synNumYears
                comm.Send(results, dest=0, tag=result_tag)

        elif (comm.size == 1) and (comm.rank == 0):
            for n, trackfile in enumerate(trackfiles):
                tracks = loadTracks(trackfile)
                self.synHist[n, :, :] = self.calculate(tracks) / \
                                        self.synNumYears

            self.calculateMeans()

    @disableOnWorkers
    def save(self):
        """Save data to file."""
        dataFile = pjoin(self.dataPath, 'density.nc')

        # Simple sanity check (should also include the synthetic data):
        if not hasattr(self, 'hist'):
            log.critical("No historical data available!")
            log.critical(("Check that data has been processed before "
                          "trying to save data"))
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
                'name': 'hist_density',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.hist),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Historical track density',
                    'units':'observations per 1-degree grid per year'
                }
            },
            1: {
                'name': 'syn_density',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistMean),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Track density - synthetic events',
                    'units':'observations per 1-degree grid per year'
                }
            },
            2: {
                'name': 'syn_density_upper',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistUpper),
                'dtype': 'f',
                'atts': {
                    'long_name': ('Track density - upper percentile '
                                  '- synthetic events'),
                    'units':' observations per 1-degree grid per year',
                    'percentile': '95'
                }
            },
            3: {
                'name': 'syn_density_lower',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistLower),
                'dtype': 'f',
                'atts': {
                    'long_name': ('Track density - lower percentile '
                                  '- synthetic events'),
                    'units': 'observations per 1-degree grid per year',
                    'percentile': '5'
                }
            }
        }

        ncSaveGrid(dataFile, dimensions, variables)


    @disableOnWorkers
    def plotTrackDensity(self):
        """Plot track density information"""

        datarange = (0, self.hist.max())
        figure = ArrayMapFigure()

        cbarlab = "TC observations/yr"
        figure.add(self.hist.T, self.X, self.Y, "Historic", datarange,
                   cbarlab, self.map_kwargs)
        figure.add(self.synHistMean.T, self.X, self.Y, "Synthetic",
                    datarange, cbarlab, self.map_kwargs)
        figure.plot()
        outputFile = pjoin(self.plotPath, 'track_density.png')
        saveFigure(figure, outputFile)

        figure2 = ArrayMapFigure()
        figure2.add(self.hist.T, self.X, self.Y, "Historic", datarange,
                    cbarlab, self.map_kwargs)
        figure2.add(self.synHist[0, :, :].T, self.X, self.Y, "Synthetic",
                    datarange, cbarlab, self.map_kwargs)
        figure2.add(self.synHist[10, :, :].T, self.X, self.Y, "Synthetic",
                    datarange, cbarlab, self.map_kwargs)
        figure2.add(self.synHist[20, :, :].T, self.X, self.Y, "Synthetic",
                    datarange, cbarlab, self.map_kwargs)
        figure2.add(self.synHist[30, :, :].T, self.X, self.Y, "Synthetic",
                    datarange, cbarlab, self.map_kwargs)
        figure2.add(self.synHist[40, :, :].T, self.X, self.Y, "Synthetic",
                    datarange, cbarlab, self.map_kwargs)
        figure2.plot()
        outputFile = pjoin(self.plotPath, 'track_density_samples.png')
        saveFigure(figure2, outputFile)


    @disableOnWorkers
    def plotTrackDensityPercentiles(self):
        """
        Plot upper and lower percentiles of track density derived from
        synthetic event sets

        """

        datarange = (0, self.hist.max())
        figure = ArrayMapFigure()

        cbarlab = "TC observations/yr"

        figure.add(self.synHistUpper.T, self.X, self.Y, "Upper percentile",
                   datarange, cbarlab, self.map_kwargs)
        figure.add(self.synHistLower.T, self.X, self.Y, "Lower percentile",
                    datarange, cbarlab, self.map_kwargs)
        figure.plot()
        outputFile = pjoin(self.plotPath, 'track_density_percentiles.png')
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

        self.plotTrackDensity()
        self.plotTrackDensityPercentiles()

        self.save()
