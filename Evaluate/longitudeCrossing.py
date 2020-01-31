"""
:mod:`LongitudeCrossing` -- calculate crossing rate of lons/lats
================================================================

.. module:: LongitudeCrossing
    :synopsis: Calculate the rate of TCs crossing lines of longitude
               and latitude, comparing historical and synthetic
               events.

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import logging

import numpy as np

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile

import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from . import interpolateTracks as interpolateTracks

from Utilities.config import ConfigParser
from Utilities.track import ncReadTrackData
from Utilities.nctools import ncSaveGrid
from Utilities.files import flProgramVersion
from Utilities.parallel import attemptParallel, disableOnWorkers
from Utilities import pathLocator
import Utilities.Intersections as Int

LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())

def loadTracks(trackfile):
    """
    Read tracks from a track .nc file and return a list of :class:`Track`
    objects.

    This calls the function `ncReadTrackData` to parse the track .nc file.

    :type  trackfile: str
    :param trackfile: the track data filename.
    """
    tracks = ncReadTrackData(trackfile)
    return tracks

class LongitudeCrossing(object):

    def __init__(self, configFile):

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

        # Determine TCRM input directory
        tcrm_dir = pathLocator.getRootDirectory()
        self.inputPath = pjoin(tcrm_dir, 'input')

        self.synNumYears = config.getint('TrackGenerator',
                                         'yearspersimulation')

        # Longitude crossing gates:
        self.gateLons = np.arange(self.lon_range.min(),
                                  self.lon_range.max() + 0.5, 10.)

        self.gateLats = np.arange(self.lat_range.min(),
                                  self.lat_range.max() + 0.5, 2.)

        # Add configuration settings to global attributes:
        self.gatts = {'history': "Longitude crossing rates for TCRM simulation",
                      'version': flProgramVersion()}

        for section in config.sections():
            for option in config.options(section):
                key = "{0}_{1}".format(section, option)
                value = config.get(section, option)
                self.gatts[key] = value


    def findCrossings(self, tracks):
        """
        Given a series of track points and a longitude, calculate
        if the tracks intersect that line of longitude.

        :param tracks: collection of :class:`Track` objects
        :return: h, ewh, weh, histograms for each line of longitude,
                 recording the rate of crossings
        """
        LOG.debug("Processing %d tracks" % (len(tracks)))
        h = np.zeros((len(self.gateLats) - 1, len(self.gateLons)))
        ewh = np.zeros((len(self.gateLats) - 1, len(self.gateLons)))
        weh = np.zeros((len(self.gateLats) - 1, len(self.gateLons)))

        for n, gLon in enumerate(self.gateLons):
            gStart = Int.Point(gLon, self.gateLats.max())
            gEnd = Int.Point(gLon, self.gateLats.min())
            lats = []
            ewlats = []
            welats = []

            for t in tracks:
                for i in range(len(t.Longitude) - 1):
                    cross = Int.Crossings()
                    start = Int.Point(t.Longitude[i], t.Latitude[i])
                    end = Int.Point(t.Longitude[i + 1], t.Latitude[i + 1])
                    r = cross.LineLine(start, end, gStart, gEnd)
                    if r.status == "Intersection":
                        lats.append(r.points[0].y)
                        startSide = Int._isLeft(gStart, gEnd, start)
                        endSide = Int._isLeft(gStart, gEnd, end)
                        if ((startSide < 0.) and (endSide >= 0.)) \
                                or ((startSide <= 0.) and (endSide > 0.)):
                            welats.append(r.points[0].y)

                        elif ((startSide > 0.) and (endSide <= 0.)) \
                                or ((startSide >= 0.) and (endSide < 0.)):
                            ewlats.append(r.points[0].y)

                    else:
                        # Track segment doesn't cross that longitude
                        continue

            # Generate the histograms to be returned:
            if len(lats) > 0:
                h[:, n], bins = np.histogram(lats, self.gateLats,
                                             density=True)
            if len(ewlats) > 0:
                ewh[:, n], bins = np.histogram(ewlats, self.gateLats,
                                               density=True)
            if len(welats) > 0:
                weh[:, n], bins = np.histogram(welats, self.gateLats,
                                               density=True)


        return h, ewh, weh

    def calcStats(self, lonCrossHist, lonCrossEW, lonCrossWE):
        """Calculate means and percentiles of synthetic event sets"""

        self.synCrossMean = np.mean(lonCrossHist, axis=0)
        self.synCrossEW = np.mean(lonCrossEW, axis=0)
        self.synCrossWE = np.mean(lonCrossWE, axis=0)

        self.synCrossUpper = percentile(lonCrossHist, per=95, axis=0)
        self.synCrossEWUpper = percentile(lonCrossEW, per=95, axis=0)
        self.synCrossWEUpper = percentile(lonCrossWE, per=95, axis=0)

        self.synCrossLower = percentile(lonCrossHist, per=5, axis=0)
        self.synCrossEWLower = percentile(lonCrossEW, per=5, axis=0)
        self.synCrossWELower = percentile(lonCrossWE, per=5, axis=0)



    @disableOnWorkers
    def historic(self):
        """Calculate historical rates of longitude crossing"""

        LOG.info("Processing historical tracks for longitude crossings")
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        source = config.get('DataProcess', 'Source')

        timestep = config.getfloat('TrackGenerator', 'Timestep')

        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.inputPath, inputFile)

        try:
            tracks = interpolateTracks.parseTracks(self.configFile,
                                                   inputFile,
                                                   source,
                                                   timestep,
                                                   interpolation_type='linear')
        except (TypeError, IOError, ValueError):
            LOG.critical("Cannot load historical track file: {0}".\
                         format(inputFile))
            raise
        else:
            self.lonCrossingHist, self.lonCrossingEWHist, \
                self.lonCrossingWEHist = self.findCrossings(tracks)

        return

    def synthetic(self):
        """Calculate synthetic rates of longitude crossing"""

        LOG.info("Processing synthetic rates of longitude crossing")

        work_tag = 0
        result_tag = 1
        filelist = os.listdir(self.trackPath)
        trackfiles = sorted([pjoin(self.trackPath, f) for f in filelist
                             if f.startswith('tracks')])

        lonCrossHist = np.zeros((len(trackfiles),
                                 len(self.gateLats) - 1,
                                 len(self.gateLons)))
        lonCrossEW = np.zeros((len(trackfiles),
                               len(self.gateLats) - 1,
                               len(self.gateLons)))
        lonCrossWE = np.zeros((len(trackfiles),
                               len(self.gateLats) - 1,
                               len(self.gateLons)))

        if (comm.rank == 0) and (comm.size > 1):

            w = 0
            n = 0
            for d in range(1, comm.size):
                comm.Send(trackfiles[w], dest=d, tag=work_tag)
                LOG.debug("Processing track file {0:d} of {1:d}".\
                          format(w + 1, len(trackfiles)))
                w += 1

            terminated = 0
            while terminated < comm.size - 1:
                results, status = comm.Recv(MPI.ANY_SOURCE, tag=result_tag,
                                               status=None)

                lonCrossHist[n, :, :], lonCrossEW[n, :, :], \
                    lonCrossWE[n, :, :] = results
                n += 1

                d = status.source

                if w < len(trackfiles):
                    comm.Send(trackfiles[w], dest=d, tag=work_tag)
                    LOG.debug("Processing track file {0:d} of {1:d}".\
                              format(w + 1, len(trackfiles)))
                    w += 1
                else:
                    comm.Send(None, dest=d, tag=work_tag)
                    terminated += 1

            self.calcStats(lonCrossHist, lonCrossEW, lonCrossWE)

        elif (comm.size > 1) and (comm.rank != 0):
            while True:
                trackfile = comm.Recv(source=0, tag=work_tag)
                if trackfile is None:
                    break

                LOG.debug("Processing {0}".format(trackfile))
                tracks = loadTracks(trackfile)
                lonCross, lonCrossEW, lonCrossWE = self.findCrossings(tracks)
                results = (lonCross, lonCrossEW, lonCrossWE)
                comm.Send(results, dest=0, tag=result_tag)

        elif (comm.size == 1) and (comm.rank == 0):
            # Assumed no Pypar - helps avoid the need to extend DummyPypar()
            for n, trackfile in enumerate(sorted(trackfiles)):
                tracks = loadTracks(trackfile)
                lonCrossHist[n, :, :], lonCrossEW[n, :, :], \
                    lonCrossWE[n, :, :] = self.findCrossings(tracks)

            self.calcStats(lonCrossHist, lonCrossEW, lonCrossWE)

    @disableOnWorkers
    def save(self):
        """Save data to file for archival and/or further processing"""

        dataFile = pjoin(self.dataPath, 'lonCrossings.nc')
        LOG.debug("Saving longitude crossing data to %s" % dataFile)

        dimensions = {
            0: {
                'name': 'lat',
                'values': self.gateLats[:-1],
                'dtype': 'f',
                'atts': {
                    'long_name':'Latitude',
                    'units':'degrees_north',
                    'axis': 'Y'
                }
            },
            1: {
                'name': 'lon',
                'values': self.gateLons,
                'dtype': 'f',
                'atts': {
                    'long_name':'Longitude',
                    'units':'degrees_east',
                    'axis': 'X'
                }
            }
        }

        variables = {
            0: {
                'name': 'hist',
                'dims': ('lat', 'lon'),
                'values' :self.lonCrossingHist,
                'dtype': 'f',
                'atts': {
                    'long_name':'Historical longitudinal crossing rate',
                    'units':'number of crossings per year'
                }
            },
            1: {
                'name': 'hist_ew',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingEWHist,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Historical longitudinal crossing rate '
                                  '- east-west crossings'),
                    'units': 'number of crossings per year'
                }
            },
            2: {
                'name': 'hist_we',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingWEHist,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Historical longitudinal crossing rate '
                                  '- west-east crossings'),
                    'units': 'number of crossings per year'
                }
            },
            3: {
                'name': 'syn_mean',
                'dims':('lat', 'lon'),
                'values': self.synCrossMean,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean synthetic longitudinal crossing rate',
                    'units': 'number of crossings per year'
                }
            },
            4: {
                'name': 'syn_mean_ew',
                'dims': ('lat', 'lon'),
                'values': self.synCrossEW,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Mean synthetic longitudinal crossing rate '
                                  '- east-west crossings'),
                    'units':'number of crossings per year'
                }
            },
            5: {
                'name': 'syn_mean_we',
                'dims': ('lat', 'lon'),
                'values': self.synCrossWE,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Mean synthetic longitudinal crossing rate '
                                  '- west-east crossings'),
                    'units': 'number of crossings per year'
                }
            },
            6: {
                'name': 'syn_upper',
                'dims': ('lat', 'lon'),
                'values': self.synCrossUpper,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Upper percentile synthetic longitudinal ',
                                  'crossing rate'),
                    'units': 'number of crossings per year',
                    'percentile': 90
                }
            },
            7: {
                'name': 'syn_upper_ew',
                'dims': ('lat', 'lon'),
                'values': self.synCrossEWUpper,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Upper percentile synthetic longitudinal '
                                  'crossing rate - east-west crossings'),
                    'units': 'number of crossings per year',
                    'percentile': 90
                }
            },
            8: {
                'name': 'syn_upper_we',
                'dims': ('lat', 'lon'),
                'values': self.synCrossWEUpper,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Upper percentile synthetic longitudinal '
                                  'crossing rate - west-east crossings'),
                    'units': 'number of crossings per year',
                    'percentile': 90
                }
            },
            9: {
                'name': 'syn_lower',
                'dims': ('lat', 'lon'),
                'values': self.synCrossLower,
                'dtype': 'f',
                'atts': {
                    'long_name':('Lower percentile synthetic longitudinal '
                                 'crossing rate'),
                    'units':'number of crossings per year',
                    'percentile': 5
                }
            },
            10: {
                'name': 'syn_lower_ew',
                'dims': ('lat', 'lon'),
                'values': self.synCrossEWLower,
                'dtype': 'f',
                'atts': {
                    'long_name':('Lower percentile synthetic longitudinal '
                                 'crossing rate - east-west crossings'),
                    'units':'number of crossings per year',
                    'percentile': 5
                }
            },
            11: {
                'name': 'syn_lower_we',
                'dims': ('lat', 'lon'),
                'values': self.synCrossWELower,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Lower percentile synthetic longitudinal '
                                  'crossing rate - west-east crossings'),
                    'units': 'number of crossings per year',
                    'percentile': 5
                }
            }
        }

        ncSaveGrid(dataFile, dimensions, variables, gatts=self.gatts)

        return

    @disableOnWorkers
    def plotCrossingRates(self):
        """Plot longitude crossing rates"""
        LOG.debug("Plotting longitude crossing rates")
        fig = Figure()
        ax1 = fig.add_subplot(2, 1, 1)
        for i in range(len(self.gateLons)):
            ax1.plot(2.* self.gateLons[i] - 100. * self.lonCrossingEWHist[:, i],
                     self.gateLats[:-1], color='r', lw=2)

            ax1.plot(2.* self.gateLons[i] - 100. * self.synCrossEW[:, i],
                     self.gateLats[:-1], color='k', lw=2)

            x1 = 2.* self.gateLons[i] - 100. * self.synCrossEWUpper[:, i]
            x2 = 2.* self.gateLons[i] - 100. * self.synCrossEWLower[:, i]
            ax1.fill_betweenx(self.gateLats[:-1], x1, x2,
                              color='0.75', alpha=0.7)

        minLonLim = 2. * (self.lon_range.min() - 10.)
        maxLonLim = 2. * (self.lon_range.max() + 10.)
        ax1.set_xlim(minLonLim, maxLonLim)
        ax1.set_xticks(2. * self.gateLons)
        ax1.set_xticklabels(self.gateLons.astype(int))
        ax1.set_xlabel("East-west crossings")
        ax1.set_ylim(self.gateLats.min(), self.gateLats[-2])
        ax1.set_ylabel('Latitude')
        ax1.grid(True)

        ax2 = fig.add_subplot(2, 1, 2)
        for i in range(len(self.gateLons)):
            ax2.plot(2.* self.gateLons[i] + 100. * self.lonCrossingWEHist[:, i],
                     self.gateLats[:-1], color='r', lw=2)

            ax2.plot(2.* self.gateLons[i] + 100. * self.synCrossWE[:, i],
                     self.gateLats[:-1], color='k', lw=2)

            x1 = 2.* self.gateLons[i] + 100. * self.synCrossWEUpper[:, i]
            x2 = 2.* self.gateLons[i] + 100. * self.synCrossWELower[:, i]
            ax2.fill_betweenx(self.gateLats[:-1], x1, x2,
                              color='0.75', alpha=0.7)

        ax2.set_xlim(minLonLim, maxLonLim)
        ax2.set_xticks(2. * self.gateLons)
        ax2.set_xticklabels(self.gateLons.astype(int))

        ax2.set_xlabel("West-east crossings")
        ax2.set_ylim(self.gateLats.min(), self.gateLats[-2])
        ax2.set_ylabel('Latitude')
        ax2.grid(True)

        fig.tight_layout()
        canvas = FigureCanvas(fig)
        canvas.print_figure(pjoin(self.plotPath, 'lon_crossing_syn.png'))

        return


    def run(self):
        """Run the longitude crossing evaluation"""
        global MPI, comm
        MPI = attemptParallel()
        comm = MPI.COMM_WORLD
        self.historic()

        comm.barrier()

        self.synthetic()

        comm.barrier()

        self.plotCrossingRates()
        self.save()
