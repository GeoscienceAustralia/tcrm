"""
:mod:`LandfallRates` -- calculate landfall rates
================================================

.. module:: LandfallRates
    :synopsis: Given a set of cyclone tracks and gates around the coastline,
               bin each coastal crossing into the gate corresponding to it's
               landfall location - then repeat for a set of synthetic events.

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import logging

import numpy as np

from os.path import join as pjoin
from scipy.stats import scoreatpercentile as percentile
from configparser import NoOptionError

from Utilities.config import ConfigParser
from Utilities.track import ncReadTrackData
from Utilities.loadData import loadTrackFile
from Utilities.parallel import attemptParallel, disableOnWorkers

from Utilities import pathLocator
import Utilities.Intersections as Int

from PlotInterface.curves import RangeCompareCurve, saveFigure

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


class LandfallRates(object):

    def __init__(self, configFile):

        config = ConfigParser()
        config.read(configFile)
        self.configFile = configFile

        outputPath = config.get('Output', 'Path')
        self.trackPath = pjoin(outputPath, 'tracks')
        self.plotPath = pjoin(outputPath, 'plots', 'stats')
        self.dataPath = pjoin(outputPath, 'process')

        # Determine TCRM input directory
        tcrm_dir = pathLocator.getRootDirectory()
        self.inputPath = pjoin(tcrm_dir, 'input')

        self.synNumYears = config.getint('TrackGenerator',
                                         'yearspersimulation')

        try:
            gateFile = config.get('Input', 'CoastlineGates')
        except NoOptionError:
            LOG.exception(("No coastline gate file specified "
                           "in configuration file"))
            raise

        gateData = np.genfromtxt(gateFile, delimiter=',')

        self.gates = Int.convert2vertex(gateData[:, 1], gateData[:, 2])
        self.coast = list(self.gates)
        self.coast.append(self.gates[0])



    def processTracks(self, tracks):
        """
        Given a collection of :class:`Track` objects and set of gate vertices,
        calculate if the tracks cross the gates in either an onshore or
        offshore direction.

        Returns the histograms for the gate counts.

        """

        landfall = []
        offshore = []

        for t in tracks:
            for i in range(1, len(t.Longitude)):
                cross = Int.Crossings()
                start = Int.Point(t.Longitude[i-1], t.Latitude[i-1])
                end = Int.Point(t.Longitude[i], t.Latitude[i])

                startOnshore = Int.inLand(start, self.coast)
                endOnshore = Int.inLand(end, self.coast)

                if not startOnshore and endOnshore:
                    # Landfall:
                    cross = Int.Crossings()
                    for j in range(1, len(self.gates) - 1):
                        r = cross.LineLine(start, end,
                                           self.gates[j-1],
                                           self.gates[j])
                        if r.status == "Intersection":
                            landfall.append(j)

                elif startOnshore and not endOnshore:
                    # Moving offshore:
                    cross = Int.Crossings()
                    for j in range(1, len(self.gates) - 1):
                        r = cross.LineLine(start, end,
                                           self.gates[j - 1],
                                           self.gates[j])
                        if r.status == "Intersection":
                            offshore.append(j)

        # Generate the histograms to be returned:
        lh, n = np.histogram(landfall, np.arange(len(self.gates)), density=True)
        oh, n = np.histogram(offshore, np.arange(len(self.gates)), density=True)
        return lh, oh

    def processResults(self, results, index):
        """
        Populate the :attr:`self.synLandfall` and :attr:`self.synOffshore`
        attributes.

        :param tuple results: tuple of arrays containing landfall and offshore
                              transition counts for gates.
        :param int index: synthetic event counter.

        """
        sLF, sOF = results
        self.synLandfall[index, :] = sLF
        self.synOffshore[index, :] = sOF

    def calculateStats(self):
        """
        Calculate mean and percentiels of landfall/offshore transition
        rates. Operates on the :attr:`self.synLandfall` and
        :attr:`self.synOffshore` attributes.

        """

        self.synMeanLandfall = np.mean(self.synLandfall, axis=0)
        self.synMeanOffshore = np.mean(self.synOffshore, axis=0)

        self.synUpperLF = percentile(self.synLandfall, per=95, axis=0)
        self.synLowerLF = percentile(self.synLandfall, per=5, axis=0)
        self.synUpperOF = percentile(self.synOffshore, per=95, axis=0)
        self.synLowerOF = percentile(self.synOffshore, per=5, axis=0)

    @disableOnWorkers
    def setOutput(self, ntracks):
        """
        Set the size of the output arrays.

        :param int ntracks: Number of track events.

        """
        self.synLandfall = np.zeros((ntracks, len(self.gates) - 1))
        self.synOffshore = np.zeros((ntracks, len(self.gates) - 1))

    @disableOnWorkers
    def historic(self):
        """Calculate historical rates of landfall"""

        LOG.info("Processing landfall rates of historical tracks")
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        source = config.get('DataProcess', 'Source')


        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.inputPath, inputFile)

        try:
            tracks = loadTrackFile(self.configFile, inputFile, source)
        except (TypeError, IOError, ValueError):
            LOG.critical("Cannot load historical track file: {0}".\
                         format(inputFile))
            raise
        else:
            self.historicLandfall, self.historicOffshore = \
                                        self.processTracks(tracks)

        return

    def synthetic(self):
        """Load synthetic data and calculate histogram"""
        LOG.info("Processing landfall rates of synthetic events")

        work_tag = 0
        result_tag = 1
        filelist = os.listdir(self.trackPath)
        trackfiles = sorted([pjoin(self.trackPath, f) for f in filelist
                             if f.startswith('tracks')])

        self.setOutput(len(trackfiles))

        if (comm.rank == 0) and (comm.size > 1):

            w = 0
            n = 0
            for d in range(1, comm.size):
                comm.Send(trackfiles[w], dest=d, tag=work_tag)
                LOG.debug("Processing track file {0:d} of {1:d}".\
                          format(w + 1, len(trackfiles)))
                w += 1

            terminated = 0
            while terminated < comm.size() - 1:
                results, status = comm.Recv(MPI.ANY_SOURCE, tag=result_tag,
                                             status=True)

                self.processResults(results, n)
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

            self.calculateStats()

        elif (comm.size > 1) and (comm.rank != 0):
            while True:
                trackfile = comm.Recv(source=0, tag=work_tag)
                if trackfile is None:
                    break

                LOG.debug("Processing %s", trackfile)
                tracks = loadTracks(trackfile)
                results = self.processTracks(tracks)
                comm.Send(results, dest=0, tag=result_tag)

        elif comm.size == 1 and comm.rank == 0:
            # Assumed no Pypar - helps avoid the need to extend DummyPypar()
            for n, trackfile in enumerate(sorted(trackfiles)):
                LOG.debug("Processing track file {0:d} of {1:d}".\
                          format(n + 1, len(trackfiles)))
                tracks = loadTracks(trackfile)
                results = self.processTracks(tracks)
                self.processResults(results, n)

            self.calculateStats()

    @disableOnWorkers
    def plot(self):
        """
        Plot the results and save to file.

        """

        figure = RangeCompareCurve()
        figure.set_size_inches(8, 3)
        xlab = "Gate number"
        ylab = "Landfall probability"

        x = np.arange(len(self.gates) - 1)

        figure.add(x, self.historicLandfall, self.synMeanLandfall,
                   self.synUpperLF, self.synLowerLF,
                   xlab, ylab, "Landfall")

        ylab = "Offshore probability"
        figure.add(x, self.historicOffshore, self.synMeanOffshore,
                   self.synUpperOF, self.synLowerOF,
                   xlab, ylab, "Offshore")

        figure.plot()
        outputFile = pjoin(self.plotPath, 'landfall_rates.png')
        saveFigure(figure, outputFile)

    def run(self):
        """Execute the analysis"""
        global MPI, comm
        MPI = attemptParallel()
        comm = MPI.COMM_WORLD
        self.historic()
        comm.barrier()
        self.synthetic()
        comm.barrier()
        self.plot()
        #self.save()
