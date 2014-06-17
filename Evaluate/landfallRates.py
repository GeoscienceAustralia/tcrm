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
from ConfigParser import NoOptionError

from functools import wraps

import interpolateTracks

from Utilities.config import ConfigParser
from Utilities.metutils import convert
from Utilities.maputils import bearing2theta
from Utilities.track import Track
from Utilities.loadData import loadTrackFile
#from Utilities.parallel import attemptParallel, disableOnWorkers

from Utilities.files import flProgramVersion
from Utilities import pathLocator
import Utilities.Intersections as Int

from PlotInterface.curves import RangeCompareCurve, saveFigure

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
    6: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[6], 'Pa'),
    7: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[7], 'Pa'),
}


def disableOnWorkers(f):
    """
    Disable function calculation on workers. Function will
    only be evaluated on the master.
    """
    @wraps(f)
    def wrap(*args, **kwargs):
        if pp.size() > 1 and pp.rank() > 0:
            return
        else:
            return f(*args, **kwargs)
    return wrap

def attemptParallel():
    """
    Attempt to load Pypar globally as `pp`.  If pypar cannot be loaded then a
    dummy `pp` is created.

    """

    global pp

    try:
        # load pypar for everyone

        import pypar as pp

    except ImportError:

        # no pypar, create a dummy one
        
        class DummyPypar(object):

            def size(self):
                return 1

            def rank(self):
                return 0

            def barrier(self):
                pass

        pp = DummyPypar()

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
    datas = []
    data = readTrackData(trackfile)
    if len(data) > 0:
        cycloneId = data['CycloneNumber']
        for i in range(1, np.max(cycloneId) + 1):
            datas.append(data[cycloneId == i])
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
            log.exception(("No coastline gate file specified "
                          "in configuration file"))
            raise
        
        gateData = np.genfromtxt(gateFile, delimiter=',')
        nGates = len(gateData)
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
        sLF, sOF = results
        self.synLandfall[index, :] = sLF
        self.synOffshore[index, :] = sOF
        
    def calculateStats(self):

        self.synMeanLandfall = np.mean(self.synLandfall, axis=0)
        self.synMeanOffshore = np.mean(self.synOffshore, axis=0)

        self.synUpperLF = percentile(self.synLandfall, per=95, axis=0)
        self.synLowerLF = percentile(self.synLandfall, per=5, axis=0)
        self.synUpperOF = percentile(self.synOffshore, per=95, axis=0)
        self.synLowerOF = percentile(self.synOffshore, per=5, axis=0)

    @disableOnWorkers
    def setOutput(self, ntracks):
        self.synLandfall = np.zeros((ntracks, len(self.gates) - 1))
        self.synOffshore = np.zeros((ntracks, len(self.gates) - 1))

    @disableOnWorkers
    def historic(self):
        """Calculate historical rates of landfall"""

        log.info("Processing landfall rates of historical tracks")
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        source = config.get('DataProcess', 'Source')
        
        timestep = config.getfloat('TrackGenerator', 'Timestep')

        if len(os.path.dirname(inputFile)) == 0:
            inputFile = pjoin(self.inputPath, inputFile)
        
        try:
            tracks = loadTrackFile(self.configFile, inputFile, source)
        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historical track file: {0}".format(inputFile))
            raise
        else:
            self.historicLandfall, self.historicOffshore = self.processTracks(tracks)

        return

    def synthetic(self):
        """Load synthetic data and calculate histogram"""
        log.info("Processing landfall rates of synthetic events")

        work_tag = 0
        result_tag = 1
        filelist = os.listdir(self.trackPath)
        trackfiles = sorted([pjoin(self.trackPath, f) for f in filelist
                             if f.startswith('tracks')])
        
        self.setOutput(len(trackfiles))
                       
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

                self.processResults(results, n)
                n += 1

                d = status.source

                if w < len(trackfiles):
                    pp.send(trackfiles[w], destination=d, tag=work_tag)
                    log.debug("Processing track file %d of %d" % (w + 1, len(trackfiles)))
                    w += 1
                else:
                    pp.send(None, destination=d, tag=work_tag)
                    terminated += 1

            self.calculateStats()
                        
        elif (pp.size() > 1) and (pp.rank() != 0):
            while(True):
                trackfile = pp.receive(source=0, tag=work_tag)
                if trackfile is None:
                    break
                
                log.debug("Processing %s" % (trackfile))
                tracks = loadTracks(trackfile)
                results = self.processTracks(tracks)
                pp.send(results, destination=0, tag=result_tag)
                
        elif pp.size() == 1 and pp.rank() == 0:
            # Assumed no Pypar - helps avoid the need to extend DummyPypar()
            for n, trackfile in enumerate(sorted(trackfiles)):
                log.debug("Processing track file %d of %d" % (w + 1, len(trackfiles)))
                tracks = loadTracks(trackfile)
                results = self.processTracks(tracks)
                self.processResults(results, n)
                
            self.calculateStats()

    @disableOnWorkers
    def plotLandfallRates(self):

        figure = RangeCompareCurve()
        figure.set_size_inches(8,3)
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
        """Run the longitude crossing evaluation"""
        
        attemptParallel()

        self.historic()

        pp.barrier()

        self.synthetic()

        pp.barrier()

        self.plotLandfallRates()
        #self.save()
