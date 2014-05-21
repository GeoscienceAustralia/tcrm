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

from matplotlib import pyplot
from mpl_toolkits.basemap import Basemap

import interpolateTracks

from Utilities.config import ConfigParser
from Utilities.metutils import convert
from Utilities.maputils import bearing2theta
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
    6: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[6], 'Pa'),
    7: lambda s: convert(float(s.strip() or 0), TRACKFILE_UNIT[7], 'Pa'),
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

def loadTracksFromFiles(trackfiles):
    for f in trackfiles:
        tracks = loadTracks(f)
        for track in tracks:
            yield track

def loadTracksFromPath(path):
    files = os.listdir(path)
    trackfiles = [pjoin(path, f) for f in files if f.startswith('tracks')]
    msg = 'Processing %d track files in %s' % (len(trackfiles), path)
    log.info(msg)
    return loadTracksFromFiles(sorted(trackfiles))
    
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
                         extend='max', pad=0.1)
    cb.ax.tick_params(direction='in')
                         
    if cb.orientation == 'horizontal':
        for t in cb.ax.get_xticklabels():
            t.set_fontsize(8)

    if clabel:
        cb.set_label(clabel)

    return

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
        self.dataPath = pjoin(outputPath, 'process','stats')

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
            lon = np.append(lon, t.Longitude)
            lat = np.append(lat, t.Latitude)  
        histogram, x, y = np.histogram2d(lon, lat,
                                         [self.lon_range,
                                          self.lat_range],
                                          normed=False)
        return histogram

    def historic(self):
        """Load historic data and calculate histogram"""
        config = ConfigParser()
        config.read(self.configFile)
        inputFile = config.get('DataProcess', 'InputFile')
        source = config.get('DataProcess', 'Source')
        
        timestep = config.getfloat('TrackGenerator', 'Timestep')
        
        path, base = os.path.split(inputFile)
        interpHistFile = pjoin(path, "interp_tracks.csv")
        try:
            tracks = interpolateTracks.parseTracks(self.configFile,
                                                   inputFile,
                                                   source,
                                                   timestep,
                                                   interpHistFile, 'linear')
        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historical track file: {0}".format(inputFile))
            raise
        else:
            #tracks = loadTrackFile(self.configFile,
            #                       inputFile,
            #                       source)
            startYr = 9999
            endYr = 0
            for t in tracks:
                startYr = min(startYr, min(t.Year))
                endYr = max(endYr, max(t.Year))
            numYears = endYr - startYr
            self.hist = self.calculate(tracks) / numYears
             


    def synthetic(self):
        """Load synthetic data and calculate histogram"""
        
        #config = ConfigParser()
        #config.read(self.configFile)
        #timestep = config.getfloat('TrackGenerator', 'Timestep')
        
        filelist = os.listdir(self.trackPath)
        trackfiles = [pjoin(self.trackPath, f) for f in filelist
                      if f.startswith('tracks')]
        self.synHist = -9999. * np.ones((len(trackfiles), 
                                         len(self.lon_range) - 1, 
                                     len(self.lat_range) - 1))
        for n, trackfile in enumerate(trackfiles):
            tracks = loadTracks(trackfile)
            self.synHist[n, :, :] = self.calculate(tracks) / self.synNumYears
            
        self.synHist = ma.masked_values(self.synHist, -9999.)
        self.synHistMean = ma.mean(self.synHist, axis=0)
        self.medSynHist = ma.median(self.synHist, axis=0)

        self.synHistUpper = percentile(ma.compressed(self.synHist), per=95, axis=0)
        self.synHistLower = percentile(ma.compressed(self.synHist), per=5, axis=0)
                    
    def save(self):
        dataFile = pjoin(self.dataPath, 'density.nc')

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
                        
    def plotHistoric(self):
        """
        Plot results
        """
        
        datarange = (0, self.hist.max())
        plotDensity(self.lon_range[:-1], self.lat_range[:-1], 
                    np.transpose(self.hist), datarange=datarange)
                    
    def plotSynthetic(self):
        """
        Plot results from synthetic events
        """
        datarange = (0, self.synHistMean.max())
        plotDensity(self.lon_range[:-1], self.lat_range[:-1], 
                    np.transpose(self.synHistMean), datarange=datarange) 
                    
    def plotTrackDensity(self):
        
        datarange = (0, self.hist.max())
        pyplot.figure()
        ax1 = pyplot.subplot(211)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1], 
                    np.transpose(self.hist), datarange=datarange,
                    clabel="TC observations/yr")
        ax1.text(self.lon_range[1], self.lat_range[1],"Historic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))
                    
        ax2 = pyplot.subplot(212)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1], 
                    np.transpose(self.synHistMean), datarange=datarange,
                    clabel="TC observations/yr") 
        ax2.text(self.lon_range[1], self.lat_range[1],"Mean synthetic",
                 bbox=dict(fc='white', ec='black', alpha=0.5))
        pyplot.savefig(pjoin(self.plotPath, 'track_density.png'))
                    
    def plotTrackDensityPercentiles(self):
        
        datarange = (0, self.hist.max())
        pyplot.figure()
        ax1 = pyplot.subplot(211)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1], 
                    np.transpose(self.synHistUpper), datarange=datarange,
                    clabel="TC observations/yr")
        ax1.text(self.lon_range[1], self.lat_range[1],"Upper percentile",
                 bbox=dict(fc='white', ec='black', alpha=0.5))
                    
        ax2 = pyplot.subplot(212)
        plotDensity(self.lon_range[:-1], self.lat_range[:-1], 
                    np.transpose(self.synHistLower), datarange=datarange,
                    clabel="TC observations/yr")
        ax2.text(self.lon_range[1], self.lat_range[1],"Lower percentile",
                 bbox=dict(fc='white', ec='black', alpha=0.5))
        pyplot.savefig(pjoin(self.plotPath, 'track_density_percentiles.png'))

                 
    def run(self):
        
        self.historic()
        self.synthetic()
        
        self.plotTrackDensity()
        self.plotTrackDensityPercentiles()
        
        #self.save()