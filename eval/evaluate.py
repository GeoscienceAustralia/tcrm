#!/usr/bin/env python
"""
 Title: evaluate.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 11/03/11 12:56:PM
 Description: perform a series of analyses to qualitatively evaluate historic
 and synthetic TC track datasets.

 Version: $Rev: 803 $

 $Id: evaluate.py 803 2013-08-13 22:10:19Z carthur $
"""

import os
from os.path import join as pjoin
import sys
import getopt
import logging as log

import numpy as np
import numpy.ma as ma
import matplotlib as mpl

from matplotlib import pyplot, cm
from matplotlib.dates import date2num
from mpl_toolkits.basemap import Basemap
from scipy.stats import scoreatpercentile as percentile
from datetime import datetime

import interpolateTracks 
from Utilities.files import flConfigFile, flStartLog
from Utilities.config import cnfGetIniValue, ConfigParser
from Utilities.loadData import loadTrackFile
from Utilities.nctools import ncSaveGrid

import Utilities.Intersections as Int
import Utilities.colours as colours

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

__version__ = "$Id: evaluate.py 803 2013-08-13 22:10:19Z carthur $"


# FIXME Set matplotlib defaults - upgrade to use matplotlibrc file
mpl.rcParams['savefig.dpi'] = 200
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'


def ShowSyntax(exit_code=0):
    """Documentation function to describe how to use this funtion"""
    if sys.stderr.isatty() and sys.stdin.isatty():
        # Ensure that STDERR can be written to
        print "{0}: ".format(sys.argv[0])
        print __version__
        print ""
        print "Perform an analysis of an historical event set and a series of "
        print "synthetic event sets."
        print ""
        print "For the synthetic event sets, the program calculates mean"
        print "values for comparison to the historic event set, and"
        print "also calculates the upper and lower percentiles (default is"
        print "the 5th and 95th percentiles) of the range of synthetic events"
        print ""
        print "Input:"
        print "Configuration file: {0}, or as specified by the -c switch".format(flConfigFile())
        print ""
        print "Output:"
        print "A series of images comparing the track density of the"
        print "historical dataset to that of the synthetic datasets"
        print "Optional log file - default is {0}".format(flConfigFile('.log'))
        print ""

    sys.exit(exit_code)


def plotDensity(x, y, data, llLon=None, llLat=None, urLon=None, urLat=None,
                res='i', dl=10., datarange=(-1.,1.), cmap='jet_r', title=None,
                xlab='Longitude', ylab='Latitude', clabel=None, maskland=False,
                maskocean=False):
    """
    Plot a grid of values using pyplot.pcolormesh() to create a pseudocolour
    plot of a 2-dimensional array, with a basemap included.

    It updates the figure object, but does not have scope to save or display it.
    That way, the function can be used to generate a number of subplots, with
    the calling program handling the saving/display requirements.
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
    if hasattr(cm, cmap):
        cmap = getattr(cm, cmap)
    else:
        cmap = colours.colourMap(cmap, 'stretched')

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
        pyplot.ylabel(ylab, fontsize=7.5)
    if xlab:
        pyplot.xlabel(xlab, fontsize=7.5)
    if title:
        pyplot.title(title)
        
    pyplot.grid(True)
    pyplot.tick_params(direction='out', right='off', top='off')

    cb = pyplot.colorbar(shrink=0.5, aspect=30,  
                         orientation='horizontal',
                         extend='max', pad=0.1)

    if cb.orientation == 'horizontal':
        for t in cb.ax.get_xticklabels():
            t.set_fontsize(8)

    if clabel:
        cb.set_label(clabel)

    return


class Evaluate:
    """
    Base class to define the input data, grids and output locations
    that are used by the other classes.
    """

    def __init__(self, **kwargs):
        self.configFile = kwargs.get('configFile')

        # Define the historical information:
        self.historicTrackFile = kwargs.get('historicTrackFile') 
        self.historicFormat = kwargs.get('historicFormat')
        self.historicNumYears = float(kwargs.get('historicNumYears'))

        # Define the synthetic information:
        self.synTrackPath = kwargs.get('synTrackPath')
        self.synFormat = kwargs.get('synFormat') 
        self.synNumYears = float(kwargs.get('synNumYears'))
        self.synNumSimulations = kwargs.get('synNumSimulations')
        self.percentile = kwargs.get('Percentile', 90)
        self.lower = (100 - self.percentile)/ 2.
        self.upper = 100 - self.lower

        # We default to Australian region:
        self.minLon = kwargs.get('MinLongitude', 60.)
        self.maxLon = kwargs.get('MaxLongitude', 180.)
        self.minLat = kwargs.get('MinLatitude', -40.)
        self.maxLat = kwargs.get('MaxLatitude', 0.)
        self.gridSize = kwargs.get('GridSize', 2.0)
        self.lonRange = np.arange(self.minLon,
                                  self.maxLon + 0.1,
                                  self.gridSize)

        self.latRange = np.arange(self.minLat,
                                  self.maxLat + 0.1,
                                  self.gridSize)

        # Longitude crossing gates:
        self.gateLons = np.arange(self.minLon, self.maxLon, 10.)
        self.gateLats = np.arange(self.minLat, self.maxLat + 0.5, 2.5)

        # location to store output:
        self.plotPath = kwargs.get('PlotPath')
        self.dataPath = kwargs.get('DataPath') 
        self.cmap = kwargs.get('ColourMap', 'hot_r')


        self.nx = len(self.lonRange) - 1
        self.ny = len(self.latRange) - 1
        self.x = self.lonRange[:-1]
        self.y = self.latRange[:-1]
        self.hist2DShape = (self.nx, self.ny)
        self.minCpRange = np.arange(800., 1020., 5.)
        self.histShape = (len(self.minCpRange) - 1,)

        self.ageBins = np.arange(0, 481, 12)
        self.histAgeDist = np.empty((len(self.ageBins) - 1))
        self.synAgeDist = np.empty((self.synNumSimulations,
                                        len(self.ageBins) - 1))

        # Create dimensions for the output netcdf files:
        self.dimensions = {
            0: {
                'name': 'lat',
                'values': self.latRange[:-1],
                'dtype': 'f',
                'atts': {
                    'long_name': 'Latitude',
                    'units': 'degrees_north'
                } 
            },
            1: {
                'name': 'lon',
                'values': self.lonRange[:-1],
                'dtype': 'f',
                'atts': {
                    'long_name': 'Longitude',
                    'units':'degrees_east'
                } 
            } 
        }

        self.lonCrossingHist = np.empty((len(self.gateLats) - 1,
                                              len(self.gateLons)))

        self.lonCrossingEWHist = np.empty((len(self.gateLats) - 1,
                                                len(self.gateLons)))

        self.lonCrossingWEHist = np.empty((len(self.gateLats) - 1,
                                                len(self.gateLons)))

        self.lonCrossingSyn = np.empty((self.synNumSimulations,
                                             len(self.gateLats) - 1,
                                             len(self.gateLons)))

        self.lonCrossingSynEW = np.empty((self.synNumSimulations,
                                               len(self.gateLats) - 1,
                                               len(self.gateLons)))

        self.lonCrossingSynWE = np.empty((self.synNumSimulations,
                                               len(self.gateLats) - 1,
                                               len(self.gateLons)))

        self.lonCrossingSynMean = np.empty((len(self.gateLats) - 1,
                                                 len(self.gateLons)))

        self.lonCrossingSynEWMean = np.empty((len(self.gateLats) - 1,
                                                   len(self.gateLons)))

        self.lonCrossingSynWEMean = np.empty((len(self.gateLats) - 1,
                                                   len(self.gateLons)))

        self.lonCrossingSynUpper = np.empty((len(self.gateLats) - 1,
                                                  len(self.gateLons)))

        self.lonCrossingSynEWUpper = np.empty((len(self.gateLats) - 1,
                                                    len(self.gateLons)))

        self.lonCrossingSynWEUpper = np.empty((len(self.gateLats) - 1,
                                                    len(self.gateLons)))

        self.lonCrossingSynLower = np.empty((len(self.gateLats) - 1,
                                                  len(self.gateLons)))

        self.lonCrossingSynEWLower = np.empty((len(self.gateLats) - 1,
                                                    len(self.gateLons)))

        self.lonCrossingSynWELower = np.empty((len(self.gateLats) - 1,
                                                    len(self.gateLons)))

        # Create global attributes for the netcdf file:
        self.gatts = {
            'historical_source': self.historicTrackFile,
            'historical_format': self.historicFormat,
            'synthetic_source': self.synTrackPath,
            'synthetic_format': self.synFormat,
            'num_synthetic_events': self.synNumSimulations,
            'program_name': sys.argv[0],
            'program_version': __version__
        }

    def calc2DHistogram(self, lon, lat):
        h, x, y = np.histogram2d(lon, lat, [self.lonRange, self.latRange],
                                   normed=False)
        return h, x, y

    def calcHistogram(self, values, bins):
        h, n = np.histogram(values, bins, normed=False)
        return h, n

class EvalPressureDistribution(Evaluate):

    def processPressureData(self, values, lon, lat):

        dataMean = np.zeros(self.hist2DShape)
        dataMin = np.zeros(self.hist2DShape)
        dataMax = np.zeros(self.hist2DShape)
        dataMed = np.zeros(self.hist2DShape)
        
        for k in xrange(self.nx):
            for l in xrange(self.ny):
                ii = np.where(((lat >= self.latRange[l]) &
                                (lat < self.latRange[l] + self.gridSize)) &
                               ((lon >= self.lonRange[k]) &
                                (lon < self.lonRange[k] + self.gridSize)))[0]
                if len(ii) > 0:
                    vv = values[ii].compress(values[ii] < sys.maxint)
                    vv = vv.compress(vv > 0.0)
                    if len(vv > 0):
                        dataMean[k, l] = np.mean(vv)
                        dataMin[k, l] = np.min(vv)
                        dataMax[k, l] = np.max(vv)
                        dataMed[k, l] = np.median(vv)

        dataMean = ma.masked_equal(dataMean, 0)
        dataMin = ma.masked_equal(dataMin, 0)
        dataMax = ma.masked_equal(dataMax, 0)
        dataMed = ma.masked_equal(dataMed, 0)
        return dataMean, dataMin, dataMax, dataMed

    def calcMinPressure(self, index, pressure):
        minpressure = 1030.
        minCP = np.zeros(index.sum())
        event = 0
        
        for i in xrange(len(index) - 1):
            if index[i + 1] == 1:
                minCP[event] = minpressure
                event += 1
                minpressure = 1030.
            if pressure[i] < minpressure:
                minpressure = pressure[i]

        minCP[-1] = minpressure
        h, n = self.calcHistogram(minCP, self.minCpRange)
        return h

    def historic(self):
        log.info("Processing historical tracks for pressure distribution")
        try:
            [i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe] = \
                interpolateTracks.parseTracks(self.configFile,
                                              self.historicTrackFile,
                                              self.historicFormat, 6.0)
        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historical track file: {0}".format(self.historicTrackFile))
            return False
        else:
            self.pHistMean, self.pHistMin, self.pHistMax, self.pHistMed = \
                self.processPressureData(p, lon, lat)
                
            self.histMinCP = self.calcMinPressure(i, p)
            return True

    def synthetic(self):
        log.info("Processing {0} synthetic events in {1}".format(
                    self.synNumSimulations, self.synTrackPath))
                    
        self.synMean = np.empty(((self.synNumSimulations,) + self.hist2DShape))
        self.synMin = np.empty(((self.synNumSimulations,) + self.hist2DShape))
        self.synMed = np.empty(((self.synNumSimulations,) + self.hist2DShape))

        self.synMeanUpper = np.empty(self.hist2DShape)
        self.synMeanLower = np.empty(self.hist2DShape)
        self.synMinUpper = np.empty(self.hist2DShape)
        self.synMinLower = np.empty(self.hist2DShape)
        self.synMinCP = np.empty((self.synNumSimulations,) + self.histShape)

        for n in xrange(self.synNumSimulations):
            trackFile = pjoin(self.synTrackPath, "tracks.%04d.csv"%(n))
            log.debug("Processing {0}".format(trackFile))
            try:
                i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe = \
                    loadTrackFile(self.configFile,
                                  trackFile,
                                  self.synFormat)
            except (TypeError, IOError, ValueError):
                log.critical("Cannot load synthetic track file: {0}".format(trackFile))
                return False
            else:
                sMean, sMin, sMax, sMed = self.processPressureData(p, lon, lat)
            finally:
                self.synMinCP[n, :] = self.calcMinPressure(i, p)
                self.synMean[n, :, :] = sMean
                self.synMin[n, :, :] = sMin
                self.synMed[n, :, :] = sMed

        return True

    def synStats(self):
        log.info("Calculating statistics...")
        msynMean = ma.masked_equal(self.synMean, 0)
        msynMin = ma.masked_equal(self.synMin, 0)
        msynMed = ma.masked_equal(self.synMed, 0)
        self.meanSynMean = np.mean(msynMean, axis=0)
        self.meanSynMin = np.mean(msynMin, axis=0)
        self.meanSynMed = np.mean(msynMed, axis=0)

        self.synMeanUpper = percentile(msynMean, per=self.upper)
        self.synMeanLower = percentile(msynMean, per=self.lower)
        self.synMinUpper = percentile(msynMin, per=self.upper)
        self.synMinLower = percentile(msynMin, per=self.lower)

        self.synMinCPUpper = [ percentile(self.synMinCP[:,i], per=self.upper) 
                                for i in xrange(len(self.minCpRange) - 1) ]
                                
        self.synMinCPLower = [ percentile(self.synMinCP[:,i], per=self.lower) 
                                for i in xrange(len(self.minCpRange) - 1) ]
                                
        self.synMinCPMean = np.mean(self.synMinCP, axis=0)

        return

    def plotPressureMaps(self):
        if not hasattr(self, 'histMinCP'):
            log.critical("No historical minimum central pressure information calculated")
            return False
        else:
            log.debug("Plotting map of mean pressures")
            outputFile = pjoin(self.plotPath, 'meanPressure.png')

            [xx, yy] = np.meshgrid(self.x, self.y)
            pyplot.figure(1)
            pyplot.clf()
            ax1 = pyplot.subplot(211)
            plotDensity(xx, yy, np.transpose(self.pHistMean),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(950., 1000.),
                        cmap='hot', clabel=None, title=None, 
                        xlab='Longitude', ylab='Latitude', 
                        maskland=False,
                        maskocean=False)

            ax2 = pyplot.subplot(212)
            plotDensity(xx, yy, np.transpose(self.meanSynMean),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(950., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', clabel='Mean central pressure (hPa)',
                        maskland=False, 
                        maskocean=False)

            pyplot.savefig(outputFile)

            log.debug("Plotting map of minimum pressures")
            pyplot.figure(2)
            pyplot.clf()
            ax1 = pyplot.subplot(211)
            plotDensity(xx, yy, np.transpose(self.pHistMin), llLon=None,
                        llLat=None, urLon=None, urLat=None, res='i',
                        dl=20., datarange=(850., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', clabel=None, 
                        maskland=False,
                        maskocean=False)
                        
            ax2 = pyplot.subplot(212)
            plotDensity(xx, yy, np.transpose(self.meanSynMin),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(850., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', 
                        clabel='Minimum central pressure (hPa)', 
                        maskland=False,
                        maskocean=False)

            minPrsImgFile = pjoin(self.plotPath, 'minPressure.png')
            pyplot.savefig(minPrsImgFile)

            pyplot.figure(3)
            ax1 = pyplot.subplot(211)
            plotDensity(xx, yy, np.transpose(self.synMeanUpper),
                        llLon=None, llLat=None, urLon=None, urLat=None, 
                        res='i', dl=20., datarange=(950., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', clabel=None, 
                        maskland=False,
                        maskocean=False)

            ax2 = pyplot.subplot(212)
            plotDensity(xx, yy, np.transpose(self.synMeanLower),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(950., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', clabel='Mean central pressure (hPa)', 
                        maskland=False,
                        maskocean=False)
                        
            pyplot.savefig(pjoin(self.plotPath, 
                                 'synMeanPressurePercentile.png'))

            pyplot.figure(4)
            ax1 = pyplot.subplot(211)
            plotDensity(xx, yy, np.transpose(self.synMinUpper),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(850., 1000.),
                        cmap='hot',title=None, xlab='Longitude',
                        ylab='Latitude', clabel=None, 
                        maskland=False,
                        maskocean=False)

            ax2 = pyplot.subplot(212)
            plotDensity(xx, yy, np.transpose(self.synMinLower),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(850., 1000.),
                        cmap='hot',title=None, xlab='Longitude',
                        ylab='Latitude', 
                        clabel='Minimum central pressure (hPa)', 
                        maskland=False,
                        maskocean=False)
                        
            pyplot.savefig(pjoin(self.plotPath, 'synMinPressurePercentile.png'))

            log.debug("Plotting map of median pressures")
            pyplot.figure(5)
            pyplot.clf()
            ax1 = pyplot.subplot(211)
            plotDensity(xx, yy, np.transpose(self.pHistMed),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20.,datarange=(950., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', clabel=None, 
                        maskland=False,
                        maskocean=False)

            ax2 = pyplot.subplot(212)
            plotDensity(xx, yy, np.transpose(self.meanSynMed),
                        llLon=None, llLat=None, urLon=None,
                        urLat=None, res='i', dl=20., datarange=(950., 1000.),
                        cmap='hot', title=None, xlab='Longitude',
                        ylab='Latitude', 
                        clabel='Median central pressure (hPa)', 
                        maskland=False,
                        maskocean=False)
                        
            pyplot.savefig(pjoin(self.plotPath, 'medPressure.png'))


            return True

    def plotPressureDistributions(self):

        if not hasattr(self, 'histMinCP'):
            log.critical("No historical minimum central pressure information calculated")
            return False
        else:
            ax1 = pyplot.subplot(111)
            ax1.plot(self.minCpRange[:-1], self.histMinCP * self.synNumYears / \
                        self.historicNumYears, color='r', lw=2)
            ax1.plot(self.minCpRange[:-1], self.synMinCPMean, color='k', lw=2)
            ax1.fill_between(self.minCpRange[:-1], self.synMinCPUpper, 
                             self.synMinCPLower, color='0.5', alpha=0.5)
                             
            pyplot.xlim(850., 1020.)
            pyplot.xticks()
            pyplot.xlabel('Minimum pressure')
            pyplot.yticks()
            pyplot.ylabel('Number')
            pyplot.grid(True)

            outputFile = pjoin(self.plotPath, 'min_pressure_dist.png')
            pyplot.savefig(outputFile)

        return True

class EvalTrackDensity(Evaluate):

    def historic(self):
        log.info("Processing historic tracks...")
        try:
            [i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe] = \
                interpolateTracks.parseTracks(self.configFile,
                                              self.historicTrackFile,
                                              self.historicFormat, 6.0)
        except (TypeError, IOError, ValueError):
            log.critical("Cannot load historic track file: {0}".\
                            format(self.historicTrackFile))
            return False
        else:
            self.hist, x, y = self.calc2DHistogram(lon, lat)

        return True

    def synthetic(self):
        log.info("Processing {0} synthetic events in {1}".\
                    format(self.synNumSimulations, self.synTrackPath))
                    
        self.synHist = np.empty(((self.synNumSimulations,) + self.hist2DShape))
        for n in xrange(self.synNumSimulations):
            trackFile = pjoin(self.synTrackPath, "tracks.%04d.csv"%(n))
            log.debug("Processing {0}".format(trackFile))
            try:
                i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe = \
                    loadTrackFile(self.configFile, 
                                  trackFile, 
                                  self.synFormat)
            except (TypeError, IOError, ValueError):
                log.critical("Cannot load synthetic track file: {0}".\
                                format(trackFile))
                return False
            else:
                self.synHist[n, :, :], x, y = self.calc2DHistogram(lon, lat)


        return True


    def synStats(self):
        """
        Calculate mean, median and percentile values for the synthetic datasets:
        """
        if not hasattr(self, 'synHist'):
            log.critical("Synthetic event sets have not been processed!")
            log.critical("Cannot calculate statistics")
            return False
        else:
            self.meanSynHist = np.mean(self.synHist, axis=0)
            self.medSynHist = np.median(self.synHist, axis=0)
            self.synHistUpper = np.empty(self.hist2DShape)
            self.synHistLower = np.empty(self.hist2DShape)
                
            """
            for k in xrange(self.nx):
                for l in xrange(self.ny):
                    self.synHistUpper[k,l] = percentile(self.synHist[:,k,l], per=self.upper)
                    self.synHistLower[k,l] = percentile(self.synHist[:,k,l], per=self.lower)
            """
            self.synHistUpper = percentile(self.synHist, per=self.upper)
            self.synHistLower = percentile(self.synHist, per=self.lower)
            
        return True

    def saveData(self):
        """
        Save data to a netcdf file for archival and/or furthur processing
        """

        dataFile = pjoin(self.dataPath, 'trackDensity.nc')

        # Simple sanity check (should also include the synthetic data):
        if not hasattr(self, 'hist'):
            log.critical("No historical data available!")
            log.critical("Check that data has been processed before trying to save data")
            return

        log.info('Saving track density data to {0}'.format(dataFile))

        # Define variables:
        variables = {
            0: {
                'name': 'hist_density',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.hist / self.historicNumYears),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Historical track density',
                    'units':'observations per 1-degree grid per year'
                } 
            },
            1: {
                'name': 'syn_density',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.meanSynHist / self.synNumYears),
                'dtype': 'f',
                'atts': {
                    'long_name': 'Track density - synthetic events',
                    'units':'observations per 1-degree grid per year'
                } 
            },
            2: {
                'name': 'syn_density_upper',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistUpper / self.synNumYears),
                'dtype': 'f',
                'atts': {
                    'long_name': ('Track density - upper percentile '
                                  '- synthetic events'),
                    'units':' observations per 1-degree grid per year',
                    'percentile': self.upper 
                } 
            },
            3: {
                'name': 'syn_density_lower',
                'dims': ('lat', 'lon'),
                'values': np.transpose(self.synHistLower / self.synNumYears),
                'dtype': 'f',
                'atts': {
                    'long_name': ('Track density - lower percentile '
                                  '- synthetic events'),
                    'units': 'observations per 1-degree grid per year',
                    'percentile': self.lower 
                } 
            } 
        }

        ncSaveGrid(dataFile, self.dimensions, variables, gatts=self.gatts)

        return

    def plotTrackDensity(self):
        """
        Plot the historic and synthetic mean track density fields
        """
        [xx, yy] = np.meshgrid(self.x, self.y)
        outputFile = pjoin(self.plotPath, 'trackDensity.png')

        log.info("Plotting data to {0}".format(self.plotPath))
        pyplot.figure()
        pyplot.clf()
        ax1 = pyplot.subplot(211)
        plotDensity(xx, yy, np.transpose(self.hist / self.historicNumYears),
                    llLon=None, llLat=None, urLon=None,
                    urLat=None, res='i', dl=20., datarange=(0., 2.),
                    cmap=self.cmap, title=None, xlab='Longitude',
                    ylab='Latitude', clabel=None,
                    maskland=False,
                    maskocean=False)

        ax2 = pyplot.subplot(212)
        plotDensity(xx, yy, np.transpose(self.meanSynHist / self.synNumYears),
                    llLon=None, llLat=None, urLon=None,
                    urLat=None, res='i', dl=20., datarange=(0., 2.),
                    cmap=self.cmap, title=None, xlab='Longitude',
                    ylab='Latitude', clabel='Track density (points/year)',
                    maskland=False,
                    maskocean=False)

        pyplot.savefig(outputFile)
        return

    def plotPercentiles(self):
        """
        Plot the upper and lower percentile values of the data
        """
        [xx, yy] = np.meshgrid(self.x, self.y)
        outputFile = pjoin(self.plotPath, 'synTrackDensityCI.png')

        log.info("Plotting data to {0}".format(self.plotPath))
        pyplot.figure()
        ax1 = pyplot.subplot(211)
        plotDensity(xx, yy, np.transpose(self.synHistUpper / self.synNumYears),
                    llLon=None, llLat=None, urLon=None,
                    urLat=None, res='i', dl=20., datarange=(0, 2.),
                    cmap=self.cmap, title=None, xlab='Longitude',
                    ylab='Latitude', clabel=None,
                    maskland=False,
                    maskocean=False)

        ax4 = pyplot.subplot(212)
        plotDensity(xx, yy, np.transpose(self.synHistLower / self.synNumYears),
                    llLon=None, llLat=None, urLon=None,
                    urLat=None, res='i', dl=20., datarange=(0, 2.),
                    cmap=self.cmap, title=None, xlab='Longitude',
                    ylab='Latitude', clabel='Track density (points/year)',
                    maskland=False,
                    maskocean=False)

        pyplot.savefig(outputFile)
        return

class EvalLongitudeCrossings(Evaluate):

    def findCrossings(self, index, lon, lat):
        """
        Given a series of track points and a longitude, calculate
        if the tracks intersect that line of longitude
        """

        h = np.zeros((len(self.gateLats) - 1, len(self.gateLons)))
        ewh = np.zeros((len(self.gateLats) - 1, len(self.gateLons)))
        weh = np.zeros((len(self.gateLats) - 1, len(self.gateLons)))

        for n, gLon in enumerate(self.gateLons):
            gStart = Int.Point(gLon, self.gateLats.max())
            gEnd = Int.Point(gLon, self.gateLats.min())
            lats = []
            ewlats = []
            welats = []

            for i in xrange(1, len(index)):
                if (index[i] == index[i - 1]) or (index[i] == 0):
                    # Continuation of a track
                    cross = Int.Crossings()
                    start = Int.Point(lon[i - 1], lat[i - 1])
                    end = Int.Point(lon[i], lat[i])
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
            h[:, n], bins = np.histogram(lats, self.gateLats, density=False)
            ewh[:, n], bins = np.histogram(ewlats, self.gateLats, density=False)
            weh[:, n], bins = np.histogram(welats, self.gateLats, density=False)

        return h, ewh, weh

    def historic(self):
        """Calculate historical rates of longitude crossing"""
        
        log.debug("Processing historical tracks for longitude crossings")
        [i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe] = \
                interpolateTracks.parseTracks(self.configFile,
                                              self.historicTrackFile,
                                              self.historicFormat, 6.0)
                                  
        self.lonCrossingHist, self.lonCrossingEWHist, \
            self.lonCrossingWEHist = self.findCrossings(i, lon, lat)

        return

    def synthetic(self):
        """Calculate synthetic rates of longitude crossing"""
        
        log.debug("Processing {0} synthetic events in {1}".\
                    format(self.synNumSimulations, self.synTrackPath))
        for n in xrange(self.synNumSimulations):
            trackFile = pjoin(self.synTrackPath, "tracks.%04d.csv"%(n))
            log.debug("Processing {0}".format(trackFile))
            try:
                i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe = \
                  loadTrackFile(self.configFile, 
                                trackFile, 
                                self.synFormat)
            except (TypeError, IOError, ValueError):
                log.critical("Cannot load synthetic track file: {0}".\
                                format(trackFile))
                return False
            else:
                self.lonCrossingSyn[n, :], self.lonCrossingSynEW[n, :], \
                  self.lonCrossingSynWE[n, :] = self.findCrossings(i, lon, lat)

        return True

    def synStats(self):
        """Calculate statistics of synthetic event sets"""
        
        log.debug("Calculating statistics for longitude crossings of synthetic events")
        if not hasattr(self, 'lonCrossingSyn'):
            log.critical("Synthetic event sets have not been processed!")
            log.critical("Cannot calculate statistics")
            return False
        else:
            self.lonCrossingSynMean = np.mean(self.lonCrossingSyn, axis=0)
            self.lonCrossingSynEWMean = np.mean(self.lonCrossingSynEW, axis=0)
            self.lonCrossingSynWEMean = np.mean(self.lonCrossingSynWE, axis=0)

            for k in xrange(len(self.gateLats) - 1):
                for l in xrange(len(self.gateLons)):
                    self.lonCrossingSynUpper[k, l] = percentile(self.lonCrossingSyn[:, k, l], per=self.upper)
                    self.lonCrossingSynEWUpper[k, l] = percentile(self.lonCrossingSynEW[:, k, l], per=self.upper)
                    self.lonCrossingSynWEUpper[k, l] = percentile(self.lonCrossingSynWE[:, k, l], per=self.upper)
                    self.lonCrossingSynLower[k, l] = percentile(self.lonCrossingSyn[:, k, l], per=self.lower)
                    self.lonCrossingSynEWLower[k, l] = percentile(self.lonCrossingSynEW[:, k, l], per=self.lower)
                    self.lonCrossingSynWELower[k, l] = percentile(self.lonCrossingSynWE[:, k, l], per=self.lower)

        return True

    def saveData(self):
        """Save data to file for archival and/or further processing"""

        log.debug("Saving longitude crossing data to file")
        dataFile = pjoin(self.dataPath, 'lonCrossings.nc')

        dimensions = {
            0: {
                'name': 'lat',
                'values': self.gateLats[:-1],
                'dtype': 'f',
                'atts': {
                    'long_name':'Latitude',
                    'units':'degrees_north'
                }
            },
            1: {
                'name': 'lon',
                'values': self.gateLons,
                'dtype': 'f',
                'atts': {
                    'long_name':'Longitude',
                    'units':'degrees_east'
                }
            }
        }

        variables = {
            0: {
                'name': 'hist',
                'dims': ('lat', 'lon'),
                'values' :self.lonCrossingHist / self.historicNumYears,
                'dtype': 'f',
                'atts': {
                    'long_name':'Historical longitudinal crossing rate',
                    'units':'number of crossings per year'
                }
            },
            1: {
                'name': 'hist_ew',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingEWHist / self.historicNumYears,
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
                'values': self.lonCrossingWEHist / self.historicNumYears,
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
                'values': self.lonCrossingSynMean / self.synNumYears,
                'dtype': 'f',
                'atts': {
                    'long_name': 'Mean synthetic longitudinal crossing rate',
                    'units': 'number of crossings per year'
                }
            },
            4: {
                'name': 'syn_mean_ew',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingSynEWMean / self.synNumYears,
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
                'values': self.lonCrossingSynWEMean / self.synNumYears,
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
                'values': self.lonCrossingSynUpper / self.synNumYears,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Upper percentile synthetic longitudinal ',
                                  'crossing rate' ),
                    'units': 'number of crossings per year',
                    'percentile': self.upper
                }
            },
            7: {
                'name': 'syn_upper_ew',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingSynEWUpper / self.synNumYears,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Upper percentile synthetic longitudinal '
                                  'crossing rate - east-west crossings'),
                    'units': 'number of crossings per year',
                    'percentile': self.upper
                }
            },
            8: {
                'name': 'syn_upper_we',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingSynWEUpper / self.synNumYears,
                'dtype': 'f',
                'atts': {
                    'long_name': ('Upper percentile synthetic longitudinal '
                                  'crossing rate - west-east crossings'),
                    'units': 'number of crossings per year',
                    'percentile': self.upper
                }
            },
            9: {
                'name': 'syn_lower',
                'dims': ('lat', 'lon'),
                'values': self.lonCrossingSynLower / self.synNumYears,
                'dtype': 'f',
                'atts': {
                    'long_name':('Lower percentile synthetic longitudinal '
                                 'crossing rate'),
                    'units':'number of crossings per year',
                    'percentile':self.lower
                }
            },
            10: {
                 'name': 'syn_lower_ew',
                 'dims': ('lat', 'lon'),
                 'values': self.lonCrossingSynEWLower / self.synNumYears,
                 'dtype': 'f',
                 'atts': {
                    'long_name':('Lower percentile synthetic longitudinal '
                                  'crossing rate - east-west crossings'),
                    'units':'number of crossings per year',
                    'percentile':self.lower
                }
            },
            11: {
                 'name': 'syn_lower_we',
                 'dims': ('lat', 'lon'),
                 'values': self.lonCrossingSynWELower / self.synNumYears,
                 'dtype': 'f',
                 'atts': {
                    'long_name': ('Lower percentile synthetic longitudinal '
                                   'crossing rate - west-east crossings'),
                    'units': 'number of crossings per year',
                    'percentile': self.lower
                 }
            }
        }

        ncSaveGrid(dataFile, dimensions, variables, gatts=self.gatts)

        return

    def plotCrossingRates(self):
        """Plot the longitude crossing rates"""

        log.debug("Plotting longitude crossing rates")
        ax1 = pyplot.subplot(211)
        for i in range(len(self.gateLons)):
            ax1.plot(2 * self.gateLons[i] - self.lonCrossingEWHist[:, i] * \
                     self.synNumYears / self.historicNumYears,
                     self.gateLats[:-1], color='r', lw=2)
                     
            ax1.plot(2 * self.gateLons[i] - self.lonCrossingSynEWMean[:, i], \
                     self.gateLats[:-1],color='k',lw=2)
                     
            x1 = 2 * self.gateLons[i] - self.lonCrossingSynEWUpper[:, i]
            x2 = 2 * self.gateLons[i] - self.lonCrossingSynEWLower[:, i]
            ax1.fill_betweenx(self.gateLats[:-1], x1, x2, 
                              color='0.5', alpha=0.5)

        minLonLim = 2 * self.minLon
        maxLonLim = 2 * self.maxLon + 20.
        pyplot.xlim(minLonLim, maxLonLim)
        pyplot.xticks(2 * self.gateLons, self.gateLons.astype(int), fontsize=8)
        pyplot.xlabel("East-west crossings")
        pyplot.ylim(self.gateLats.min(), self.gateLats[-2])
        pyplot.yticks(fontsize=8)
        pyplot.ylabel('Latitude')
        ax1.tick_params(direction='out', top='off', right='off')
        pyplot.grid(True)

        ax2 = pyplot.subplot(212)
        for i in range(len(self.gateLons)):
            ax2.plot(2 * self.gateLons[i] + self.lonCrossingWEHist[:, i] * \
                     self.synNumYears / self.historicNumYears,
                     self.gateLats[:-1], color='r', lw=2)
                     
            ax2.plot(2 * self.gateLons[i] + self.lonCrossingSynWEMean[:, i],
                     self.gateLats[:-1], color='k', lw=2)
                     
            x1 = 2 * self.gateLons[i] + self.lonCrossingSynWEUpper[:, i]
            x2 = 2 * self.gateLons[i] + self.lonCrossingSynWELower[:, i]
            ax2.fill_betweenx(self.gateLats[:-1], x1, x2, 
                              color='0.5', alpha=0.5)

        pyplot.xlim(minLonLim, maxLonLim)  
        pyplot.xticks(2 * self.gateLons, self.gateLons.astype(int), fontsize=8)

        pyplot.xlabel("West-east crossings")
        pyplot.ylim(self.gateLats.min(), self.gateLats[-2])
        pyplot.yticks(fontsize=8)
        pyplot.ylabel('Latitude')
        pyplot.grid(True)
        ax2.tick_params(direction='out', top='off', right='off')
        pyplot.savefig(pjoin(self.plotPath,'lon_crossing_syn.png'))

        return

class EvalAgeDistribution(Evaluate):
    def calculateAge(self, index, yr, mon, day, hr, mn):
        """Calculate the age of all TC events in the input dataset"""

        start = np.where(index == 1)[0]
        end = np.empty(len(start))
        end[:-1] = start[1:] - 1
        end[-1] = len(index) - 1
        startDT = [datetime(yr[i], mon[i], day[i], hr[i], mn[i]) for i in start ]
        endDT = [datetime(yr[i], mon[i], day[i], hr[i], mn[i]) for i in end ]

        startDN = date2num(startDT)
        endDN = date2num(endDT)

        age = (endDN - startDN) * 24.

        hist, bins = np.histogram(age, self.ageBins, density=False)

        return hist

    def historic(self):
        """Calculate historical rates of longitude crossing"""
        log.debug("Processing historical tracks for age distribution")
        [i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe] = \
                interpolateTracks.parseTracks(self.configFile,
                                              self.historicTrackFile,
                                              self.historicFormat, 6.0)
                                              
        self.histAgeDist = self.calculateAge(i, y, m, d, h, mn)
        return True

    def synthetic(self):
        log.debug("Processing {0} synthetic events in {1}".\
                    format(self.synNumSimulations, self.synTrackPath))

        for n in xrange(self.synNumSimulations):
            trackFile = pjoin(self.synTrackPath, "tracks.%04d.csv"%(n))
            log.debug("Processing {0}".format(trackFile))
            try:
                i, y, m, d, h, mn, lon, lat, p, s, b, w, r, pe = \
                    loadTrackFile(self.configFile,
                                  trackFile,
                                  self.synFormat)
            except (TypeError, IOError, ValueError):
                log.critical("Cannot load synthetic track file: {0}".\
                                format(trackFile))
                return False
            else:
                self.synAgeDist[n, :] = self.calculateAge(i, y, m, d, h, mn)


        return True

    def synStats(self):
        if not hasattr(self, 'synAgeDist'):
            log.critical("Synthetic event sets have not been processed!")
            log.critical("Cannot calculate statistics")
            return False
        else:
            self.synAgeUpper = [ percentile(self.synAgeDist[:, i], per=self.upper)
                                    for i in xrange(len(self.ageBins) - 1) ]
            self.synAgeLower = [ percentile(self.synAgeDist[:, i], per=self.lower)
                                    for i in xrange(len(self.ageBins) - 1) ]
            self.synAgeMean = np.mean(self.synAgeDist, axis=0)

            return True

    def plotAgeDistribution(self):
        x = np.arange(len(self.ageBins)-1)
        pyplot.figure(1)
        ax1 = pyplot.subplot(111)
        ax1.plot(x, self.histAgeDist * self.synNumYears / self.historicNumYears,
                 color='r', lw=2)
                 
        ax1.plot(x, self.synAgeMean, color='k', lw=2)
        ax1.fill_between(x, self.synAgeLower, self.synAgeUpper, 
                         color='0.5', alpha=0.5)
                         
        pyplot.xticks(x[0:-1:4], self.ageBins[0:-1:4])
        pyplot.xlabel('Age (hours)')
        pyplot.yticks()
        pyplot.ylabel('Count')
        pyplot.grid(True)

        pyplot.savefig(pjoin(self.plotPath, 'tcAgeDistribution.png'))

def run(config_file):
    """Run the process"""
    

    args = dict(configFile        = config_file,
                historicTrackFile = cnfGetIniValue(config_file, 'Input', 'HistoricTrackFile'),
                historicFormat    = cnfGetIniValue(config_file, 'Input', 'HistoricSource'),
                historicNumYears  = cnfGetIniValue(config_file, 'Input', 'HistoricNumYears'),
                synTrackPath      = cnfGetIniValue(config_file, 'Input', 'SyntheticTrackPath'),
                synFormat         = cnfGetIniValue(config_file, 'Input', 'SyntheticSource'),
                synNumSimulations = cnfGetIniValue(config_file, 'Input', 'NumSyntheticEvents', 100),
                synNumYears       = cnfGetIniValue(config_file, 'Input', 'SyntheticNumYears', 1),
                MinLongitude      = cnfGetIniValue(config_file, 'Region', 'MinimumLongitude', 60.),
                MaxLongitude      = cnfGetIniValue(config_file, 'Region', 'MaximumLongitude', 180.),
                MinLatitude       = cnfGetIniValue(config_file, 'Region', 'MinimumLatitude', -40.),
                MaxLatitude       = cnfGetIniValue(config_file, 'Region', 'MaximumLatitude', 0.),
                GridSize          = cnfGetIniValue(config_file, 'Region', 'GridSize', 1.0),
                PlotPath          = cnfGetIniValue(config_file, 'Output', 'PlotPath', os.getcwd()),
                DataPath          = cnfGetIniValue(config_file, 'Output', 'DataPath', os.getcwd()),
                ColourMap         = cnfGetIniValue(config_file, 'Output', 'ColourMap','hot_r'))

    log.info("Processing track density information")
    tD = EvalTrackDensity(**args)
    tD.historic()
    tD.synthetic()
    tD.synStats()
    tD.plotTrackDensity()
    tD.plotPercentiles()
    tD.saveData()

    log.info("Processing pressure distribution information")
    pD = EvalPressureDistribution(**args)
    pD.historic()
    pD.synthetic()
    pD.synStats()
    pD.plotPressureMaps()
    pD.plotPressureDistributions()


    log.info("Processing longitude crossing information")
    lc = EvalLongitudeCrossings(**args)
    lc.historic()
    lc.synthetic()
    lc.synStats()
    lc.plotCrossingRates()
    lc.saveData()

    log.info("Processing age distribution information")
    aD = EvalAgeDistribution(**args)
    aD.historic()
    aD.synthetic()
    aD.synStats()
    aD.plotAgeDistribution()

def process_args(argv):
    """Main process to read the config settings and process the data"""

    verbose = False
    gConfigFile = flConfigFile()
    try:
        opts, args = getopt.getopt(argv, 'c:hv', ['config=', 'help', 'verbose'])
    except getopt.GetoptError:
        ShowSyntax(2)
    except IndexError:
        ShowSyntax(2)
    else:
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                ShowSyntax()
            elif opt in ("-c", "--config"):
                gConfigFile = arg
            elif opt in ("-v", "--verbose"):
                verbose = True

    
    log = flStartLog(cnfGetIniValue(gConfigFile, 'Logging', 'LogFile', flConfigFile('.log')),
                     cnfGetIniValue(gConfigFile, 'Logging', 'LogLevel', 'DEBUG'),
                     cnfGetIniValue(gConfigFile, 'Logging', 'Verbose', True),
                     cnfGetIniValue(gConfigFile, 'Logging', 'DateStamp', False)) 
           
    run(gConfigFile)

    log.info("Completed {0}".format(sys.argv[0]))

if __name__ == '__main__':
    process_args(sys.argv[1:])
