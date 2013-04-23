#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Title: trackAnalysis.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 11/13/07 11:09:AM
Description: Read in a track file and plot a series of diagnostics to
allow comparison of synthetic and historic event sets

Version :$Rev: 810 $

$Id: trackAnalysis.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

import my_tool as myutils
import metutils
import map as maputils
import stats
import numpy
import pylab
from files import flConfigFile, flLoadFile
from config import cnfGetIniValue, gConfigFile

global gConfigFile

try:
    gConfigFile = sys.argv[1]
except IndexError:
    gConfigFile = flConfigFile()


inputFile = cnfGetIniValue(gConfigFile, 'Input', 'File')
data = flLoadFile(inputFile, ",")
cNum = data[:, 0]
cAge = data[:, 1]
cLon = data[:, 2]
cLat = data[:, 3]
cVfm = data[:, 4]
cTheta = data[:, 5]
cPrs = data[:, 6]
cPe = data[:, 7]
cRm = data[:, 8]

cInd = ones(len(cNum))
ind_ = diff(cNum)
cInd[1:] = ind_

tempTracks = transpose([cInd, cNum, cAge, cLon, cLat, cVfm, cTheta, cPrs,
                        cPe, cRm])

cTracks = separateTracks(tempTracks)


def separateTracks(cTrack):
    """separateTracks(cTrack):
    Reads in the tracks of cyclones and returns a dictionary with
    separate entries for each individual cyclone.
    This makes plotting the intensity easier.
    """
    tracks = {}
    n = 0

    for (index, num, age, lon, lat, vFm, thetaFm, pCentre,
         pEnv, rMax) in cTrack:
        if int(index) == 1 and n > 0:
            tracks[n] = {"Number":numbers,
                         "Age":ages,
                         "Lon":lons,
                         "Lat":lats,
                         "vFm":vFms,
                         "thetaFm":thetaFms,
                         "pCentre":pCentres,
                         "pEnv":pEnvs,
                         "rMax":rMaxs}

        if int(index) == 1:
            n += 1
            numbers = ages = lons = lats = vFms = thetaFms = pCentres = \
            pEnvs = rMaxs = []

        ages.append(float(age))
        numbers.append(int(num))
        lons.append(float(lon))
        lats.append(float(lat))
        vFms.append(float(vFm))
        thetaFms.append(float(thetaFm))
        pCentres.append(float(pCentre))
        pEnvs.append(float(pEnv))
        rMaxs.append(float(rMax))
        return tracks

"""def histSize(data):
    bins = numpy.arange(0,200,5)
    n, bins,patches = pylab.hist(data,bins)
    return n,bins
def histPressure(data):
    bins = numpy.arange(800,10,1010)
    n,bins,patches = pylab.hist(data,bins)
    return n,bins
def histSpeed(data):
    bins = numpy.arange(0,50,5)
    n,bins,patches = pylab.hist(data,bins)
    return n,bins
def histLat(data,gridLimit,gridSpace):
    bins = numpy.arange(gridLimit['yMin'],gridLimit['yMax']+gridSpace['y'],gridSpace['y'])
    n,bins,patches = pylab.hist(data,bins)
    return n,bins
def histLon(data):
    bins = numpy.arange(gridLimit['xMin'],gridLimit['xMax']+gridSpace['y'],gridSpace['x'])
    n,bins,patches = pylab.hist(data,bins)
    return n,bins
def histBearing(histogram):
    n,bins,patches = pylab.hist(histogram.data,histogram.bins)
    return n


def calcMean(data, ax=None):
    mean = numpy.mean(data,ax)
    std = numpy.std(data,ax)
    var = numpy.var(data,ax)
    return mean, std, var

def histogram(h):
    n,b,p = pylab.hist(h.data, h.bins)
    return n, b, p
"""
class histogram:
    """histogram:

    Description:

    Parameters:
    Members:
    Methods:
    Internal methods:
    """
    def __init__(self, bins):
        self.bins = bins
        self.data = []
        self.mean = []
        self.stdev = []
    def freq(self):
        n, b, p = pylab.hist(self.data, self.bins)
        return n

class Analysis:
    """Analysis:

    Description:

    Parameters:
    Members:
    Methods:
    Internal methods:
    """

    def __init__(self, cfgFile):
        self.Plot = PlotData()
        self.figNum = 0
        self.cfg = myutils.loadConfig(open(cfgFile))
        self.gridLimit = eval(self.cfg['gridLimit'])
        self.gridSpace = eval(self.cfg['gridSpace'])
        self.cols = numpy.array(self.cfg['Columns'].split(","))
        self.dataPath = myutils.checkPath(self.cfg['Data.Path'])
        # Output path should not be the same as the data path!
        self.outputPath = myutils.checkPath(self.cfg['Output.Path'])

        self.fileList = []
        dirList = os.listdir(self.dataPath)
        # At this time, assume only the track files to be analysed reside
        # in the input data path
        self.fileList = dirList  # [file for file in dirList if os.path.isfile(file)]
        self.data = []
        m = len(self.fileList)
        for i in range(m):
            self.data.append(loadData(self.dataPath+self.fileList[i]))

        bearing = histogram(numpy.arange(-22.5, 382.5, 45.))

        for i in range(m):
            n,bins,patches = pylab.hist(data[i][:, 4])
            h.speed.append(n)
            n,bins = histBear(data[i][:, 5])
            h.bear.append(n)
            n,bins = histPressure(data[i][:, 6])
            h.pressure.append(n)

        for param in h:


    def loadData(self, filename):
        """loadData(fileName):
        Load data into an array
        """
        data = myutils.load(filename, delimiter=",")
        return data




class PlotData:
    """PlotData:

    Description:

    Parameters:
    Members:
    Methods:
    Internal methods:
    """

    def __init__(self):

    def plotHist(bins, data):
        """plotHist(bins, data):
        Plot a histogram of data given a set of bins
        """
        self.figNum += 1
        pylab.figure(self.figNum)
        pylab.hist(data, bins)
        pylab.show()

    def plotErrBar(data, xerr=None, yerr=None):


