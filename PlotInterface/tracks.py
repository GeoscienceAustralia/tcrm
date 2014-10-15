"""
:mod:`tracks` -- plot TC tracks
===================================

.. module:: tracks
    :synopsis: Plot TC tracks on a map.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""
import os
import sys
import logging

import numpy as np

import Utilities.shptools as shptools
from Utilities.config import ConfigParser
from Utilities.files import flStartLog
from Utilities.stats import between

from maps import MapFigure, saveFigure

class TrackMapFigure(MapFigure):

    def add(self, tracks, xgrid, ygrid, title, map_kwargs):
        self.subfigures.append((tracks, xgrid, ygrid, title, map_kwargs))

    def subplot(self, axes, subfigure):
        tracks, xgrid, ygrid, title, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        for t in tracks:
            mlon, mlat = mapobj(t.Longitude, t.Latitude)
            mapobj.plot(mlon, mlat, 'k-', linewidth=1.5)
        axes.set_title(title)
        self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.fillContinents(mapobj)
        self.addMapScale(mapobj)

class SingleTrackMap(TrackMapFigure):

    def plot(self, tracks, xgrid, ygrid, title, map_kwargs):
        self.add(tracks, xgrid, ygrid, title, map_kwargs)
        super(SingleTrackMap, self).plot()


def saveTrackMap(tracks, xgrid, ygrid, title, map_kwargs, filename):
    fig = SingleTrackMap()
    fig.plot(tracks, xgrid, ygrid, title, map_kwargs)
    saveFigure(fig, filename)

def main(configFile):
    from Utilities.loadData import loadTrackFile
    from Utilities.config import ConfigParser
    from os.path import join as pjoin, normpath, dirname
    baseDir = normpath(pjoin(dirname(__file__), '..'))
    inputPath = pjoin(baseDir, 'input')
    config = ConfigParser()
    config.read(configFile)

    inputFile = config.get('DataProcess', 'InputFile')
    source = config.get('DataProcess', 'Source')

    gridLimit = config.geteval('Region', 'gridLimit')

    xx = np.arange(gridLimit['xMin'], gridLimit['xMax'] + .1, 0.1)
    yy = np.arange(gridLimit['yMin'], gridLimit['yMax'] + .1, 0.1)

    xgrid, ygrid = np.meshgrid(xx, yy)

    if len(dirname(inputFile)) == 0:
        inputFile = pjoin(inputPath, inputFile)

    try:
        tracks = loadTrackFile(configFile, inputFile, source)
    except (TypeError, IOError, ValueError):
        log.critical("Cannot load historical track file: {0}".format(inputFile))
        raise

    title = source
    outputPath = config.get('Output', 'Path')
    outputPath = pjoin(outputPath, 'plots','stats')
    outputFile = pjoin(outputPath, 'tctracks.png')

    map_kwargs = dict(llcrnrlon=xgrid.min(),
                      llcrnrlat=ygrid.min(),
                      urcrnrlon=xgrid.max(),
                      urcrnrlat=ygrid.max(),
                      projection='merc',
                      resolution='i')

    figure = TrackMapFigure()
    figure.add(tracks, xgrid, ygrid, title, map_kwargs)
    figure.plot()
    saveFigure(figure, outputFile)

if __name__ == "__main__":
    configFile = sys.argv[1]
    main(configFile)
