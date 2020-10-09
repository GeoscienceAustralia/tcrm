"""
:mod:`tracks` -- plot TC tracks
===================================

.. module:: tracks
    :synopsis: Plot TC tracks on a map.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

Note: This uses the Australian Tropical Cyclone Intensity Scale
      for colourizing the track segments, based on maximum
      10-minute wind speeds.

"""
import sys
import logging as log

import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

from PlotInterface.maps import MapFigure, saveFigure

def makeSegments(xx, yy):
    """
    Create a list of line segments from x,y coordinates, in the
    correct format for LineCollection.

    :param x: :class:`numpy.ndarray` of x-coordinates.
    :param y: :class:`numpy.ndarray` of y-coordinates.
    """

    points = np.array([xx, yy]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return segments



class TrackMapFigure(MapFigure):
    """
    Base class for plotting track maps.
    """

    def colorline(self, xdata, ydata, zdata=None, alpha=0.9,
                  colours=['0.75', '#0FABF6', '#0000FF', 
                           '#00FF00', '#FF8100', '#ff0000'],
                  intervals=[0, 17.5, 24.5, 32.5, 44.2, 55.5, 1000]):
        """
        Create and add line collections to an axes instance, using
        an optional magnitude value to colourize the line segments.

        :param x: :class:`numpy.ndarray` of x-coordinates for lines.
        :param y: :class:`numpy.ndarray` of y-coordinates for lines.
        :param z: (Optional) :class:`numpy.ndarray` of magnitudes to
                  colourize the line segments.
        :param float linewidth: Line width of the line segments to plot.
        :param float alpha: Transparency level of the line segments.
        :param list colours: List of HTML colour codes to use for colourizing
                             the line segments.
        :param list intervals: List of break points for colourizing the 
                               line segments.

        """

        if zdata is None:
            zdata = np.linspace(0.0, 1.0, len(xdata))

        if not hasattr(zdata, '__iter__'):
            zdata = np.array([zdata])

        zdata = np.asarray(zdata)

        segments = makeSegments(xdata, ydata)
        cmap = ListedColormap(colours)
        norm = BoundaryNorm(intervals, cmap.N)
        lc = LineCollection(segments, array=zdata, cmap=cmap,
                            norm=norm, alpha=alpha)

        labels = ['No data', 'Category 1', 'Category 2',
                  'Category 3', 'Category 4', 'Category 5']
        handles = []
        for c, l in zip(cmap.colors, labels):
            handles.append(Line2D([0], [0], color=c, label=l))

        ax = self.gca()
        ax.add_collection(lc)
        ax.legend(handles, labels, loc=2, frameon=True)

    def add(self, tracks, xgrid, ygrid, title, map_kwargs):
        self.subfigures.append((tracks, xgrid, ygrid, title, map_kwargs))

    def subplot(self, axes, subfigure):
        tracks, xgrid, ygrid, title, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        for track in tracks:
            mlon, mlat = mapobj(track.Longitude, track.Latitude)
            self.colorline(mlon, mlat, track.WindSpeed, alpha=0.75)
                           
        axes.set_title(title)
        #self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.fillContinents(mapobj)
        self.addMapScale(mapobj)

class SingleTrackMap(TrackMapFigure):
    """
    Plot single TC track on a map.
    """

    def plot(self, tracks, xgrid, ygrid, title, map_kwargs):
        self.add(tracks, xgrid, ygrid, title, map_kwargs)
        super(SingleTrackMap, self).plot()


def saveTrackMap(tracks, xgrid, ygrid, title, map_kwargs, filename):
    """
    Create a track map and save to file.

    :param tracks: collection of :class:`Track` objects
    :param xgrid: :class:`numpy.ndarray` of longitude points defining the
                  domain.
    :param ygrid: :class:`numpy.ndarray` of latitude points defining the
                  domain.
    :param str title: Title string for the plot.
    :param dict map_kwargs: Keyword args that will define the
                            :class:`GeoAxes` instance.
    :param str filename: Path to the file to save the image.

    """

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
    outputPath = pjoin(outputPath, 'plots', 'stats')
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
    CONFIG = sys.argv[1]
    main(CONFIG)
