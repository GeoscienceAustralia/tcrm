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


 Title: plotTracks.py

 Author: Geoff Xu
 Email:

 CreationDate: 2006-02-16

 Description: Plots the cyclone tracks using the data in 'cycloneTracks'.

 Version:
 ModifiedBy: Geoff Xu
 ModifiedDate: 2006-02-17
 Modification:

 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate:2006-10-25
 Modification: Added descriptive headers and metadata

 ModifiedBy: Nariman Habili, nariman.habili@ga.gov.au
 ModifiedDate: 2006-11-21
 Modification: Conformance with new style guide.
               Using ptools instead of passing object tool in constructor

 Version 412
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-02-01
 Modification: Use Basemap toolkit to generate coastlines

 Version: 634
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2007-03-29
 Modification: Tracks are colourised based on estimated pressure.

 Version: $Rev$
 ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
 ModifiedDate: 2009-04-07
 Modification: Can now read in tracks from a shapefile


 SeeAlso: (related programs)
 Constraints:

 $Id: PlotTracks.py 634 2007-12-16 21:16:58Z carthur $
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
    
    
class PlotTracks(object):
    """PlotTracks: plots tropical cyclone tracks
    Parameters
    ----------
    cyclone_tracks : string (file name including path) or array
        Latitude and longitude data for cyclone tracks
    gridLimit: dictionary
        Contains bounds for the region of interest
    gridSpace: dictionary
        Contains spacing used for gridded assessment of data


    Members
    -------
    map: Basemap object

    cycloneTracks : string (file name including path) or array
        Latitude and longitude data for cyclone tracks
    gridLimit: dictionary
        Contains bounds for the region of interest
    gridSpace: dictionary
        Contains spacing used for gridded assessment of data


    Methods
    -------
    plotTracks()
        Plot cyclone tracks with Australian coast line displayed
    """

    def __init__(self, configFile, cycloneTracks, gridLimit, gridSpace ):
        """
        Initialize the data needed for the plots including
        gridLimit, gridSpace and cyclone track data.
        """
        self.configFile = configFile
        config = ConfigParser()
        config.read(configFile)

        if gridSpace['x'] <= 0. or gridSpace['y'] <= 0.:
            raise InvalidArguments, 'Invalid input on grid spacing: grid spacing cannot be negative or zero'

        ar = abs((gridLimit['yMax'] - gridLimit['yMin']) /
                 (gridLimit['xMax'] - gridLimit['xMin']))

        prj = config.get('PlotInterface', 'Projection')
        res = config.get('PlotInterface', 'Resolution')

        self.map = Basemap(projection=prj, resolution=res,
                            llcrnrlon=gridLimit['xMin'],
                            urcrnrlon=gridLimit['xMax'],
                            llcrnrlat=gridLimit['yMin'],
                            urcrnrlat=gridLimit['yMax'])

        self.gridSpace = gridSpace
        self.gridLimit = gridLimit
        self.meridian = (range(gridLimit['xMin'],
                               gridLimit['xMax']+gridSpace['x'],
                               gridSpace['x']))
        self.parallel = (range(gridLimit['yMin'],
                               gridLimit['yMax']+gridSpace['y'],
                               gridSpace['y']))

        if type(cycloneTracks) is str:
            if cycloneTracks.endswith('shp'):
                # We have a shapefile. Use Basemap's inbuilt features to handle this.
                self.map.readshapefile(cycloneTracks, 'tctrack',
                                       drawbounds=False)
            else:
                # Assume the track file is a CSV file, and the fields are in a specific order...
                self.cycloneTracks = files.flLoadFile(cycloneTracks, '%', ',')
        elif type(cycloneTracks) is np.ndarray:
            self.cycloneTracks = cycloneTracks
        elif type(cycloneTracks) is array:
            self.cycloneTracks = cycloneTracks
        else:
            raise TypeError, "Unknown type for cycloneTracks"



    def __doc__(self):
        """
        documentation on what this class does
        """
        return "Plot a set of cyclone tracks \
        The tracks can be in one of several formats: \
        "

    def plotTracks(self, tracks, titlestr=None, colours=False):
        """plotTracks([titlestr]):
        Plot cyclone tracks with coastline displayed. Title string is
        optional
        """
        counter = 0
        for t in tracks:
            counter += 1
            track_length = len(t.Longitude)
            self.map.plot(t.Longitude[0], t.Latitude[0],'k.')
            if track_length > 1:
                for i in range(track_length - 1):
                    seg_lon = [t.Longitude[i], t.Longitude[i + 1]]
                    seg_lat = [t.Latitude[i], t.Latitude[i + 1]]
                    if t.CentralPressure[i] == sys.maxint:
                        col = 'k'
                    elif between(t.CentralPressure[i], 985., sys.maxint):
                        col = 'b'
                    elif between(t.CentralPressure[i], 970., 985.):
                        col = 'g'
                    elif between(t.CentralPressure[i], 955., 970.):
                        col = 'y'
                    elif between(t.CentralPressure[i], 930., 955.):
                        col = 'r'
                    else:
                        col = 'm'
                    self.map.plot(seg_lon, seg_lat, col+'-')

        self.map.drawcoastlines()
        self.map.drawparallels(self.parallel, labels=[1,0,0,1], fontsize=9)
        self.map.drawmeridians(self.meridian, labels=[1,0,0,1], fontsize=9)
        self.map.fillcontinents()

        if titlestr is None:
            pyplot.title(str(counter) + ' Cyclone Tracks')
        else:
            pyplot.title(titlestr)
        print "Plotted %d cyclone tracks" % counter

    def plotTracksShp(self, titlestr=None, colours=False):
        """
        Plot tracks from a shapefile. The data is stored in the map
        instance using the readshapefile() function
        """
        for shapedict,shape in zip(self.map.tctrack_info,self.map.tctrack):
            print shapedict['Index']
            cp = shapedict['Pressure']
            xx, yy = zip(*shape)
            # show part of track where storm > Cat 4 as thick red.
            if colours:
                if cp > 985.:
                    self.map.plot(xx, yy, color='b')
                elif cp > 970. and cp <= 985.:
                    self.map.plot(xx, yy, color='g')
                elif cp > 955. and cp <= 970.:
                    self.map.plot(xx, yy, color='y')
                elif cp > 930. and cp <= 955.:
                    self.map.plot(xx, yy, linewidth=1.5, color='r')
                else:
                    self.map.plot(xx, yy, linewidth=1.5, color='m')
            else:
                self.map.plot(xx,yy,linewidth=1.5,color='k')
        self.map.drawcoastlines()
        self.map.drawparallels(self.parallel, labels=[1,0,0,1], fontsize=9)
        self.map.drawmeridians(self.meridian, labels=[1,0,0,1], fontsize=9)
        self.map.fillcontinents()
        if titlestr:
            pyplot.title(titlestr)

    def _separateTracks(self):
        """_separateTracks():
        Reads in the tracks of cyclones and returns a dictionary with
        separate entries for each individual cyclone.
        This makes plotting the intensity easier.
        """
        self.cyclone = {}
        n = 0

        for (index, lon, lat, pressure) in self.cycloneTracks:
            if index == 1 and n > 0:
                self.cyclone[n] = {"lon":lons, "lat":lats, "pressure":cps}

            if index == 1:
                n += 1
                lons = []
                lats = []
                cps = []

            lons.append(lon)
            lats.append(lat)
            cps.append(pressure)


