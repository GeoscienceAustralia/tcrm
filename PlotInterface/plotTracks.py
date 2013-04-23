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
import os, sys, pdb, logging

import pylab, numpy
from mpl_toolkits.basemap import Basemap

import Utilities.shptools as shptools
import Utilities.config as config
import Utilities.files as files
import Utilities.grid as grid
import Utilities.my_tool as myutils

class PlotTracks:
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
        
        if gridSpace['x'] <= 0. or gridSpace['y'] <= 0.:
            raise InvalidArguments, 'Invalid input on grid spacing: grid spacing cannot be negative or zero'

        ar = abs((gridLimit['yMax'] - gridLimit['yMin']) /
                 (gridLimit['xMax'] - gridLimit['xMin']))

        prj = config.cnfGetIniValue(self.configFile, 'PlotInterface',
                                    'Projection', 'cyl')
        res = config.cnfGetIniValue(self.configFile, 'PlotInterface',
                                    'Resolution', 'l')

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
        elif type(cycloneTracks) is numpy.ndarray:
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

    def plotTracks(self, titlestr=None, colours=False):
        """plotTracks([titlestr]):
        Plot cyclone tracks with coastline displayed. Title string is
        optional
        """
        self._separateTracks()
        counter = 0
        for n in self.cyclone:
            counter += 1
            self.map.plot([self.cyclone[n]['lon'][0]],
                          [self.cyclone[n]['lat'][0]], 'k.')
            for i in range(1, len(self.cyclone[n]['lon'])):
                if colours:
                    if self.cyclone[n]['pressure'][i-1] == sys.maxint:
                        col = 'k'
                    elif ((self.cyclone[n]['pressure'][i-1] >= 985.) and
                          (self.cyclone[n]['pressure'][i-1] < sys.maxint)):
                        col = 'b'
                    elif ((self.cyclone[n]['pressure'][i-1] >= 970.) and
                          (self.cyclone[n]['pressure'][i-1] < 985.)):
                        col = 'g'
                    elif ((self.cyclone[n]['pressure'][i-1] >= 955.) and
                          (self.cyclone[n]['pressure'][i-1] < 970.)):
                        col = 'y'
                    elif ((self.cyclone[n]['pressure'][i-1] >= 930.) and
                          (self.cyclone[n]['pressure'][i-1] < 955.)):
                        col = 'r'
                    else:
                        col = 'm'
                else:
                    col = 'k'
                self.map.plot([self.cyclone[n]['lon'][i-1], self.cyclone[n]['lon'][i]],
                              [self.cyclone[n]['lat'][i-1], self.cyclone[n]['lat'][i]],
                              col+'-')
        self.map.drawcoastlines()
        self.map.drawparallels(self.parallel, labels=[1,0,0,1], fontsize=9)
        self.map.drawmeridians(self.meridian, labels=[1,0,0,1], fontsize=9)
        self.map.fillcontinents()

        if titlestr is None:
            pylab.title(str(counter) + ' Cyclone Tracks')
        else:
            pylab.title(titlestr)
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
            pylab.title(titlestr)

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



if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        configFile = 'plotTracks.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python plotTracks.py {config filename}.ini"
            raise IOError, error_msg
    # If specified config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError, error_msg
        
    #cycloneTracks = cfgFiles['Output.cyclone_tracks']
    #cycloneData = myutils.textread('DataProcess\\dataprocess.ini')
    #indicator = numpy.array(cycloneData['indicator'], int)
    #lat = numpy.array(cycloneData['lat'], float)
    #lon = numpy.array(cycloneData['lon'], float)
    #pressure = numpy.array(cycloneData['pressure'], float)

    #cycloneTracks = numpy.concatenate(([indicator], [lon], [lat], [pressure]), axis=0)
    #cycloneTracks = numpy.transpose(cycloneTracks)

    #data = files.flLoadFile("C:\\WorkSpace\\atcram\\data\\output\\post1960\\syn_tracks.0001.txt", '%', ',')
    #ind = numpy.ones(numpy.size(data[:,0]))
    #ind[1:] = numpy.diff(data[:,0])
    #cycloneTracks = numpy.concatenate(([ind], [data[:,2]], [data[:,3]], [data[:,6]]), axis=0)
    #cycloneTracks = numpy.transpose(cycloneTracks)
    cycloneTracks = 'C:/WorkSpace/data/PAGASA/tracks/tracks_0001.shp'
    gridLimit = {'xMin':110, 'xMax':140, 'yMin':5, 'yMax':30} #eval(config.cnfGetIniValue(configFile, 'Map', 'gridLimit') #{'xMin':90, 'xMax':180, 'yMin':-40, 'yMax':0}
    gridSpace = {'x':1, 'y':1} #eval(config.cnfGetIniValue(configFile, 'Map', 'gridSpace') #{'x':5, 'y':5}
    pT = PlotTracks(configFile, cycloneTracks, gridLimit, gridSpace)
    pT.plotTracksShp()
    pylab.savefig('C:/WorkSpace/data/PAGASA/tracks/tracks_0001.png')
    pylab.show()
