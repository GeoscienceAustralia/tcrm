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
"""
#-------------------------------------------------------------------------------
#** Title: grid_interface.py
#**
#** Author: Geoff Xu
#** CreationDate: 2005-12-14
#**
#** Description: Wrapper class that allows a user to interact with the n by n grid
#**              related classes for representing historical cyclone data.
#**      This is a revision of the Python code written by Geoff Xu (2006).
#*
#* ModifiedBy: 2006-06-02 - Geoff Xu.
#*             2006-10-24 - C. Arthur - Added descriptive headers and metadata
#*
#* SeeAlso:
#* Constraints: See list of imports
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Imports:
#-------------------------------------------------------------------------------
from pylab import *
from scipy import gplt
from figure_number import *
from plot_origin_gui import PlotOriginGUI
from plot_tracks import PlotTracks
from plot_bearing import PlotBearing
from plot_speed import PlotSpeed
from plot_pressure import PlotPressure
from plot_cell import PlotCell
from grid_tool import *
from warning import *

#-------------------------------------------------------------------------------
# 'GridInterface' class:
#-------------------------------------------------------------------------------
class GridInterface:
    """
    Parameters
    ----------
    aus_coast_line : string (file name including path)
        latitude and longitude data for Australian coast line
    init_lon_lat : string (file name including path)
        latitude and longitude data for cyclone origins
    all_lon_lat : string (file name including path)
        latitude and longitude data for all cyclone records
    init_bearing : string (file name including path)
        initial bearings of cyclones
    bearing_no_init : string (file name including path)
        bearings of cyclones with no initial bearing
    all_bearing : string (file name including path)
        all bearings of cyclones
    init_speed : string (file name including path)
        initial speed of cyclones from init_speed file
    speed_no_init : string (file name including path)
        speed of cyclones with no initial speed from speed_no_init file
    all_speed : string (file name including path)
        all speed of cyclones from all_speed file
    init_pressure : string (file name including path)
        initial pressures of cyclones
    pressure_no_init : string (file name including path)
        pressures of cyclones with no initial pressure
    all_pressure : string (file name including path)
        all pressures of cyclones
    a : int
        longitude a of a X b grid in a cell
    b : int
        latitude b of a X b grid in a cell
    Nan : int
        integer representation for not a number
    parameter_type : int
        represents the type of cyclone parameter that is going to be plotted

    Members
    -------
    aus_coast_line : string (file name including path)
        latitude and longitude data for Australian coast line
    init_lon_lat : string (file name including path)
        latitude and longitude data for cyclone origins
    all_lon_lat : string (file name including path)
        latitude and longitude data for all cyclone records
    init_bearing : string (file name including path)
        initial bearings of cyclones
    bearing_no_init : string (file name including path)
        bearings of cyclones with no initial bearing
    all_bearing : string (file name including path)
        all bearings of cyclones
    init_speed : string (file name including path)
        initial speeds of cyclones
    speed_no_init : string (file name including path)
        speeds of cyclones with no initial speed
    all_speed : string (file name including path)
        all speeds of cyclones
    init_pressure : string (file name including path)
        initial pressures of cyclones
    pressure_no_init : string (file name including path)
        pressures of cyclones with no initial pressure
    all_pressure : string (file name including path)
        all pressures of cyclones
    a : int
        longitude a of a X b grid in a cell
    b : int
        latitude b of a X b grid in a cell
    parameter_type : int
        represents the type of cyclone parameter that is going to be plotted
    tool : GridTool object
        a class instance of GridTool

    Methods
    -------
    plotOrigin()
        generate plots relating to the origins of cyclones
    plotBearing()
        generate plots relating to the bearings of cyclones
    plotSpeed()
        generate plots relating to the speeds of cyclones
    plotPressure()
        generate plots relating to the pressures of cyclones
    plotCell(cellNo)
        generate cyclone parameter plots relating to a particular cell

    """
    def __init__(self, aus_coast_line, init_lon_lat, all_lon_lat, \
                init_bearing, bearing_no_init, all_bearing,\
                init_speed, speed_no_init, all_speed,\
                init_pressure, pressure_no_init, all_pressure,\
                a=5, b=5, Nan=-1, parameter_type=1):
        """
        initialize the data needed for the interface including
        australia coast line data,
        longitudes & latitudes of cyclone origins and all cyclone data,
        initial cyclone bearings,
        cyclone bearings excluding initial ones at origin,
        all cyclone bearings,
        and a & b for a X b grid in a cell
        """
        if a <= 0 or b <= 0:
            raise InvalidArguments, 'invalid input on a or b: a or b cannot be negative or zero'
        if parameter_type < 1 or parameter_type > 2:
            raise InvalidArguments, 'invalid input on parameter_type: type can only be 1 or 2'
        self.aus_coast_line = aus_coast_line
        self.init_lon_lat = init_lon_lat
        self.all_lon_lat = all_lon_lat
        self.init_bearing = init_bearing
        self.bearing_no_init = bearing_no_init
        self.all_bearing = all_bearing
        self.init_speed = init_speed
        self.speed_no_init = speed_no_init
        self.all_speed = all_speed
        self.init_pressure = init_pressure
        self.pressure_no_init = pressure_no_init
        self.all_pressure = all_pressure
        self.a = a
        self.b = b
        self.parameter_type = parameter_type
        self.tool = GridTool(Nan)

    def __doc__(self):
        """
        documentation on what this class does
        """
        return 'A wrapper class that allows user to interact with all \
               the n by n grid related classes for representing historical \
               data of cyclones'

    def plotOrigin(self, fig, oricounts, lon_bin, lat_bin):
        """
        generate plots relating to the origins of cyclones
        """
        print 'generating plots relating to the origins of cyclones'
        print 'reading data from ' + self.aus_coast_line
        if self.parameter_type == 1:
            print 'reading data from ' + self.init_lon_lat
            print 'reading data from ' + self.init_bearing
            print 'reading data from ' + self.init_speed
            print 'reading data from ' + self.init_pressure
            p = PlotOriginGUI(self.aus_coast_line, self.init_lon_lat, self.a,
                              self.b, self.init_bearing, self.init_speed,
                              self.init_pressure, self.tool)
        else:
            print 'reading data from ' + self.all_lon_lat
            print 'reading data from ' + self.all_bearing
            print 'reading data from ' + self.all_speed
            print 'reading data from ' + self.all_pressure
            p = PlotOriginGUI(self.aus_coast_line, self.all_lon_lat, self.a,
                              self.b, self.all_bearing, self.all_speed,
                              self.all_pressure, self.tool)
        print 'outputting data into ' + oricounts
        print 'outputting data into ' + lon_bin
        print 'outputting data into ' + lat_bin
        p.generatePlots(fig, oricounts, lon_bin, lat_bin)

    def plotTracks(self, fig, cyclone_tracks):
        """
        generate plots relating to cyclone tracks
        """
        print 'generating plots relating to cyclone tracks'
        print 'reading data from ' + self.aus_coast_line
        print 'reading data from ' + cyclone_tracks
        p = PlotTracks(self.aus_coast_line, cyclone_tracks, self.a, self.b,
                       self.tool)
        figure(fig.getNum())
        p.plotTracks()
        fig.next()

    def plotBearing(self, fig):
        """
        generate plots relating to the bearings of cyclones
        """
        print 'generating plots relating to the bearings of cyclones'
        print 'reading data from ' + self.init_bearing
        print 'reading data from ' + self.bearing_no_init
        print 'reading data from ' + self.all_bearing
        if self.parameter_type == 1:
            print 'reading data from ' + self.init_lon_lat
            p = PlotBearing(self.init_lon_lat, self.init_bearing, self.a,
                            self.b, self.tool)
        else:
            print 'reading data from ' + self.all_lon_lat
            p = PlotBearing(self.all_lon_lat, self.all_bearing, self.a,
                            self.b, self.tool)
        figure(fig.getNum())
        p.plotHistGrid()
        fig.next()
        figure(fig.getNum())
        p.plotHistGridTopHalf()
        fig.next()
        figure(fig.getNum())
        p.plotHistParameter(self.init_bearing, self.bearing_no_init,
                            self.all_bearing)
        fig.next()

    def plotSpeed(self, fig):
        """
        generate plots relating to the speeds of cyclones
        """
        print 'generating plots relating to the speeds of cyclones'
        print 'reading data from ' + self.init_speed
        print 'reading data from ' + self.speed_no_init
        print 'reading data from ' + self.all_speed
        if self.parameter_type == 1:
            print 'reading data from ' + self.init_lon_lat
            p = PlotSpeed(self.init_lon_lat, self.init_speed, self.a,
                          self.b, self.tool)
        else:
            print 'reading data from ' + self.all_lon_lat
            p = PlotSpeed(self.all_lon_lat, self.all_speed, self.a, self.b,
                          self.tool)
        figure(fig.getNum())
        p.plotHistGrid()
        fig.next()
        figure(fig.getNum())
        p.plotHistGridTopHalf()
        fig.next()
        figure(fig.getNum())
        p.plotHistParameter(self.init_speed, self.speed_no_init,
                            self.all_speed)
        fig.next()

    def plotPressure(self, fig):
        """
        generate plots relating to the pressures of cyclones
        """
        print 'generating plots relating to the pressures of cyclones'
        print 'reading data from ' + self.init_pressure
        print 'reading data from ' + self.pressure_no_init
        print 'reading data from ' + self.all_pressure
        if self.parameter_type == 1:
            print 'reading data from ' + self.init_lon_lat
            p = PlotPressure(self.init_lon_lat, self.init_pressure, self.a,
                             self.b, self.tool)
        else:
            print 'reading data from ' + self.all_lon_lat
            p = PlotPressure(self.all_lon_lat, self.all_pressure, self.a,
                             self.b, self.tool)
        figure(fig.getNum())
        p.plotHistGrid()
        fig.next()
        figure(fig.getNum())
        p.plotHistGridTopHalf()
        fig.next()
        figure(fig.getNum())
        p.plotHistParameter(self.init_pressure, self.pressure_no_init,
                            self.all_pressure)
        fig.next()

    def plotCell(self, fig, cellNo):
        """
        generate cyclone parameter plots relating to a particular cell
        """
        print 'generating cyclone parameter plots relating to a particular cell'
        if self.parameter_type == 1:
            print 'reading data from ' + self.init_lon_lat
            print 'reading data from ' + self.init_bearing
            print 'reading data from ' + self.init_speed
            print 'reading data from ' + self.init_pressure
            p = PlotCell(self.init_lon_lat, self.init_bearing,
                         self.init_speed, self.init_pressure, self.a, self.b,
                         self.tool)
        else:
            print 'reading data from ' + self.all_lon_lat
            print 'reading data from ' + self.all_bearing
            print 'reading data from ' + self.all_speed
            print 'reading data from ' + self.all_pressure
            p = PlotCell(self.all_lon_lat, self.all_bearing, self.all_speed,
                         self.all_pressure, self.a, self.b, self.tool)

        p.generateParameter(cellNo)
        figure(fig.getNum())
        subplot(1, 3, 1)
        p.plotBearing()
        subplot(1, 3, 2)
        p.plotSpeed()
        subplot(1, 3, 3)
        p.plotPressure()
        fig.next()
