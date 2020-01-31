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
#** Title: PlotOriginGUI.py
#**
#** Author: Geoff Xu
#** Email:
#**
#** CreationDate: 2006-01-19
#**
#** Description: Define the class for PlotOriginGUI. This is a Python
#**              implementation
#**              implementation of a graphical user interface for the
#**              PlotOrigin class.  This is a revision of the Python
#**              code written by Geoff Xu (2006).
#**
#** ModifiedBy: 2006-02-22 - G. Xu
#**             2006-10-24 - C. Arthur - Added descriptive headers and metadata
#**     2006-11-21 - N. Habili - Conformance with style guide
#*
#*  SeeAlso: (related programs)
#*  Constraints:
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Imports:
#-------------------------------------------------------------------------------
import PlotOrigin
import PlotCell
import pylab

#-------------------------------------------------------------------------------
# 'PlotOriginGUI' class:
#-------------------------------------------------------------------------------
class PlotOriginGUI:
    """
    Parameters
    ----------
    aus_coast_line : string (file name including path)
        latitude and longitude data for Australian coast line
    lon_lat : string (file name including path)
        latitude and longitude data for cyclone origins
    a : int
        longitude a of a X b grid in a cell
    b : int
        latitude b of a X b grid in a cell
    init_bearing : string (file name including path)
        initial cyclone bearings
    init_speed : string (file name including path)
        initial cyclone speeds
    init_pressure : string (file name including path)
        initial cyclone pressures

    Members
    -------
    p : PlotOrigin object
        the object instance of PlotOrigin class
    pc : PlotCell object
        the object instance of PlotCell class

    External Methods
    ----------------
    generatePlots(fig,oricounts,lon_bin,lat_bin)
        a method that generates cyclone origin plots

    Internal Methods(Slots)
    ----------------
    on_click(event)
        plot the individual cell plots when mouse is clicked

    Internal Methods
    ----------------
    plotCell(cellNo)
        generate cyclone parameter plots relating to a particular cell
    """
    def __init__(self, aus_coast_line, lon_lat, origin_counts, longitude_bin,
                 latitude_bin, a, b, init_bearing, init_speed, init_pressure):
        """
        initialize the data needed including
        longitudes & latitudes of cyclone origins(lon_lat),
        longitude & latitude dimensions for the grid(a & b),
        initial cyclone bearings,
        initial cyclone speeds,
        initial cyclone pressures,
        as well as plotOrigin object
        """
        self.p = PlotOrigin.PlotOrigin(aus_coast_line, lon_lat, origin_counts,
                                       longitude_bin, latitude_bin, a, b)
        self.pc = PlotCell.PlotCell(lon_lat, init_bearing, init_speed,
                                    init_pressure, a, b)

        self.cellFigNo = 0

    def generatePlots(self, fig):
        """
        a method that generates cyclone origin plots
        """
        pylab.figure(fig.getNum())
        self.p.plotOrigin()
        next(fig)
        pylab.connect('button_press_event', self.on_click)
        pylab.figure(fig.getNum())
        self.p.plotOriginCounts()
        next(fig)
        pylab.connect('button_press_event', self.on_click)
        pylab.figure(fig.getNum())
        self.p.plotOriginCells()
        next(fig)
        pylab.connect('button_press_event', self.on_click)
        pylab.figure(fig.getNum())
        self.p.plotOriginCellsHist()
        next(fig)
        self.cellFigNo = fig.getNum()
        next(fig)
        self.p.plot3DOriginCounts()

    def on_click(self, event):
        """
        plot the individual cell plots when mouse is clicked
        """
        if event.button == 1 and event.inaxes is not None:
            cellNo = self.p.getCellNo(event.xdata, event.ydata)
            self.plotCell(cellNo)

    def plotCell(self, cellNo):
        """
        generate cyclone parameter plots relating to a particular cell
        """
        self.pc.generateParameter(cellNo)
        pylab.figure(self.cellFigNo)
        pylab.subplot(1, 3, 1)
        self.pc.plotBearing()
        pylab.hold(False)
        pylab.subplot(1, 3, 2)
        self.pc.plotSpeed()
        pylab.hold(False)
        pylab.subplot(1, 3, 3)
        self.pc.plotPressure()
        pylab.hold(False)

if __name__ == "__main__":
    import utils.figure_number as FigureNumber
    import utils.my_tool as mtools

    cfgMain = mtools.loadConfig(open('main.ini'))
    cfgFiles = mtools.loadConfig(open(cfgMain['Files.configFile']))

    aus_coast_line = cfgFiles['Input.aus_coast_line']
    init_lon_lat = cfgFiles['Output.init_lon_lat']
    longitude_bin = cfgFiles['Output.longitude_bin']
    latitude_bin = cfgFiles['Output.latitude_bin']
    origin_counts = cfgFiles['Output.origin_counts']
    init_bearing = cfgFiles['Output.init_bearing']
    init_speed = cfgFiles['Output.init_speed']
    init_pressure = cfgFiles['Output.init_pressure']

    fig = FigureNumber.FigureNumber()

    lonDim = 5
    latDim = 5

    plotOrgGUI = PlotOriginGUI(aus_coast_line, init_lon_lat, origin_counts,
                               longitude_bin, latitude_bin, lonDim, latDim,
                               init_bearing, init_speed, init_pressure)
    plotOrgGUI.generatePlots(fig)
    #pylab.show()
