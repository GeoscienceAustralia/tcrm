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


 Title:plotCell.py
 Author: Geoff Xu
 Email:

 CreationDate: 2005-12-23
 ModifiedBy: 2006-04-25 - Geoff Xu
             2006-10-24 - C. Arthur - Added descriptive headers and metadata
     2006-11-20 - N. Habili - Removed usage of tool class
                  _validCellNo and _findLonLat are now private methods.
                  Uses removeNum instead of removeNan
                  Removed for loops in _validCellNo and _findLonLat
                  Update to numpy
 Description: Defines the class for plotting cyclone parameters for a
              particular cell including bearing, speed and pressure.
              This is a revision of the Python code written by Geoff Xu
              (2006).
 SeeAlso: plot_parameters.py
 Constraints: This is a child of PlotParameters

 Version: $Rev$

 $Id$
"""

#-------------------------------------------------------------------------------
# Imports:
#-------------------------------------------------------------------------------
import os
import sys
import pylab
from matplotlib import cm
import numpy as np

import Utilities.stats as stats
import Utilities.files as files
import Utilities.my_tool as myutils

#-------------------------------------------------------------------------------
# 'PlotCell' class:
#-------------------------------------------------------------------------------
class PlotCell:
    """
    Parameters
    ----------
    lon_lat : string (file name including path)
        latitude and longitude data for cyclone origins
    bearing_list : string (file name including path)
        list of cyclone bearings
    speed_list : string (file name including path)
        list of cyclone speeds
    pressure_list : string (file name including path)
        list of cyclone pressures
    gridLimit : dictionary of limits of the grid
    gridSpace : dictionary of grid spacing
    cellNum : cell number to plot

    Members
    -------
    bearing_list : 1D array of float
        initial bearing of cyclones from bearing_list file
    speed_list : 1D array of float
        initial speed of cyclones from speed_list file
    pressure_list : 1D array of float
        initial pressure of cyclones from pressure_list file
    lonlat : 2D array of float
        latitude and longitude data for cyclone origins
    gridLimit : dictionary of limits of the grid
    gridSpace : dictionary of grid spacing
    indij : 1D list of int
        a list of index of cyclone occurences in the cell
    cellNum : int
        a particular cell number


    Methods
    -------
    generateParameter(cellNum)
        generates the cyclone parameter data for a particular cell
    plotBearing()
        plots the cyclone bearing historgram of a particular cell
    plotSpeed()
        plots the cyclone speed historgram of a particular cell
    plotPressure()
        plots the cyclone pressure historgram of a particular cell

    Internal Methods
    ----------------


    """

    def __init__(self, lon_lat, bearing_list, speed_list, pressure_list,
                 bearing_rate_list, speed_rate_list, pressure_rate_list,
                 cellNum, gridLimit, gridSpace):
        """
        initialize the data needed for the plots including
        longitudes & latitudes of cyclone origins,
        cyclone bearings,
        cyclone speeds,
        cyclone pressures,
        as well as calculating and storing the
        cyclone index counts for a specific inputted cell number
        """
        self.gridLimit = gridLimit
        self.gridSpace = gridSpace
        self.lonlat = files.flLoadFile(lon_lat, '%', ',')

        self.bearing_list = files.flLoadFile(bearing_list)
        self.speed_list = files.flLoadFile(speed_list)
        self.pressure_list = files.flLoadFile(pressure_list)

        self.bearingRate_list = files.flLoadFile(bearing_rate_list)
        self.speedRate_list = files.flLoadFile(speed_rate_list)
        self.pressureRate_list = files.flLoadFile(pressure_rate_list)
        self.cellNum = cellNum
        self.extractParameters(cellNum)

    def extractParameters(self, cellNum):
        """
        Extracts cyclone parameters for a given cell
        """
        if not stats.validCellNum(cellNum, self.gridLimit, self.gridSpace):
            raise ValueError('Invalid input on cellNum: cell number is out of range')
        lon = self.lonlat[:,0]
        lat = self.lonlat[:,1]
        wLon, nLat = stats.getCellLonLat(cellNum, self.gridLimit,
                                         self.gridSpace)
        eLon = wLon + self.gridSpace['x']
        sLat = nLat - self.gridSpace['y']

        indij = np.where(((lat >= sLat) & (lat < nLat)) & \
                         (lon >= wLon) & (lon < eLon))
        self.bearing = myutils.removeNum(self.bearing_list[indij])
        self.speed = myutils.removeNum(self.speed_list[indij])
        self.pressure = myutils.removeNum(self.pressure_list[indij])
        self.bearingRate = myutils.removeNum(self.bearingRate_list[indij])
        self.speedRate = myutils.removeNum(self.speedRate_list[indij])
        self.pressureRate = myutils.removeNum(self.pressureRate_list[indij])

    def plotBearing(self):
        """
        Plot historgram of TC bearing as a polar diagram
        """
        fig = pylab.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
        theta = np.arange(0.0, 2*np.pi, np.pi/6.)

        bearing = (np.pi/2. - self.bearing * np.pi/180.)
        ii = np.where(bearing < 0)
        bearing[ii] += 2 * np.pi
        rads, bins = np.histogram(bearing, bins=12, range=(0., 2.*np.pi))

        bars = ax.bar(theta, rads, width=np.pi/6., bottom=0.0)

        [b.set_facecolor(cm.RdYlBu(float(r)/float(max(rads))))
         for r, b in zip(rads, bars)]

        pylab.title('Cyclone Bearing', fontsize=10)



    def plotBearingRate(self):
        """
        Plot histogram of rate of change of bearing
        """
        pylab.figure()
        pylab.hist(self.bearingRate, list(range(-90, 91, 10)))
        pylab.title('Bearing rate', fontsize=10)
        pylab.ylabel('Counts', fontsize=8)
        pylab.xlim(-90, 90)
        pylab.xlabel('Change of bearing (degrees/hour)', fontsize=8)

    def plotSpeed(self):
        """
        Plot histogram of TC speeds
        """
        pylab.figure()
        pylab.hist(self.speed, list(range(0, 101, 5)))
        #pylab.title('Cell Number ' + str(self.cellNum) + '\nTotal Cyclone Occurence ' + str(len(self.indij)) + '\nCyclone Speeds',fontsize=10)
        pylab.ylabel('Counts', fontsize=8)
        pylab.xlim(0, 100)
        pylab.xlabel('Cyclone speed (km/h)', fontsize=8)
        pylab.grid()

    def plotSpeedRate(self):
        """
        Plot histogram of rate of change of speed
        """
        pylab.figure()
        pylab.hist(self.speedRate, list(range(-20, 21, 1)))
        #pylab.title('Cell Number ' + str(self.cellNum) + '\nTotal Cyclone Occurence ' + str(len(self.indij)) + '\nSpeed rate',fontsize=10)
        pylab.ylabel('Counts', fontsize=8)
        pylab.xlim(-20, 20)
        pylab.xlabel('Change of speed (km/h/h)', fontsize=8)
        pylab.grid()

    def plotPressure(self):
        """
        Plot histogram of TC central pressure
        """
        pylab.figure()
        pylab.hist(self.pressure, list(range(900, 1011, 5)))
        pylab.title('Central Pressure (hPa)', fontsize=10)
        pylab.ylabel('Counts', fontsize=8)
        pylab.xlim(900, 1010)
        pylab.xlabel('Pressure (hPa)', fontsize=8)
        pylab.grid()

    def plotPressureRate(self):
        """
        Plot histogram of pressure rate of change
        """
        pylab.figure()
        pylab.hist(self.pressureRate, list(range(-10, 10, 1)))
        pylab.title('Pressure Rate (hPa/h)', fontsize=10)
        pylab.ylabel('Counts', fontsize=8)
        pylab.xlim(-10, 10)
        pylab.xlabel('Change of pressure (hPa/h)', fontsize=8)
        pylab.grid()

if __name__ == "__main__":
    try:
        configFile = sys.argv[1]
    except IndexError:
        # Try loading config file with same name as python script
        configFile = __file__.rstrip('.py') + '.ini'
        # If no filename is specified and default filename doesn't exist => raise error
        if not os.path.exists(configFile):
            error_msg = "No configuration file specified, please type: python main.py {config filename}.ini"
            raise IOError(error_msg)
    # If config file doesn't exist => raise error
    if not os.path.exists(configFile):
        error_msg = "Configuration file '" + configFile +"' not found"
        raise IOError(error_msg)

    #dataPath = config.cnfGetIniValue(configFile,'Input','Path', os.getcwd())
    dataPath = '/home/carthur/atcram/data/output/australia/run1'

    lon_lat = os.path.join(dataPath, 'all_lon_lat')
    bearing = os.path.join(dataPath, 'all_bearing')
    speed = os.path.join(dataPath, 'all_speed')
    pressure = os.path.join(dataPath, 'all_pressure')
    bearing_rate = os.path.join(dataPath, 'bearing_rate')
    speed_rate = os.path.join(dataPath, 'speed_rate')
    pressure_rate = os.path.join(dataPath, 'pressure_rate')

    #gridSpace = eval(config.cnfGetIniValue(configFile, 'Parameters', 'gridSpace')
    #gridLimit = eval(config.cnfGetIniValue(configFile, 'Parameters', 'gridLimit')
    gridLimit = {'xMin':90, 'xMax':180, 'yMin':-40, 'yMax':0}
    gridSpace = {'x':5, 'y':5}

    cellNum = stats.getCellNum(115., -15., gridLimit, gridSpace)

    p = PlotCell(lon_lat, bearing, speed, pressure, bearing_rate, speed_rate,
                 pressure_rate, cellNum, gridLimit, gridSpace)
    p.plotBearing()
    p.plotSpeed()
    p.plotSpeedRate()
    p.plotPressure()
    p.plotPressureRate()
    pylab.show()
