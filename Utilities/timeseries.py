#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

 Title: timeseries.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 10/02/08 1:20:PM
 Description: Extract station timeseries from each timestep of a simulation.
 This samples the regional wind speed, not the site-specific wind speed.
 To include site-specific effects, you will first need to include the multiplier
 values for each site in the station file, then run tsmultipliers.py
 to apply said multipliers to the output.

 Version :$Rev: 528 $

 $Id: timeseries.py 528 2011-11-23 21:53:18Z carthur $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

from config import cnfGetIniValue
from files import flLoadFile, flSaveFile
import numpy
from matplotlib.dates import num2date

__version__ = '$Id: timeseries.py 528 2011-11-23 21:53:18Z carthur $'

ISO_FORMAT = "%Y-%m-%d %H:%M"

class timeseries:
    """timeseries:

    Description:

    Parameters:
    Members:
    Methods:
    Internal methods:
    """

    def __init__(self, configFile, stnFile=None, outputPath=None):
        self.configFile = configFile
        self.logger = logging.getLogger()
        self.meta = False
        if stnFile is None:
            stnFile = cnfGetIniValue(self.configFile, 'Timeseries', 'StationFile')

        if outputPath is None:
            self.outputPath = os.path.join(cnfGetIniValue(self.configFile, 'Output', 'Path'), 'timeseries')
        else:
            self.outputPath = outputPath
        self.logger.debug("Loading stations from %s"%stnFile)
        self.logger.debug("Timeseries data will be written into %s"%self.outputPath)
        stndata = flLoadFile(stnFile, delimiter=',')
        # If there is more than 3 columns, save the additional columns as 'metadata'
        if stndata.shape[1] > 3:
            self.metadata = stndata[:,3:]
            self.meta = True
        self.stnid = stndata[:,0]
        self.stnlon = stndata[:,1].astype(float)
        self.stnlat = stndata[:,2].astype(float)
        self.data = {}
        for stn in self.stnid:
            self.data[stn] = []

    def extract(self, spd, uu, vv, prs, gridx, gridy, tstep):
        """
        Extract data from the grid at the given locations.
        Data is stored in a dictionary, with keys as the station id's.
        """
        dt = num2date(tstep,tz=None)
        dtfmt = dt.isoformat(' ')
        if spd is None:
            # Set these values to zero:
            for i in range(len(self.stnid)):
                if self.meta:
                    self.data[self.stnid[i]].append(list(numpy.concatenate(([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs],self.metadata[i,:]))))
                else:
                    self.data[self.stnid[i]].append([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs])
        else:
            self.logger.debug("Extracting timeseries data")
            for i in range(len(self.stnid)):
                if (float(self.stnlon[i]) < gridx.min() or \
                    float(self.stnlon[i]) > gridx.max() or \
                    float(self.stnlat[i]) < gridy.min() or \
                    float(self.stnlat[i]) > gridy.max()):
                    self.logger.debug("Station %s (%8.4f, %8.4f) is outside of windfield grid" % \
                                       (repr(self.stnid[i]), float(self.stnlon[i]), float(self.stnlat[i])))
                    if self.meta:
                        self.data[self.stnid[i]].append(list(numpy.concatenate(([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs[0,0]],self.metadata[i,:]))))
                    else:
                        self.data[self.stnid[i]].append([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs[0,0]])
                else:
                    self.logger.debug("Sampling data for station %s (%8.4f, %8.4f)" % \
                                       (repr(self.stnid[i]), float(self.stnlon[i]), float(self.stnlat[i])))
                    x = gridx.searchsorted(float(self.stnlon[i]))
                    y = gridy.searchsorted(float(self.stnlat[i]))
                    s = spd[y,x]
                    u = uu[y,x]
                    v = vv[y,x]
                    b = numpy.mod((180./numpy.pi)*numpy.arctan2(-u,-v), 360.)
                    p = prs[y,x]
                    if self.meta:
                        self.data[self.stnid[i]].append(list(numpy.concatenate(([dt,self.stnlon[i], self.stnlat[i],s,u,v,b,p],self.metadata[i,:]))))
                    else:
                        self.data[self.stnid[i]].append([dt,self.stnlon[i],self.stnlat[i],s,u,v,b,p])

    def shutdown(self):
        """
        Write the data to file, each station to a separate file.
        """
        header = 'Time,Longitude,Latitude,Speed,UU,VV,Bearing,Pressure'
        for stn in self.stnid:
            if numpy.any(numpy.array(self.data[stn])[:,3]>0.0):
                tmpdata = numpy.array(self.data[stn])
                fname = os.path.join(self.outputPath, 'ts.%s.csv'%str(stn))
                tmpdata[:,0] = numpy.array([tmpdata[i,0].strftime(ISO_FORMAT) for i in xrange(tmpdata[:,0].size)])
                flSaveFile(fname, numpy.array(tmpdata), header, ',', '%s')
        self.logger.info("Station data written to file")
