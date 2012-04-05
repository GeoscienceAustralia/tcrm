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

 Title: timeseries.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 10/02/08 1:20:PM
 Description: Extract station timeseries from each timestep of a simulation.
 This samples the regional wind speed, not the site-specific wind speed.
 To include site-specific effects, you will first need to include the multiplier
 values for each site in the station file, then run tsmultipliers.py
 to apply said multipliers to the output.

 Version :$Rev: 686 $

 $Id: timeseries.py 686 2012-03-29 04:24:59Z carthur $
"""
import os, sys, pdb, logging, traceback
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

import datetime
import numpy

from config import cnfGetIniValue
from files import flLoadFile, flSaveFile, flLogFatalError
from matplotlib.dates import num2date
from shptools import *

__version__ = '$Id: timeseries.py 686 2012-03-29 04:24:59Z carthur $'

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
        if stnFile.endswith("shp"):
            vertices = shpGetVertices(stnFile, keyName=cnfGetIniValue( self.configFile, 'Timeseries', 'StationID', None ) )
            self.stnid = vertices.keys()
            self.stnlon = []
            self.stnlat = []
            for stn in self.stnid:
                self.stnlon.append( vertices[stn][0] )
                self.stnlat.append( vertices[stn][1] )

        else:
            stndata = flLoadFile(stnFile, delimiter=',')
            # If there is more than 3 columns, save the additional columns as 'metadata'
            if stndata.shape[1] > 3:
                self.metadata = stndata[:,3:]
                self.meta = True
            self.stnid = stndata[:,0]
            self.stnlon = stndata[:,1].astype(float)
            self.stnlat = stndata[:,2].astype(float)

        self.data = {}
        self.maxdata = {}
        self.mindata = {}
        for stn in self.stnid:
            self.data[stn] = []
        for i in xrange( len( self.stnid ) ):
            # For storing the maximum wind speed:
            self.maxdata[self.stnid[i]] = [self.stnid[i],self.stnlon[i],self.stnlat[i],'',0.0,0.0,0.0,0.0,0.0]
            # For storing the minimum pressure:
            self.mindata[self.stnid[i]] = [self.stnid[i],self.stnlon[i],self.stnlat[i],'',0.0,0.0,0.0,0.0,9999999.]

    def _write(self, filename, data, header=None, delimiter=',', fmt='%.18e'):
        """
        Write the data to a csv file. To get around the problem of writing
        strings and floats to the one file, we change the way 'fmt' is
        used. In this function, it refers to the whole line, not just each
        element as in flSaveFile()
        """
        try:
            directory, fname = os.path.split( filename )
        except AttributeError:
            self.logger.exception( 'Input filename is not a string' )
            flLogFatalError( traceback.format_exc( ).splitlines( ) )

        if not os.path.isdir( directory ):
            try:
                os.makedirs( directory )
            except:
                self.logger.exception( 'Cannot build path: %s'%( directory ) )
                flLogFatalError( traceback.format_exc( ).splitlines( ) )

        self.logger.debug( 'Saving data to %s'%( filename ) )

        if type( filename ) == str:
            if fname.endswith( '.gz' ):
                import gzip
                fh = gzip.open( filename, 'wb' )
            else:
                fh = open( filename, 'w' )

        elif hasattr( filename, 'seek' ):
            fh = filename

        else:
            self.logger.error( 'Filename must be a string or file handle' )
            raise IOError( 'Filename must be a string or file handle' )

        if header:
            fh.write( '%' + header + '\n' )

        X = numpy.asarray( data )
        origShape = None

        if len( X.shape ) == 1:
            origShape = X.shape
            X.shape = len( X ), 1

        for row in X:
            try:
                if type(fmt) == list:
                    fh.write( delimiter.join( [f%v for f,v in zip(fmt,row)] ) + '\n' )
                elif type(fmt) == str:
                    fh.write( delimiter.join( [fmt%val for val in row] ) + '\n' )
                else:
                    self.logger.exception( "Mismatch between format string and values in _write" )
                    raise TypeError, "Mismatch between format string and values in _write"

            except ValueError:
                self.logger.exception( "Cannont write data to file" )
                raise

        fh.close( )

        if origShape is not None:
            X.shape = origShape


    def extract(self, spd, uu, vv, prs, gridx, gridy, tstep):
        """
        Extract data from the grid at the given locations.
        Data is stored in a dictionary, with keys as the station id's.
        """
        dt = num2date(tstep,tz=None)
        dtfmt = dt.isoformat(' ')
        if spd is None:
            # Set these values to zero:
            for i in xrange(len(self.stnid)):
                if self.meta:
                    self.data[self.stnid[i]].append(list(numpy.concatenate(([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs],self.metadata[i,:]))))
                else:
                    self.data[self.stnid[i]].append([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs])

        else:
            #self.logger.debug("Extracting timeseries data")
            for i in xrange(len(self.stnid)):
                if (float(self.stnlon[i]) < gridx.min() or \
                    float(self.stnlon[i]) > gridx.max() or \
                    float(self.stnlat[i]) < gridy.min() or \
                    float(self.stnlat[i]) > gridy.max()):
                    #self.logger.debug("Station %s (%8.4f, %8.4f) is outside of windfield grid" % \
                    #                   (repr(self.stnid[i]), float(self.stnlon[i]), float(self.stnlat[i])))
                    if self.meta:
                        self.data[self.stnid[i]].append(list(numpy.concatenate(([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs[0,0]],self.metadata[i,:]))))
                    else:
                        self.data[self.stnid[i]].append([dt,self.stnlon[i],self.stnlat[i],0.0,0.0,0.0,0.0,prs[0,0]])
                else:
                    #self.logger.debug("Sampling data for station %s (%8.4f, %8.4f)" % \
                    #                   (repr(self.stnid[i]), float(self.stnlon[i]), float(self.stnlat[i])))
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
                    if s > self.maxdata[self.stnid[i]][4]:
                        self.maxdata[self.stnid[i]][3:] = [dt,s,u,v,b,p]
                    if p < self.mindata[self.stnid[i]][8]:
                        self.mindata[self.stnid[i]][3:] = [dt,s,u,v,b,p]
                        

    def shutdown(self):
        """
        Write the data to file, each station to a separate file.
        """

        header = 'Time,Longitude,Latitude,Speed,UU,VV,Bearing,Pressure'
        maxheader = 'Station,Longitude,Latitude,Time,Speed,UU,VV,Bearing,Pressure'
        for stn in self.stnid:
            if numpy.any(numpy.array(self.data[stn])[:,3]>0.0):
                tmpdata = numpy.array(self.data[stn])
                fname = os.path.join(self.outputPath, 'ts.%s.csv'%str(stn))
                tmpdata[:,0] = numpy.array([tmpdata[i,0].strftime(ISO_FORMAT) for i in xrange(tmpdata[:,0].size)])
                self._write(fname, numpy.array(tmpdata), header, ',', 
                            ['%s','%7.3f','%7.3f','%6.2f','%6.2f','%6.2f','%6.2f','%7.2f'])
        
        maxfname = os.path.join(self.outputPath, 'maxwind.csv' )
        mindata = []
        maxdata = []
        for stn in self.stnid:
            if type(self.maxdata[stn][3]) == datetime.datetime:
                self.maxdata[stn][3] = self.maxdata[stn][3].strftime(ISO_FORMAT)
                self.mindata[stn][3] = self.mindata[stn][3].strftime(ISO_FORMAT)
                self.maxdata[stn][0] = str(int(stn))
                self.mindata[stn][0] = str(int(stn))
                maxdata.append(self.maxdata[stn])
                mindata.append(self.mindata[stn])
            else:
                pass

        self._write( maxfname, numpy.array(maxdata), maxheader, ',', fmt='%s' )
                     #['%s','%7.3f','%7.3f','%s','%6.2f','%6.2f','%6.2f','%6.2f','%7.2f'] )
        minfname = os.path.join(self.outputPath,'minpressure.csv')

        self._write( minfname, numpy.array(mindata), maxheader, ',',fmt='%s') 
                     #['%s','%7.3f','%7.3f','%s','%6.2f','%6.2f','%6.2f','%6.2f','%7.2f'] )
        self.logger.info("Station data written to file")

