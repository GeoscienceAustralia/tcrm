"""
:mod:`timeseries` - Extract timeseries from each timestep of a simulation
=========================================================================

Extract station timeseries from each timestep of a simulation.
This samples the regional wind speed, not the site-specific wind speed.
To include site-specific effects, you will first need to include the multiplier
values for each site in the station file, then run tsmultipliers.py
to apply said multipliers to the output.

"""

import logging

import datetime
import numpy as np

from os.path import join as pjoin

from Utilities.config import ConfigParser
from Utilities.files import flLoadFile
from Utilities.maputils import find_index
from matplotlib.dates import num2date
from shptools import shpGetVertices


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

ISO_FORMAT = "%Y-%m-%d %H:%M"

OUTPUT_NAMES = ('Time', 'Longitude', 'Latitude',
                'Speed', 'UU', 'VV', 'Bearing',
                'Pressure')

CONFIG_DEFAULTS = """
[Timeseries]
StationID=None
"""

class Timeseries(object):
    """Timeseries:

    Description:

    Parameters:
    Members:
    Methods:
    Internal methods:
    """

    def __init__(self, configFile):
        """
        Read configuration settings, load station data and set up
        output recarrays.
        
        :param str configFile: path to a configuration file.
        """

        config = ConfigParser()
        config.read(configFile)

        self.meta = False

        stnFile = config.get('Timeseries', 'StationFile')
        self.outputPath = pjoin(config.get('Output', 'Path'), 
                                    'process', 'timeseries')

        log.debug("Loading stations from %s"%stnFile)
        log.debug("Timeseries data will be written into %s"%self.outputPath)
        
        if stnFile.endswith("shp"):
            try:
                key_name = config.get('Timeseries', 'StationID')
            except ConfigParser.NoOptionError:
                key_name = None
                
            vertices = shpGetVertices(stnFile, key_name=key_name)
            self.stnid = vertices.keys()
            self.stnlon = []
            self.stnlat = []
            for stn in self.stnid:
                self.stnlon.append(vertices[stn][0])
                self.stnlat.append(vertices[stn][1])

        else:
            stndata = flLoadFile(stnFile, delimiter=',')
            # If there are more than 3 columns, save the additional 
            # columns as 'metadata'
            if stndata.shape[1] > 3:
                self.metadata = stndata[:, 3:]
                self.meta = True
            self.stnid = stndata[:, 0]
            self.stnlon = stndata[:, 1].astype(float)
            self.stnlat = stndata[:, 2].astype(float)

        self.data = {}
        self.maxdata = {}
        self.mindata = {}
            
        self.data.fromkeys(self.stnid, [])
        
        for i in xrange(len(self.stnid)):
            # For storing the maximum wind speed:
            self.maxdata[self.stnid[i]] = [self.stnid[i],
                                           self.stnlon[i],
                                           self.stnlat[i],
                                           '', 0.0, 0.0, 0.0, 0.0, 0.0]
            # For storing the minimum pressure:
            self.mindata[self.stnid[i]] = [self.stnid[i],
                                           self.stnlon[i],
                                           self.stnlat[i],
                                           '', 0.0, 0.0, 0.0, 0.0, 9999999.]
                                           
    def sample(self, lon, lat, spd, uu, vv, prs, gridx, gridy):
        x = find_index(gridx, float(lon))
        y = find_index(gridy, float(lat))
        s = spd[y, x]
        u = uu[y, x]
        v = vv[y, x]
        b = np.mod((180. / np.pi) * np.arctan2(-u, -v), 360.)
        p = prs[y, x]
        
        return (s, u, v, b, p)
        

    def extract(self, tstep, spd, uu, vv, prs, gridx, gridy):
        """
        Extract data from the grid at the given locations.
        Data is stored in a dictionary, with keys as the station id's.
        
        
        :param float tstep: time step being evaluated, as a float (output
                            from matplotlib.num2date)
        """
        dt = num2date(tstep, tz=None)
        if spd is None:
            # Set these values to zero:
            for i in xrange(len(self.stnid)):
                if self.meta:
                    self.data[self.stnid[i]].append(list(np.concatenate(([dt,
                             self.stnlon[i], self.stnlat[i], 0.0, 0.0, 0.0, 
                             0.0, prs], self.metadata[i, :]))))
                else:
                    self.data[self.stnid[i]].append([dt, self.stnlon[i],
                              self.stnlat[i], 0.0, 0.0, 0.0, 0.0, prs])

        else:
            # Inside the grid:
            for i, (stnid, lon, lat) in enumerate(zip(self.stnid, 
                                                      self.stnlon, 
                                                      self.stnlat)):
                if (float(lon) < gridx.min() or float(lon) > gridx.max() or \
                    float(lat) < gridy.min() or float(lat) > gridy.max()):
                    # Outside the grid:
                    if self.meta:
                        self.data[stnid].append(list(np.concatenate(([dt, lon, 
                                        lat, 0.0, 0.0, 0.0, 0.0, prs[0, 0]],
                                        self.metadata[i, :]))))
                    else:
                        self.data[stnid].append([dt, lon, lat, 0.0, 0.0, 0.0, 
                                                    0.0, prs[0, 0]])
                else:
                    result = self.sample(lon, lat, spd, uu, vv, prs, 
                                                  gridx, gridy)
                    s, u, v, b, p = result
                    if self.meta:
                        self.data[stnid].append(list(np.concatenate(([dt, lon, 
                                lat, s, u, v, b, p], self.metadata[i, :]))))
                    else:
                        self.data[stnid].append([dt, lon, lat, s, u, v, b, p])
                    if s > self.maxdata[stnid][4]:
                        self.maxdata[stnid][3:] = [dt, s, u, v, b, p]
                    if p < self.mindata[stnid][8]:
                        self.mindata[stnid][3:] = [dt, s, u, v, b, p]
                        

    def shutdown(self):
        """
        Write the data to file, each station to a separate file.
        """

        header = 'Time,Longitude,Latitude,Speed,UU,VV,Bearing,Pressure'
        maxheader = ('Station,Longitude,Latitude,Time,Speed,'
                        'UU,VV,Bearing,Pressure')
        fmt = ['%s', '%7.3f', '%7.3f', '%6.2f', '%6.2f', 
               '%6.2f', '%6.2f', '%7.2f']
        for stn in self.stnid:
            if np.any(np.array(self.data[stn])[:, 3] > 0.0):
                tmpdata = np.array(self.data[stn])
                fname = pjoin(self.outputPath, 'ts.%s.csv' % str(stn))
                tmpdata[:, 0] = np.array([tmpdata[i, 0].strftime(ISO_FORMAT) 
                                        for i in xrange(tmpdata[:, 0].size)])
                                        
                np.savetxt(fname, np.array(tmpdata), fmt=fmt,
                           delimiter=',', header=header)
        
        maxfname = pjoin(self.outputPath, 'maxwind.csv')
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

        np.savetxt(maxfname, np.array(maxdata), fmt='%s', 
                   header=maxheader, delimiter=',')
                     #['%s','%7.3f','%7.3f','%s','%6.2f','%6.2f',
                     #  '%6.2f','%6.2f','%7.2f'] )
        minfname = pjoin(self.outputPath, 'minpressure.csv')

        np.savetxt(minfname, np.array(mindata), fmt='%s', 
                   header=maxheader, delimiter=',')             
                   #['%s','%7.3f','%7.3f','%s','%6.2f','%6.2f',
                   #  '%6.2f','%6.2f','%7.2f'] )
        log.info("Station data written to file")

