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
import pdb

import datetime
import numpy as np

from os.path import join as pjoin

from Utilities.config import ConfigParser
from Utilities.files import flLoadFile
from Utilities.maputils import find_index
from matplotlib.dates import num2date
from shptools import shpGetVertices

#from config import NoOptionError
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

ISO_FORMAT = "%Y-%m-%d %H:%M"

OUTPUT_NAMES = ('Time', 'Longitude', 'Latitude',
                'Speed', 'UU', 'VV', 'Bearing',
                'Pressure')
OUTPUT_FMT = ['%s', '%7.3f', '%7.3f', 
              '%6.2f', '%6.2f', '%6.2f', '%6.2f', 
              '%7.2f']

CONFIG_DEFAULTS = """
[Timeseries]
StationID=None
"""

class Station(object):
    def __init__(self, station_id, longitude, latitude):
        
        self.id = station_id
        self.lon = longitude
        self.lat = latitude
        self.data = []

    def __getattr__(self, key):
        """
        Get the `key` from the `data` object.

        :type  key: str
        :param key: the key to lookup in the `data` object.
        """
        if key.startswith('__') and key.endswith('__'):
            return super(Track, self).__getattr__(key)
        return self.data[key]
        
    def insideGrid(self, gridx, gridy):
        if (float(self.lon) >= gridx.min() and float(self.lon) <= gridx.max() and \
            float(self.lat) >= gridy.min() and float(self.lat) <= gridy.max()):
            return True
        else:
            return False

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
        self.stations = []
        if stnFile.endswith("shp"):
            try:
                key_name = config.get('Timeseries', 'StationID')
            except NoOptionError:
                key_name = None
                
            vertices = shpGetVertices(stnFile, key_name=key_name)

            for stn in vertices.keys():
                self.stations.append(Station(stn, vertices[stn][0][0], 
                                                  vertices[stn][0][1]))

        
        else:
            stndata = flLoadFile(stnFile, delimiter=',')
            # If there are more than 3 columns, save the additional 
            # columns as 'metadata'
            if stndata.shape[1] > 3:
                self.metadata = stndata[:, 3:]
                self.meta = True
            stnid = stndata[:, 0]
            stnlon = stndata[:, 1].astype(float)
            stnlat = stndata[:, 2].astype(float)
            for id, lon, lat in zip(stnid, stnlon, stnlat):
                self.stations.append(Station(id, lon, lat))
        
    def sample(self, lon, lat, spd, uu, vv, prs, gridx, gridy):
        """
        Extract values from 2-dimensional grids at the given lat/lon.
        
        :param float lon: Longitude of the point to extract.
        :param float lat: Latitude of the point to extract.
        :param spd: :class:`numpy.ndarray` of speed values.
        :param uu: :class:`numpy.ndarray` of eastward wind speed values.
        :param vv: :class:`numpy.ndarray` of northward wind speed values.    
        :param prs: :class:`numpy.ndarray` of pressure values.
        :param gridx: :class:`numpy.ndarray` of grid longitudes.
        :param gridy: :class:`numpy.ndarray` of grid latitudes.

        :return: speed, esatward and northward wind components, and pressure
                 values at the given location
        :rtype: tuple
        """
        x = find_index(gridx, float(lon))
        y = find_index(gridy, float(lat))
        s = spd[y, x]
        u = uu[y, x]
        v = vv[y, x]
        b = np.mod((180. / np.pi) * np.arctan2(-u, -v), 360.)
        p = prs[y, x]
        
        return (s, u, v, b, p)
        

    def extract(self, dt, spd, uu, vv, prs, gridx, gridy):
        """
        Extract data from the grid at the given locations.
        Data is stored in a dictionary, with keys as the station id's.
        
        :param float tstep: time step being evaluated, as a float (output
                            from matplotlib.num2date)
        :param spd: :class:`numpy.ndarray` of speed values.
        :param uu: :class:`numpy.ndarray` of eastward wind speed values.
        :param vv: :class:`numpy.ndarray` of northward wind speed values.    
        :param prs: :class:`numpy.ndarray` of pressure values.
        :param gridx: :class:`numpy.ndarray` of grid longitudes.
        :param gridy: :class:`numpy.ndarray` of grid latitudes.
        
        """

        for stn in self.stations:
            if stn.insideGrid(gridx, gridy):
                result = self.sample(stn.lon, stn.lat, spd, uu, vv, prs,
                                      gridx, gridy)
                s, u, v, b, p = result
                stn.data.append([dt, stn.lon, stn.lat, s, u, v, b, p])

            else:
                stn.data.append([dt, stn.lon, stn.lat, 0.0, 0.0,  
                                          0.0, 0.0, prs[0, 0]])
                    

    def shutdown(self):
        """
        Write the data to file, each station to a separate file.
        """

        header = 'Time,Longitude,Latitude,Speed,UU,VV,Bearing,Pressure'
        #maxheader = ('Station,Longitude,Latitude,Time,Speed,'
        #                'UU,VV,Bearing,Pressure')
                
        for stn in self.stations:
            
            if np.any(np.array(stn.data)[:, 3] > 0.0):
            
                tmpdata = np.array(stn.data)
                fname = pjoin(self.outputPath, 'ts.%s.csv' % str(stn.id))
                tmpdata[:, 0] = np.array([tmpdata[i, 0].strftime(ISO_FORMAT) 
                                    for i in xrange(tmpdata[:, 0].size)])

                np.savetxt(fname, np.array(tmpdata), fmt=OUTPUT_FMT,
                           delimiter=',', header=header)
        

        """
        for stn in self.stations:
            if type(self.maxdata[stn.id][3]) == datetime.datetime:
                self.maxdata[stn.id][3] = self.maxdata[stn.id][3].strftime(ISO_FORMAT)
                self.mindata[stn.id][3] = self.mindata[stn.id][3].strftime(ISO_FORMAT)
                self.maxdata[stn.id][0] = str(int(stn.id))
                self.mindata[stn.id][0] = str(int(stn.id))
                maxdata.append(self.maxdata[stn.id])
                mindata.append(self.mindata[stn.id]) 
        
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
        """
        log.info("Station data written to file")

