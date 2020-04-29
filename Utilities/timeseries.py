"""
:mod:`timeseries` - Extract timeseries from each timestep of a simulation
=========================================================================

Extract station timeseries from each timestep of a simulation.
This samples the regional wind speed, not the site-specific wind speed.
To include site-specific effects, you will first need to include the
multiplier values for each site in the station file, then run
tsmultipliers.py to apply said multipliers to the output.

"""

import logging
from os.path import join as pjoin
from configparser import NoOptionError

import numpy as np

from Utilities.config import ConfigParser
from Utilities.files import flLoadFile
from Utilities.maputils import find_index
from Utilities.dynarray import DynamicRecArray
from .shptools import shpGetVertices

#from config import NoOptionError
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

ISO_FORMAT = "%Y-%m-%d %H:%M"

OUTPUT_NAMES = ('Station', 'Time', 'Longitude', 'Latitude',
                'Speed', 'UU', 'VV', 'Bearing',
                'Pressure')
OUTPUT_TYPES = ['|U16', '|U16', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
OUTPUT_FMT = ['%s', '%s', '%9.5f', '%9.5f',
              '%6.2f', '%6.2f', '%6.2f', '%6.2f',
              '%7.2f']

MINMAX_NAMES = ('Station', 'Time', 'Longitude', 'Latitude',
                'Speed', 'UU', 'VV', 'Bearing', 'Pressure')
MINMAX_TYPES = ['|U16', '|U16', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8']
MINMAX_FMT = ['%s', '%s', '%9.5f', '%9.5f',
              '%6.2f', '%6.2f', '%6.2f', '%6.2f',
              '%7.2f']

CONFIG_DEFAULTS = """
[Timeseries]
StationID=None
"""

class Station(object):
    """Station:

    Description: An object to represent a location for which time series 
                 data will be extracted

    Members:
    `id`: Unique id string for the station
    `lon`: Longitude of the station (geographic coordinates)
    `lat`: Latitude of the station (geographic coordinates)
    `data`: A `DynamicRecArray` to hold the time series data

    Methods:
    `insideGrid`: Determine if the station is inside the simulation domain.
    """

    def __init__(self, stationid, longitude, latitude):

        self.id = stationid
        self.lon = longitude
        self.lat = latitude
        self.data = DynamicRecArray(dtype={'names': OUTPUT_NAMES,
                                           'formats':OUTPUT_TYPES})

    def __getattr__(self, key):
        """
        Get the `key` from the `data` object.

        :param str key: the key to lookup in the `data` object.
        """
        if key.startswith('__') and key.endswith('__'):
            return super(Station, self).__getattr__(key)
        return self.data.data[key]

    def insideGrid(self, gridx, gridy):
        """
        Determine if a point is within the defined grid

        """
        if (float(self.lon) >= gridx.min() \
            and float(self.lon) <= gridx.max() and \
            float(self.lat) >= gridy.min() and \
            float(self.lat) <= gridy.max()):
            return True
        else:
            return False

class Timeseries(object):
    """Timeseries:

    Description: Extract data at a set of :class:`Station`s 

    Parameters:

    :param str configFile: Path to a TCRM configuration file

    Members:
    `meta`: Boolean whether additional metadata is attached to the `Station`s
    `outputPath`: Directory where extracted data will be stored in csv-format files
    `minfile`: Name of the file where minima for all `Station`s will be stored. 
               This will be the `outputPath` folder
    `maxfile`: As above, but for maxima (e.g. maximum wind speeds)
    `stations`: A list of `Station` objects, read from a file containing details of the stations
    
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

        stnFile = config.get('Timeseries', 'LocationFile')
        self.outputPath = pjoin(config.get('Output', 'Path'),
                                'process', 'timeseries')

        self.maxfile = pjoin(config.get('Output', 'Path'),
                             'process', 'maxima.csv')
        self.minfile = pjoin(config.get('Output', 'Path'),
                             'process', 'minima.csv')


        log.info(f"Loading timeseries stations from {stnFile}")
        log.debug(f"Timeseries data will be written into {self.outputPath}")
        self.stations = []
        if stnFile.endswith("shp"):
            try:
                key_name = config.get('Timeseries', 'StationID')
            except NoOptionError:
                key_name = None

            vertices = shpGetVertices(stnFile, key_name=key_name)

            for stn in list(vertices.keys()):
                lat = vertices[stn][0][1]
                lon = vertices[stn][0][0]
                lon = np.where(lon < 0., lon + 360., lon)
                self.stations.append(Station(stn, lon, lat))


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
            for sid, lon, lat in zip(stnid, stnlon, stnlat):
                self.stations.append(Station(sid, lon, lat))
        log.info(f"There are {len(self.stations)} stations that will collect timeseries data")

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
        xx = find_index(gridx, float(lon))
        yy = find_index(gridy, float(lat))
        ss = spd[yy, xx]
        ux = uu[yy, xx]
        vy = vv[yy, xx]
        bb = np.mod((180. / np.pi) * np.arctan2(-ux, -vy), 360.)
        pp = prs[yy, xx]

        return (ss, ux, vy, bb, pp)


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
        stns = 0
        for stn in self.stations:
            if stn.insideGrid(gridx, gridy):
                stns += 1
                result = self.sample(stn.lon, stn.lat, spd, uu, vv, prs,
                                     gridx, gridy)
                ss, ux, vy, bb, pp = result
                stn.data.append((str(stn.id), dt, stn.lon, stn.lat, ss,
                                 ux, vy, bb, pp))

            else:
                stn.data.append((str(stn.id), dt, stn.lon, stn.lat, 0.0, 0.0,
                                 0.0, 0.0, prs[0, 0]))
        log.debug("Extracted data for {0} stations".format(stns))

    def shutdown(self):
        """
        Write the data to file, each station to a separate file.
        """

        header = 'Station,Time,Longitude,Latitude,Speed,UU,VV,Bearing,Pressure'
        maxheader = ('Station,Time,Longitude,Latitude,Speed,'
                     'UU,VV,Bearing,Pressure')

        max_data = DynamicRecArray(dtype={'names': MINMAX_NAMES,
                                          'formats':MINMAX_TYPES})

        min_data = DynamicRecArray(dtype={'names': MINMAX_NAMES,
                                          'formats':MINMAX_TYPES})

        for stn in self.stations:

            if np.any(stn.data.data['Speed'] > 0.0):
                fname = pjoin(self.outputPath, 'ts.%s.csv' % str(stn.id))
                log.debug("Saving time series data to {0}".format(fname))
                with open(fname, 'wb') as fh:
                    np.savetxt(fh, np.array(stn.data.data), fmt=OUTPUT_FMT,
                               delimiter=',', header=header, comments='', encoding='ascii')

                max_step = np.argmax(stn.data.data['Speed'])
                min_step = np.argmin(stn.data.data['Pressure'])
                max_data.append(tuple(stn.data.data[max_step]))
                min_data.append(tuple(stn.data.data[min_step]))


        np.savetxt(self.maxfile, max_data.data, fmt=MINMAX_FMT, delimiter=',',
                   header=maxheader, comments='')
        np.savetxt(self.minfile, min_data.data, fmt=MINMAX_FMT, delimiter=',',
                   header=maxheader, comments='')
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
