"""
Track-related attributes
"""

import numpy as np
from datetime import datetime

trackFields = ('Indicator', 'CycloneNumber', 'Year', 'Month', 
               'Day', 'Hour', 'Minute', 'TimeElapsed', 'Datetime','Longitude',
               'Latitude', 'Speed', 'Bearing', 'CentralPressure',
               'WindSpeed', 'rMax', 'EnvPressure')

trackTypes = ('i', 'i', 'i', 'i', 
              'i', 'i', 'i', 'f', datetime,
              'f', 'f', 'f', 'f', 'f', 
              'f', 'f', 'f')

trackFormats = ('%i, %i, %i, %i,' 
                '%i, %i, %i, %5.1f,' '%s',
                '%8.3f, %8.3f, %6.2f, %6.2f, %7.2f,'
                '%6.2f, %6.2f, %7.2f')
                


class Track(object):

    """
    A single tropical cyclone track.

    The object exposes the track data through the object attributes.
    For example, If `data` contains the tropical cyclone track data
    (`numpy.array`) loaded with the :meth:`readTrackData` function,
    then the central pressure column can be printed out with the
    code::

        t = Track(data)
        print(t.CentralPressure)

    :type  data: numpy.ndarray
    :param data: the tropical cyclone track data.
    """

    def __init__(self, data):
        """
        :type  data: numpy.ndarray
        :param data: the tropical cyclone track data.
        """
        self.data = data
        self.trackId = None
        self.trackfile = None
        self.trackMinPressure = None
        self.trackMaxWind = None


    def __getattr__(self, key):
        """
        Get the `key` from the `data` object.

        :type  key: str
        :param key: the key to lookup in the `data` object.
        """
        if key.startswith('__') and key.endswith('__'):
            return super(Track, self).__getattr__(key)

        return self.data[key]

    def inRegion(self, gridLimit):
        """
        Check if the tropical cyclone track falls within a region.

        :type  gridLimit: :class:`dict`
        :param gridLimit: the region to check.
                          The :class:`dict` should contain the keys
                          :attr:`xMin`, :attr:`xMax`, :attr:`yMin` and
                          :attr:`yMax`. The *x* variable bounds the
                          latitude and the *y* variable bounds the
                          longitude.

        """
        xMin = gridLimit['xMin']
        xMax = gridLimit['xMax']
        yMin = gridLimit['yMin']
        yMax = gridLimit['yMax']

        return ((xMin <= np.min(self.Longitude)) and
                (np.max(self.Latitude) <= xMax) and
                (yMin <= np.min(self.Latitude)) and
                (np.max(self.Latitude) <= yMax))
                

