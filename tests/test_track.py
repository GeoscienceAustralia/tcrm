"""
:mod:`test_track` -- test suite for :class:`Utilities.track` module
===================================================================

"""

import os
import sys
from os.path import join as pjoin
import unittest
import NumpyTestCase
import numpy as np
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime, timedelta

try:
    import pathLocate
except:
    from unittests import pathLocate

unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

from Utilities import track

class TestTrack(NumpyTestCase.NumpyTestCase):
    """
    Test track class methods
    """

    def setUp(self):
        self.test_track_file = pjoin(unittest_dir, 'test_data', 'tracks.01.csv')
        self.nc_track_file = pjoin(unittest_dir, 'test_data', 'tracks.01.nc')

        data = np.loadtxt(self.test_track_file, comments='%', delimiter=',',
                          skiprows=1, dtype={'names':track.TCRM_COLS,
                                             'formats':track.TCRM_FMTS},
                          converters=track.TCRM_CNVT)
        datas = []
        self.tracks = []
        cycId = data['CycloneNumber']
        for i in range(1, np.max(cycId) + 1):
            datas.append(data[cycId] == i)

        n = len(datas)
        for i, data in enumerate(datas):
            mytrack = track.Track(data)
            mytrack.trackfile = self.test_track_file
            mytrack.trackId = (i, n)
            self.tracks.append(mytrack)

