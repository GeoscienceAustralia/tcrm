"""
:mod:`test_track` -- test suite for :class:`Utilities.track` module
===================================================================

"""

import os
import sys
from os.path import join as pjoin
import unittest
from . import NumpyTestCase
import numpy as np
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime, timedelta

try:
    from . import pathLocate
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
        pass
        """
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
        """

class TestLoadTracks(NumpyTestCase.NumpyTestCase):
    """
    Test functions to load tracks from files
    """

    def setUp(self):
        self.test_track_file_nc = pjoin(unittest_dir, 'test_data',
                                        'tracks.01.nc')
        self.dummy_track_file = pjoin(unittest_dir, 'test_data',
                                      'tracks.02.nc')
        self.dummy_path = pjoin(unittest_dir, 'test_path')
        self.test_track_file_csv = pjoin(unittest_dir, 'test_data',
                                         'tracks.01.csv')

    def test_NoPathExists(self):
        """loadTracksFromPath raises IOError if path doesn't exist"""
        self.assertRaises(IOError, track.loadTracksFromPath,
                          self.dummy_path)

    def test_NotListInput(self):
        """loadTracksFromFiles raises TypeError if passed non-list"""
        with self.assertRaises(TypeError):
            list(track.loadTracksFromFiles(self.dummy_track_file))

    def test_NotFileInput(self):
        """loadTracksFromFiles raises TypeError if list contents not files"""
        with self.assertRaises(TypeError):
            list(track.loadTracksFromFiles([0, 1, 2, 3, 4]))

        with self.assertRaises(IOError):
            list(track.loadTracksFromFiles(['string01', 'string02']))


class TestncReadTrackData(NumpyTestCase.NumpyTestCase):
    """
    Test functions to load tracks from netCDF files
    """
    def setUp(self):
        self.test_track_file = pjoin(unittest_dir, 'test_data',
                                     'tracks.01.nc')
        self.dummy_track_file = pjoin(unittest_dir, 'test_data',
                                      'tracks.02.nc')
    def test_NotFileInput(self):
        """ncReadTrackData raises IOError if passed non-file"""
        self.assertRaises(IOError, track.ncReadTrackData,
                          self.dummy_track_file)

if __name__ == "__main__":
    unittest.main()
