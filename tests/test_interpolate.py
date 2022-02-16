"""
:mod:`test_interpolate` -- test suite for :class:`Evaluate.interpolateTracks` module
===================================================================

"""

import os

import sys
from os.path import join as pjoin
import unittest
from tests import NumpyTestCase
import numpy as np

from datetime import datetime, timedelta

from cftime import date2num as cfdate2num, num2date as cfnum2date
import matplotlib.dates as mdates

from Evaluate import interpolateTracks

class TestDateInterpolation(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        base = datetime(2000, 1, 1)
        self.dates = np.array([base + timedelta(hours=i) for i in range(24)])
        self.absurddates = np.array([datetime(4500, 1, 1) + 
                                     timedelta(hours=i) for i in range(24)])

if __name__ == "__main__":
    unittest.main()