import unittest
import numpy as np
from numpy.testing import assert_equal

from os.path import join as pjoin
from . import pathLocate

unittest_dir = pathLocate.getUnitTestDirectory()

import Utilities.columns as columns
from Utilities import config

class TestReadCSV(unittest.TestCase):

    def setUp(self):
        config.reset()
        self.configFile = pjoin(unittest_dir, 'test_data','test_cols.ini')
        self.inputFile = pjoin(unittest_dir, 'test_data','test_cols.csv')
        self.source = 'TEST'
        self.badinput = pjoin(unittest_dir, 'test_data','test.csv')
        self.nonfileinput = 3

        self.test_data = np.array([(1, 2, 3), (4, 5, 6), (7, 8, 9)],
                                  dtype=[('A', float),
                                         ('B', float),
                                         ('C', int)])
    def tearDown(self):
        config.reset()

    def test_badinput(self):
        self.assertRaises(IOError, columns.colReadCSV,
                          self.configFile, self.badinput, self.source)
        self.assertRaises(TypeError, columns.colReadCSV,
                          self.configFile, self.nonfileinput, self.source)

    def test_readCSV(self):
        loadedData = columns.colReadCSV(self.configFile, self.inputFile, self.source)
        assert_equal(loadedData, self.test_data)
        
if __name__ == "__main__":
    unittest.main()
