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


 Title: TestStatut.py
 Author: N. Habili, nariman.habili@ga.gov.au
 CreationDate: 2006-12-20
 Description: Unit testing module for statutils.py

 Version: $Rev: 563 $

 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: TestStat.py 563 2007-10-24 02:52:40Z carthur $
"""
import os, sys
import unittest
import pickle
from scipy import array, zeros
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
import Utilities.stats as statutils
from Utilities.files import flStartLog


class TestStats(NumpyTestCase.NumpyTestCase):

    gridLimit = {'xMin':70, 'xMax':180, 'yMin':-40, 'yMax':0}
    gridSpace = {'x':5, 'y':5}

    def test_cdf(self):
        """Testing cdf"""
        pdf = array([0.33350065, 0.71365127, 0.42428029,
                     0.99204143, 0.01738811])
        cdf = array([0.13442936, 0.42209201, 0.59311334,
                     0.9929911 , 1.0])
        y = array([0, 1, 2, 3, 4])

        self.numpyAssertAlmostEqual(statutils.cdf(y, pdf), cdf)
        
    def test_cdfzeros(self):
        """Test cdf returns zero array for zero input"""
        x = array([0, 1, 2, 3, 4])
        y = zeros(len(x))
        self.numpyAssertAlmostEqual(statutils.cdf(y, x),
                                    zeros(len(x)))
    def test_cdf2d(self):
        """Testing cdf2d"""
        pdf = array([[0.07383899, 0.66509486, 0.21795259, 0.349258],
                    [0.72284544, 0.56474223, 0.16157475, 0.71899118],
                    [0.73459315, 0.29610227, 0.94534262, 0.4801483 ],
                    [0.58482537, 0.25089279, 0.40459775, 0.65046793]])

        cdf = array([[0.0094408, 0.0944775, 0.12234415, 0.16699906],
                    [0.10186128, 0.25910395, 0.30762902, 0.44421163],
                    [0.19578378, 0.39088508, 0.5602783 , 0.75825101],
                    [0.27055752, 0.49773705, 0.71886075, 1.0]])

        x = array([0, 1, 2, 3])
        y = array([0, 1, 2, 3])

        self.numpyAssertAlmostEqual(statutils.cdf2d(x, y, pdf), cdf)

    def test_GetCellNum(self):
        """Testing getCellNum"""
        #valid values
        knownValues = ((173, -39, 174), (173, -30, 152), (133, -30, 144),
                       (70, 0, 0), (173.5, -0.5, 20))
        for lon, lat, cell in knownValues:
            result = statutils.getCellNum(lon, lat, self.gridLimit,
                                          self.gridSpace)
            self.assertEqual(cell, result)

        #outside of limits
        invalidLatLongs = ((60, -20), (190, -20), (90, -50), (90, 20),
                           (180, 0), (70, -40))
        for lon, lat in invalidLatLongs:
            self.assertRaises(ValueError, statutils.getCellNum,
                              lon, lat, self.gridLimit, self.gridSpace)

    def test_GetCellLonLat(self):
        """Testing getCellLonLat"""
        #valid values
        knownValues = ((0, 70, 0), (50, 100, -10), (150, 160, -30),
                       (175, 175, -35))
        for cell, lon, lat in knownValues:
            (lon_, lat_) = statutils.getCellLonLat(cell, self.gridLimit,
                                                   self.gridSpace)
            self.assertEqual((lon_,lat_), (lon,lat))

        invalidCells = (-10, 176, 2000)
        for cell in invalidCells:
            self.assertRaises(IndexError, statutils.getCellLonLat,
                              cell, self.gridLimit, self.gridSpace)

    def test_ValidCellNum(self):
        """Testing validCellNum"""
        knownValues = ((-10, False), (0, True), (100, True),
                       (175, True), (176, False), (2000, False))
        for cell, valid in knownValues:
            result = statutils.validCellNum(cell, self.gridLimit,
                                            self.gridSpace)
            self.assertEqual(result, valid)

    def test_MaxCellNum(self):
        """Testing maxCellNum"""
        maxCellNum = 175
        self.assertEqual(maxCellNum, statutils.maxCellNum(self.gridLimit,
                                                          self.gridSpace))

    def test_GetOccurence(self):
        """Testing getOccurance"""
        occurList = [-10, 20, 30, 40]

        indList = [0, 2]
        self.numpyAssertEqual(array([-10, 30]),
                              statutils.getOccurence(occurList, indList))

        indList = []
        self.numpyAssertEqual(array([]),
                              statutils.getOccurence(occurList, indList))

        indList = [-10, 2]
        self.assertRaises(IndexError, statutils.getOccurence,
                          occurList, indList)

        indList = [0, 20]
        self.assertRaises(IndexError, statutils.getOccurence,
                          occurList, indList)

        occurList = []
        self.assertRaises(IndexError, statutils.getOccurence,
                          occurList, indList)

    def test_GetMaxRange(self):
        """Testing getMaxRange"""
        knownValues = ((10, 20, 3, 22), (10, 200, 9, 208))
        for minVal, maxVal, step, maxRange in knownValues:
            result = statutils.statMaxRange(minVal, maxVal, step)
            self.assertEqual(maxRange, result)

    def test_GetMinRange(self):
        """Testing getMinRange"""
        knownValues = ((10, 20, 9, 2), (-123, 2002, 3, -125))
        for minVal, maxVal, step, minRange in knownValues:
            result = statutils.statMinRange(minVal, maxVal, step)
            self.assertEqual(minRange, result)
            
    def test_BandWidth(self):
        """Testing bandwidth function"""
        self.assertRaises(TypeError, statutils.bandwidth, 'string')
        self.assertRaises(TypeError, statutils.bandwidth, 0.5)
        data = pickle.load(open(os.path.join(unittest_dir, 'test_data', 
                                             'kde_parameters_pressure_rate.pkl'), 'rb'))
        result = statutils.bandwidth(data)
        self.assertAlmostEqual(result, 0.116510487, places=7)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestStats,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
