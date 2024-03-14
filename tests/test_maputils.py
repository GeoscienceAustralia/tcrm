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


 Title: testMapUtils.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2006-11-03
 Description: Test output from utils/map.py
 The functions included in this test were derived from MATLAB code. The MATLAB code
 is used to generate the results against which output from the Python functions are
 compared.

 SeeAlso: (related programs)
 Constraints:
"""

import os, sys
from scipy import array, arange, pi
import numpy as np
import numpy.ma as ma
import unittest
from tests import NumpyTestCase
try:
    from . import pathLocate
except:
    from tests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import maputils
from Utilities.files import flStartLog

class TestMapUtils(NumpyTestCase.NumpyTestCase):
    X = np.array([3, 5, 9, 11, 65])
    Y = np.array([4, 12, 40, 60, 72])
    Z = np.array([5, 13, 41, 61, 97])
    lat = np.arange(0, -21, -1, 'd')
    lon = np.arange(130, 151, 1, 'd')
    theta = np.arange(-2.*np.pi, 2.*np.pi, np.pi/4, 'd')

    findpts = [ 131.5, 135, 140., 140.9 ]
    indices = [ 1, 5, 10, 11]
    nearvals = [131.0, 135.0, 140.0, 141.0]

    XN = np.array([-111.1251, -111.1251, -111.1251, -111.1251, -111.1251,
                   -111.1251, -111.1251, -111.1251, -111.1251, -111.1251, -111.1251,
                   -111.1251, -111.1251, -111.1251, -111.1251, -111.1251, -111.1251,
                   -111.1251, -111.1251, -111.1251])

    YE = np.array([111.1082, 111.0574, 110.9728, 110.8544, 110.7022, 110.5164,
                   110.2968, 110.0437, 109.7570, 109.4369, 109.0834, 108.6968, 108.2770,
                   107.8242, 107.3386, 106.8203, 106.2695, 105.6863, 105.0709, 104.4234])

    AZI = np.array([134.81195977, 134.82952992, 134.85594569, 134.89121489,
                    134.93535056, 134.98837098, 135.05029967, 135.12116539,
                    135.20100215, 135.2898492 , 135.38775107, 135.4947575 ,
                    135.61092353, 135.73630944, 135.87098074, 136.01500822,
                    136.16846789, 136.33144099, 136.50401397, 136.68627847])

    Dist = np.array([156.89956829, 156.8761494 , 156.82932912, 156.75914243,
                     156.66564185, 156.54889743, 156.40899678, 156.24604514,
                     156.0601654 , 155.85149814, 155.62020174, 155.36645239,
                     155.09044421, 154.79238931, 154.4725179 , 154.13107836,
                     153.76833737, 153.38457998, 152.98010978, 152.55524897])

    def test_ValuesXY2r(self):
        """Test xy2R function"""
        result = maputils.xy2r(self.X, self.Y)
        for r, zz in zip(result, self.Z):
            self.assertEqual(r, zz)

    def test_ValueslatLon2XY(self):
        """Test latLon2XY function"""
        [X,Y] = maputils.latLon2XY(0, 0, self.lat, self.lon)
        for xx, xn in zip(X, self.XN):
            self.assertAlmostEqual(xx, xn, 2)
        for yy, ye in zip(Y, self.YE):
            self.assertAlmostEqual(yy, ye, 2)

    def test_ValueslatLon2Azi(self):
        """Test latLon2Azi function"""
        [A,D] = maputils.latLon2Azi(self.lat, self.lon)
        for aa, azi in zip(A, self.AZI):
            self.assertAlmostEqual(aa, azi, 2)
        for dd, dist in zip(D, self.Dist):
            self.assertAlmostEqual(dd, dist, 2)

    def  test_GridLatLonDist(self):
        """Test gridLatLonDist function"""
        cLon = 0.5
        cLat = 1.0
        lonArray = np.array([0.59297447, 0.20873497, 0.44271653, 0.36579662, 0.06680392])
        latArray = np.array([0.5019297, 0.42174226, 0.23712093, 0.02745615, 0.13316245])
        expected = ma.array([[56.30400807, 64.11584285, 55.7129098, 57.32175097, 73.35094702],
                             [65.08411729, 71.94898913, 64.57343322, 65.9665517 , 80.28821237],
                             [85.4022051, 90.74293694, 85.0136491, 86.07661618, 97.4877433],
                             [108.56672733, 112.8162381 , 108.2613338, 109.09805046, 118.30941319],
                             [96.87985127, 101.61920664, 96.537498, 97.47489149, 107.6850083]])

        dist = maputils.gridLatLonDist(cLon, cLat, lonArray, latArray)
        self.numpyAssertAlmostEqual(dist, expected)

    def test_GridLatLonBear(self):
        """Test gridLatLon function"""
        cLon = 0.5
        cLat = 1.0
        lonArray = np.array([0.59297447, 0.20873497, 0.44271653, 0.36579662, 0.06680392])
        latArray = np.array([0.5019297, 0.42174226, 0.23712093, 0.02745615, 0.13316245])
        expected = np.array([[2.95705149, -2.61243611, -3.0270878 , -2.87840196, -2.42573357],
                             [2.98217456, -2.67499093, -3.04285356, -2.91354889, -2.4986278 ],
                             [3.02031491, -2.77686502, -3.06664317, -2.96745336, -2.62513214],
                             [3.04627842, -2.85059019, -3.08275714, -3.0044598 , -2.72252399],
                             [3.03474017, -2.81742223, -3.07560296, -2.9879869 , -2.67812704]])

        bear = maputils.gridLatLonBear(cLon, cLat, lonArray, latArray)
        self.numpyAssertAlmostEqual(bear, expected)

    def test_Bearing(self):
        """Test conversion from bearing to theta and back again"""
        for th in self.theta:
            bearing = maputils.theta2bearing(th)
            result = maputils.bearing2theta(bearing)
            if th < 0.0:
                self.assertAlmostEqual(th+2.*pi, result)
            else:
                self.assertAlmostEqual(th, result)

    def test_findindex_Err(self):
        """Test that find_index raises ValueError if second arg is an array"""
        self.assertRaises(ValueError, maputils.find_index, self.lon, self.findpts)

    def test_findindex(self):
        """Test find_index function"""
        for pt,idx in zip(self.findpts, self.indices):
            i = maputils.find_index(self.lon, pt)
            self.assertEqual(i,idx)

    def test_findnearest(self):
        """Test find_nearest function"""
        for pt,val in zip(self.findpts, self.nearvals):
            v = maputils.find_nearest(self.lon, pt)
            self.assertEqual(v,val)

    def test_findnearest_Err(self):
        """Test that find_nearest raises ValueError if second arg is an array"""
        self.assertRaises(ValueError, maputils.find_nearest, self.lon, self.findpts)

#class TestInput(unittest.TestCase):
#   xx=[1, 3, 5, 9, 11]
#   yy=[1, 4, 12, 40, 60]
#   lat=np.arange(-20,0)
#   lon=np.arange(120,140)

#   def testInputxy2r(self):
#       self.assertRaises(maputils.ArrayMismatch, maputils.xy2r, xx, yy[0:len(y)-1])

#   def testInputll2xy(self):
#       self.assertRaises(maputils.ArrayMismatch, maputils.latLon2XY, lat, lon[0:len(lon)-1])

#   def testInputll2azi(self):
#       self.assertRaises(maputils.ArrayMismatch, maputils.latLon2Azi, lat, lon[0:len(lon)-1])

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestMapUtils,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
