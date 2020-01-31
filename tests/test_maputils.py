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
import numpy
import numpy.ma as ma
import unittest
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import maputils
from Utilities.files import flStartLog

class TestMapUtils(NumpyTestCase.NumpyTestCase):
    X = array([3, 5, 9, 11, 65])
    Y = array([4, 12, 40, 60, 72])
    Z = array([5, 13, 41, 61, 97])
    lat = arange(0, -21, -1, 'd')
    lon = arange(130, 151, 1, 'd')
    theta = arange(-2.*pi,2.*pi,pi/4, 'd')

    findpts = [ 131.5, 135, 140., 140.9 ]
    indices = [ 1, 5, 10, 11]
    nearvals = [131.0, 135.0, 140.0, 141.0]

    XN = array([-111.1251, -111.1251, -111.1251, -111.1251, -111.1251,
                -111.1251, -111.1251, -111.1251, -111.1251, -111.1251, -111.1251,
                -111.1251, -111.1251, -111.1251, -111.1251, -111.1251, -111.1251,
                -111.1251, -111.1251, -111.1251])

    YE = array([111.1082, 111.0574, 110.9728, 110.8544, 110.7022, 110.5164,
                110.2968, 110.0437, 109.7570, 109.4369, 109.0834, 108.6968, 108.2770,
                107.8242, 107.3386, 106.8203, 106.2695, 105.6863, 105.0709, 104.4234])

    AZI = array([135.0044, 135.0175, 135.0393, 135.0699, 135.1092, 135.1574,
                135.2143, 135.2802, 135.3549, 135.4385, 135.5312, 135.6329, 135.7437,
                135.8637, 135.9930, 136.1315, 136.2795, 136.4370, 136.6041, 136.7808])

    Dist = array([157.1427, 157.1068, 157.0470, 156.9633, 156.8559, 156.7248,
                156.5700, 156.3918, 156.1902, 155.9654, 155.7176, 155.4470, 155.1538,
                154.8382, 154.5004, 154.1407, 153.7595, 153.3570, 152.9336, 152.4895])

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
        lonArray = array([0.59297447, 0.20873497, 0.44271653, 0.36579662, 0.06680392])
        latArray = array([0.5019297, 0.42174226, 0.23712093, 0.02745615, 0.13316245])
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
        lonArray = array([0.59297447, 0.20873497, 0.44271653, 0.36579662, 0.06680392])
        latArray = array([0.5019297, 0.42174226, 0.23712093, 0.02745615, 0.13316245])
        expected = array([[2.95705149, -2.61243611, -3.0270878 , -2.87840196, -2.42573357],
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
#   lat=range(-20,0)
#   lon=range(120,140)

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
