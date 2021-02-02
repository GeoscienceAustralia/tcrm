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


Title: testMetUtils.py

Author: Craig Arthur
Email: craig.arthur@ga.gov.au
CreationDate: 2006-11-14
Description: Unit test for utils/met.py

TODO: Update to test all functions in :mod:`metutils`

"""
import os, sys
import unittest
from numpy import array, arange, pi
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import metutils
from Utilities.files import flStartLog

class TestConvert(NumpyTestCase.NumpyTestCase):
    lat = array(arange(0,-21,-1, 'd'))

    f = array([  0.00000000e+00,  -2.53834962e-06,  -5.07592604e-06,
                -7.61195628e-06,  -1.01456678e-05,  -1.26762889e-05,
                -1.52030487e-05,  -1.77251775e-05,  -2.02419070e-05,
                -2.27524707e-05,  -2.52561037e-05,  -2.77520434e-05,
                -3.02395297e-05,  -3.27178046e-05,  -3.51861134e-05,
                -3.76437042e-05,  -4.00898283e-05,  -4.25237407e-05,
                -4.49447000e-05,  -4.73519686e-05,  -4.97448134e-05])

    def test_mps2kmphr(self):
        """Convert from m/s to km/h"""
        self.assertEqual(metutils.convert(0, "mps", "kph"), 0.)
        self.assertEqual(metutils.convert(1, "mps", "kph"), 3.6)
        self.assertEqual(metutils.convert(5, "mps", "kph"), 18.)
        self.assertEqual(metutils.convert(15, "mps", "kph"), 54.)
        self.assertEqual(metutils.convert(100, "mps", "kph"), 360.)

    def test_km2deg(self):
        """Convert distance in km to distance in degrees"""
        self.assertEqual(metutils.convert(0, "km", "deg"), 0)
        self.assertAlmostEqual(metutils.convert(1, "km", "deg"), 360/(2*pi*6367), 3)
        self.assertAlmostEqual(metutils.convert(2, "km", "deg"), 720/(2*pi*6367), 3)
        self.assertAlmostEqual(metutils.convert(10, "km", "deg"), 3600/(2*pi*6367), 3)

    def test_deg2km(self):
        """Convert distance in degrees to distance in km"""
        self.assertEqual(metutils.convert(0, "deg", "km"), 0)
        self.assertAlmostEqual(metutils.convert(1, "deg", "km"), (1/(360/(2*pi*6367))), 3)
        self.assertAlmostEqual(metutils.convert(2, "deg", "km"), (2/(360/(2*pi*6367))), 3)
        self.assertAlmostEqual(metutils.convert(10, "deg", "km"), (10/(360/(2*pi*6367))), 3)
        
    def test_km2m(self):
        """Convert distance in km to distance in m"""
        self.assertEqual(metutils.convert(0, "km", "m"), 0)
        self.assertAlmostEqual(metutils.convert(1, "km", "m"), 1000)

    def test_m2km(self):
        """Convert distance in m to distance in km"""
        self.assertEqual(metutils.convert(0, "m", "km"), 0)
        self.assertAlmostEqual(metutils.convert(1000., "m", "km"), 1.)
        self.assertAlmostEqual(metutils.convert(10000., "m", "km"), 10.)

    def test_m2nm(self):
        self.assertEqual(metutils.convert(0, "m", "nm"), 0)
        self.assertAlmostEqual(metutils.convert(1000., "m", "nm"), 0.539957)


    def test_hPa2Pa(self):
        """Convert pressure from hPa to Pa"""
        self.assertEqual(metutils.convert(0, "hPa", "Pa"), 0.0)
        self.assertEqual(metutils.convert(1, "hPa", "Pa"), 100.0)
        self.assertEqual(metutils.convert(10, "hPa", "Pa"), 1000.0)
        self.assertEqual(metutils.convert(15, "hPa", "Pa"), 1500.0)
        self.assertEqual(metutils.convert(600, "hPa", "Pa"), 60000.0)

    def test_kgmetre2hPa(self):
        """Convert from Pa to hPa"""
        self.assertEqual(metutils.convert(0, "Pa", "hPa"), 0.0)
        self.assertEqual(metutils.convert(1, "Pa", "hPa"), 0.01)
        self.assertEqual(metutils.convert(100, "Pa", "hPa"), 1.0)
        self.assertEqual(metutils.convert(200, "Pa", "hPa"), 2.0)
        self.assertEqual(metutils.convert(600, "Pa", "hPa"), 6.0)

    def test_celcius2F(self):
        """Convert temperatures in Celcius to Farenheit"""
        self.assertEqual(metutils.convert(0, "C", "F"), 32.0)
        self.assertAlmostEqual(metutils.convert(37.78, "C", "F"), 100, 2)
        self.assertAlmostEqual(metutils.convert(-40, "C", "F"), -40, 2)

    def test_farenheit2C(self):
        """Convert temperatures in Farenheit to Celcius"""
        self.assertEqual(metutils.convert(32, "F", "C"), 0.0)
        self.assertAlmostEqual(metutils.convert(100, "F", "C"), 37.78, 2)
        self.assertAlmostEqual(metutils.convert(-40, "F", "C"), -40, 2)


class TestCoriolis(NumpyTestCase.NumpyTestCase):
    lat = array(arange(0,-21,-1,'f'))

    f = array([  0.00000000e+00,  -2.53834962e-06,  -5.07592604e-06,
                -7.61195628e-06,  -1.01456678e-05,  -1.26762889e-05,
                -1.52030487e-05,  -1.77251775e-05,  -2.02419070e-05,
                -2.27524707e-05,  -2.52561037e-05,  -2.77520434e-05,
                -3.02395297e-05,  -3.27178046e-05,  -3.51861134e-05,
                -3.76437042e-05,  -4.00898283e-05,  -4.25237407e-05,
                -4.49447000e-05,  -4.73519686e-05,  -4.97448134e-05], dtype='float32')


    def test_coriolisArray(self):
        """Test Coriolis with an array of latitudes"""
        self.numpyAssertAlmostEqual(metutils.coriolis(self.lat), self.f)

    def test_coriolisScalar(self):
        """Test Coriolis with a scalar value for latitude"""
        self.assertAlmostEqual(metutils.coriolis(self.lat[-1]), self.f[-1], 3)

class TestWetBulb(NumpyTestCase.NumpyTestCase):
    """
    TODO:
    Test other wet bulb calculations

    """

    def setUp(self):
        self.T = [-5, 0, 5, 10, 15, 20, 25, 30, 35, 40]
        self.Td =  [-10., -5, 0, 5, 10, 15, 20, 20, 20, 20]
        self.prs = [1023.5,]*10
        self.Tw = [-6.4, -1.7, 2.9, 7.57, 12.2, 16.89, 21.59, 23.07, 24.48, 25.8]

    def test_dp2wb_input_error(self):
        """Test ValueError is raised when inputs are non-sensical"""
        self.assertRaises(ValueError, metutils.dewPointToWetBulb, 20, 25, 1010)

    def test_dp2wb(self):
        """Test dewPointToWetBulb conversion"""
        for T, Td, prs, Tw in zip(self.T, self.Td, self.prs, self.Tw):
            self.assertAlmostEqual(metutils.dewPointToWetBulb(T, Td, prs), Tw)


if __name__ == "__main__":
    unittest.main()

