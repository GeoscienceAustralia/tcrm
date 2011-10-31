"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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
"""

#Test of redfearns formula. Tests can be verified at
#
#http://www.cellspark.com/UTM.html
#http://www.ga.gov.au/nmd/geodesy/datums/redfearn_geo_to_grid.jsp

import os, sys
import unittest
from numpy import allclose

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))
from Utilities.lat_long_UTM_conversion import *
from Utilities.files import flStartLog
from redfearn import degminsec2decimal_degrees, decimal_degrees2degminsec

#-------------------------------------------------------------

class TestCase(unittest.TestCase):

    def test_UTM_1(self):
        #latitude:  -37 39' 10.15610" 
        #Longitude: 143 55' 35.38390" 
        #Site Name:    GDA-MGA: (UTM with GRS80 ellipsoid) 
        #Zone:   54    
        #Easting:  758173.797  Northing: 5828674.340 
        #Latitude:   -37  39 ' 10.15610 ''  Longitude: 143  55 ' 35.38390 '' 
        #Grid Convergence:  1  47 ' 19.36 ''  Point Scale: 1.00042107 

        lat = degminsec2decimal_degrees(-37,39,10.15610)
        lon = degminsec2decimal_degrees(143,55,35.38390) 
        assert allclose(lat, -37.65282114)
        assert allclose(lon, 143.9264955)


        zone, easting, northing = LLtoUTM(lat,lon)

        assert zone == 54
        assert allclose(easting, 758173.797)
        assert allclose(northing, 5828674.340)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert allclose(lat,  lat_calced)
        assert allclose(lon, long_calced)


    def test_UTM_2(self):
        #TEST 2

        #Latitude:  -37 57 03.7203
        #Longitude: 144 25 29.5244
        #Zone:   55    
        #Easting:  273741.297  Northing: 5796489.777 
        #Latitude:   -37  57 ' 3.72030 ''  Longitude: 144  25 ' 29.52440 '' 
        #Grid Convergence:  -1  35 ' 3.65 ''  Point Scale: 1.00023056 

        lat = degminsec2decimal_degrees(-37,57,03.7203)
        lon = degminsec2decimal_degrees(144,25,29.5244) 
        #print lat, lon

        zone, easting, northing = LLtoUTM(lat,lon)

        
        assert zone == 55
        assert allclose(easting, 273741.297)
        assert allclose(northing, 5796489.777)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert allclose(lat,  lat_calced)
        assert allclose(lon, long_calced)
        
        
    def test_UTM_3(self):
        #Test 3
        lat = degminsec2decimal_degrees(-60,0,0)
        lon = degminsec2decimal_degrees(130,0,0) 

        zone, easting, northing = LLtoUTM(lat,lon)
        #print zone, easting, northing

        assert zone == 52
        assert allclose(easting, 555776.267)
        assert allclose(northing, 3348167.264)

        Lat, Long = UTMtoLL(northing, easting, zone)

    def test_UTM_4(self):
        #Test 4 (Kobenhavn, Northern hemisphere)
        lat = 55.70248
        dd,mm,ss = decimal_degrees2degminsec(lat)

        lon = 12.58364
        dd,mm,ss = decimal_degrees2degminsec(lon)

        zone, easting, northing = LLtoUTM(lat,lon)
        
        assert zone == 33
        assert allclose(easting, 348157.631)
        assert allclose(northing, 6175612.993) 

        lat_calced, long_calced = UTMtoLL(northing, easting, zone,
                                          isSouthernHemisphere=False) 
        assert allclose(lat,  lat_calced)
        assert allclose(lon, long_calced)

    def test_UTM_5(self):
        #Test 5 (Wollongong)

        lat = degminsec2decimal_degrees(-34,30,0.)
        lon = degminsec2decimal_degrees(150,55,0.) 
        
        zone, easting, northing = LLtoUTM(lat,lon)

        #print zone, easting, northing

        assert zone == 56
        assert allclose(easting, 308728.009)
        assert allclose(northing, 6180432.601)

        lat_calced, long_calced = UTMtoLL(northing, easting, zone) 
        assert allclose(lat,  lat_calced)
        assert allclose(lon, long_calced)
#-------------------------------------------------------------

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestCase,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
