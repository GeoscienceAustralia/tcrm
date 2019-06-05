"""

 Version: $Rev: 733 $

 $Id: test_Intersections.py 733 2012-11-14 02:52:01Z carthur $
"""

import os
import sys
import unittest
import logging

import numpy
from . import NumpyTestCase

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))
from Utilities.files import flStartLog

from Utilities import Intersections

class IntersectionsTest(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        """Set up the test environment"""

        xverts = [100., 100., 110., 110., 100.]
        yverts = [-20., -30., -30., -20., -20.]
        self.point1 = Intersections.Point(xverts[0],yverts[0])
        self.point2 = Intersections.Point(xverts[1],yverts[1])
        self.leftPoint = Intersections.Point( 105., -25. )
        self.rightPoint = Intersections.Point( 95., -25. )
        self.circleCentre = Intersections.Point( 105., -25. )
        self.circleRadius0 = 0.0
        self.circleRadius1 = 5.0
        self.circleRadius2 = 7.5
        self.circleRadius3 = 7.0
        self.Crossings = Intersections.Crossings()
        self.vertices = Intersections.convert2vertex( xverts, yverts )


    def test_isLeft(self):
        """Test _isLeft() returns +ve value for a point to the left of the line"""
        result = Intersections._isLeft(self.point1, self.point2, self.leftPoint)
        self.assertTrue( result > 0 )

    def test_isRight(self):
        """Test _isLeft() returns -ve value for a point to the right of the line"""
        result = Intersections._isLeft(self.point1, self.point2, self.rightPoint)
        self.assertTrue( result < 0 )

    def test_isOnLine(self):
        """Test _isLeft() returns zero for a point lying on the line"""
        result = Intersections._isLeft(self.point1, self.point2, self.point2)
        self.assertTrue( result==0 )

    def test_inLand(self):
        """Test inLand() returns correct value for points in/outside vertices"""
        self.assertTrue( Intersections.inLand( self.leftPoint, self.vertices ) )
        self.assertFalse( Intersections.inLand( self.rightPoint, self.vertices ) )

    def test_CircleLine(self):
        """Test Crossings.CircleLine()"""

        v = self.Crossings.CircleLine(self.circleCentre,
                                      self.circleRadius0,
                                      self.point1,
                                      self.point2 )
        self.assertEqual(v.status, "Outside" )

        v = self.Crossings.CircleLine(self.circleCentre,
                                      self.circleRadius1,
                                      self.point1,
                                      self.point2 )
        self.assertEqual(v.status, "Tangent" )

        v = self.Crossings.CircleLine(self.circleCentre,
                                      self.circleRadius2,
                                      self.point1,
                                      self.point2 )
        self.assertEqual(v.status, "Inside" )

        v = self.Crossings.CircleLine(self.circleCentre,
                                      self.circleRadius3,
                                      self.point1,
                                      self.point2 )
        self.assertEqual(v.status, "Intersection" )
        self.assertEqual(len(v.points), 2)
        self.assertAlmostEqual(v.points[0].x, 100.0)
        self.assertAlmostEqual(v.points[0].y,-29.89897948)

    #def test_LineLine(self):
    #    """Test Crossings.LineLine()"""
    #    v = self.Crossings.LineLine(self.Intersections.Point(100.,-20.),
    #                                self.Intersections.Point(100.,-40.),
    #                                self.Intersections.Point(


if __name__ == '__main__':
    flStartLog('','CRITICAL',False)
    testSuite = unittest.TestLoader().loadTestsFromTestCase(IntersectionsTest)
    unittest.TextTestRunner(verbosity=2).run(testSuite)

