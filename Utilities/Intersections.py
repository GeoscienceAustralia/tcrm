#!/usr/bin/env python
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

 Title: Intersections.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2007-03-13
 Description: Functions to determine
 Reference:
 SeeAlso:
 Constraints:

 Version: $Rev: 512 $
 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: Intersections.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, math, numpy, pdb
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
__version__ = '$Id: Intersections.py 512 2011-10-31 07:20:38Z nsummons $'
def convert2vertex(a1, a2):
    """
    Converts 2 1D arrays into a list of Points
    """
    result = []
    for i in range(len(a1)):
        result.append(Point(a1[i], a2[i]))
    return result

def inLand(P, V):
    """
    test to see if a point is within the list of vertices
    """
    return _cnPnPoly(P, V) and _wnPnPoly(P, V)

def _cnPnPoly(P, V):
    """
    _cnPnPoly(): crossing number test for a point in a polygon
        Input:   P = a point,
                V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
        Return:  False if point outside, True if inside
    Copyright 2001, softSurfer (www.softsurfer.com)
    This code may be freely used and modified for any purpose
    providing that this copyright notice is included with it.
    SoftSurfer makes no warranty for this code, and cannot be held
    liable for any real or imagined damage resulting from its use.
    Users of this code must verify correctness for their application.
    """
    cn = 0    # the crossing number counter
    n = len(V) - 1
    # loop through all edges of the polygon
    for i in range(n):    # edge from V[i] to V[i+1]
        if (((V[i].y <= P.y) and (V[i+1].y > P.y)) \
            or ((V[i].y > P.y) and (V[i+1].y <= P.y))):   # a downward crossing
            # compute the actual edge-ray intersect x-coordinate
            vt = (P.y - V[i].y) / (V[i+1].y - V[i].y)
            if (P.x < V[i].x + vt * (V[i+1].x - V[i].x)):  # P.x < intersect
                cn += 1   # a valid crossing of y=P.y right of P.x
    return (cn % 2 == 1)  # 0 if even (out), and 1 if odd (in)

def _wnPnPoly(P, V):
    """
    _wnPnPoly(): winding number test for a point in a polygon
        Input:   P = a point,
                V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
        Return:  False if point outside, True if inside
                (the winding number, wn, =0 only if P is outside V[])
    Copyright 2001, softSurfer (www.softsurfer.com)
    This code may be freely used and modified for any purpose
    providing that this copyright notice is included with it.
    SoftSurfer makes no warranty for this code, and cannot be held
    liable for any real or imagined damage resulting from its use.
    Users of this code must verify correctness for their application.
    """
    wn = 0    # the winding number counter
    n = len(V) - 1
    # loop through all edges of the polygon
    for i in range(n):    # edge from V[i] to V[i+1]
        if (V[i].y <= P.y):          # start y <= P.y
            if (V[i+1].y > P.y):      # an upward crossing
                if (_isLeft(V[i], V[i+1], P) > 0):  # P left of edge
                    wn += 1            # have a valid up intersect
            elif (V[i+1].y <= P.y):     # a downward crossing
                if (_isLeft(V[i], V[i+1], P) < 0):  # P right of edge
                    wn -= 1            # have a valid down intersect
    return (wn != 0)

def _isLeft(P0, P1, P2):
    """
    "isLeft(): tests if a point is Left|On|Right of an infinite line.
    Input:  three points P0, P1, and P2
    Return: >0 for P2 left of the line through P0 and P1
               =0 for P2 on the line
               <0 for P2 right of the line
    Copyright 2001, softSurfer (www.softsurfer.com)
    This code may be freely used and modified for any purpose
    providing that this copyright notice is included with it.
    SoftSurfer makes no warranty for this code, and cannot be held
    liable for any real or imagined damage resulting from its use.
    Users of this code must verify correctness for their application.
    """
    return (P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y)


class Intersection:
    """
    Description:
    Parameters:
    Members:
    Methods:
    Internal Methods:

    """

    def __init__(self, state = None):
        """
        Initialise required fields
        """
        self.status = state
        self.points = []


class Crossings:
    """
    Description:
    Parameters:
    Members:
    Methods:
    Internal Methods:
    """

    def __init__(self):
        """
        Initialise required fields
        """

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return ' '

    def CircleLine(self, c, r, a1, a2):
        """
        Short description:
        """

        a = (a2.x - a1.x) * (a2.x - a1.x) + (a2.y - a1.y) * (a2.y - a1.y)
        b = 2. * ((a2.x - a1.x) * (a1.x - c.x) + (a2.y - a1.y) * (a1.y - c.y))
        cc = c.x*c.x + c.y*c.y + a1.x*a1.x + a1.y*a1.y - 2. * \
            (c.x * a1.x + c.y * a1.y) - r*r
        deter = b*b - 4.*a*cc

        if (deter < 0.):
            result = Intersection("Outside")
        elif (deter == 0.):
            result = Intersection("Tangent")
        else:
            e = math.sqrt(deter)
            u1 = (-b + e) / (2.*a)
            u2 = (-b - e) / (2.*a)
            if ((u1 < 0. or u1 > 1.) and (u2 < 0. or u2 > 1.)):
                if ((u1 < 0. and u2 < 0.) or (u1 > 1. and u2 > 1.)):
                    result = Intersection("Outside")
                else:
                    result = Intersection("Inside")
            else:
                result = Intersection("Intersection")
                if (0. <= u1 and u1 <= 1.):
                    p = self.lerp(a1, a2, u1)
                    result.points.append(Point(p[0], p[1]))
                if (0. <= u2 and u2 <= 1.):
                    p = self.lerp(a1, a2, u2)
                    result.points.append(Point(p[0], p[1]))

        return result

    def CirclePolygon(self, c, r, points):
        """CirclePolygon(c,r):
        Determine the intersection points of a circle and a polygon
        """

        result = Intersection("No Intersection")
        for i in range(len(points)-1):
            a1 = points[i]
            a2 = points[i+1]
            inter = self.CircleLine(c, r, a1, a2)
            if inter.status == "Intersection":
                result.points.append(inter.points)
                result.status = "Intersection"
        #if ( len(result.points) > 0. ):


        return result

    def lerp(self, a1, a2, u):
        """lerp(a1, a2, u)
        Linear interpolation between two points:
        """
        return (a1.x + (a2.x - a1.x) * u, a1.y + (a2.y - a1.y) * u)

    def LineLine(self, a1, a2, b1, b2):
        """LineLine():
        Determine if two lines intersect
        """

        ua_t = (b2.x - b1.x) * (a1.y - b1.y) - (b2.y - b1.y) * (a1.x - b1.x)
        ub_t = (a2.x - a1.x) * (a1.y - b1.y) - (a2.y - a1.y) * (a1.x - b1.x)
        u_b = (b2.y - b1.y) * (a2.x - a1.x) - (b2.x - b1.x) * (a2.y - a1.y)

        if (u_b != 0):
            ua = ua_t / u_b
            ub = ub_t / u_b

            if (0 <= ua and ua <= 1 and 0 <= ub and ub <= 1):
                result = Intersection("Intersection")
                result.points.append(Point(a1.x + ua * (a2.x - a1.x),
                                           a1.y + ua * (a2.y - a1.y)))
            else:
                result = Intersection("No Intersection")
        else:
            if (ua_t == 0 or ub_t == 0):
                result = Intersection("Coincident")
            else:
                result = Intersection("Parallel")

        return result

    def LinePolygon(self, a1, a2, points):
        """LinePolygon():
        Determine if a line intersects a polygon
        """
        result = Intersection("No Intersection")

        for i in range(len(points)-1):
            b1 = points[i]
            b2 = points[i+1]
            inter = self.LineLine(a1, a2, b1, b2)
            result.points.append(inter.points)

        if (len(result.points) > 0):
            result.status = "Intersection"

        return result


class Point:
    """
    Class representation for a Point
    Contains x and y members
    Created by Geoff Xu, 2006
    """
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def getX():
        return self.x

    def getY():
        return self.y

