#!/usr/bin/env python
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

Title: trackutils.py - tools to help with track generation

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-12-06
Description:


Version: $Rev: 686 $
ModifiedBy:
ModifiedDate: yyyy-mm-dd
SeeAlso: (related programs)
Constraints:

$Id: trackutils.py 686 2012-03-29 04:24:59Z carthur $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
import metutils
import maputils

__version__ = '$Id: trackutils.py 686 2012-03-29 04:24:59Z carthur $'

def convert2vertex(a1, a2):
    """
    converts 2 1D arrays into a list of Points
    """
    result = []
    for i in xrange(len(a1)):
        result.append(Point(a1[i], a2[i]))
    return result

def inLand(P, V):
    """
    Test to see if a point is within the list of vertices
    """
    return _cnPnPoly(P, V) and _wnPnPoly(P, V)

def _cnPnPoly(P, V):
    """
    cn_PnPoly(): crossing number test for a point in a polygon
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
    for i in xrange(n):    # edge from V[i] to V[i+1]
        if (((V[i].y <= P.y) and (V[i+1].y > P.y)) \
            or ((V[i].y > P.y) and (V[i+1].y <= P.y))):   # a downward crossing
            # compute the actual edge-ray intersect x-coordinate
            vt = (P.y - V[i].y) / (V[i+1].y - V[i].y)
            if (P.x < V[i].x + vt * (V[i+1].x - V[i].x)): # P.x < intersect
                cn += 1   # a valid crossing of y=P.y right of P.x
    return (cn % 2 == 1)    # 0 if even (out), and 1 if odd (in)

def _wnPnPoly(P, V):
    """
    wn_PnPoly(): winding number test for a point in a polygon
    Input: P = a point,
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
    for i in xrange(n):    # edge from V[i] to V[i+1]
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
    isLeft(): tests if a point is Left|On|Right of an infinite line.
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

class Point:
    """
    Class representation for a Point
    contains x and y members
    created by Geoff Xu, 2006
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def getX():
        return self.x

    def getY():
        return self.y

