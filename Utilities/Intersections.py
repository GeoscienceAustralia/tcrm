"""
:mod:`Intersections` -- determine line intersections for polygons and lines
===========================================================================

.. module:: Intersections
    :synopsis: Functions to determine line intersections for polygons and
               line features.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""
import math


def convert2vertex(a1, a2):
    """
    Converts 2 1D arrays into a list of :class:`Points`.

    :param a1: ordinate
    :param a2: abscissa

    :returns: List of :class:`Point` objects.

    """

    result = [Point(x, y) for x, y in zip(a1, a2)]
    return result


def inLand(P, V):
    """
    Test to see if a point is within the list of vertices.

    :param P: :class:`Point` object.
    :param V: List of :class:`Point` objects that represent vertices of
              the shape to be tested. Must be a closed set.

    :returns: ``True`` if the point lies inside the vertices, ``False``
              otherwise.
    :rtype: boolean

    """
    return _cnPnPoly(P, V) and _wnPnPoly(P, V)


def _cnPnPoly(P, V):
    """
    Crossing number test for a point in a polygon.

    :param P: :class:`Point` object.
    :param V: List of :class:`Point` objects that represent vertices of
              the shape to be tested. Must be a closed set.

    :returns: 0 if outside, 1 if inside the polygon.

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
        if (V[i].y <= P.y and V[i + 1].y > P.y) \
                or (V[i].y > P.y and V[i + 1].y <= P.y):   # a downward crossing
            # compute the actual edge-ray intersect x-coordinate
            vt = (P.y - V[i].y) / (V[i + 1].y - V[i].y)
            if P.x < V[i].x + vt * (V[i + 1].x - V[i].x):  # P.x < intersect
                cn += 1   # a valid crossing of y=P.y right of P.x
    return cn % 2 == 1  # 0 if even (out), and 1 if odd (in)


def _wnPnPoly(P, V):
    """
    Winding number test for a point in a polygon. The winding number
    is zero only if the point is outside the polygon.

    :param P: :class:`Point` object.
    :param V: List of :class:`Point` objects that represent vertices of
              the shape to be tested. Must be a closed set.

    :returns: 0 if outside, 1 if inside the polygon.

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
        if V[i].y <= P.y:          # start y <= P.y
            if V[i + 1].y > P.y:      # an upward crossing
                if _isLeft(V[i], V[i + 1], P) > 0:  # P left of edge
                    wn += 1            # have a valid up intersect
            elif V[i + 1].y <= P.y:     # a downward crossing
                if _isLeft(V[i], V[i + 1], P) < 0:  # P right of edge
                    wn -= 1            # have a valid down intersect
    return wn != 0


def _isLeft(P0, P1, P2):
    """
    Tests if a point is Left|On|Right of an infinite line.

    :param P0: One :clas:`Point` on the line.
    :param P1: Second :class:`Point` on the line.
    :param P2: :class:`Point` to be tested.

    :returns: >0 for ``P2`` left of the line through ``P0`` and ``P1``
              =0 for ``P2`` on the line
              <0 for ``P2`` right of the line

    Copyright 2001, softSurfer (www.softsurfer.com)
    This code may be freely used and modified for any purpose
    providing that this copyright notice is included with it.
    SoftSurfer makes no warranty for this code, and cannot be held
    liable for any real or imagined damage resulting from its use.
    Users of this code must verify correctness for their application.
    """
    return (P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y)


class Intersection(object):

    """
    An Intersection object.

    :param str state: Initial state of the :class:`Intersection` (default ``None``).

    Parameters: None
    Members: status (string) and points (list)
    Methods: None
    Internal Methods: None

    """

    def __init__(self, state=None):
        """
        Initialise the members of the object
        """
        self.status = state
        self.points = []


class Crossings(object):

    """
    Determine if a a line intersects some other geometric feature
    (another line, a circle, a polygon).

    """

    def __init__(self):
        """
        Initialise required fields - none required
        """

    def CircleLine(self, c, r, a1, a2):
        """
        Determine the intersection points of a circle and a line segment.

        :param c: :class:`Point` object describing the centre of the circle.
        :param r: Eadius of the circle in map units.
        :param a1: :class:`Point` object describing the start of the line segment.
        :param a2: :class:`Point` object describing the end of the line segment.

        :returns: :class:`Intersection` object with the
                  :attr:`Intersection.status` updated. "Inside" means the line
                  segment is wholly contained within the circle; "Outside"
                  means the line segment is wholly outside the circle;
                  "Intersection" means the line intersects the circle. In this
                  final case, the :attr:`points` attribute of the
                  :class:`Intersection` object is populated with a
                  :class:`Point` object of the location of the intersection.

        """

        a = (a2.x - a1.x) * (a2.x - a1.x) + (a2.y - a1.y) * (a2.y - a1.y)
        b = 2. * ((a2.x - a1.x) * (a1.x - c.x) + (a2.y - a1.y) * (a1.y - c.y))
        cc = c.x * c.x + c.y * c.y + a1.x * a1.x + a1.y * a1.y - 2. * \
            (c.x * a1.x + c.y * a1.y) - r * r
        deter = b * b - 4. * a * cc

        if deter < 0.:
            result = Intersection("Outside")
        elif deter == 0.:
            result = Intersection("Tangent")
        else:
            e = math.sqrt(deter)
            u1 = (-b + e) / (2. * a)
            u2 = (-b - e) / (2. * a)
            if (u1 < 0. or u1 > 1.) and (u2 < 0. or u2 > 1.):
                if (u1 < 0. and u2 < 0.) or (u1 > 1. and u2 > 1.):
                    result = Intersection("Outside")
                else:
                    result = Intersection("Inside")
            else:
                result = Intersection("Intersection")
                if 0. <= u1 and u1 <= 1.:
                    p = self.lerp(a1, a2, u1)
                    result.points.append(Point(p[0], p[1]))
                if 0. <= u2 and u2 <= 1.:
                    p = self.lerp(a1, a2, u2)
                    result.points.append(Point(p[0], p[1]))

        return result

    def CirclePolygon(self, c, r, points):
        """
        Determine the intersection points of a circle and a polygon.
        Uses :func:`Crossings.CircleLine` on each line segment in the
        polygon to determine if the two features intersect.

        :param c: :class:`Point` representing the centre of the circle.
        :param r: Radius of the circle in map units.
        :paran points: Array of :class:`Point` objects representing the
                       polygon.

        :returns: :class:`Intersection` object with updated status and points
                  attributes.
        """

        result = Intersection("No Intersection")
        if isinstance(points, list):
            for p in points:
                for i in range(len(p) - 1):
                    a1 = p[i]
                    a2 = p[i - 1]
                    inter = self.CircleLine(c, r, a1, a2)
                    if inter.status == "Intersection":
                        result.points.append(inter.points)
                        result.status = "Intersection"
        else:
            for i in range(len(points) - 1):
                a1 = points[i]
                a2 = points[i + 1]
                inter = self.CircleLine(c, r, a1, a2)
                if inter.status == "Intersection":
                    result.points.append(inter.points)
                    result.status = "Intersection"
        # if ( len(result.points) > 0. ):

        return result

    def lerp(self, a1, a2, u):
        """
        Linear interpolation between two points:

        :param a1: First :class:`Point`.
        :param a2: Second :class:`Point`.
        :param u: Fractional distance along the line between the two points.

        :return: Coordinates of the interpolated point as a tuple.
        """
        return (a1.x + (a2.x - a1.x) * u, a1.y + (a2.y - a1.y) * u)

    def LineLine(self, a1, a2, b1, b2):
        """
        Determine if two line segments intersect.

        :param a1: Starting :class:`Point` of line 1.
        :param a2: Ending :class:`Point` of line 1.
        :param b1: Starting :class:`Point` of line 2.
        :param b2: Ending :class:`Point` of line 2.

        :returns: :class:`Intersection` object with :attr:`status` set to
                  ``Intersection`` if the lines intersect, ``Coincident``
                  if the lines overlap, ``Parallel`` if the lines are
                  parallel or ``No Intersection`` if the lines do not
                  intersect. If the lines intersect, then the
                  :attr:`points` attribute is set to the location of the
                  intersection.

        """

        ua_t = (b2.x - b1.x) * (a1.y - b1.y) - (b2.y - b1.y) * (a1.x - b1.x)
        ub_t = (a2.x - a1.x) * (a1.y - b1.y) - (a2.y - a1.y) * (a1.x - b1.x)
        u_b = (b2.y - b1.y) * (a2.x - a1.x) - (b2.x - b1.x) * (a2.y - a1.y)

        if u_b != 0:
            ua = ua_t / u_b
            ub = ub_t / u_b

            if 0 <= ua and ua <= 1 and 0 <= ub and ub <= 1:
                result = Intersection("Intersection")
                result.points.append(Point(a1.x + ua * (a2.x - a1.x),
                                           a1.y + ua * (a2.y - a1.y)))
            else:
                result = Intersection("No Intersection")
        else:
            if ua_t == 0 or ub_t == 0:
                result = Intersection("Coincident")
            else:
                result = Intersection("Parallel")

        return result

    def LinePolygon(self, a1, a2, points):
        """
        Determine if a line intersects a polygon.

        :param a1: Starting :class:`Point` of the line.
        :param a2: Ending :class:`Point` of the line.
        :param points: Collection of :class:`Point` objects that
                       represent the vertices of a polygon.


        """
        result = Intersection("No Intersection")
        if isinstance(points, list):
            for p in points:
                for i in range(len(p) - 1):
                    b1 = p[i]
                    b2 = p[i + 1]
                    inter = self.LineLine(a1, a2, b1, b2)
                    result.points.append(inter.points)
        else:
            for i in range(len(points) - 1):
                b1 = points[i]
                b2 = points[i + 1]
                inter = self.LineLine(a1, a2, b1, b2)
                result.points.append(inter.points)

        if len(result.points) > 0:
            result.status = "Intersection"

        return result


class Point(object):

    """
    Class representation for a Point
    Contains x and y members
    Created by Geoff Xu, 2006
    """

    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def getX(self):
        return self.x

    def getY(self):
        return self.y
