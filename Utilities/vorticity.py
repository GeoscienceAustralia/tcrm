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

 Title: vorticity.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 07/09/08 4:40:PM
 Description: Calculates the relative or absolute vorticity of a wind
 field described by U and V components of wind.
 Assumes the array of components is defined on a lat/lon grid.

 Version :$Rev: 686 $

 $Id: vorticity.py 686 2012-03-29 04:24:59Z carthur $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
from numpy import *
import metutils
__version__ = '$Id: vorticity.py 686 2012-03-29 04:24:59Z carthur $'

def relative(u, v, lon, lat):
    """relative(u, v, lon, lat):
    Calculates the relative vorticity (curl(u,v)) of a wind field using
    a basic centred difference scheme. The centred differencing means
    the returned array is reduced in size by 2 elements in each
    dimension.

    zeta((len(lat)-2,len(lon)-2)) = vorticity.relative(u,v,lon,lat)

    """
    dx = numpy.zeros((len(lat), len(lon)-2))
    dy = numpy.zeros((len(lat)-2, len(lon)))
    du = numpy.zeros((len(lat)-2, len(lon)))
    dv = numpy.zeros((len(lat), len(lon)-2))
    zeta = numpy.zeros((len(lat)-2, len(lon)-2))

    for i in xrange(1, len(lon)-1):
        for j in xrange(0, len(lat)):
            dx[j,i-1] = metutils.convert((lon[i+1]-lon[i-1])*numpy.cos(pi*lat[j]/180.), "deg", "m")
            dv[j,i-1] = v[i+1,j]-v[i-1,j]

    for i in xrange(0, len(lon)):
        for j in xrange(1,len(lat)-1):
            dy[j-1,i] = metutils.convert((lat[j+1]-lat[j-1]), "deg", "m")
            du[j-1,i] = u[i,j+1]-u[i,j-1]

    for i in xrange(len(lat)-2):
        for j in xrange(len(lon)-2):
            zeta[i,j] = dv[i,j]/dx[i,j] - du[i,j]/dy[i,j]
    return zeta

def absolute(u, v, lon, lat):
    """absolute(u,v,lon,lat):
    Calculates the absolute vorticity (curl(u,v)) of a wind field using
    a basic centred difference scheme. The centred differencing means
    the returned array is reduced in size by 2 elements in each
    dimension.

    zeta((len(lat)-2,len(lon)-2)) = vorticity.absolute(u,v,lon,lat)
    """
    dx = numpy.zeros((len(lat), len(lon)-2))
    dy = numpy.zeros((len(lat)-2, len(lon)))
    du = numpy.zeros((len(lat)-2, len(lon)))
    dv = numpy.zeros((len(lat), len(lon)-2))
    zeta = numpy.zeros((len(lat)-2, len(lon)-2))

    for i in xrange(1, len(lon)-1):
        for j in xrange(0, len(lat)):
            dx[j,i-1] = metutils.convert((lon[i+1]-lon[i-1])*numpy.cos(pi*lat[j]/180.), "deg", "m")
            dv[j,i-1] = v[i+1,j]-v[i-1,j]

    for i in xrange(0, len(lon)):
        for j in xrange(1, len(lat)-1):
            dy[j-1,i] = metutils.convert((lat[j+1]-lat[j-1]), "deg", "m")
            du[j-1,i] = u[i,j+1]-u[i,j-1]

    for i in xrange(len(lat)-2):
        for j in xrange(len(lon)-2):
            zeta[i,j] = dv[i,j]/dx[i,j] - du[i,j]/dy[i,j] + metutils.coriolis(lat[i+1])

    return zeta
