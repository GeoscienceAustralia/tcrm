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

 Title: convolve.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2009-10-13 1:02:PM
 Description: Provides functions to convolve two 2D fields
 References: From: http://scipy.org/Cookbook/SignalSmooth

 Version :$Rev: 642 $

 $Id: convolve.py 642 2012-02-21 07:54:04Z nsummons $
"""
import os, sys, pdb, logging

import numpy
from scipy import signal

def getKernel(d,m,r=25.,h=5.):
    """
    Define an appropriate kernel for smoothing the data:

    Input: d - direction - one of ['N','NE','E','SE','S','SW','W','NW']
           m - multiplier type - either 'terrain' or 'shield'
           r - resolution of the dataset (in metres) (default is 25)
           h - nominal height of buildings, for use in evaluating the
               shielding multiplier (default is 5)

    Output: 2-d kernel to apply to the raw multiplier values to generate
            directional multiplier values.

    Example: g = getKernel( 'NE', 'terrain', 25., 5. )
    """

    if d not in ['N','NE','E','SE','S','SW','W','NW']:
        raise ValueError, "No valid direction provided"

    if m=="terrain":
        # We assume that the boundary layer develops over a 1 km distance:
        n = int(1000./r)
        g = numpy.zeros((2*n+1,2*n+1))
        x = numpy.arange(-n,n+1)
        y = numpy.arange(n,-1*n-1,-1)
        xx,yy = numpy.meshgrid(x,numpy.transpose(y))
        rr = r*numpy.sqrt(xx**2+yy**2)
        bear = 90.-(180./numpy.pi)*numpy.arctan2(yy,xx)
        i = numpy.where(bear>180.)
        bear[i] -= 360.

        if d=="S":
            ii = numpy.where((bear==0.) & (rr<=1000.))
        elif d=="SW":
            ii = numpy.where((bear==45.) & (rr<=1000.))
        elif d=="W":
            ii = numpy.where((bear==90.) & (rr<=1000.))
        elif d=="NW":
            ii = numpy.where((bear==135.) & (rr<=1000.))
        elif d=="N":
            ii = numpy.where(((bear==180.) | (bear==-180.)) & (rr<=1000.))
        elif d=="NE":
            ii = numpy.where((bear==-135.) & (rr<=1000.))
        elif d=="E":
            ii = numpy.where((bear==90.) & (rr<=1000.))
        elif d=="SE":
            ii = numpy.where((bear==-45.) & (rr<=1000.))

    elif m=="shield":
        # AS/NZS 1170.2 specifies the shielding distance at 20 times the
        # nominal height of the buildings - assume a 5 m height for
        # residential buildings, hence a 100 m radius.
        n = int(20.*h/r)
        g = numpy.zeros((2*n+1,2*n+1))
        x = numpy.arange(-n,n+1)
        y = numpy.arange(n,-1*n-1,-1)
        xx,yy = numpy.meshgrid(x,numpy.transpose(y))
        rr = r*numpy.sqrt(xx**2+yy**2)
        bear = 90.-(180./numpy.pi)*numpy.arctan2(yy,xx)
        i = numpy.where(bear>180.)
        bear[i] -= 360.
        if d=="S":
            ii = numpy.where((bear>=-22.5) & (bear<=22.5) & (rr<=20.*h))
        elif d=="SW":
            ii = numpy.where((bear>=22.5) & (bear<=67.5) & (rr<=20.*h))
        elif d=="W":
            ii = numpy.where((bear>=67.5) & (bear<=112.5) & (rr<=20.*h))
        elif d=="NW":
            ii = numpy.where((bear>=112.5) & (bear<=157.5) & (rr<=20.*h))
        elif d=="N":
            ii = numpy.where((bear>=157.5) | (bear<=-157.5) & (rr<=20.*h))
        elif d=="NE":
            ii = numpy.where((bear>=-157.5) & (bear<=-112.5) & (rr<=20.*h))
        elif d=="E":
            ii = numpy.where((bear>=-112.5) & (bear<=-67.5) & (rr<=20.*h))
        elif d=="SE":
            ii = numpy.where((bear>=-67.5) & (bear<=-22.5) & (rr<=20.*h))

    g[ii] = 1.
    return g/g.sum()



def convolve(im, dir, m="terrain", res=25.,height=5.):
    """
    Smooth a 2D array im by convolving with a kernel of size n
    Input:
    im (2D array): Array of values to be smoothed
    dir (string): one of 'N','NE','E','SE','S','SW','W','NW' to define the
        direction of the site multiplier to evaluate
    m (string): model type = either "terrain" or "shield"
    res (float): resolution of the input grid dataset
    height (float): nominal height of the buildings to be used in evaluating
        the shielding multiplier.

    Output: 2D array of convolved data (convolved with the appropriate kernel)
            The output array is the same size as the input array, with 
            boundary values set to be filled to a value of 1.0

    Example: arr = convolve(im, 'NW', 'terrain', 25., 5. )
    """
    g = getKernel(dir,m,res,height)
    improc = signal.convolve2d(im, g, mode='same', boundary='fill',fillvalue=1.0)
    return(improc)
