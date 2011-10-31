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

 Title: smooth.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2009-10-13 1:02:PM
 Description: Provides functions to smooth a 2D field using a Gaussian
              filter
 References: From: http://scipy.org/Cookbook/SignalSmooth

 Version :$Rev: 512 $

 $Id: smooth.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

import numpy
from scipy import signal



def gaussKern(size):
    """
    Calculate a normalised Gaussian kernel to apply as a smoothing
    function.
    Input:
    size (int): the size of the kernel to use (how many points will be
                used in the smoothing operation).
    Output:
    g (array(size,size)): Normalised 2D kernel array for use in
                          convolutions
    """
    size = int(size)
    x,y = numpy.mgrid[-size:size+1,-size:size+1]
    g = numpy.exp(-(x**2/float(size)+y**2/float(size)))
    return g / g.sum()
def smooth(im, n=15):
    """
    Smooth a 2D array im by convolving with a Gaussian kernel of size n
    Input:
    im (2D array): Array of values to be smoothed
    n (int) : number of points to include in the smoothing
    Output:
    improc(2D array): smoothed array (same dimensions as the input array)
    """
    g = gaussKern(n)
    improc = signal.convolve2d(im, g, mode='same', boundary='symm')
    return(improc)
