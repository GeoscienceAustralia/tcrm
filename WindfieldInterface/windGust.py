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


Title: windGust.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2007-02-02
Description: Determine a suitable wind gust factor to apply to
wind fields. Uses the results of Paulsen & Schroeder (2005).
Reference: Paulsen, B.M. & J.L. Schroeder (2005) An examination of
Tropical and Extratropical Gust Factors and the Associated Wind Speed
Histograms, J. Appl. Meteorol., 44, pp. 270-280

SeeAlso:
Constraints:

Version: $Rev: 649 $
ModifiedBy:
ModifiedDate:
Modification:

$Id: windGust.py 649 2011-10-31 05:35:00Z nsummons $
"""

import os, sys, pdb, logging

import numpy

def windGust(size, m=1.39, s=0.24):
    """windGust():
    Return an estimated gust factor. If the random gust factor is
    less than 1.18 (the lowest observation in Paulsen & Schroeder),
    the factor is set to 1.18.
    """
    fac = numpy.random.lognormal(mean=numpy.log(m), sigma=s, size=size)
    #mtools.rnorm(n=1, m=1.59, s=0.24)
    #if fac < 1.18:
    #	fac = 1.18
    return fac
