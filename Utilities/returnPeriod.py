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

Title: returnPeriod.py - calculate return period based on environmental
       parameters.
Author: Craig Arthur, craig.arthur@ga.gov.au
Created: 2008-01-19
Description: Calculate a return period based on evironmental parameters
             of genesis potential index and maximum potential intensity.
References: Emanuel and Nolan (2004)
        Emanuel (2000)
        Emanuel (1999)
SeeAlso:

Version: $Rev: 642 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-02-26
Modification: Changed f3 to have non-constant value to ensure some
              changes in the frequency of low-probability events.

$Id: returnPeriod.py 642 2012-02-21 07:54:04Z nsummons $
"""

import os, pdb, gc
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

import scipy
import numpy
import my_tool as myutils
#import nctools
import met
import map
__version__ = '$Id: returnPeriod.py 642 2012-02-21 07:54:04Z nsummons $'
def intensityDistribution(v):
    """intensityDistribution(v):
    The distribution of intensity, normalised by the maximum potential
    intensity.  This is based on the statistical analysis by Emanuel
    (2000).
    The final tail of the distribution is chosen to ensure there is some
    change in frequency of low probability events.
    """
    f1 = 237.*(0.75-v)
    f2 = 120.*(0.92-v)
    f3 = 20.0*(1.00-v)
    f = numpy.max([f1, f2, f3], axis=0)
    # Normalise the function:
    f = 1. - f/max(f)
    return f

def returnPeriod(cdf, frequency):
    """returnPeriod(cdf, frequency):
    Calculate the return period based on an annual frequency and a CDF
    of intensity.
    """
    rp = 1.0/((1-cdf[:-1])*frequency)
    return rp


def calculateRetunPeriod(vmax, frequency, years=None):
    """calculateRetunPeriod(vmax, frequency):
    Given a maximum intensity and a frequency of occurrence, calculate
    the return period wind speed at set intervals. The return periods
    can be chosen by the user, or a default set are used.
    """
    if years is None:
        years = [20., 50., 100., 200., 500., 1000., 2000., 5000.]
    vmin = 17./vmax
    v = numpy.arange(vmin, 1, 0.0001)
    cy = intensityDistribution(v)
    rp = returnPeriod(cy, frequency)
    Rp = scipy.interp(years, rp, vmax*v[:-1])
    return Rp
