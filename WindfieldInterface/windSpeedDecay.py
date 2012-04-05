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


Title: windSpeedDecay.py - determine the filling rate of a cyclone after
landfall

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-12-01
Description: Calculate the filling rate (decay rate) of the maximum wind
speed in a cyclone following landfall. Based on the paper by Kaplan and
DeMaria (1995). This determines the filling rate of the cyclone, with
the maximum speed converted to a pressure difference by inversion of the
vMax relations.

SeeAlso: (related programs)
Constraints:
Version: $Rev: 810 $

ModifiedBy:
ModifiedDate: yyyy-mm-dd

References:
Kaplan, J. and M. DeMaria, 1995:
    A Simple Empirical Model for Predicting the Decay of Tropical
    Cyclone Winds after Landfall.
    J. Appl. Met., 34, 2499-2512

$Id: windSpeedDecay.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

from numpy import *
import metutils

def windSpeedDecay(vMaxLandfall, inlandDist, inlandTime):
    """
    From Kaplan and DeMaria (1995): A Simple Empirical Model for
    Predicting the Decay of Tropical Cyclone Winds after Landfall,
    J. Appl. Met., 34, 2499-2512
    """
    #pdb.set_trace()
    vMax = metutils.convert(vMaxLandfall, "mps", "kts")
    vBackground = 26.7
    alpha = 0.095
    R = 0.9
    c1 = 0.0109
    d1 = -0.0503
    t0 = 50.

    if inlandDist > 1.0:
        m = c1*inlandTime*(t0 - inlandTime)
        b = d1*inlandTime*(t0 - inlandTime)
        C = m*log(inlandDist) + b
        inlandVmax = vBackground + \
                     (R*vMax - vBackground)*exp(-alpha*inlandTime) - C
    else:
        inlandVmax = vBackground + \
                     (R*vMax - vBackground)*exp(-alpha*inlandTime)

    inlandVmax = metutils.convert(inlandVmax, "kts", "mps")

    return inlandVmax

