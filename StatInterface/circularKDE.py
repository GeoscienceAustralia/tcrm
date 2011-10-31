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


Title: circularKDE.py
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2007-07-24
Description: Create a probability distribution function for radial-type
data, e.g. bearings. Returns the grid on which the PDF is defined, the
PDF itself and the corresponding CDF. By default, the grid is specified
on [0,2*pi)
Reference:

SeeAlso:

Constraints:

Version: $Rev: 644 $
ModifiedBy:
ModifiedDate:
Modification:

$Id: circularKDE.py 644 2011-10-31 05:32:50Z nsummons $
"""

import numpy
import KPDF
from scipy.special import i0
from math import pi
import stats

def circularKDE(parameters, kdeStep=pi/16.):
    """circularKDE(parameters, kdeStep=pi/16.)
    Create a probability distribution function for radial-type data,
    e.g. bearings. Returns the grid on which the PDF is defined, the
    PDF itself and the corresponding CDF

    By default, the grid is specified on [0,2*pi)

    $Id: circularKDE.py 644 2011-10-31 05:32:50Z nsummons $
    """
    bw = KPDF.UPDFOptimumBandwidth(parameters)
    grid = numpy.arange(0, 2*pi+kdeStep, kdeStep)
    pdf = numpy.empty(len(grid), 'float')
    chi = 1./(2*pi*i0(bw))
    for k in parameters:
        kH = chi*numpy.exp(bw*numpy.cos(grid-k))
        pdf += kH/kH.sum()

    pdf = pdf/len(pdf)
    cy = stats.cdf(grid, pdf)

    return grid, pdf, cy

