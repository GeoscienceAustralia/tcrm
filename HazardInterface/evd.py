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


Title: evd.py - calculate extreme value distributions
Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2010-04-23
Description: Calculate the parameters for the GEV distribution, using
             the method of L-moments to estimate the parameters. The
             function will not only fit the distribution parameters,
             but also calculate the return period values for specified
             return periods.

References:
Hosking, J. R. M., 1990: L-moments: Analysis and Estimation of Distributions
using Linear Combinations of Order Statistics. Journal of the Royal
Statistical Society, 52, 1, 105-124.

Version: 335
ModifiedBy: Nicholas Summons
ModifiedDate: 14 July 2010
Modification: Code simplified to return GEVD estimation for single vector of pre-sorted gust wind speeds.
              This allows for clean integration with the hazard interface of TCRM.

Version: $Rev: 758 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2011-03-11 8:56:AM
Modification: Minimum number of records required for calculation made a
              a kwarg

$Id: evd.py 758 2011-11-24 07:30:15Z nsummons $
"""

import os, sys, pdb, logging, traceback
import numpy
__version__ = "$Id: evd.py 758 2011-11-24 07:30:15Z nsummons $"

logger = logging.getLogger()

try:
    # Test if LMOMENTS package is installed
    import lmoments as lmom
except:
    # If not, load ported python version of required functions
    # (note: this runs about half the speed of the LMOMENTS fortran package)
    import Utilities.lmomentFit as lmom
    logger.debug('LMOMENTS package not found - reverting to slower python version of code')

def estimate_EVD(v, years, missingValue=-9999.,minRecords=50,yrspersim=10):

    # Convert to float to prevent integer division & ensure consistent data types for output variables
    yrspersim = numpy.array(yrspersim, dtype='float32')
    missingValue = numpy.array(missingValue, dtype='float32')
    years = numpy.array(years, dtype='float32')

    # Initialise variables:
    loc,scale,shp = [missingValue, missingValue, missingValue]
    w = missingValue*numpy.ones(len(years), dtype='float32')

    if (v.max() > 0.):
        ii = numpy.flatnonzero(v)
        # Only calculate l-moments for those grid points
        # where the values are not all equal, and where
        # there are 50 or more valid values.
        if (v[ii].min() != v[ii].max()) and (len(ii)>=minRecords):
            l1,l2,l3 = lmom.samlmu(v[ii], 3)
            t3 = l3/l2
            if (l2<=0.) or (numpy.abs(t3) >= 1.):
                # Reject points where the second l-moment is negative
                # or the ratio of the third to second is > 1.
                logger.debug("Invalid l-moments")
            else:
                # Parameter estimation returns the location,
                # scale and shape parameters
                xmom = [l1,l2,t3]
                loc,scale,shp = numpy.array(lmom.pelgev(xmom), dtype='float32')
                # We only store the values if the first parameter is
                # finite (i.e. the location parameter is finite)
                if not numpy.isfinite(loc):
                    loc,scale,shp = [missingValue, missingValue, missingValue]

    for i,t in enumerate(years):
        if shp == -9999:
            w[i] = missingValue
        else:
            w[i] = numpy.transpose(loc + (scale/shp)*(1.-numpy.power(-1.*numpy.log(1.-(yrspersim/t)),shp)))
            # Replace any non-finite numbers with the missing value:
            if not numpy.isfinite(w[i]):
                w[i] = missingValue

    return w, loc, scale, shp
