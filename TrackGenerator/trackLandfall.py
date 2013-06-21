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


Title: trackLandfall.py
Author: C. Arthur, craig.arthur@ga.gov.au
CreationDate: 2007-07-04
Description: Determines central pressure as a function of time over
land. Based on Vickery and Twisdale (1995).
References: Vickery, P. J. and L. A. Twisdale (1995). Wind-Field and
            Filling Models for Hurricane Wind-Speed Predictions.
            Journal of Structural Engineering 121(11): 1700-1709.
            Powell, M., G. Soukup, S. Cocke, S. Gulati,
            N. Morriseau-Leroy, S. Hamid, N. Dorst and L. Axe (2005).
            State of Florida hurricane loss projection model:
            Atmospheric science component. Journal of Wind Engineering
            and Industrial Aerodynamics, 93(8): 651-674.


SeeAlso:

Constraints:
Version: 95
ModifiedBy:
ModifiedDate:
Modification:

Version: $Rev: 810 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-01-13
Modification: Documentation improved
$Id: trackLandfall.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

import numpy
from Utilities.grid import SampleGrid
from Utilities import pathLocator
from config import ConfigParser

__version__ = '$Id: trackLandfall.py 810 2012-02-21 07:52:50Z nsummons $'

class LandfallDecay:
    """
    Description: Calculates the decay rate of a tropical cyclone after
    it has made landfall.  Based on the work of on the model of Vickery
    & Twisdale (1995) and employed in the Florida Hurricane Loss model
    of Powell et al. (2005). The model is tuned for Florida conditions,
    and so may require some further tuning for the region of interest.

    Parameters:
    dt - time step of the generated cyclone tracks

    Members:
    dt - time step of the generated cyclone tracks
    tol - time the cyclone has been over land (resets to zero if the
          cyclone moves offshore)
    landLon, landLat - array of longitude & latitude of the landmask
                       data
    landMask - 2D array of 0/1 indicating ocean/land respectively (any
               positive non-zero value can be used to indicate land
               points)

    Methods:
    onLand - Determine if a cyclone centred at (cLon, cLat) is over land
             or not.
    pChange - If the cyclone centre is over land, then this function
              determines the filling as a function of the time elapsed
              since the storm made landfall.

    Internal Methods:
    None

    """
    def __init__(self, configFile, dt):
        """
        Initialise required fields
        """
        self.configFile = configFile

        config = ConfigParser()
        config.read(configFile)

        landMaskFile = config.get('Input', 'LandMask')

        self.logger = logging.getLogger()
        
        self.landMask = SampleGrid(landMaskFile)
        self.tol = 0 # Time over land
        self.dt = dt

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Determines the decay rate of a cyclone after it makes landfall.\
        This function is based on the model of Vickery & Twisdale (1995) and employed\
        in the Florida Hurricane Loss model of Powell et al. (2005). The model is tuned\
        for Florida conditions, and so may require some further tuning for the region of interest.'

    def onLand(self, cLon, cLat):
        """
        Determine if a cyclone centred at (cLon, cLat) is over land or not.
        """

        if self.landMask.sampleGrid(cLon,cLat) > 0.0:
            self.tol += self.dt
            self.logger.debug("Storm centre: %6.2f, %6.2f"%(cLon, cLat))
            self.logger.debug("Time over land: %d hours"%self.tol)
            return True
        else:
            self.tol = 0
            return False

    def pChange(self, pCentre, pEnv):
        """
        If the cyclone centre is over land, then this function
        determines the filling as a function of the time elapsed since
        the storm made landfall.  a0 and a1 control the filling rate of
        the cyclone. Larger values increase the filling rate (i.e.
        cyclones decay faster). An additional random perturbation is
        applied to induce some additional variability.
        Original (V&T1995 values) for a0=0.006, a1=0.00046,
        eps=random.normal(0.0,0.0025)
        """
        deltaP = pEnv - pCentre
        a0 = 0.008
        a1 = 0.0008
        epsilon = numpy.random.normal(0.0, 0.001)
        alpha = a0 + a1*deltaP + epsilon
        deltaPnew = deltaP*numpy.exp(-alpha*self.tol)
        pCentreNew = pEnv - deltaPnew
        return pCentreNew
