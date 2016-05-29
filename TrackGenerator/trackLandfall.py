"""
:mod:`trackLandfall` -- calculate central pressure for over-land TCs
====================================================================

Determines the decay rate of a cyclone after it makes landfall.
Based on the central pressure deficit at landfall and the time over
land, calculate the revised pressure deficit.

References:

* Vickery, P. J. and L. A. Twisdale (1995). Wind-Field and Filling
  Models for Hurricane Wind-Speed Predictions. Journal of Structural
  Engineering 121(11): 1700-1709.
* Powell, M., G. Soukup, S. Cocke, S. Gulati, N. Morriseau-Leroy,
  S. Hamid, N. Dorst and L. Axe (2005). State of Florida hurricane
  loss projection model: Atmospheric science component. Journal of
  Wind Engineering and Industrial Aerodynamics, 93(8): 651-674.

"""

import logging

import numpy as np

from Utilities.grid import SampleGrid
from Utilities.config import ConfigParser

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

class LandfallDecay:
    """
    Description: Calculates the decay rate of a tropical cyclone after
    it has made landfall.  Based on the work of on the model of Vickery
    & Twisdale (1995) and employed in the Florida Hurricane Loss model
    of Powell et al. (2005). The model is tuned for Florida conditions,
    and so may require some further tuning for the region of interest.

    :param str configFile: Configuration file
    :param float dt: time step of the generated cyclone tracks

    Members:
    dt - time step of the generated cyclone tracks
    landMask - 2D array of 0/1 indicating ocean/land respectively (any
    positive non-zero value can be used to indicate land points)

    Methods:
    onLand - Determine if a cyclone centred at (cLon, cLat) is over land
    or not.
    pChange - If the cyclone centre is over land, then this function
    determines the filling as a function of the time elapsed
    since the storm made landfall.

    """

    def __init__(self, configFile, dt):
        """
        Initialise required fields

        """

        self.configFile = configFile

        config = ConfigParser()
        config.read(configFile)

        landMaskFile = config.get('Input', 'LandMask')

        self.landMask = SampleGrid(landMaskFile)
        self.tol = 0 # Time over land
        self.dt = dt

    def onLand(self, cLon, cLat):
        """
        Determine if a cyclone centred at (cLon, cLat) is over land or not.

        :param float cLon: TC longitude
        :param float cLat: TC latitude

        :rtype: boolean
        :return: True if TC is over land, False otherwise

        """

        if self.landMask.sampleGrid(cLon, cLat) > 0.0:
            self.tol += self.dt
            log.debug("Storm centre: %6.2f, %6.2f"%(cLon, cLat))
            log.debug("Time over land: %d hours"%self.tol)
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

        :param float pCentre: TC central pressure (hPa)
        :param float pEnv: TC environmental pressure (hPa)

        :rtype: float
        :return: revised central pressure

        """
        deltaP = pEnv - pCentre
        a0 = 0.008
        a1 = 0.0008
        epsilon = np.random.normal(0.0, 0.001)
        alpha = a0 + a1 * deltaP + epsilon
        deltaPnew = deltaP * np.exp(-alpha * self.tol)
        pCentreNew = pEnv - deltaPnew
        return pCentreNew
