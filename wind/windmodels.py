"""
This provides classes and methods that calculate the gradient level
and surface (10-m above ground) wind speed, based on a number of
parametric profiles and boundary layer models. These classes and
methods provide the wind field at a single point in time. The complete
wind swath is evaluated in the calling classes.

The :class:`windmodels.WindSpeedModel` classes define the
wind-pressure relations that define the maximum wind speed for a given
central pressure deficit.

The :class:`windmodels.WindProfileModel` classes define parametric
radial profiles of gradient level wind speed around the primary TC
vortex. These do not account for TC motion or interaction with the
surface.

The :class:`windmodels.WindFieldModel` classes define the boundary
layer models implemented, which relate the gradient level vortex to
the surface wind speed, incorporating the wavenumber-1 assymetry due
to storm forward motion, and the effects of (uniform) surface
roughness.

Wind speeds are assumed to represent a 1-minute mean wind
speed. Conversion to other averaging periods is performed by
application of gust factors in the calling classes.

:Note: Not all models are fully implemented - some cannot be fully
       implemented, as the mathematical formulation results in
       discontinuous profiles (e.g. Rankine vortex), for which a
       mathematical representation of the vorticity cannot be
       defined. The vorticity is required for the
       :class:`windmodels.KepertWindFieldModel`. Users should
       carefully select the combinations of wind profiles and wind
       fields to ensure sensible results.

"""

import numpy as np
from math import exp, sqrt
import Utilities.metutils as metutils
from Utilities._vorticity import relative
import logging

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

class WindSpeedModel(object):

    """
    Abstract wind speed model.
    """

    def __init__(self, windProfileModel):
        """
        Abstract wind speed model.
        """
        self.profile = windProfileModel

    @property
    def eP(self):
        """
        Environment pressure.
        """
        eP = self.profile.eP
        if eP < 10000:
            eP = metutils.convert(eP, 'hPa', 'Pa')
        return eP

    @property
    def cP(self):
        """
        Current pressure.
        """
        cP = self.profile.cP
        if cP < 10000:
            cP = metutils.convert(cP, 'hPa', 'Pa')
        return cP

    @property
    def dP(self):
        """
        Pressure difference.
        """
        return self.eP - self.cP

    def maximum(self):
        """
        Maximum wind speed.
        """
        raise NotImplementedError


class WilloughbyWindSpeed(WindSpeedModel):

    """
    Willoughby & Rahn (2004), Parametric Representation of the Primary
    Hurricane Vortex. Part I: Observations and Evaluation of the
    Holland (1980) Model.  Mon. Wea. Rev., 132, 3033-3048
    """

    def maximum(self):
        return 0.6252 * sqrt(self.dP)


class HollandWindSpeed(WindSpeedModel):

    """
    .. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

    Holland (1980), An Analytic Model of the Wind and Pressure Profiles
    in Hurricanes. Mon. Wea. Rev, 108, 1212-1218 Density of air is
    assumed to be 1.15 kg/m^3.  |beta| is assumed to be 1.3. Other values
    can be specified.  Gradient level wind (assumed maximum).

    """

    def maximum(self):
        beta = self.profile.beta
        rho = 1.15
        return np.sqrt(beta * self.dP / (np.exp(1) * rho))

class NewHollandWindSpeed(WindSpeedModel):
    """"
    .. |beta|    unicode:: U+003B2 .. GREEK SMALL LETTER BETA
    
    Holland et al. (2010)
    """

    def maximum(self):
        beta = self.profile.beta
        rho = 1.15
        return np.sqrt(beta * self.dP/(rho * np.exp(1)))

class AtkinsonWindSpeed(WindSpeedModel):

    """
    Atkinson and Holliday (1977), *Tropical Cyclone Minimum Sea
    Level Pressure / Maximum Sustained Wind Relationship for
    the Western North Pacific* . Mon. Wea. Rev., **105**, 421-427
    Maximum 10m, 1-minute wind speed. Uses ``pEnv`` as 1010 hPa.
    """

    def maximum(self):
        cP = metutils.convert(self.cP, 'Pa', 'hPa')
        return 3.04 * pow(1010.0 - cP, 0.644)


class WindProfileModel(object):

    """
    The base wind profile model.

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax, windSpeedModel):
        self.rho = 1.15  # density of air
        self.grid = grid
        self.lat = grid.cLat
        self.lon = grid.cLon
        self.eP = eP
        self.cP = cP
        self.rMax = rMax
        self.speed = windSpeedModel(self)
        self.f = metutils.coriolis(grid.cLat)
        self.vMax_ = None

    @property
    def dP(self):
        """
        Pressure difference.
        """
        return self.eP - self.cP

    @property
    def vMax(self):
        """
        Maximum wind speed.
        """
        if self.vMax_:
            return self.vMax_
        else:
            return self.speed.maximum()

    @vMax.setter
    def vMax(self, value):
        """
        Set the maximum wind speed.

        :param float value: The maximum wind speed value to set.
        """
        self.vMax_ = value

    def velocity(self, grid):
        """
        Calculate velocity as a function of radial distance `R`.
        Represents the velocity of teh gradient level vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre.

        :returns: Array of gradient level wind speed.
        :rtype: :class:`numpy.ndarray`

        """
        raise NotImplementedError

    def vorticity(self, grid):
        """
        Calculate the vorticity associated with the (gradient level)
        vortex at radius `R`.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre.

        :returns: Array of gradient level (relative) vorticity.
        :rtype: :class:`numpy.ndarray`

        """
        raise NotImplementedError


class JelesnianskiWindProfile(WindProfileModel):

    """
    Jelesnianski model of the wind profile.

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (m).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax,
                 windSpeedModel=WilloughbyWindSpeed):
        rMax = rMax
        WindProfileModel.__init__(self, grid, eP, cP, rMax,
                                  windSpeedModel)

    def velocity(self, grid):
        """
        Calculate velocity as a function of radial distance.
        Represents the velocity of teh gradient level vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre (metres).

        :returns: Array of gradient level wind speed.
        :rtype: :class:`numpy.ndarray`

        """

        V = 2 * self.vMax * self.rMax * grid.R / (self.rMax ** 2 + grid.R ** 2)
        V = np.sign(self.f) * V
        return V

    def vorticity(self, grid):
        """
        Calculate the vorticity associated with the (gradient level)
        vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre (metres).

        :returns: Array of gradient level (relative) vorticity.
        :rtype: :class:`numpy.ndarray`

        """
        Z = (np.sign(self.f) * 2 * self.vMax * self.rMax / (self.rMax
             ** 2 + grid.R ** 2) + np.sign(self.f) * 2 * self.vMax *
             self.rMax * (self.rMax ** 2 - grid.R ** 2) /
             (self.rMax ** 2 + grid.R ** 2) ** 2)
        return Z


class HollandWindProfile(WindProfileModel):

    """
    .. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

    Holland profile. For `r < rMax`, we reset the wind field to a
    cubic profile to avoid the barotropic instability mentioned in
    Kepert & Wang (2001).

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (m).
    :param float beta: |beta| parameter.
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax, beta,
                 windSpeedModel=HollandWindSpeed):
        WindProfileModel.__init__(self, grid, eP, cP, rMax,
                                  windSpeedModel)
        self.beta = beta

    def secondDerivative(self):
        """
        Second derivative of profile at rMax.
        """

        beta = self.beta
        dP = self.dP
        rho = self.rho
        f = self.f
        rMax = self.rMax

        E = exp(1)
        
        d2Vm = ((beta * dP * (-4 * beta ** 3 * dP / rho -
                (-2 + beta ** 2) * E * (np.abs(f) * rMax) ** 2)) /
                (E * rho * sqrt((4 * beta * dP) / (E * rho)
                 + (f * rMax) ** 2) * (4 * beta * dP * rMax ** 2 / rho
                 + E * (f * rMax ** 2) ** 2)))

        try:
            assert d2Vm < 0.0
        except AssertionError:
            log.critical("Pressure deficit: %f, RMW: %f" % (dP, rMax))
            raise

        return d2Vm

    def firstDerivative(self):
        """
        First derivative of profile at rMax
        """
        beta = self.beta
        dP = self.dP
        rho = self.rho
        f = self.f
        rMax = self.rMax

        E = exp(1)
        
        dVm = (-np.abs(f)/2 + (E*(f**2)*rMax*np.sqrt((4*beta*dP/rho)/E + \
                                                     (f*rMax)**2))/ \
                  (2*(4*beta*dP/rho + E*(f*rMax)**2)))
        return dVm

    def velocity(self, grid):
        """
        Calculate velocity as a function of radial distance.
        Represents the velocity of teh gradient level vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre (metres).

        :returns: Array of gradient level wind speed.
        :rtype: :class:`numpy.ndarray`

        """

        d2Vm = self.secondDerivative()
        dVm = self.firstDerivative()
        aa = ((d2Vm / 2. - (dVm - self.vMax /self.rMax) 
               / self.rMax) / self.rMax)
        bb = (d2Vm - 6 * aa * self.rMax) / 2.
        cc = dVm - 3 * aa * self.rMax ** 2 - 2 * bb * self.rMax
        delta = (self.rMax / grid.R) ** self.beta
        edelta = np.exp(-delta)

        V = (np.sqrt((self.dP * self.beta / self.rho)
             * delta * edelta + (grid.R * self.f / 2.) ** 2) - grid.R *
             np.abs(self.f) / 2.)

        icore = np.where(grid.R <= self.rMax)
        V[icore] = (grid.R[icore] * (grid.R[icore] * (grid.R[icore] * aa + bb) + cc))
        V = np.sign(self.f) * V
        return V

    def vorticity(self, grid):
        """
        Calculate the vorticity associated with the (gradient level)
        vortex. 
        
        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre (metres).

        :returns: Array of gradient level (relative) vorticity.
        :rtype: :class:`numpy.ndarray`

        """

        beta = self.beta
        delta = (self.rMax / grid.R) ** beta
        edelta = np.exp(-delta)

        Z = np.abs(self.f) +\
            (beta**2*self.dP*(delta**2)*edelta/\
             (2*self.rho*grid.R) - beta**2*self.dP*delta*edelta/\
             (2*self.rho*grid.R) + grid.R*self.f**2/4) / \
            np.sqrt(beta*self.dP*delta*edelta/self.rho + \
                    (grid.R*self.f/2)**2) + \
            (np.sqrt(beta*self.dP*delta*edelta/self.rho + \
                    (grid.R*self.f/2)**2))/grid.R

        # Calculate first and second derivatives at R = Rmax:
        d2Vm = self.secondDerivative()
        dVm = self.firstDerivative()
        aa = ((d2Vm / 2 - (dVm - self.vMax /
              self.rMax) / self.rMax) / self.rMax)
        bb = (d2Vm - 6 * aa * self.rMax) / 2
        cc = dVm - 3 * aa * self.rMax ** 2 - 2 * bb * self.rMax

        icore = np.where(grid.R <= self.rMax)
        Z[icore] = grid.R[icore] * (grid.R[icore] * 4 * aa + 3 * bb) + 2 * cc
        Z = np.sign(self.f) * Z
        return Z


class WilloughbyWindProfile(HollandWindProfile):

    """
    .. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

    The Willoughby & Rahn (2004) relation, which makes |beta| a function
    of Vmax, rMax and latitude. We use Willoughby & Rahn's (2004)
    relation for Vmax *only*.  This determines the |beta| parameter then
    calls Holland (which means the profile is cubic within Rmax) to
    calculate the wind profile.  The |beta| term calculation is based on
    Atlantic and Eastern Pacific cyclone data, not Australian data.

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax,
                 windSpeedModel=WilloughbyWindSpeed):
        HollandWindProfile.__init__(self, grid, eP, cP, rMax, 1.0,
                                    windSpeedModel)
        self.beta = (1.0036 + 0.0173 * self.vMax - 0.313 * np.log(rMax/1000.)
                     + 0.0087 * np.abs(grid.cLat))
        self.speed = HollandWindSpeed(self)


class RankineWindProfile(WindProfileModel):

    """
    Rankine vortex profile. Vmax determined by dp using, by default,
    the Willoughby & Rahn method.

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax,
                 windSpeedModel=WilloughbyWindSpeed):
        WindProfileModel.__init__(self, grid, eP, cP, rMax,
                                  windSpeedModel)
        self.alpha = 0.5

    def velocity(self, grid):
        """
        Calculate velocity as a function of radial distance.
        Represents the velocity of teh gradient level vortex.

        This includes an assumption about the shape of the profile
        outside Rmax. The literature indicates 0.4 < alpha < 0.6
        (e.g. see Holland, 1980).

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre.

        :returns: Array of gradient level wind speed.
        :rtype: :class:`numpy.ndarray`

        """

        V = self.vMax * (self.rMax / grid.R) ** self.alpha
        icore = np.where(grid.R <= self.rMax)
        V[icore] = self.vMax * (grid.R[icore] / self.rMax)
        V = np.sign(self.f) * V
        return V

    def vorticity(self, grid):
        """
        Calculate the vorticity associated with the (gradient level)
        vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre.

        :returns: Array of gradient level (relative) vorticity.
        :rtype: :class:`numpy.ndarray`

        """
        V = self.velocity(grid)
        Z = V/grid.R - self.alpha*self.vMax*self.rMax/(grid.R**(self.alpha+1))
        icore = np.where(grid.R <= self.rMax)
        Z[icore] = np.sign(self.f) * ((self.vMax * (grid.R[icore] /
                    self.rMax) + self.vMax / self.rMax))

        return Z


class SchloemerWindProfile(HollandWindProfile):

    """
    .. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

    Schloemer's (1954) is the same as the Holland relation with
    |beta| = 1

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax):
        HollandWindProfile.__init__(self, grid, eP, cP, rMax, 1.0)


class DoubleHollandWindProfile(WindProfileModel):

    """
    .. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA

    McConochie *et al*'s double Holland vortex model (based on Cardone *et
    al*, 1994).  This application is the Coral Sea adaptation of the
    double vortex model (it can also be used for concentric eye-wall
    configurations).

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param float beta1: |beta| parameter for the primary vortex.
    :param float beta2: |beta| parameter for the secondary vortex.
    :param float rmax2: Optional. Radius to the secondary wind
                        maximum (default=250 km).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax, beta1, beta2, rMax2=150.,
                 windSpeedModel=HollandWindSpeed):
        WindProfileModel.__init__(self, grid, eP, cP, rMax,
                                  windSpeedModel)

        # Scale dp2 if dP is less than 800 Pa

        if self.dP < 1500.:
            self.dp2 = ((self.dP / 1500.) * (800. + (self.dP - 800.) /
                        2000.))
        else:
            self.dp2 = 800. + (self.dP - 800.) / 2000.

        self.dp1 = self.dP - self.dp2

        self.beta = beta1
        self.beta1 = beta1
        self.beta2 = beta2
        self.rMax2 = rMax2

        if self.beta1 is None:
            self.beta1 = 7.3 - self.cP / 16000.

        if self.beta2 is None:
            self.beta2 = 7.2 - self.cP / 16000.

    def secondDerivative(self):
        """
        Second derivative of the profile.
        """
        beta1 = self.beta1
        beta2 = self.beta2
        rMax1 = self.rMax
        rMax2 = self.rMax2
        rho = self.rho
        dp1 = self.dp1
        dp2 = self.dp2
        f = self.f

        E = exp(1)
        nu = pow((rMax2 / rMax1), beta2)

        d2Vm = (-1 /
                (8 *
                 (4 * beta1 * dp1 / (rho * E) +
                  (4 * beta2 * dp2 / rho) * nu * exp(-nu) +
                  (rMax1 * f) ** 2) ** 1.5)
                * (-(4 * (beta1 ** 2) * dp1 / (rho * rMax1 * E)) +
                    (4 * (beta1 ** 2) * dp1 / (rho * rMax1 * E)) -
                    (4 * (beta2 ** 2) * dp2 / rho) *
                    (nu / rMax1) * exp(-nu)
                    + (4 * (beta2 ** 2) * dp2 / rho) *
                    ((nu ** 2) / rMax1) * exp(-nu)
                    + 2 * rMax1 * f ** 2) ** 2
                + 1 / (4 * sqrt((4 * beta1 * dp1 / (rho * E)) +
                                (4 * beta2 * dp2 / rho) * nu * 2 +
                                exp(-nu) + (rMax1 * f) ** 2))
                * ((4 * (beta1 ** 3) * dp1 / (rho * (rMax1 ** 2) * E))
                   + (4 * (beta1 ** 2) * dp1 / (rho * (rMax1 ** 2) * E))
                   - (12 * (beta1 ** 3) * dp1 / (rho * (rMax1 ** 2) * E))
                   - (4 * (beta1 ** 2) * dp1 / (rho * (rMax1 ** 2) * E))
                   + (4 * (beta1 ** 3) * dp1 / (rho * (rMax1 ** 2) * E))
                   + (4 * (beta2 ** 3) * dp2 / rho) *
                     (nu / (rMax1 ** 2)) * exp(-nu)
                   + (4 * (beta2 ** 2) * dp2 / rho) *
                     (nu / (rMax1 ** 2)) * exp(-nu)
                   - (12 * (beta2 ** 3) * dp2 / rho) *
                     (nu ** 2) / (rMax1 ** 2) * exp(-nu)
                   - (4 * (beta2 ** 2) * dp2 / rho) *
                     (nu ** 2) / (rMax1 ** 2) * exp(-nu)
                   + (4 * (beta2 ** 3) * dp2 / rho) *
                     (nu ** 3) / (rMax1 ** 2) * exp(-nu)
                   + 2 * f ** 2))

        assert d2Vm < 0.0

        return d2Vm

    def velocity(self, grid):
        """
        Calculate velocity as a function of radial distance.
        Represents the velocity of teh gradient level vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre.

        :returns: Array of gradient level wind speed.
        :rtype: :class:`numpy.ndarray`

        """
        rMax = self.rMax
        rMax2 = self.rMax2

        # Scale dp2 if dP is less than 800 Pa

        if self.dP < 1500.:
            dp2 = (self.dP / 1500.) * (800. + (self.dP - 800.) / 2000.)
        else:
            dp2 = 800. + (self.dP - 800.) / 2000.

        dp1 = self.dP - dp2

        # The two gradient wind components

        mu = (rMax / grid.R) ** self.beta1
        nu = (rMax2 / grid.R) ** self.beta2
        emu = np.exp(-mu)
        enu = np.exp(-nu)

        gradientV1 = (self.beta1 * dp1 / self.rho) * mu * emu
        gradientV2 = (self.beta2 * dp2 / self.rho) * nu * enu

        V = (np.sign(self.f) * np.sqrt(gradientV1 + gradientV2 + (grid.R *
             self.f / 2.) ** 2) - grid.R * np.abs(self.f) / 2.)

        vMax = np.abs(V).max()

        d2Vm = self.secondDerivative()
        aa = (d2Vm / 2. - (-vMax / rMax) / rMax) / rMax
        bb = (d2Vm - 6 * aa * rMax) / 2.
        cc = -3 * aa * rMax ** 2 - 2 * bb * rMax

        # Replace all values within rMax of the storm centre with the
        # cubic profile to eliminate barotropic instability

        if self.dP >= 1500.:
            icore = np.where(grid.R <= rMax)
            V[icore] = (np.sign(self.f) * grid.R[icore] * (grid.R[icore] *
                        (grid.R[icore] * aa + bb) + cc))

        return V

    def vorticity(self, grid):
        """
        Calculate the vorticity associated with the (gradient level)
        vortex.

        :param R: :class:`numpy.ndarray` of distance of grid from
                  the TC centre.

        :returns: Array of gradient level (relative) vorticity.
        :rtype: :class:`numpy.ndarray`

        """

        rm = self.rMax
        rm2 = self.rMax2

        # Scale dp2 if dP is less than 1500 Pa:
        if self.dP < 1500.:
            dp2 = (self.dP / 1500.) * (800. + (self.dP - 800.) / 2000.)
        else:
            dp2 = 800. + (self.dP - 800.) / 2000.

        dp1 = self.dP - dp2

        chi = self.beta1 * dp1 / self.rho
        psi = self.beta2 * dp2 / self.rho

        delta = (self.rMax / grid.R) ** self.beta1
        gamma = (self.rMax2 / grid.R) ** self.beta2
        edelta = np.exp(-delta)
        egamma = np.exp(-gamma)

        # Derivatives:

        ddelta = (-self.beta1 * (self.rMax ** self.beta1) / (grid.R **
                  (self.beta1 + 1)))
        dgamma = (-self.beta2 * (self.rMax2 ** self.beta2) / (grid.R **
                  (self.beta2 + 1)))

        Z = (np.sign(self.f) * np.sqrt(chi * delta * edelta + psi *
             gamma * egamma + (self.f * grid.R / 2) ** 2) / grid.R -
             np.abs(self.f) + (1 / 2) *
             (chi * ddelta * edelta * (1 - delta) +
              psi * dgamma * egamma * (1 - gamma) +
              grid.R * self.f ** 2) /
             np.sqrt(chi * delta * edelta + psi * gamma *
                     egamma + (self.f * grid.R / 2) ** 2))

        d2Vm = self.secondDerivative()
        aa = ((d2Vm / 2.0 - (-1.0 * np.sign(self.f) * self.vMax /
                             self.rMax) / self.rMax) / self.rMax)
        bb = (d2Vm - 6.0 * aa * self.rMax) / 2.0
        cc = -3.0 * aa * self.rMax ** 2.0 - 2.0 * bb * self.rMax

        if self.dP >= 1500.:
            icore = np.where(grid.R <= self.rMax)
            Z[icore] = (grid.R[icore] * (grid.R[icore] * 4.0 * aa + 3.0 * bb) +
                        2.0 * cc)

        return Z


class PowellWindProfile(HollandWindProfile):

    """
    Powell *et al*, 2005. Another definition of the beta parameter
    inserted into the Holland model.  Unlike Willoughby and Rahn's
    model, there is no reliance on vMax.  Powell *et al*. also included a
    small random term, but since the beta value is also used in the
    vorticity calculation, we need to ensure the values used in this
    function and the corresponding vorticity function match.

    :param float lat: Latitude of TC centre.
    :param float lon: Longitude of TC centre.
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param windSpeedModel: A maximum wind speed model to apply.
    :type  windSpeedModel: :class:`windmodels.WindSpeedModel` instance.

    """

    def __init__(self, grid, eP, cP, rMax):
        beta = 1.881093 - 0.010917 * np.abs(grid.cLat) - 0.005567 * rMax/1000.
        if beta < 0.8:
            beta = 0.8
        if beta > 2.2:
            beta = 2.2
        HollandWindProfile.__init__(self, grid, eP, cP, rMax, beta,
                                    windSpeedModel=HollandWindSpeed)


class NewHollandWindProfile(WindProfileModel):

    """
    Holland et al. 2010.  In this version, the exponent is allowed to
    vary linearly outside the radius of maximum wind. i.e. rather than
    take the sqare root, the exponent varies around 0.5.  

    This uses a vorticity field calculated using teh curl of the vector
    wind field, rather than an analytic solution.

    :param grid: :class:`ModelGrid` instance defining the simulation grid.
    :type grid: :class:`ModelGrid`
    :param float eP: environmental pressure (hPa).
    :param float cP: centrral pressure of the TC (hPa).
    :param float rMax: Radius to maximum wind (km).
    :param float vFm: Forward speed of the cyclone (m/s).
    :param float dcP: Rate of change of central pressure (Pa/hr).
    :param float rGale: Radius of gale force winds (default=150 km).

    """

    def __init__(self, grid, eP, cP, rMax, vFm, dcP, rGale=150.):
        WindProfileModel.__init__(self, grid, eP, cP, rMax,
                                  windSpeedModel=NewHollandWindSpeed)
        self.rGale = rGale
        self.dcP = dcP
        self.vFm = vFm

    def velocity(self, grid):
        """
        Calculate velocity as a function of radial distance.
        Represents the velocity of the gradient level vortex.

        :param grid: :class:`ModelGrid` instance defining the simulation grid.
        :returns: Array of gradient level wind speed.
        :rtype: :class:`numpy.ndarray`

        NOTE: dP and dcP attributes are nominally held in units of 
        Pa, but the relation for `Bs` uses hPa.

        """

        Bs = (-4.4e-5 * (self.dP/100.) ** 2. + 0.01 * (self.dP/100.) 
              - 0.014 * np.abs(self.lat) + 1.0
              + 0.15 * self.vFm ** (0.6 * (1 - self.dP / 21500.))
              + 0.00015 * self.dcP)
        self.beta = Bs

        deltag = np.power(self.rMax / self.rGale, Bs)
        edeltag = np.exp(1 - deltag)
        rgterm = deltag * edeltag
        xn = np.log(17./np.abs(self.vMax)) / np.log(rgterm)
        xx = 0.5 * np.ones(grid.R.shape)

        i = np.where(grid.R > self.rMax)
        if xn <= 0.5:
            xx[i] = 0.5
        else:
            xx[i] = (0.5 + (grid.R[i] - self.rMax) * (xn - 0.5) / (self.rGale -
                                                              self.rMax))

        delta = (self.rMax / grid.R) ** Bs
        edelta = np.exp(1-delta)
        rterm = delta * edelta
        rterm = np.where(rterm > 1, 1, rterm)
        V = self.vMax * np.power(rterm, xx)

        #V = (np.sign(self.f) * np.power((self.dP * Bs / self.rho) *
        #     delta * edelta, xx))

        return V

    def vorticity(self, grid): 
        """
        Calculate the vorticity associated with the (gradient level)
        vortex.

        :param grid: :class:`ModelGrid` instance defining the simulation grid.

        :returns: Array of gradient level (relative) vorticity.
        :rtype: :class:`numpy.ndarray`

        """
        V = self.velocity(grid)
        UU = V * np.sin(-grid.theta)
        VV = V * np.cos(-grid.theta)
        
        # Need to determine lon/lat of the grid, based on 
        # storm centre, grid resolution and margin:
        Z = relative(UU, VV, grid.lonGrid, grid.latGrid)
        return Z
        


class WindFieldModel(object):

    """
    Wind field (boundary layer) models. These define the boundary
    layer models implemented, which relate the gradient level vortex
    to the surface wind speed, incorporating the wavenumber-1
    assymetry due to storm forward motion, and the effects of
    (uniform) surface roughness.

    :param windProfileModel: A `wind.WindProfileModel` instance.

    """

    def __init__(self, windProfileModel):
        self.profile = windProfileModel
        self.V = None
        self.Z = None

    @property
    def rMax(self):
        """
        Helper property to return the maximum radius from the
        wind profile.
        """
        return self.profile.rMax

    @property
    def f(self):
        """
        Helper property to return the coriolis force from the
        wind profile.
        """
        return self.profile.f

    def velocity(self, grid):
        """
        Helper property to return the wind velocity at radiuses `R`
        from the wind profile or the precalculated attribute.
        """
        if self.V is None:
            return self.profile.velocity(grid)
        else:
            return self.V

    def vorticity(self, grid):
        """
        Helper property to return the wind vorticity at radiuses `R`
        from the wind profile or the precalculated attribute.
        """
        if self.Z is None:
            return self.profile.vorticity(grid)
        else:
            return self.Z

    def field(self, grid, vFm, thetaFm, thetaMax=0.):
        """
        The wind field.
        """
        raise NotImplementedError


class HubbertWindField(WindFieldModel):

    """
    Hubbert, G.D., G.J. Holland, L.M. Leslie and M.J. Manton, 1991:
    A Real-Time System for Forecasting Tropical Cyclone Storm Surges.
    *Weather and Forecasting*, **6**, 86-97

    :param R: Distance from the storm centre to the grid (km).
    :type  R: :class:`numpy.ndarray`
    :param lam: Direction (geographic bearing, positive clockwise)
                from storm centre to the grid.
    :type  lam: :class:`numpy.ndarray`
    :param float vFm: Foward speed of the storm (m/s).
    :param float thetaFm: Forward direction of the storm (geographic
                          bearing, positive clockwise).
    :param float thetaMax: Bearing of the location of the maximum
                           wind speed, relative to the direction of
                           motion.
    """

    def field(self, grid, vFm, thetaFm, thetaMax=0.):
        V = self.velocity(grid)
        thetaFm = np.pi/2.0 - thetaFm
        Km = .70
        inflow = -np.sign(self.f) * 25. * np.ones(np.shape(grid.R))
        core = np.where(grid.R < self.rMax)
        inflow[core] = 0
        inflow = inflow * np.pi / 180

        thetaMaxAbsolute = np.array(thetaFm) + thetaMax * np.pi / 180.
        asym = vFm * np.cos(thetaMaxAbsolute - grid.theta + np.pi)
        Vsf = Km * V + asym
        phi = inflow - grid.theta

        # Additional factor to convert to 1-munite sustained wind speed:
        Ux = Vsf * np.sin(phi) * 1.069
        Vy = Vsf * np.cos(phi) * 1.069

        return Ux, Vy


class McConochieWindField(WindFieldModel):
    """
    McConochie, J.D., T.A. Hardy and L.B. Mason, 2004:
    Modelling tropical cyclone over-water wind and pressure fields.
    Ocean Engineering, 31, 1757-1782.

    It appears this paper used 10-minute mean wind speeds for calibration.
    Therefore, we convert to return a 1-minute sustained wind speed equivalent.

    """

    def field(self, grid, vFm, thetaFm, thetaMax=0.):
        """
        :param R: Distance from the storm centre to the grid (km).
        :type  R: :class:`numpy.ndarray`
        :param lam: Direction (geographic bearing, positive clockwise)
                    from storm centre to the grid.
        :type  lam: :class:`numpy.ndarray`
        :param float vFm: Foward speed of the storm (m/s).
        :param float thetaFm: Forward direction of the storm (geographic
                              bearing, positive clockwise).
        :param float thetaMax: Bearing of the location of the maximum
                               wind speed, relative to the direction of
                               motion.
        """
        V = self.velocity(grid)
        thetaFm = np.pi/2.0 - thetaFm
        inflow = 25. * np.ones(np.shape(grid.R))
        mid = np.where(grid.R < 1.2 * self.rMax)
        inflow[mid] = 10. + 75. * (grid.R[mid] / self.rMax - 1.)
        inner = np.where(grid.R < self.rMax)
        inflow[inner] = 10. * grid.R[inner] / self.rMax
        inflow = -np.sign(self.f) * inflow * np.pi / 180.

        thetaMaxAbsolute = np.array(thetaFm) + thetaMax * np.pi / 180.
        phi = inflow - grid.theta

        asym = (0.5 * (1. + np.cos(thetaMaxAbsolute - grid.theta)) * 
                vFm * (V / np.abs(V).max()))
        Vsf = V + asym

        # Surface wind reduction factor:
        swrf = 0.81 * np.ones(np.shape(Vsf))
        low = np.where(np.abs(Vsf) >= 6)
        med = np.where(np.abs(Vsf) >= 19.5)
        high = np.where(np.abs(Vsf) >= 45)
        swrf[low] = 0.81 - (2.93 * (np.abs(Vsf[low]) - 6.) / 1000.)
        swrf[med] = 0.77 - (4.31 * (np.abs(Vsf[med]) - 19.5) / 1000.)
        swrf[high] = 0.66

        # Additional factor to convert to 1-munite sustained wind speed:
        Ux = swrf * Vsf * np.sin(phi) * 1.069
        Vy = swrf * Vsf * np.cos(phi) * 1.069

        return Ux, Vy


class KepertWindField(WindFieldModel):

    """
    Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the
    Tropical Cyclone Core. Part I: Linear Theory.  J. Atmos. Sci., 58,
    2469-2484

    """

    def field(self, grid, vFm, thetaFm, thetaMax=0.):
        """
        :param R: Distance from the storm centre to the grid (km).
        :type  R: :class:`numpy.ndarray`
        :param lam: Direction (geographic bearing, positive clockwise)
                    from storm centre to the grid.
        :type  lam: :class:`numpy.ndarray`
        :param float vFm: Foward speed of the storm (m/s).
        :param float thetaFm: Forward direction of the storm (geographic
                              bearing, positive clockwise, radians).
        :param float thetaMax: Bearing of the location of the maximum
                               wind speed, relative to the direction of
                               motion.

        """

        V = self.velocity(grid)
        Z = self.vorticity(grid)
        K = 50.  # Diffusivity
        Cd = 0.002  # Constant drag coefficient
        thetaFm = np.pi/2. - thetaFm
        Vm = np.abs(V).max()
        if (vFm > 0) and (Vm/vFm < 5.):
            Umod = vFm * (1. - (1. - Vm/vFm)/5.)
        else:
            Umod = vFm
        Vt = Umod * np.ones(V.shape)
        
        core = np.where(grid.R > 2. * self.rMax)
        Vt[core] = Umod * np.exp(-((grid.R[core] / (2.*self.rMax)) - 1.) ** 2. )

        al = ((2. * V / grid.R ) + self.f) / (2. * K)
        be = (self.f + Z) / (2. * K)
        gam = (V / (2. * K * grid.R))
        
        albe = np.sqrt(al / be)

        ind = np.where(np.abs(gam) > np.sqrt(al * be))
        chi = np.abs((Cd / K) * V / np.sqrt(np.sqrt(al * be)))
        eta = np.abs((Cd / K) * V / np.sqrt(np.sqrt(al * be) + np.abs(gam)))
        psi = np.abs((Cd / K) * V / np.sqrt(np.abs(np.sqrt(al * be)
                                                   - np.abs(gam))))

        i = complex(0., 1.)
        A0 = -(chi * (1 + i * (1 + chi))*V) / (2*chi**2 + 3*chi +2)
        
        # Symmetric surface wind component
        u0s = np.sign(self.f) * albe * A0.real
        v0s =                          A0.imag

        Am = -(psi * (1 + 2 * albe + (1 + i) * (1 + albe) * eta) * Umod) / \
             (albe * ((2 + 2 * i) * (1 + eta * psi) + 3 * psi + 3* i * eta))
        AmIII = -(psi * (1 + 2 * albe + (1 + i) * (1 + albe) * eta) * Umod) / \
                (albe * ((2 - 2 * i + 3 * (eta + psi) + (2 + 2 * i)*eta * psi)))
        Am[ind] = AmIII[ind]

        # First asymmetric surface component
        ums =            albe * (Am * np.exp(-i * grid.theta *
                                             np.sign(self.f))).real
        vms = np.sign(self.f) * (Am * np.exp(-i * grid.theta *
                                             np.sign(self.f))).imag

        Ap = -(eta * (1 - 2 * albe + (1 + i)*(1 - albe) * psi) * Umod) / \
             (albe * ((2 + 2*i)*(1 + eta * psi) + 3*eta + 3*i*psi))
        ApIII = -(eta * (1 - 2 * albe + (1 - i)*(1 - albe)*psi)*Umod) / \
                (albe * (2 + 2 * i + 3 * (eta + psi) + (2 -2 * i)*eta*psi))
        Ap[ind] = ApIII[ind]

        # Second asymmetric surface component
        ups =            albe * (Ap * np.exp(i * grid.theta *
                                             np.sign(self.f))).real
        vps = np.sign(self.f) * (Ap * np.exp(i * grid.theta *
                                             np.sign(self.f))).imag

        # Total surface wind in (moving coordinate system)
        us =     u0s + ups + ums
        vs = V + v0s + vps + vms

        usf = us + Umod * np.cos(grid.theta - thetaFm)
        vsf = vs - Umod * np.sin(grid.theta - thetaFm)
        phi = np.arctan2(usf, vsf)

        # Surface winds, cartesian coordinates
        Ux = np.sqrt(usf ** 2. + vsf ** 2.) * np.sin(phi - grid.theta)
        Vy = np.sqrt(usf ** 2. + vsf ** 2.) * np.cos(phi - grid.theta)

        return Ux, Vy


# Automatic discovery of models and required parameters


def allSubclasses(cls):
    """
    Recursively find all subclasses of a given class.
    """
    return cls.__subclasses__() + \
        [g for s in cls.__subclasses__() for g in allSubclasses(s)]


def profile(name):
    """
    Helper function to return the appropriate wind profile
    model given a `name`.
    """
    return PROFILES[name]


def profileParams(name):
    """
    List of additional parameters required for a wind profile model.
    """
    from inspect import getargspec
    std = getargspec(WindProfileModel.__init__)[0]
    new = getargspec(profile(name).__init__)[0]
    params = [p for p in new if p not in std]
    return params


def field(name):
    """
    Helper function to return the appropriate wind field
    model given a `name`.
    """
    return FIELDS[name]


def fieldParams(name):
    """
    List of additional parameters required for a wind field model.
    """
    from inspect import getargspec
    std = getargspec(WindFieldModel.__init__)[0]
    new = getargspec(field(name).__init__)[0]
    params = [p for p in new if p not in std]
    return params


PROFILES = dict([(k.__name__.replace('WindProfile', '').lower(), k)
                 for k in allSubclasses(vars()['WindProfileModel'])])

FIELDS = dict([(k.__name__.replace('WindField', '').lower(), k)
               for k in allSubclasses(vars()['WindFieldModel'])])
