"""
Wind Models

"""

import numpy as np
from math import exp, sqrt
import Utilities.metutils as metutils


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
    Holland (1980), An Analytic Model of the Wind and Pressure Profiles
    in Hurricanes. Mon. Wea. Rev, 108, 1212-1218 Density of air is
    assumed to be 1.15 kg/m^3.  beta is assumed to be 1.3. Other values
    can be specified.  Gradient level wind (assumed maximum).
    """

    def maximum(self):
        beta = self.profile.beta
        rho = 1.15
        return sqrt(beta * self.dP / (exp(1) * rho))


class AtkinsonWindSpeed(WindSpeedModel):

    """
    Atkinson and Holliday (1977), Tropical Cyclone Minimum Sea
    Level Pressure / Maximum Sustained Wind Relationship for
    the Western North Pacific. Mon. Wea. Rev., 105, 421-427
    Maximum 10m, 1-minute wind speed. Uses pEnv as 1010 hPa
    """

    def maximum(self):
        cP = metutils.convert(self.cP, 'Pa', 'hPa')
        return 3.04 * pow(1010.0 - cP, 0.644)


class WindProfileModel(object):

    """
    Wind profile model.
    """

    def __init__(self, lat, lon, eP, cP, rMax, windSpeedModel):
        self.rho = 1.15  # density of air
        self.lat = lat
        self.lon = lon
        self.eP = eP
        self.cP = cP
        self.rMax = rMax
        self.speed = windSpeedModel(self)
        self.f = metutils.coriolis(lat)
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
        """
        self.vMax_ = value

    def velocity(self, R):
        """
        Wind velocity at radiuses `R`.
        """
        raise NotImplementedError

    def vorticity(self, R):
        """
        Wind vorticity at radiuses `R`.
        """
        raise NotImplementedError


class JelesnianskiWindProfile(WindProfileModel):

    """
    Jelesnianski model of the wind profile
    """

    def __init__(self, lat, lon, eP, cP, rMax,
                 windSpeedModel=WilloughbyWindSpeed):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax,
                                  windSpeedModel)

    def velocity(self, R):
        V = 2 * self.vMax * self.rMax * R / (self.rMax ** 2 + R ** 2)
        V = np.sign(self.f) * V
        return V

    def vorticity(self, R):
        Z = (np.sign(self.f) * 2 * self.vMax * self.rMax / (self.rMax
             ** 2 + R ** 2) + np.sign(self.f) * 2 * self.vMax *
             self.rMax * (self.rMax ** 2 - R ** 2) /
             (self.rMax ** 2 + R ** 2) ** 2)
        return Z


class HollandWindProfile(WindProfileModel):

    """
    Holland profile. For r < rMax, we reset the wind field to a
    cubic profile to avoid the barotropic instability mentioned in
    Kepert & Wang (2001).
    """

    def __init__(self, lat, lon, eP, cP, rMax, beta,
                 windSpeedModel=HollandWindSpeed):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax,
                                  windSpeedModel)
        self.beta = beta

    def secondDerivative(self):
        """
        Second derivative of profile at rMax
        """

        beta = self.beta
        dP = self.dP
        rho = self.rho
        f = self.f
        rMax = self.rMax

        E = exp(1)
        d2Vm = ((beta * dP * (-4 * beta ** 3 * dP / rho -
                (-2 + beta ** 2) * E * (f * rMax) ** 2)) /
                (E * rho * sqrt((4 * beta * dP) / (E * rho)
                 + (f * rMax) ** 2) * (4 * beta * dP * rMax ** 2 / rho
                 + E * (f * rMax ** 2) ** 2)))

        assert d2Vm < 0.0

        return d2Vm

    def velocity(self, R):
        d2Vm = self.secondDerivative()
        aa = ((d2Vm / 2. - (-self.vMax / self.rMax) / self.rMax) /
              self.rMax)
        bb = (d2Vm - 6 * aa * self.rMax) / 2.
        cc = -3 * aa * self.rMax ** 2 - 2 * bb * self.rMax
        delta = (self.rMax / R) ** self.beta
        edelta = np.exp(-delta)

        V = (np.sign(self.f) * np.sqrt((self.dP * self.beta / self.rho)
             * delta * edelta + (R * self.f / 2.) ** 2) - R *
             np.abs(self.f) / 2.)

        icore = np.where(R <= self.rMax)
        V[icore] = (np.sign(self.f) * R[icore] * (R[icore] * (R[icore]
                    * aa + bb) + cc))

        return V

    def vorticity(self, R):
        beta = self.beta
        delta = (self.rMax / R) ** beta
        edelta = np.exp(-delta)

        Z = (np.sign(self.f) * (np.sqrt((self.dP * beta / self.rho) *
             delta * edelta + (R * self.f / 2.) ** 2)) / R -
             np.abs(self.f) + edelta *
             (2 * (beta ** 2) * self.dP * (delta - 1) * delta +
              self.rho * edelta * (self.f * R) ** 2) /
             (2 * self.rho * R *
              np.sqrt(4 * (beta * self.dP / self.rho) * delta * edelta
                      + (self.f * R) ** 2)))

        # Calculate first and second derivatives at R = Rmax:
        d2Vm = self.secondDerivative()
        aa = ((d2Vm / 2 - (-1.0 * np.sign(self.f) * self.vMax /
              self.rMax) / self.rMax) / self.rMax)
        bb = (d2Vm - 6 * aa * self.rMax) / 2
        cc = -3 * aa * self.rMax ** 2 - 2 * bb * self.rMax

        icore = np.where(R <= self.rMax)
        Z[icore] = R[icore] * (R[icore] * 4 * aa + 3 * bb) + 2 * cc

        return Z


class WilloughbyWindProfile(HollandWindProfile):

    """
    The Willoughby & Rahn (2004) relation, which makes beta a function
    of Vmax, rMax and latitude. We use Willoughby & Rahn's (2004)
    relation for Vmax *only*.  This determines the beta parameter then
    calls Holland (which means the profile is cubic within Rmax) to
    calculate the wind profile.  The beta term calculation is based on
    Atlantic and Eastern Pacific cyclone data, not Australian data.
    """

    def __init__(self, lat, lon, eP, cP, rMax,
                 windSpeedModel=WilloughbyWindSpeed):
        HollandWindProfile.__init__(self, lat, lon, eP, cP, rMax, 1.0,
                                    windSpeedModel)
        self.beta = (1.0036 + 0.0173 * self.vMax - 0.313 * np.log(rMax)
                     + 0.0087 * np.abs(lat))
        self.speed = HollandWindSpeed(self)


class RankineWindProfile(WindProfileModel):

    """
    Rankine vortex profile. Vmax determined by dp using, by default,
    the Willoughby & Rahn method.
    """

    def __init__(self, lat, lon, eP, cP, rMax,
                 windSpeedModel=WilloughbyWindSpeed):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax,
                                  windSpeedModel)
        self.alpha = 0.5

    def velocity(self, R):
        # An assumption about the shape of the profile outside Rmax.
        # The literature indicates 0.4 < alpha < 0.6 (e.g. see Holland,
        # 1980)
        V = self.vMax * (self.rMax / R) ** self.alpha
        icore = np.where(R <= self.rMax)
        V[icore] = self.vMax * (R[icore] / self.rMax)
        V = np.sign(self.f) * V
        return V

    def vorticity(self, R):
        Z = (np.sign(self.f) * self.vMax * ((self.rMax / R) **
             self.alpha) / R - self.alpha * self.vMax * (self.rMax **
             self.alpha) / (R ** self.alpha))
        icore = np.where(R <= self.rMax)
        Z[icore] = (np.sign(self.f) * (self.vMax * (R[icore] /
                    self.rMax) + self.vMax / self.rMax))
        return Z


class SchloemerWindProfile(HollandWindProfile):

    """
    Schloemer's (1954) is the same as the Holland relation with
    beta = 1
    """

    def __init__(self, lat, lon, eP, cP, rMax):
        HollandWindProfile.__init__(self, lat, lon, eP, cP, rMax, 1.0)


class DoubleHollandWindProfile(WindProfileModel):

    """
    McConochie et al's double Holland vortex model (based on Cardone et
    al, 1994).  This application is the Coral Sea adaptation of the
    double vortex model (it can also be used for concentric eye-wall
    configurations).
    """

    def __init__(self, lat, lon, eP, cP, rMax, beta1, beta2, rMax2):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax,
                                  WindSpeedModel)

        # Scale dp2 if dP is less than 800 Pa

        if self.dP < 1500.:
            self.dp2 = ((self.dP / 1500.) * (800. + (self.dP - 800.) /
                        2000.))
        else:
            self.dp2 = 800. + (self.dP - 800.) / 2000.

        self.dp1 = self.dP - self.dp2

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

    def velocity(self, R):
        rMax = self.rMax
        rMax2 = self.rMax2

        # Scale dp2 if dP is less than 800 Pa

        if self.dP < 1500.:
            dp2 = (self.dP / 1500.) * (800. + (self.dP - 800.) / 2000.)
        else:
            dp2 = 800. + (self.dP - 800.) / 2000.

        dp1 = self.dP - dp2

        # The two gradient wind components

        mu = (rMax / R) ** self.beta1
        nu = (rMax2 / R) ** self.beta2
        emu = np.exp(-mu)
        enu = np.exp(-nu)

        gradientV1 = (self.beta1 * dp1 / self.rho) * mu * emu
        gradientV2 = (self.beta2 * dp2 / self.rho) * nu * enu

        V = (np.sign(self.f) * np.sqrt(gradientV1 + gradientV2 + (R *
             self.f / 2.) ** 2) - R * np.abs(self.f) / 2.)

        vMax = np.abs(V).max()

        d2Vm = self.secondDerivative()
        aa = (d2Vm / 2. - (-vMax / rMax) / rMax) / rMax
        bb = (d2Vm - 6 * aa * rMax) / 2.
        cc = -3 * aa * rMax ** 2 - 2 * bb * rMax

        # Replace all values within rMax of the storm centre with the
        # cubic profile to eliminate barotropic instability

        if self.dP >= 1500.:
            icore = np.where(R <= rMax)
            V[icore] = (np.sign(self.f) * R[icore] * (R[icore] *
                        (R[icore] * aa + bb) + cc))

        return V

    def vorticity(self, R):

        # Scale dp2 if dP is less than 1500 Pa:
        if self.dP < 1500.:
            dp2 = (self.dP / 1500.) * (800. + (self.dP - 800.) / 2000.)
        else:
            dp2 = 800. + (self.dP - 800.) / 2000.

        dp1 = self.dP - dp2

        chi = self.beta1 * dp1 / self.rho
        psi = self.beta2 * dp2 / self.rho

        delta = (self.rMax / R) ** self.beta1
        gamma = (self.rMax2 / R) ** self.beta2
        edelta = np.exp(-delta)
        egamma = np.exp(-gamma)

        # Derivatives:

        ddelta = (-self.beta1 * (self.rMax ** self.beta1) / (R **
                  (self.beta1 + 1)))
        dgamma = (-self.beta2 * (self.rMax2 ** self.beta2) / (R **
                  (self.beta2 + 1)))

        Z = (np.sign(self.f) * np.sqrt(chi * delta * edelta + psi *
             gamma * egamma + (self.f * R / 2) ** 2) / R -
             np.abs(self.f) + (1 / 2) *
             (chi * ddelta * edelta * (1 - delta) +
              psi * dgamma * egamma * (1 - gamma) +
              R * self.f ** 2) /
             np.sqrt(chi * delta * edelta + psi * gamma *
                     egamma + (self.f * R / 2) ** 2))

        d2Vm = self.secondDerivative()
        aa = ((d2Vm / 2.0 - (-1.0 * np.sign(self.f) * self.vMax /
              self.rMax) / self.rMax) / self.rMax)
        bb = (d2Vm - 6.0 * aa * self.rMax) / 2.0
        cc = -3.0 * aa * self.rMax ** 2.0 - 2.0 * bb * self.rMax

        if self.dP >= 1500.:
            icore = np.where(R <= self.rMax)
            Z[icore] = (R[icore] * (R[icore] * 4.0 * aa + 3.0 * bb) +
                        2.0 * cc)

        return Z


class PowellWindProfile(HollandWindProfile):

    """
    Powell et al, 2005. Another definition of the beta parameter
    inserted into the Holland model.  Unlike Willoughby and Rahn's
    model, there is no reliance on vMax.  Powell et al. also included a
    small random term, but since the beta value is also used in the
    vorticity calculation, we need to ensure the values used in this
    function and the corresponding vorticity function match.
    """

    def __init__(self, lat, lon, eP, cP, rMax):
        beta = 1.881093 - 0.010917 * np.abs(lat) - 0.005567 * rMax
        if beta < 0.8:
            beta = 0.8
        if beta > 2.2:
            beta = 2.2

        HollandWindProfile.__init__(self, lat, lon, eP, cP, rMax, beta)


class NewHollandWindProfile(WindProfileModel):

    """
    Holland et al. 2010.  In this version, the exponent is allowed to
    vary linearly outside the radius of maximum wind. i.e. rather than
    take the sqare root, the exponent varies around 0.5.  Currently
    this version does not have a corresponding vorticity profile set up
    in windVorticity, so it cannot be applied in wind field modelling.
    """

    def __init__(self, lat, lon, eP, cP, rMax, rGale=150.):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax)
        self.rGale = rGale

    def velocity(self, R):
        # In this incarnation, we are assuming the pressure rate of
        # change and forward velocity is zero, so there is no
        # requirement for a first-pass guess at x

        Bs = (-0.000044 * (self.dP / 100.) ** 2. + 0.01 * (self.dP /
              100.) - 0.014 * np.abs(self.lat) + 1.0)

        deltag = np.power(self.rMax / self.rGale, Bs)
        edeltag = np.exp(-1. * deltag)
        rgterm = Bs * self.dP * deltag * edeltag / self.rho
        xn = np.log(17.) / np.log(rgterm)
        xx = 0.5 * np.ones(R.shape)

        i = np.where(R > self.rMax)
        xx[i] = (0.5 + (R[i] - self.rMax) * (xn - 0.5) / (self.rGale -
                 self.rMax))

        delta = (self.rMax / R) ** Bs
        edelta = np.exp(-delta)

        V = (np.sign(self.f) * np.power((self.dP * Bs / self.rho) *
             delta * edelta, xx))

        return V

    def vorticity(self, R):
        raise Exception


class WindFieldModel(object):

    """
    Wind field model.
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

    def velocity(self, R):
        """
        Helper property to return the wind velocity at radiuses `R`
        from the wind profile or the precalculated attribute.
        """
        if self.V is None:
            return self.profile.velocity(R)
        else:
            return self.V

    def vorticity(self, R):
        """
        Helper property to return the wind vorticity at radiuses `R`
        from the wind profile or the precalculated attribute.
        """
        if self.Z is None:
            return self.profile.vorticity(R)
        else:
            return self.Z

    def field(self, R, lam, vFm, thetaFm, thetaMax=0.):
        """
        The wind field.
        """
        raise NotImplementedError


class HubbertWindField(WindFieldModel):

    """
    Hubbert, G.D., G.J. Holland, L.M. Leslie and M.J. Manton, 1991:
    A Real-Time System for Forecasting Tropical Cyclone Storm Surges.
    Weather and Forecasting, 6, 86-97
    """

    def field(self, R, lam, vFm, thetaFm, thetaMax=0.):
        V = self.velocity(R)

        Km = .70
        inflow = 25. * np.ones(np.shape(R))
        core = np.where(R < self.rMax)
        inflow[core] = 0
        inflow = inflow * np.pi / 180

        thetaMaxAbsolute = thetaFm + thetaMax
        asym = vFm * np.cos(thetaMaxAbsolute - lam + np.pi)
        Vsf = Km * V + asym
        phi = inflow - lam

        Ux = Vsf * np.sin(phi)
        Vy = Vsf * np.cos(phi)

        return Ux, Vy


class McConochieWindField(WindFieldModel):

    """
    McConochie, J.D., T.A. Hardy and L.B. Mason, 2004:
    Modelling tropical cyclone over-water wind and pressure fields.
    Ocean Engineering, 31, 1757-1782
    """

    def field(self, R, lam, vFm, thetaFm, thetaMax=0.):
        V = self.velocity(R)

        inflow = 25. * np.ones(np.shape(R))
        mid = np.where(R < 1.2 * self.rMax)
        inflow[mid] = 10. + 75. * (R[mid] / self.rMax - 1.)
        inner = np.where(R < self.rMax)
        inflow[inner] = 10. * R[inner] / self.rMax
        inflow = inflow * np.pi / 180.

        thetaMaxAbsolute = thetaFm + thetaMax
        phi = inflow - lam

        asym = (0.5 * (1. + np.cos(thetaMaxAbsolute - lam)) * vFm * (V
                / np.abs(V).max()))
        Vsf = V + asym

        # Surface wind reduction factor:
        swrf = 0.81 * np.ones(np.shape(Vsf))
        low = np.where(Vsf >= 6)
        med = np.where(Vsf >= 19.5)
        high = np.where(Vsf >= 45)
        swrf[low] = 0.81 - (2.93 * (Vsf[low] - 6.) / 1000.)
        swrf[med] = 0.77 - (4.31 * (Vsf[med] - 19.5) / 1000.)
        swrf[high] = 0.66

        Ux = swrf * Vsf * np.sin(phi)
        Vy = swrf * Vsf * np.cos(phi)

        return Ux, Vy


class KepertWindField(WindFieldModel):

    """
    Kepert, J., 2001: The Dynamics of Boundary Layer Jets within the
    Tropical Cyclone Core. Part I: Linear Theory.  J. Atmos. Sci., 58,
    2469-2484
    """

    def field(self, R, lam, vFm, thetaFm, thetaMax=0.):
        V = self.velocity(R)
        Z = self.vorticity(R)
        K = 50.  # Diffusivity
        Cd = 0.002  # Constant drag coefficient

        al = (2. * V / R + self.f) / (2. * K)
        be = (self.f + Z) / (2. * K)
        gam = V / (2. * K * R)
        albe = np.sqrt(al / be)

        ind = np.where(np.abs(gam) > np.sqrt(al * be))
        chi = (Cd / K) * V / np.sqrt(np.sqrt(al * be))
        eta = (Cd / K) * V / np.sqrt(np.sqrt(al * be) + np.abs(gam))
        psi = (Cd / K) * V / np.sqrt(np.abs(np.sqrt(al * be) - gam))

        i = complex(0., 1.)
        A0 = (-chi * V * (1. + i * (1. + chi)) / (2. * chi ** 2. + 3.
              * chi + 2.))

        # Symmetric surface wind component

        u0s = albe * A0.real
        v0s = A0.imag

        Am = (-((1. + (1. + i) * eta) / albe + (2. + (1. + i) * eta))
              * psi * vFm /
              ((2. + 2. * i) * (1 + eta * psi) + 3. * psi + 3. * i * eta))

        Am[ind] = (-((1. + (1. + i) * eta[ind]) / albe[ind] +
                   (2. + (1. + i) * eta[ind])) * psi[ind] * vFm /
                   (2. - 2. * i + 3. * (psi[ind] + eta[ind])
                    + (2. + 2. * i) * eta[ind] * psi[ind]))

        # First asymmetric surface component

        ums = albe * (Am * np.exp(-i * lam)).real
        vms = (Am * np.exp(-i * lam)).imag

        Ap = (-((1. + (1. + i) * psi) / albe - (2. + (1. + i) * psi)) *
              eta * vFm /
              ((2. + 2. * i) * (1. + eta * psi) + 3. * eta +
               3. * i * psi))

        Ap[ind] = (-((1. + (1. - i) * psi[ind]) / albe[ind] -
                     (2. + (1. - i) * psi[ind])) *
                   eta[ind] * vFm /
                   (2. + 2. * i + 3. * (eta[ind] + psi[ind]) +
                    (2. - 2. * i) * eta[ind] * psi[ind]))

        # Second asymmetric surface component

        ups = albe * (Ap * np.exp(i * lam)).real
        vps = (Ap * np.exp(i * lam)).imag

        # Total surface wind in (moving coordinate system)

        us = u0s + ups + ums
        vs = V + v0s + vps + vms

        usf = us + vFm * np.cos(lam - thetaFm)
        vsf = vs - vFm * np.sin(lam - thetaFm)
        phi = np.arctan2(usf, vsf)

        # Surface winds, cartesian coordinates

        Ux = (np.sqrt(usf ** 2. + vsf ** 2.) * np.sin(phi - lam))
        Vy = (np.sqrt(usf ** 2. + vsf ** 2.) * np.cos(phi - lam))

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
    return WIND_PROFILES[name]


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
    return WIND_FIELDS[name]


def fieldParams(name):
    """
    List of additional parameters required for a wind field model.
    """
    from inspect import getargspec
    std = getargspec(WindFieldModel.__init__)[0]
    new = getargspec(field(name).__init__)[0]
    params = [p for p in new if p not in std]
    return params


WIND_PROFILES = {k.__name__.replace('WindProfile', '').lower(): k
                 for k in allSubclasses(vars()['WindProfileModel'])}

WIND_FIELDS = {k.__name__.replace('WindField', '').lower(): k
               for k in allSubclasses(vars()['WindFieldModel'])}
