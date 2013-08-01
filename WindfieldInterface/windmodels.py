import numpy as np
from math import exp, sqrt, log


class WindSpeedModel(object):

    def __init__(self, windProfileModel):

        # Convert from hPa to Pa if necessary

        if windProfileModel.cP < 10000:
            self.cP = metutils.convert(pCentre, 'hPa', 'Pa')

        if windProfileModel.eP < 10000:
            self.eP = metutils.convert(pEnv, 'hPa', 'Pa')

        self.dP = eP - cP

    def maximum(self):
        raise NotImplementedError


class WilloughbyWindSpeed(WindSpeedModel):

    """
    Willoughby & Rahn (2004), Parametric Representation of the
    Primary Hurricane Vortex. Part I: Observations and
    Evaluation of the Holland (1980) Model.
    Mon. Wea. Rev., 132, 3033-3048
    """

    def maximum(self):
        return 0.6252 * sqrt(self.dP)


class HollandWindSpeed(WindSpeedModel):

    """
    Holland (1980), An Analytic Model of the Wind and Pressure
    Profiles in Hurricanes. Mon. Wea. Rev, 108, 1212-1218
    Density of air is assumed to be 1.15 kg/m^3.
    beta is assumed to be 1.3. Other values can be specified.
    Gradient level wind (assumed maximum).
    """

    def __init__(self, windProfileModel):
        super(HollandWindSpeed, self).__init__(windProfileModel)
        self.beta = windProfileModel.beta
        self.rho = 1.15

    def maximum(self):
        return sqrt(self.beta * self.dP / (exp(1) * self.rho))


class AtkinsonWindSpeed(WindSpeedModel):

    """
    Atkinson and Holliday (1977), Tropical Cyclone Minimum Sea
    Level Pressure / Maximum Sustained Wind Relationship for
    the Western North Pacific. Mon. Wea. Rev., 105, 421-427
    Maximum 10m, 1-minute wind speed. Uses pEnv as 1010 hPa
    """

    def maxSpeed(self):
        return 3.04 * pow(1010.0 - metutils.convert(self.cP, 'Pa', 'hPa'), 0.644)


class WindProfileModel(object):

    def __init__(self, lat, lon, eP, cP, rMax):
        self.rho = 1.15  # density of air
        self.lat = lat
        self.lon = lon
        self.eP = eP
        self.cP = cP
        self.rMax = rMax
        self.f = metutils.coriolis(lat)
        self.dP = eP - cP

    def atRadiuses(self, R):
        raise NotImplementedError


class JelesnianskiWindProfile(WindProfileModel):

    """
    Jelesnianski model of the wind profile
    """

    def __init__(self, lat, lon, eP, cP, rMax, windSpeed=WilloughbyWindSpeed):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax)
        self.vMax = windSpeed(self).maximum()

    def atRadiuses(self, R):
        V = 2 * self.vMax * self.rMax * R / (self.rMax ** 2 + R ** 2)
        V = np.sign(f) * V
        return V


class HollandWindProfile(WindProfileModel):

    """
    Holland profile. For r < rMax, we reset the wind field to a
    cubic profile to avoid the barotropic instability mentioned in
    Kepert & Wang (2001).
    """

    def __init__(self, lat, lon, eP, cP, rMax, beta, windSpeed=HollandWindSpeed):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax)
        self.vMax = windSpeed(self).maximum()
        self.beta = beta

    def atRadiuses(self, R):
        dP = self.dP

        # Calculate first and second derivatives at R = Rmax

        d2Vm = derivatives.holland(self.f, self.rMax, self.beta, dP, self.rho)
        aa = (d2Vm / 2. - (-self.vMax / self.rMax) / self.rMax) / self.rMax
        bb = (d2Vm - 6 * aa * self.rMax) / 2.
        cc = -3 * aa * self.rMax ** 2 - 2 * bb * self.rMax

        # The profile as originally stated

        delta = (self.rMax / R) ** self.beta
        edelta = np.exp(-delta)

        V = np.sign(self.f) * np.sqrt((dP * self.beta / self.rho) *
                                      delta * edelta + (R * self.f / 2.) ** 2)
                            - R * np.abs(self.f) / 2.

        # Replace all values within rMax of the storm centre with the
        # cubic profile to eliminate barotropic instability:

        icore = np.where(R <= self.rMax)
        V[icore] = np.sign(self.f) * R[icore] \
                                   * (R[icore] * (R[icore] * aa + bb) + cc)

        return V


class WilloughbyWindProfile(HollandWindProfile):

    """
    The Willoughby & Rahn (2004) relation, which makes beta a
    function of Vmax, rMax and latitude. We use Willoughby & Rahn's
    (2004) relation for Vmax *only*.
    This determines the beta parameter then calls Holland (which
    means the profile is cubic within Rmax) to calculate the wind
    profile.  The beta term calculation is based on Atlantic and
    Eastern Pacific cyclone data, not Australian data.
    """

    def __init__(self, lat, lon, eP, cP, rMax, windSpeed=WilloughbyWindSpeed):
        #NOTE: beta depends on vMax so we set a dummy beta and redefine later
        HollandProfileModel.__init__(self, lat, lon, eP, cP, rMax, 1.0, windSpeed)
        self.beta = 1.0036 + 0.0173 * self.vMax - 0.313 * log(rMax) + 0.0087 * abs(lat)
 

class RankineWindProfile(WilloughbyWindProfile):

    """
    Rankine vortex profile. Vmax determined by dp using, by default,
    the Willoughby & Rahn method.
    """

    def __init__(self, lat, lon, eP, cP, rMax, windSpeed=WilloughbyWindSpeed):
        WilloughbyWindProfile.__init__(self, lat, lon, eP, cP, rMax, windSpeed)
        self.alpha = 0.5

    def atRadiuses(self, R):
        # An assumption about the shape of the profile outside Rmax.
        # The literature indicates 0.4 < alpha < 0.6 (e.g. see Holland,
        # 1980)
        V = self.vMax * (self.rMax / R) ** self.alpha
        icore = np.where(R <= self.rMax)
        V[icore] = self.vMax * (R[icore] / self.rMax)
        V = np.sign(self.f) * V
        return V


class SchloemerWindProfile(HollandWindProfile):

    """
    Schloemer's (1954) is the same as the Holland relation with
    beta = 1
    """

    def __init__(self, lat, lon, eP, cP, rMax):
        HollandWindProfile.__init__(self, lat, lon, eP, cP, rMax, 1.0)


class DoubleHollandWindProfile(WindProfileModel):

    """
    McConochie et al's double Holland vortex model (based on Cardone
    et al, 1994).  This application is the Coral Sea adaptation of
    the double vortex model (it can also be used for concentric
    eye-wall configurations).  The tunable parameters in this
    relation are 'dp1', 'dp2', 'b1', 'b2' and 'rMax2'
    """

    def __init__(self, lat, lon, eP, cP, rMax, beta1, beta2, dp1, dp2, rMax2):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax)
        self.beta1 = beta1
        self.beta2 = beta2
        self.dp1 = dp1
        self.dp2 = dp2
        self.rMax2 = rMax2

    def atRadiuses(self, R):
        # Scale dp2 if dP is less than 800 Pa:
        if self.dP < 1500.:
            dp2 = (self.dP / 1500.) * (800. + (self.dP - 800.) / 2000.)
        else:
            dp2 = 800. + (self.dP - 800.) / 2000.
        dp1 = self.dP - dp2
        if self.beta1 is None:
            self.beta1 = 7.3 - self.cP / 16000.
        if self.beta2 is None:
            self.beta2 = 7.2 - self.cP / 16000.

        # The two gradient wind components:
        mu = (rMax / R) ** self.beta1
        nu = (rMax2 / R) ** self.beta2
        emu = np.exp(-mu)
        enu = np.exp(-nu)

        gradientV1 = (self.beta1 * dp1 / self.rho) * mu * emu
        gradientV2 = (self.beta2 * dp2 / self.rho) * nu * enu

        V = np.sign(self.f) * np.sqrt(gradientV1 + gradientV2 + (R * self.f / 2.) ** 2) \
            - R * np.abs(self.f) / 2.

        vMax = np.abs(V).max()
        # aa, bb, cc = doubleHollandCoefficient(vMax, rMax1, rMax2, dp1, dp2,
        # beta1, beta2, f, rho)

        # delta = (rMax/R)**self.beta1
        # gamma = (rMax2/R)**self.beta2

        # Calculate first and second derivatives at R = Rmax:
        d2Vm = derivatives.doubleHolland(self.f, rMax, rMax2,
                                         self.beta1, self.beta2, dp1, dp2,
                                         self.rho)
        aa = (d2Vm / 2. - (-vMax / rMax) / rMax) / rMax
        bb = (d2Vm - 6 * aa * rMax) / 2.
        cc = -3 * aa * rMax ** 2 - 2 * bb * rMax

        # Replace all values within rMax of the storm centre with the cubic
        # profile to eliminate barotropic instability:
        if self.dP >= 1500.:
            icore = np.where(R <= rMax)
            V[icore] = np.sign(self.f) * R[icore] \
                * (R[icore] * (R[icore] * aa + bb) + cc)

        return V


class PowellWindProfile(HollandWindProfile):

    """
    Powell et al, 2005
    Another definition of the B parameter inserted into the Holland
    model.  Unlike Willoughby and Rahn's model, there is no reliance
    on vMax.  Powell et al. also included a small random term, but
    since the beta value is also used in the vorticity calculation,
    we need to ensure the values used in this function and the
    corresponding vorticity function match.
    """

    def __init__(self, lat, lon, eP, cP, rMax):
        beta = 1.881093 - 0.010917 * np.abs(lat) - 0.005567 * rMax
        if beta < 0.8:
            beta = 0.8
        if beta > 2.2:
            beta = 2.2
        HollandProfileModel.__init__(self, lat, lon, eP, cP, rMax, beta)


class NewHollandWindModel(WindProfileModel):

    """
    Holland et al. 2010
    In this version, the exponent is allowed to vary linearly
    outside the radius of maximum wind. i.e. rather than take the
    sqare root, the exponent varies around 0.5.
    Currently this version does not have a corresponding vorticity
    profile set up in windVorticity, so it cannot be applied in
    wind field modelling.
    """

    def __init__(self, lat, lon, eP, cP, rMax, rGale=150.):
        WindProfileModel.__init__(self, lat, lon, eP, cP, rMax)
        self.rGale = rGale

     def atRadiuses(self, R):
        # In this incarnation, we are assuming the pressure rate of change and
        # forward velocity is zero, so there is no requirement for a first-pass
        # guess at x
        Bs = -0.000044 * (self.dP / 100.) ** 2. + 0.01 * ( self.dP / 100.) -
        0.014 * np.abs(lat) + 1.0 deltag = np.power(rMax / rGale, Bs) edeltag =
        np.exp(-1. * deltag) rgterm = Bs * self.dP * deltag * edeltag /
        self.rho xn = np.log(17.) / np.log(rgterm) xx = 0.5 * np.ones(R.shape)
        i = np.where(R > rMax) xx[i] = 0.5 + (R[i] - rMax) * (xn - 0.5) /
        (rGale - rMax)

        delta = (rMax / R) ** Bs
        edelta = np.exp(-delta)

        V = np.sign(self.f) * np.power(
            (self.dP * Bs / self.rho) * delta * edelta, xx)

        return V
