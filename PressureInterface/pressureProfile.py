"""
:mod:`pressureProfile` -- radial pressure profile for a TC
==========================================================

.. module:: pressureProfile
    :synopsis: Returns the radial pressure field for a range
               of parametric radial wind models.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

The available wind profiles are:
Rankine vortex - undefined!
Jelesnianski - undefined!
Holland (with cubic core)
Schloemer (just the Holland with beta = 1)
Willoughby and Rahn (Holland with beta = beta(vmax, rmax, latitude))
McConochie et al. (double exponential profile)
Powell et al. (Holland with beta = beta(rmax, latitude))

SeeAlso:
Constraints:
Version: $Rev: 810 $

References:
Holland, G.J., 1980:
An Analytic model of the Wind and Pressure Profiles in Hurricanes.
Mon. Wea. Rev., 108, 1212-1218.
Jelesnianski, C.P., 1966:
Numerical Computations of Storm Surges without Bottom Stress.
Mon. Wea. Rev., 94(6), 379-394
McConochie, J.D., T.A. Hardy and L.B. Mason, 2004:
Modelling tropical cyclone over-water wind and pressure fields.
Ocean Engineering, 31, 1757-1782
Powell, M., G. Soukup, S. Cocke, S. Gulati, N. Morrisuea-Leroy,
S. Hamid, N. Dorst and L. Axe, 2005:
State of Florida hurricane loss projection model: Atmospheric science component.
Journal of Wind Engineering and Industrial Aerodynamics, 93 (8), 651-674
Schloemer, R.W., 1954:
Analysis and synthesis of hurricane wind patterns over Lake Okeechobee.
NOAA Hydromet. Rep. 31, 49 pp.
Willoughby, H.E. and M.E. Rahn, 2004:
Parametric Representation of the Primary Hurricane Vortex. Part I:
Observations and Evaluation of the Holland (1980) Model.
Mon. Wea. Rev., 132, 3033-3048

$Id: pressureProfile.py 810 2012-02-21 07:52:50Z nsummons $
"""

import logging

import Utilities.metutils as metutils

import numpy
import wind.vmax as vmax
import time


class PrsProfile:
    """
    Description: Define the radial wind profiles used in tropical
    cyclone modelling. These are radial profiles only and do not include
    asymmetries that arise due to the forward motion of the storm.

    Parameters:
        R: grid of distances from the storm centre (distances in km)
        pEnv: Environmental pressure (Pa)
        pCentre: Central pressure of storm (Pa)
        rMax: Radius of maximum winds (km)
        cLat: Latitude of storm centre
        cLon: Longitude of storm centre
        beta: Holland beta parameter
    Members:
        R: grid of distances from the storm centre (distances in km)
        pEnv: Environmental pressure (Pa)
        pCentre: Central pressure of storm (Pa)
        rMax: Radius of maximum winds (km)
        cLat: Latitude of storm centre
        cLon: Longitude of storm centre
        beta: Holland beta parameter

    Methods:
        (rankine: Rankine vortex)
        (jelesnianski: Jelesnianski's storm surge model wind field)
        holland: Holland's radial wind field
        willoughby: Holland profile with beta a function of vMax, rMax
                    and cLat
        schloemer: Holland profile with beta==1
        doubleHolland: McConochie's double vortex model

    Internal Methods:
        None
    """

    def __init__(self, R, pEnv, pCentre, rMax, cLat, cLon, beta=1.3,
                 rMax2=250., beta1=None, beta2=None ):
        """
        Initialise required fields
        """
        self.R = R
        self.cLon = cLon
        self.cLat = cLat
        self.rMax = rMax
        self.dP = pEnv-pCentre
        self.pCentre = pCentre
        self.pEnv = pEnv
        # Density of air:
        self.rho = 1.15
        self.f = metutils.coriolis(cLat)
        self.beta = beta
        self.rMax2 = rMax2
        self.rMax2 = rMax2
        self.beta1 = beta1
        self.beta2 = beta2
        self.logger = logging.getLogger()
        self.logger.debug("Storm centre: %3f %3f" %(self.cLon, self.cLat))
        self.logger.debug("Coriolis parameter: %3f" % self.f)

#    def rankine(self, vMaxType="willoughby"):
#        """
#        Rankine vortex profile. Vmax determined by dp using, by default, the
#        Willoughby & Rahn method.
#        """
#        t0=time.time()
#        vMax=vmax.vmax(self.pCentre, self.pEnv, vMaxType)
        # An assumption about the shape of the profile outside Rmax.
        # The literature indicates 0.4 < alpha < 0.6 (e.g. see Holland, 1980)
#        alpha=0.5
#        V = vMax*(self.rMax/self.R)**alpha
#        icore = where(self.R <= self.rMax)
#        V[icore] = vMax*(self.R[icore]/self.rMax)
#        V = sign(self.f)*V
#        self.logger.debug( "Timing for rankine wind profile calculation: %.3f" %(time.time()-t0) )
#        return V

#    def jelesnianski(self, vMaxType="willoughby"):
#        """
#        Jelesnianski model of the wind profile
#        """
#        t0=time.time()
#        vMax=vmax.vmax(self.pCentre, self.pEnv, vMaxType)
#        V = 2*vMax*self.rMax*self.R/(self.rMax**2 + self.R**2)
#        V = sign(self.f)*V
#        self.logger.debug( "Timing for jelesnianski wind profile calculation: %.3f" %(time.time()-t0) )
#        return V

    def holland(self, beta=None):
        """
        Holland profile.
        """
        if beta == None:
            beta = self.beta
        t0 = time.time()
        P = numpy.zeros(self.R.shape)
        P = self.pCentre + self.dP*numpy.exp(-(self.rMax/self.R)**beta)
        self.logger.debug("Timing for holland wind profile calculation: %.3f"
                           % (time.time()-t0))
        return P

    def willoughby(self):
        """
        The Willoughby & Rahn (2004) relation, which makes beta a function of
        Vmax, rMax and latitude. We use Willoughby & Rahn's (2004) relation
        for Vmax *only*.
        This determines the beta parameter then calls Holland (which means the
        profile is cubic within Rmax) to calculate the wind profile.
        The beta term calculation is based on Atlantic and Eastern Pacific cyclone
        data, not Australian data.
        """
        vMax = vmax.vmax(self.pCentre, self.pEnv, type="willoughby")
        beta = 1.0036 + 0.0173*vMax  - 0.313*numpy.log(self.rMax) \
               + 0.0087*numpy.abs(self.cLat)
        P = self.holland(beta)
        return P

    def schloemer(self):
        """
        Schloemer's (1954) is the same as the Holland relation with
        beta = 1
        """
        beta = 1.
        P = self.holland(beta)
        return P

    def doubleHolland(self, rMax2=250.):
        """
        McConochie et al's double Holland vortex model (based on Cardone
        et al, 1994).  This application is the Coral Sea adaptation of
        the double vortex model (it can also be used for concentric
        eye-wall configurations).
        The tunable parameters in this relation are 'dp1', 'dp2', 'b1',
        'b2' and 'rMax2'
        """
        t0 = time.time()
        # Scale dp2 if dP is less than 800 Pa:
        if self.dP < 1500.:
            dp2 = (self.dP/1500.)*(800. + (self.dP - 800.)/2000.)
        else:
            dp2 = 800. + (self.dP - 800.)/2000.
        dp1 = self.dP - dp2
        if self.beta1 is None:
            self.beta1 = 7.3 - self.pCentre/16000.
        if self.beta2 is None:
            self.beta2 = 7.2 - self.pCentre/16000.

        # The two gradient wind components:
        mu = (self.rMax/self.R)**self.beta1
        nu = (self.rMax2/self.R)**self.beta2
        emu = numpy.exp(-mu)
        enu = numpy.exp(-nu)
        P = self.pCentre + dp1*emu +dp2*enu
#        gradientV1 = (self.beta1*dp1/self.rho)*mu*emu
#        gradientV2 = (self.beta2*dp2/self.rho)*nu*enu
#
#        P = sign(self.f)*sqrt(gradientV1+gradientV2+(self.R*self.f/2)**2)-self.R*abs(self.f)/2

#        vMax=abs(P).max()
        #aa, bb, cc = doubleHollandCoefficient(vMax, rMax1, rMax2, dp1, dp2, beta1, beta2, f, rho)

        #delta = (self.rMax/self.R)**self.beta1
        #gamma = (self.rMax2/self.R)**self.beta2

        # Calculate first and second derivatives at R = Rmax:
#        d2Vm = derivatives.doubleHolland(self.f, self.rMax, self.rMax2, self.beta1, self.beta2, dp1, dp2, self.rho)
#        aa = (d2Vm/2 - (-vMax/self.rMax)/self.rMax) / self.rMax
#        bb = (d2Vm - 6*aa*self.rMax) / 2
#        cc = -3*aa*self.rMax**2 - 2*bb*self.rMax

        # Replace all values within rMax of the storm centre with the cubic
        # profile to eliminate barotropic instability:
#        if self.dP >= 1500.:
#            icore = where(self.R <= self.rMax)
#            P[icore] = sign(self.f)*self.R[icore]*(self.R[icore]*(self.R[icore]*aa + bb) + cc)
        self.logger.debug("Timing for doubleHolland wind profile calculation: %.3f" % (time.time()-t0))
        return P

    def powell(self):
        """
        Powell et al, 2005
        Another definition of the B parameter inserted into the Holland
        model.  Unlike Willoughby and Rahn's model, there is no reliance
        on vMax.  Powell et al. also included a small random term, but
        since the beta value is also used in the vorticity calculation,
        we need to ensure the values used in this function and the
        corresponding vorticity function match.
        """

        beta = 1.881093 - 0.010917*abs(self.cLat) - 0.005567*self.rMax

        # Include the censoring of beta to lie in the interval 0.8 - 2.2:
        if beta < 0.8:
            beta = 0.8
        elif beta > 2.2:
            beta = 2.2

        P = self.holland(beta)
        return P
