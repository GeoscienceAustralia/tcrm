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


Title: windProfile.py - radial wind profile for a given instance of a
       cyclone

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-11-20
Description: Return the radial velocity field for a range of wind
             profiles.
The available wind profiles are:
Rankine vortex
Jelesnianski
Holland (with cubic core)
Schloemer (just the Holland with beta = 1)
Willoughby and Rahn (Holland with beta = beta(vmax, rmax, latitude))
McConochie et al. (double exponential profile)
Powell et al. (Holland with beta = beta(rmax, latitude))

For the Holland and double Holland, the region inside rMax is modified
to a cubic function of R to avoid the barotropic instability noted by
Kepert (2001).

SeeAlso: windVorticity.py, derivative.py
Constraints:
Version: $Rev: 649 $

Version: 573
ModifiedBy: Craig Arthur
ModifiedDate: 2006-12-04
Modification: Now calculates V correctly, accounting for the change of
sign in V due to hemisphere [note the abs(self.f) references].

Version: 573
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-04-04
Modification: dVm variable removed - dVm is zero at rMax

Version :$Rev: 649 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-01-06
Modification: Added Powell et al (2005) relation for beta

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
    State of Florida hurricane loss projection model: Atmospheric
    science component.
    Journal of Wind Engineering and Industrial Aerodynamics, 93 (8),
    651-674
Schloemer, R.W., 1954:
    Analysis and synthesis of hurricane wind patterns over Lake
    Okeechobee.
    NOAA Hydromet. Rep. 31, 49 pp.
Willoughby, H.E. and M.E. Rahn, 2004:
    Parametric Representation of the Primary Hurricane Vortex. Part I:
    Observations and Evaluation of the Holland (1980) Model.
    Mon. Wea. Rev., 132, 3033-3048

$Id: windProfile.py 649 2011-10-31 05:35:00Z nsummons $
"""
#-----------------------------------------------------------------------
# Imports:
#-----------------------------------------------------------------------
import os, sys, pdb, logging

import Utilities.metutils as metutils
import derivatives
import numpy
import vmax
import time

class WindProfile:
    """
    Description: Define the radial wind profiles used in tropical
                 cyclone modelling. These are radial profiles only and
                 do not include asymmetries that arise due to the
                 forward motion of the storm.

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
        rankine: Rankine vortex
        jelesnianski: Jelesnianski's storm surge model wind field
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
        self.logger.debug("Coriolis parameter: %3f" %self.f)


    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Generate the radial wind profile for a given instance of a \
                tropical cyclone. \
                Profiles available are: \
                Rankine vortex \
                Jelesnianski \
                Holland \
                Schloemer (a simplification of the Holland profile) \
                Willoughby & Rahn (a more complex version of the Holland profile) \
                McConochie (double Holland vortex)\
                For the first two, the maximum wind speed is required - this can \
                be calculated using vmax.py '

    def rankine(self, vMaxType="willoughby"):
        """
        Rankine vortex profile. Vmax determined by dp using, by default,
        the Willoughby & Rahn method.
        """
        t0 = time.time()
        vMax = vmax.vmax(self.pCentre, self.pEnv, vMaxType)
        # An assumption about the shape of the profile outside Rmax.
        # The literature indicates 0.4 < alpha < 0.6 (e.g. see Holland,
        # 1980)
        alpha = 0.5
        V = vMax*(self.rMax/self.R)**alpha
        icore = numpy.where(self.R <= self.rMax)
        V[icore] = vMax*(self.R[icore]/self.rMax)
        V = numpy.sign(self.f)*V
        self.logger.debug("Timing for rankine wind profile calculation: %.3f" %(time.time()-t0))
        return V

    def jelesnianski(self, vMaxType="willoughby"):
        """
        Jelesnianski model of the wind profile
        """
        t0 = time.time()
        vMax = vmax.vmax(self.pCentre, self.pEnv, vMaxType)
        V = 2*vMax*self.rMax*self.R/(self.rMax**2 + self.R**2)
        V = numpy.sign(self.f)*V
        self.logger.debug("Timing for jelesnianski wind profile calculation: %.3f" %(time.time()-t0))
        return V

    def holland(self, beta=None):
        """
        Holland profile. For r < rMax, we reset the wind field to a
        cubic profile to avoid the barotropic instability mentioned in
        Kepert & Wang (2001).
        Air density (rho) is assumed to be 1.15 kg/m**3
        """
        if beta==None:
            beta=self.beta
        t0 = time.time()
        V = numpy.zeros(self.R.shape)
        vMax = vmax.vmax(self.pCentre, self.pEnv, "holland", beta)

        # Calculate first and second derivatives at R = Rmax:
        d2Vm = derivatives.holland(self.f, self.rMax, beta,
                                   self.dP, self.rho)
        aa = (d2Vm/2. - (-vMax/self.rMax)/self.rMax) / self.rMax
        bb = (d2Vm - 6*aa*self.rMax) / 2.
        cc = -3*aa*self.rMax**2 - 2*bb*self.rMax

        # The profile as originally stated:
        delta = (self.rMax/self.R)**beta
        edelta = numpy.exp(-delta)

        V = numpy.sign(self.f)*numpy.sqrt((self.dP*beta/self.rho)*delta*edelta + (self.R*self.f/2.)**2) \
            - self.R*numpy.abs(self.f)/2.

        # Replace all values within rMax of the storm centre with the
        # cubic profile to eliminate barotropic instability:
        icore = numpy.where(self.R <= self.rMax)
        V[icore] = numpy.sign(self.f) * self.R[icore] \
                   * (self.R[icore]*(self.R[icore]*aa + bb) + cc)
        self.logger.debug("Timing for holland wind profile calculation: %.3f" %(time.time()-t0))
        return V

    def willoughby(self):
        """
        The Willoughby & Rahn (2004) relation, which makes beta a
        function of Vmax, rMax and latitude. We use Willoughby & Rahn's
        (2004) relation for Vmax *only*.
        This determines the beta parameter then calls Holland (which
        means the profile is cubic within Rmax) to calculate the wind
        profile.  The beta term calculation is based on Atlantic and
        Eastern Pacific cyclone data, not Australian data.
        """
        vMax = vmax.vmax(self.pCentre, self.pEnv, type="willoughby")
        beta = 1.0036 + 0.0173*vMax - 0.313*numpy.log(self.rMax) \
               + 0.0087*numpy.abs(self.cLat)
        V = self.holland(beta)
        return V

    def schloemer(self):
        """
        Schloemer's (1954) is the same as the Holland relation with
        beta = 1
        """
        beta = 1.
        V = self.holland(beta)
        return V

    def doubleHolland(self, rMax2=250.):
        """
        McConochie et al's double Holland vortex model (based on Cardone
        et al, 1994).  This application is the Coral Sea adaptation of
        the double vortex model (it can also be used for concentric
        eye-wall configurations).  The tunable parameters in this
        relation are 'dp1', 'dp2', 'b1', 'b2' and 'rMax2'
        """
        t0=time.time()
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

        rMax1 = self.rMax

        # The two gradient wind components:
        mu = (self.rMax/self.R)**self.beta1
        nu = (self.rMax2/self.R)**self.beta2
        emu = numpy.exp(-mu)
        enu = numpy.exp(-nu)

        gradientV1 = (self.beta1*dp1/self.rho)*mu*emu
        gradientV2 = (self.beta2*dp2/self.rho)*nu*enu

        V = numpy.sign(self.f)*numpy.sqrt(gradientV1+gradientV2+(self.R*self.f/2.)**2) \
            - self.R*numpy.abs(self.f)/2.

        vMax = numpy.abs(V).max()
        #aa, bb, cc = doubleHollandCoefficient(vMax, rMax1, rMax2, dp1, dp2, beta1, beta2, f, rho)

        #delta = (self.rMax/self.R)**self.beta1
        #gamma = (self.rMax2/self.R)**self.beta2

        # Calculate first and second derivatives at R = Rmax:
        d2Vm = derivatives.doubleHolland(self.f, self.rMax, self.rMax2,
                                         self.beta1, self.beta2, dp1, dp2,
                                         self.rho)
        aa = (d2Vm/2. - (-vMax/self.rMax)/self.rMax) / self.rMax
        bb = (d2Vm - 6*aa*self.rMax) / 2.
        cc = -3*aa*self.rMax**2 - 2*bb*self.rMax

        # Replace all values within rMax of the storm centre with the cubic
        # profile to eliminate barotropic instability:
        if self.dP >= 1500.:
            icore = numpy.where(self.R <= self.rMax)
            V[icore] = numpy.sign(self.f) * self.R[icore] \
                       * (self.R[icore]*(self.R[icore]*aa + bb) + cc)
        self.logger.debug("Timing for doubleHolland wind profile calculation: %.3f" %(time.time()-t0))
        return V

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

        beta = 1.881093 - 0.010917*numpy.abs(self.cLat) - 0.005567*self.rMax

        # Include the censoring of beta to lie in the interval 0.8 - 2.2:
        if beta < 0.8:
            beta = 0.8
        if beta > 2.2:
            beta = 2.2
        self.logger.debug("Beta parameter: %4.2f"%beta)
        V = self.holland(beta)
        return V

    def newHolland(self,r_gale=150.):
        """
        Holland et al. 2010
        In this version, the exponent is allowed to vary linearly
        outside the radius of maximum wind. i.e. rather than take the
        sqare root, the exponent varies around 0.5.
        Currently this version does not have a corresponding vorticity
        profile set up in windVorticity, so it cannot be applied in
        wind field modelling.
        """
        # In this incarnation, we are assuming the pressure rate of change and forward velocity
        # is zero, so there is no requirement for a first-pass
        # guess at x
        Bs = -0.000044*(self.dP/100.)**2. + 0.01*(self.dP/100.) - 0.014*numpy.abs(self.cLat) + 1.0
        deltag = numpy.power(self.rMax/r_gale,Bs)
        edeltag = numpy.exp(-1.*deltag)
        rgterm = Bs*self.dP*deltag*edeltag/self.rho
        xn = numpy.log(17.)/numpy.log(rgterm)
        xx = 0.5*numpy.ones(self.R.shape)
        i = numpy.where(self.R>self.rMax)
        xx[i] = 0.5 + (self.R[i] - self.rMax)*(xn - 0.5)/(r_gale - self.rMax)

        delta = (self.rMax/self.R)**Bs
        edelta = numpy.exp(-delta)
        pdb.set_trace()
        V = numpy.sign(self.f)*numpy.power((self.dP*Bs/self.rho)*delta*edelta,xx)
        return V
