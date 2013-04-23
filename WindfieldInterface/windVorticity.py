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


Title: windVorticity.py - the radial vorticity field for a given
       instance of a cyclone.
Author: Craig Arthur
Email: craig.arthur@ga.gov.au
CreationDate: 2006-11-20
Description: Return the radial vorticity field for a range of wind
profiles.  The Rankine and Jelesnianski profiles are functions of vMax,
for consistency, vMax is passed to all vorticity profiles.
Timing information for the Willoughby and Schloemer profiles is reported
as Holland, as these two profiles call the Holland profile with a
different beta value. This prevents the timing being reported twice.

SeeAlso: windProfile.py, derivative.py
Constraints:


Version: 88
ModifiedBy: Craig Arthur
ModifiedDate: 2006-12-04
Modification: Now calculates V/r + dV/dr correctly, accounting for the
change of sign in V due to hemisphere [note the abs(self.f) references].

Version :$Rev: 810 $
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
    Journal of Wind Engineering and Industrial Aerodynamics, 93 (8), 651-674
Schloemer, R.W., 1954:
    Analysis and synthesis of hurricane wind patterns over Lake
    Okeechobee.
    NOAA Hydromet. Rep. 31, 49 pp.
Willoughby, H.E. and M.E. Rahn, 2004:
    Parametric Representation of the Primary Hurricane Vortex. Part I:
    Observations and Evaluation of the Holland (1980) Model.
    Mon. Wea. Rev., 132, 3033-3048

$Id: windVorticity.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

import time
import numpy
import Utilities.metutils as metutils
import derivatives

class WindVorticity:
    """
    Description: Return the radial vorticity field for a range of wind profiles.
        This corresponds to wind_profile.py
        The Rankine and Jelesnianski profiles are functions of vMax,
        for consistency, vMax is passed to all vorticity profiles.

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
                 rMax2=250., beta1=None, beta2=None):
        """
        Initialise required fields:
        """
        self.R = R
        self.rMax= rMax
        self.cLon = cLon
        self.cLat = cLat
        self.dP = pEnv-pCentre
        self.pCentre = pCentre
        self.pEnv = pEnv
        # Density of air:
        self.rho = 1.15
        self.f = metutils.coriolis(cLat)
        self.beta = beta
        self.rMax2 = rMax2
        self.beta1 = beta1
        self.beta2 = beta2
        self.logger = logging.getLogger()

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Generate the radial vorticity profile for a given instance of a \
                tropical cyclone. \
                Profiles available are: \
                Rankine vortex \
                Jelesnianski \
                Holland \
                Schloemer (a simplification of the Holland profile) \
                Willoughby & Rahn (a more complex version of the Holland profile)\
                McConochie (double Holland vortex)\
                For the first two, the maximum wind speed is required - this can \
                be calculated using vmax.py '

    def rankine(self, vMax):
        """
        Rankine vortex profile. Vmax determined by dp, using the
        default Willoughby & Rahn method.
        """
        t0=time.time()
        alpha = 0.5
        Z = numpy.sign(self.f)*vMax*((self.rMax/self.R)**alpha)/self.R - \
            alpha*vMax*(self.rMax**alpha)/(self.R**alpha)

        icore = numpy.where(self.R <= self.rMax)
        Z[icore] = numpy.sign(self.f) \
                   * (vMax*(self.R[icore]/self.rMax) + vMax/self.rMax)
        self.logger.debug("Timing for rankine vorticity calculation: %.3f" %(time.time()-t0))
        return Z

    def jelesnianski(self, vMax):
        """
        Jelesnianski model of the wind profile:
        """
        t0=time.time()
        Z = numpy.sign(self.f)*2*vMax*self.rMax/(self.rMax**2 + self.R**2) + \
            numpy.sign(self.f)*2*vMax*self.rMax*(self.rMax**2 - self.R**2) \
            / (self.rMax**2 + self.R**2)**2
        self.logger.debug("Timing for jelesnianski vorticity calculation: %.3f" %(time.time()-t0))
        return Z

    def holland(self, vMax, beta=None):
        """
        Holland profile. For r < rMax, we reset the wind field to a
        cubic profile to avoid the barotropic instability mentioned in
        Kepert & Wang (2001).:
        """
        if beta==None:
            beta=self.beta
        t0=time.time()
        delta = (self.rMax/self.R)**beta
        edelta = numpy.exp(-delta)

        Z = numpy.sign(self.f)*(numpy.sqrt((self.dP*beta/self.rho)*delta*edelta + (self.R*self.f/2.)**2))/self.R \
            - numpy.abs(self.f) \
            + edelta*(2*(beta**2)*self.dP*(delta-1)*delta + self.rho*edelta*(self.f*self.R)**2) \
            / (2*self.rho*self.R*numpy.sqrt(4*(beta*self.dP/self.rho)*delta*edelta + (self.f*self.R)**2))

        # Calculate first and second derivatives at R = Rmax:
        d2Vm = derivatives.holland(self.f, self.rMax, beta,
                                   self.dP, self.rho)
        aa = (d2Vm/2 - (-1.0*numpy.sign(self.f)*vMax/self.rMax)/self.rMax) \
             / self.rMax
        bb = (d2Vm - 6*aa*self.rMax) / 2
        cc = -3*aa*self.rMax**2 - 2*bb*self.rMax

        icore = numpy.where(self.R <= self.rMax)
        Z[icore] = self.R[icore] * (self.R[icore] * 4*aa + 3*bb) + 2*cc
        self.logger.debug( "Timing for holland vorticity calculation: %.3f" %(time.time()-t0) )
        return Z

    def willoughby(self,vMax):
        """
        The Willoughby & Rahn (2004) relation, which makes beta a
        function of Vmax, rMax and latitude. We use Willoughby & Rahn's
        (2004) relation for Vmax *only*.
        This determines the beta parameter then calls holland to
        calculate the vorticity profile.
        """

        beta = 1.0036 + 0.0173*vMax  - 0.313*numpy.log(self.rMax) \
               + 0.0087*numpy.abs(self.cLat)
        Z = self.holland(vMax, beta=beta)
        return Z

    def schloemer(self,vMax):
        """
        Schloemer's (1954) is the same as the Holland relation with
        beta = 1
        """
        beta = 1.
        Z = self.holland(vMax, beta=beta)
        return Z

    def doubleHolland(self, vMax):
        """
        McConochie et al's double Holland vortex model (based on Cardone
        et al, 1994).  This application is the Coral Sea adaptation of
        the double vortex model (it can also be used for concentric
        eye-wall configurations).  The tunable parameters in this
        relation are 'dp1', 'dp2', 'b1', 'b2' and 'rMax2'
        """
        #pdb.set_trace()
        t0=time.time()

        # Scale dp2 if dP is less than 1500 Pa:
        if self.dP < 1500.:
            dp2 = (self.dP/1500.)*(800. + (self.dP - 800.)/2000.)
        else:
            dp2 = 800. + (self.dP - 800.)/2000.
        dp1 = self.dP - dp2
        if self.beta1 is None:
            self.beta1 = 7.3 - self.pCentre/16000.

        if self.beta2 is None:
            self.beta2 = 7.2 - self.pCentre/16000.

        chi_ = self.beta1*dp1/self.rho
        psi_ = self.beta2*dp2/self.rho

        delta_ = (self.rMax/self.R)**self.beta1
        gamma_ = (self.rMax2/self.R)**self.beta2
        edelta_ = numpy.exp(-delta_)
        egamma_ = numpy.exp(-gamma_)

        # Derivatives:
        ddelta_ = -self.beta1*(self.rMax**self.beta1)/(self.R**(self.beta1+1))
        dgamma_ = -self.beta2*(self.rMax2**self.beta2)/(self.R**(self.beta2+1))

        """
Z = (numpy.sqrt((b1*dp1/rho)*delta*edelta + (b2*dp2/rho)*gamma*egamma + (r*f/2)**2) - r*f/2)/r \
    -f/2 + (-4*(b1**2)*(dp1/rho)*(delta/r)*edelta + 4*(b1**2)*(dp1/rho)*((delta**2)/r)*edelta - \
    4*(b2**2)*(dp2/rho)*(gamma/r)*egamma + 4*(b2**2)*(dp2/rho)*((gamma**2)/r)*egamma+2*r*f**2)/ \
    (4*numpy.sqrt(4*(b1*dp1/rho)*delta*edelta + 4*(b2*dp2/rho)*gamma*egamma+(r*f)**2))
    """

        # Vorticity:

        Z = numpy.sign(self.f)*numpy.sqrt(chi_*delta_*edelta_ + psi_*gamma_*egamma_ + (self.f*self.R/2)**2)/self.R \
            - numpy.abs(self.f) + \
            (1/2)*(chi_*ddelta_*edelta_*(1-delta_) + psi_*dgamma_*egamma_*(1-gamma_) + self.R*self.f**2) \
             / numpy.sqrt(chi_*delta_*edelta_ + psi_*gamma_*egamma_ + (self.f*self.R/2)**2)
        """
        Z = (numpy.sqrt((beta1*dp1/self.rho)*delta*edelta + \
            (beta2*dp2/self.rho)*gamma*egamma+(self.R*self.f/2)**2)-self.R*numpy.abs(self.f)/2)/self.R \
            - numpy.abs(self.f)/2 + (-4*(beta1**2)*(dp1/(self.rho*self.R))*delta*edelta + \
            4*(beta1**2)*(dp1/(self.rho*self.R))*(delta**2)*edelta - \
            4*(beta2**2)*(dp2/(self.rho*self.R))*gamma*egamma + \
            4*(beta2**2)*(dp2/(self.rho*self.R))*(gamma**2)*egamma + 2*self.R*(self.f**2))/ \
            (4*numpy.sqrt(4*(beta1*dp1/self.rho)*delta*edelta + \
            4*(beta2*dp2/self.rho)*gamma*egamma+(self.R*self.f)**2))
        """
        
        # Calculate first and second derivatives at R = Rmax:
        d2Vm = derivatives.doubleHolland(self.f, self.rMax, self.rMax2,
                                         self.beta1, self.beta2, dp1,dp2,
                                         self.rho)
        aa = (d2Vm/2.0 - (-1.0*numpy.sign(self.f)*vMax/self.rMax)/self.rMax) \
             / self.rMax
        bb = (d2Vm - 6.0*aa*self.rMax) / 2.0
        cc = -3.0*aa*self.rMax**2.0 - 2.0*bb*self.rMax
        if self.dP >= 1500.:
            icore = numpy.where( self.R <= self.rMax)
            Z[icore] = self.R[icore] * (self.R[icore] \
                       * 4.0*aa + 3.0*bb) + 2.0*cc
        self.logger.debug("Timing for doubleHolland vorticity calculation: %.3f" %(time.time()-t0))
        return Z

    def powell(self,vMax):
        """
        Powell et al, 2005
        Another definition of the B parameter inserted into the Holland
        model.  Unlike Willoughby and Rahn's model, there is no reliance
        on vMax.  Powell et al. also included a small random term, but
        since the beta value is also used in the radial profile
        calculation, we need to ensure the values used in this function
        and the corresponding profile function match.
        """

        beta = 1.881093 - 0.010917*numpy.abs(self.cLat) - 0.005567*self.rMax

        # Include the censoring of beta to lie in the interval 0.8 - 2.2:
        if beta < 0.8:
            beta = 0.8
        if beta > 2.2:
            beta = 2.2

        Z = self.holland(vMax, beta=beta)
        return Z
