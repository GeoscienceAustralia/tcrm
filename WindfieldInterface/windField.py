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


Title: windField.py - determine the asymmetric wind field around a
       cyclone

Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2006-11-30
Description:
Determine the wind field for a given cyclone position, intensity, speed,
bearing and size. This applies the asymmetry corrections to a radial
wind profile and returns x- and y-components (east and north components)
 of the wind field.
Available options:
Hubbert
McConochie
Kepert
The Hubbert algorithm includes an arbitrary surface wind reduction
factor (SWRF).  The McConochie SWRF is described in Harper et al. 2001.

SeeAlso: (related programs)
Constraints:
Version: $Rev: 810 $

Version:
ModifiedBy:
ModifiedDate: yyyy-mm-dd

Reference:
Harper, B.A, T.A. Hardy, L.B. Mason, L. Bode, I.R. Young and P. Nielsen, 2001:
    Queensland climate change and community vulnerability to tropical
    cyclones, ocean hazards assessment. Stage 1 report.
    Dept. of Natural Resources and Mines, Queensland, Brisbane,
    Australia, 368 pp.
Hubbert, G.D., G.J. Holland, L.M. Leslie and M.J. Manton, 1991:
    A Real-Time System for Forecasting Tropical Cyclone Storm Surges.
    Weather and Forecasting, 6, 86-97
Kepert, J., 2001:
    The Dynamics of Boundary Layer Jets within the Tropical Cyclone
    Core. Part I: Linear Theory.
    J. Atmos. Sci., 58, 2469-2484
McConochie, J.D., T.A. Hardy and L.B. Mason, 2004:
    Modelling tropical cyclone over-water wind and pressure fields.
    Ocean Engineering, 31, 1757-1782

$Id: windField.py 810 2012-02-21 07:52:50Z nsummons $
"""

import os, sys, pdb, logging

from numpy import *
import time

logger = logging.getLogger()

class WindField:
    """
    Description: Determine the wind field for a given cyclone position,
                 intensity, speed, bearing and size. This applies the
                 asymmetry corrections to a radial wind profile and
                 returns x- and y-components (east and north components)
                 of the wind field.

    Parameters:
        R: grid of distances from the storm centre (distances in km)
        lam: grid of bearings from the storm centre
        rMax: Radius of maximum winds (km)
        f: Coriolis parameter (1/s)
        V: grid of radial wind speed (m/s)
        Z: grid of radial vorticity (1/s)
        vFm: Forward speed of storm (m/s)
        thetaFm: Forward direction of storm
        thetaMax: angle from direction of storm to maximum wind

    Members:
        R: grid of distances from the storm centre (distances in km)
        lam: grid of bearings from the storm centre
        rMax: Radius of maximum winds (km)
        f: Coriolis parameter (1/s)
        V: grid of radial wind speed (m/s)
        Z: grid of radial vorticity (1/s)
        vFm: Forward speed of storm (m/s)
        thetaFm: Forward direction of storm
        thetaMax: angle from direction of storm to maximum wind

    Methods:
        hubbert:
        mcconochie:
        kepert:

    Internal Methods:

    """

    def __init__(self, R, lam, rMax, f, V, Z, vFm, thetaFm, thetaMax=None):
        """
        Initialise required fields:
        """
        self.R = R
        self.lam = lam
        self.rMax = rMax
        self.f = f
        self.V = V
        self.Z = Z
        self.vFm = vFm
        self.thetaFm = thetaFm
        self.thetaMax = thetaMax


    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return 'Determine the wind field for a given cyclone position, intensity, speed, bearing\
            and size. This applies the asymmetry corrections to a radial wind profile and\
            returns x- and y-components (east and north components) of the wind field.'


    def hubbert(self):
        """
        Hubbert, G.D., G.J. Holland, L.M. Leslie and M.J. Manton, 1991:
        A Real-Time System for Forecasting Tropical Cyclone Storm Surges.
        Weather and Forecasting, 6, 86-97
        """
        t0 = time.time()
        Km = .70
        inflow = 25.*ones(shape(self.R))
        core = where(self.R < self.rMax)
        inflow[core] = 0
        inflow = inflow*pi/180

        thetaMaxAbsolute = self.thetaFm + self.thetaMax
        asym = self.vFm*cos(thetaMaxAbsolute - self.lam+pi)
        Vsf = Km*self.V + asym
        phi = inflow - self.lam

        Ux = Vsf*sin(phi)
        Vy = Vsf*cos(phi)
        logger.debug("Timing for wind field: %.3f" %(time.time()-t0))
        return Ux, Vy

    def mcconochie(self):
        """
        McConochie, J.D., T.A. Hardy and L.B. Mason, 2004:
        Modelling tropical cyclone over-water wind and pressure fields.
        Ocean Engineering, 31, 1757-1782
        """
        t0=time.time()

        inflow = 25.*ones(shape(self.R))
        mid = where(self.R < 1.2*self.rMax)
        inflow[mid] = 10. + 75.*(self.R[mid]/self.rMax-1.)
        inner = where(self.R < self.rMax)
        inflow[inner] = 10.*self.R[inner]/self.rMax
        inflow = inflow*pi/180.

        thetaMaxAbsolute = self.thetaFm + self.thetaMax
        phi = inflow - self.lam

        asym = 0.5*(1.+cos(thetaMaxAbsolute - self.lam))*self.vFm*(self.V/abs(self.V).max())
        Vsf = self.V + asym

        # Surface wind reduction factor:
        swrf = 0.81*ones(shape(Vsf))
        low = where(Vsf >= 6)
        med = where(Vsf >= 19.5)
        high = where(Vsf >= 45)
        swrf[low] = 0.81-(2.93*(Vsf[low] - 6.)/1000.)
        swrf[med] = 0.77-(4.31*(Vsf[med] - 19.5)/1000.)
        swrf[high] = 0.66

        Ux = swrf*Vsf*sin(phi)
        Vy = swrf*Vsf*cos(phi)
        logger.debug("Timing for wind field: %.3f" %(time.time()-t0))

        return Ux, Vy

    def kepert(self):
        """
        Kepert, J., 2001:
        The Dynamics of Boundary Layer Jets within the Tropical Cyclone
        Core. Part I: Linear Theory.
        J. Atmos. Sci., 58, 2469-2484
        """
        t0 = time.time()
        K = 50.    # Diffusivity

        # Constant drag coefficient.
        Cd = 0.002

        # Large and Pond definition of Cd
        #Cd = (0.61 + 0.063*self.V)/1000

        # Capped version of Large and Pond:
        # I arbitarily impose the cap at 25 m/s (CA)
        #i = where(self.V>25.)
        #Cd = (0.61 + 0.063*self.V)/1000.
        #Cd[i] = (0.61 + 0.063*25.)/1000.

        al = (2. * self.V / self.R + self.f) / (2.*K)
        be = (self.f + self.Z) / (2.*K)
        gam = self.V / (2.*K*self.R)
        albe = sqrt(al/be)

        ind = where(abs(gam) > sqrt(al*be))
        chi = (Cd/K)*self.V / sqrt(sqrt(al*be))
        eta = (Cd/K)*self.V / sqrt(sqrt(al*be) + abs(gam))
        psi = (Cd/K)*self.V / sqrt(abs(sqrt(al*be) - (gam)))
        i = complex(0.,1.)
        A0 =  -chi*self.V*(1. + i*(1. + chi)) / (2.*chi**2. + 3.*chi + 2.)

        # Symmetric surface wind component:
        u0s = albe * A0.real
        v0s =        A0.imag

        Am = -((1. + (1.+i)*eta)/albe   + (2. + (1.+i)*eta))*psi*self.vFm \
             / ((2. + 2.*i)*(1 + eta*psi) + 3.*psi + 3.*i*eta)

        Am[ind] = -((1. + (1.+i)*eta[ind])/albe[ind] + \
                    (2. + (1.+i)*eta[ind]))*psi[ind]*self.vFm \
                  /(2. - 2.*i + 3.*(psi[ind] + eta[ind]) + \
                    (2. + 2.*i)*eta[ind]*psi[ind])
        # First asymmetric surface component:
        ums = albe * (Am * exp(-i*self.lam)).real
        vms =        (Am * exp(-i*self.lam)).imag

        Ap = -((1. + (1.+i)*psi)/albe - (2. + (1.+i)*psi))*eta*self.vFm \
             / ((2. + 2.*i)*(1. + eta*psi) + 3.*eta + 3.*i*psi)
        Ap[ind] = -((1. + (1.-i)*psi[ind])/albe[ind] - \
                    (2. + (1.-i)*psi[ind]))*eta[ind]*self.vFm \
                  /(2. + 2.*i + 3.*(eta[ind] + psi[ind]) + \
                    (2. - 2.*i)*eta[ind]*psi[ind])
        # Second asymmetric surface component:
        ups = albe * (Ap * exp(i*self.lam)).real
        vps =        (Ap * exp(i*self.lam)).imag

        # Total surface wind in (moving coordinate system):
        us =          u0s + ups + ums
        vs = self.V + v0s + vps + vms

        usf = us + self.vFm*cos(self.lam-self.thetaFm)
        vsf = vs - self.vFm*sin(self.lam-self.thetaFm)
        #Uf =   + vFm*cos(thetaFm)
        #Vf = V - vFm*sin(thetaFm)
        phi = arctan2(usf,vsf)

        # Surface winds, cartesian coordinates:
        Ux = (sqrt(usf**2.+vsf**2.)*sin(phi - self.lam))
        Vy = (sqrt(usf**2.+vsf**2.)*cos(phi - self.lam))

        # Gradient level winds, cartesian coordinates:
        #Ugx = (Uf**2+Vf**2)*cos(lam-phi)
        #Vgx = (Uf**2+Vf**2)*sin(lam-phi)
        logger.debug("Timing for wind field: %.3f" %(time.time()-t0))

        return Ux,Vy
