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


Title: derivatives.py

Author: Craig Arthur
Email: craig.arthur@ga.gov.au
CreationDate: 2006-11-29
Description: Calculates required 1st and 2nd derivatives at r = rMax.
Allows for a smooth transition between the actual profile and the cubic
core.


Version: $Rev: 646 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-04-04
Modification: Weaved doubleHolland routine. Also eliminated need to
              calculate first derivative, as it's (by definition) zero
              at rMax.

SeeAlso: (related programs)
Constraints:

$Id: derivatives.py 646 2008-02-07 22:58:51Z nhabili $
"""

import os, sys, pdb, logging

from scipy import exp, sqrt, empty
#from scipy import weave
#from scipy.weave import converters

def holland(f, rMax, beta, dP, rho):
    """
    Return the first and second derivatives of the Holland profile at
    r = rMax
    """
#    E = exp(1)
#    dVm = -abs(f)/2 + (E*(f**2)*rMax*sqrt((4*beta*dP/rho)/E + (f*rMax)**2))/ \
#           (2*(4*beta*dP/rho + E*(f*rMax)**2))
#    d2Vm = (beta*dP*(-4*beta**3*dP/rho - (-2 + beta**2)*E*(f*rMax)**2)) / \
#           (E*rho*sqrt((4*beta*dP)/(E*rho) + (f*rMax)**2) * \
#           (4*beta*dP*rMax**2/rho + E*(f*rMax**2)**2))
#
#    return dVm, d2Vm
  
#    Note: C Code disabled as offering only negligible overall speed improvement
#    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#    d = empty(1, 'd');
#
#    code = """
#
#    double f_ = f;
#    double rMax_ = rMax;
#    double beta_ = beta;
#    double dP_ = dP;
#    double rho_ = rho;
#    double E = exp(1);
#
#    d(0) = (beta_*dP_*(-4*pow(beta_, 3)*dP_/rho_ - (-2 + pow(beta_, 2))*E*pow(f_*rMax_, 2))) /
#           (E*rho_*sqrt((4*beta_*dP_)/(E*rho_) + pow(f_*rMax_, 2)) *
#           (4*beta_*dP_*pow(rMax_, 2)/rho_ + E*pow(f_*pow(rMax_, 2), 2)));
#    """
#    err = weave.inline(code,
#                      ['f', 'rMax', 'beta', 'dP', 'rho', 'd'],
#                      type_converters=converters.blitz,
#                      compiler = 'gcc')
#    #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                      
    E = exp(1)
    d2Vm = (beta*dP*(-4*beta**3*dP/rho - (-2 + beta**2)*E*(f*rMax)**2)) / \
           (E*rho*sqrt((4*beta*dP)/(E*rho) + (f*rMax)**2) * \
           (4*beta*dP*rMax**2/rho + E*(f*rMax**2)**2))
                      
    if d2Vm >= 0.0:
        raise ValueError
    else:
        return d2Vm

def doubleHolland(f, rMax1, rMax2, beta1, beta2, dp1, dp2, rho):
    """
    Return the first and second derivatives of the double Holland profile
    (as described by McConochie et al) at r = rMax1
    """

    """
    delta = (rMax1/self.R)**beta1
    gamma = (rMax2/self.R)**beta2
    edelta = exp(-delta)
    egamma = exp(-gamma)
    The wind profile:
    V = sqrt(beta1*dp1/rho*delta*exp(-delta)+beta2*dp2/rho*gamma*exp(-gamma)+(self.R*self.f/2)**2)-self.R*self.f/2

    First derivative:
    dV = -self.f/2 + \
        (-4*beta1**2*dp1/self.rho*delta/self.R*exp(-delta) + \
        (4*beta1**2*dp1/self.rho)*((delta**2)/self.R)*exp(-delta) - \
        (4*beta2**2*dp2/self.rho)*gamma/self.R*exp(-gamma) + \
        (4*beta2**2*dp2/self.rho)*((gamma**2)/self.R)*exp(-gamma)+2*self.R*self.f**2)/ \
        (4*sqrt(4*beta1*dp1/self.rho*delta*exp(-delta) + \
        (4*beta2*dp2/self.rho)*gamma*exp(-gamma)+(self.R*self.f)**2))
    dVm = -self.f/2 + \
            (-4*beta1**2*dp1/(self.rho*rMax1*E) + \
            4*beta1**2*dp1/(self.rho*rMax1*E) - \
            (4*beta2**2*dp2/self.rho)*(nu/rMax1)*exp(-nu) + \
            (4*beta2**2*dp2/self.rho)*((nu**2)/rMax1)*exp(-nu)+2*rMax1*self.f**2)/ \
            (4*sqrt((4*beta1*dp1/self.rho*E) + \
            (4*beta2*dp2/self.rho)*nu*exp(-nu)+(rMax1*self.f)**2))


    Second derivative:
    d2V = -1/(8*((4*beta1*dp1/rho)*((rMax1/self.R)**beta1)*exp(-(rMax1/self.R)**beta1) + \
            (4*beta2*dp2/rho)*((rMax2/self.R)**beta2)*exp(-(rMax2/self.R)**beta2) + self.R**2*self.f**2)**(1.5)) * \
            (-4*(beta1**2)*dp1/rho*((rMax1/self.R)**beta1/self.R)*exp(-(rMax1/self.R)**beta1) + \
            (4*(beta1**2)*dp1/rho)*((((rMax1/self.R)**beta1)**2)/self.R)*exp(-(rMax1/self.R)**beta1) - \
            (4*(beta2**2)*dp2/rho)*(((rMax2/self.R)**beta2)/self.R)*exp(-(rMax2/self.R)**beta2) + \
            (4*(beta2**2)*dp2/rho)*((((rMax2/self.R)**beta2)**2)/self.R)*exp(-(rMax2/self.R)**beta2) + 2*self.R*self.f**2)**2 + \
            1/(4*sqrt((4*beta1*dp1/rho)*(rMax1/self.R)**beta1*exp(-(rMax1/self.R)**beta1) + \
            (4*b2*dp2/rho)*(rMax2/self.R)**b2*2+exp(-(rMax2/self.R)**b2) + (self.R*self.f)**2)) * \
            ((4*(beta1**3)*dp1/rho)*(((rMax1/self.R)**beta1)/(self.R**2))*exp(-(rMax1/self.R)**beta1) + \
            (4*(beta1**2)*dp1/rho)*(((rMax1/self.R)**beta1)/(self.R**2))*exp(-(rMax1/self.R)**beta1) - \
            (12*(beta1**3)*dp1/rho)*((((rMax1/self.R)**beta1)**2)/(self.R**2))*exp(-(rMax1/self.R)**beta1) - \
            (4*(beta1**2)*dp1/rho)*((((rMax1/self.R)**beta1)**2)/(self.R**2))*exp(-(rMax1/self.R)**beta1) + \
            (4*(beta1**3)*dp1/rho)*((((rMax1/self.R)**beta1)**3)/(self.R**2))*exp(-(rMax1/self.R)**beta1) + \
            (4*(beta2**3)*dp2/rho)*(((rMax2/self.R)**beta2)/(self.R**2))*exp(-(rMax2/self.R)**beta2) + \
            (4*(beta2**2)*dp2/rho)*(((rMax2/self.R)**beta2)/(self.R**2))*exp(-(rMax2/self.R)**beta2) - \
            (12*(beta2**3)*dp2/rho)*((((rMax2/self.R)**b2)**2)/(self.R**2))*exp(-(rMax2/self.R)**beta2) - \
            (4*(beta2**2)*dp2/rho)*((((rMax2/self.R)**beta2)**2)/(self.R**2))*exp(-(rMax2/self.R)**beta2) + \
            (4*(beta2**3)*dp2/rho)*((((rMax2/self.R)**beta2)**3)/(self.R**2))*exp(-(rMax2/self.R)**beta2)+2*self.f**2)

    # Second derivative at r = rMax1:
    d2Vm = -1/(8*(4*beta1*dp1/(self.rho*E) + \
            (4*beta2*dp2/self.rho)*(nu)*exp(-nu) + (rMax1*self.f)**2)**(3/2)) * \
            (-(4*(beta1**2)*dp1/(self.rho*rMax1*E)) + \
            (4*(beta1**2)*dp1/(self.rho*rMax1*E)) - \
            (4*(beta2**2)*dp2/self.rho)*(nu/rMax1)*exp(-nu) + \
            (4*(beta2**2)*dp2/self.rho)*((nu**2)/rMax1)*exp(-nu) + 2*rMax1*self.f**2)**2 + \
            1/(4*sqrt((4*beta1*dp1/(self.rho*E)) + \
            (4*beta2*dp2/self.rho)*nu*2+exp(-nu) + (rMax1*self.f)**2)) * \
            ((4*(beta1**3)*dp1/(self.rho*(rMax1**2)*E)) + \
            (4*(beta1**2)*dp1/(self.rho*(rMax1**2)*E)) - \
            (12*(beta1**3)*dp1/(self.rho*(rMax1**2)*E)) - \
            (4*(beta1**2)*dp1/(self.rho*(rMax1**2)*E)) + \
            (4*(beta1**3)*dp1/(self.rho*(rMax1**2)*E)) + \
            (4*(beta2**3)*dp2/self.rho)*(nu/(rMax1**2))*exp(-nu) + \
            (4*(beta2**2)*dp2/self.rho)*(nu/(rMax1**2))*exp(-nu) - \
            (12*(beta2**3)*dp2/self.rho)*((nu**2)/(rMax1**2))*exp(-nu) - \
            (4*(beta2**2)*dp2/self.rho)*((nu**2)/(rMax1**2))*exp(-nu) + \
            (4*(beta2**3)*dp2/self.rho)*((nu**3)/(rMax1**2))*exp(-nu)+2*self.f**2)
    """

    #dVm = -abs(f)/2 + (-4*(beta1**2)*dp1/(rho*rMax1*E) + \
    #        4*(beta1**2)*dp1/(rho*rMax1*E) - \
    #        (4*(beta2**2)*dp2/rho)*(nu/rMax1)*exp(-nu) + \
    #        (4*(beta2**2)*dp2/rho)*((nu**2)/rMax1)*exp(-nu)+2*rMax1*f**2)/ \
    #        (4*sqrt((4*beta1*dp1/(rho*E)) + \
    #        (4*beta2*dp2/rho)*nu*exp(-nu)+(rMax1*f)**2))
    #
    #d2Vm = -1/(8*(4*beta1*dp1/(rho*E) + \
    #        (4*beta2*dp2/rho)*(nu)*exp(-nu) + (rMax1*f)**2)**(1.5)) * \
    #        (-(4*(beta1**2)*dp1/(rho*rMax1*E)) + \
    #        (4*(beta1**2)*dp1/(rho*rMax1*E)) - \
    #        (4*(beta2**2)*dp2/rho)*(nu/rMax1)*exp(-nu) + \
    #        (4*(beta2**2)*dp2/rho)*((nu**2)/rMax1)*exp(-nu) + 2*rMax1*f**2)**2 + \
    #        1/(4*sqrt((4*beta1*dp1/(rho*E)) + \
    #        (4*beta2*dp2/rho)*nu*2+exp(-nu) + (rMax1*f)**2)) * \
    #        ((4*(beta1**3)*dp1/(rho*(rMax1**2)*E)) + \
    #        (4*(beta1**2)*dp1/(rho*(rMax1**2)*E)) - \
    #        (12*(beta1**3)*dp1/(rho*(rMax1**2)*E)) - \
    #        (4*(beta1**2)*dp1/(rho*(rMax1**2)*E)) + \
    #        (4*(beta1**3)*dp1/(rho*(rMax1**2)*E)) + \
    #        (4*(beta2**3)*dp2/rho)*(nu/(rMax1**2))*exp(-nu) + \
    #        (4*(beta2**2)*dp2/rho)*(nu/(rMax1**2))*exp(-nu) - \
    #        (12*(beta2**3)*dp2/rho)*((nu**2)/(rMax1**2))*exp(-nu) - \
    #        (4*(beta2**2)*dp2/rho)*((nu**2)/(rMax1**2))*exp(-nu) + \
    #        (4*(beta2**3)*dp2/rho)*((nu**3)/(rMax1**2))*exp(-nu)+2*f**2)

    #return dVm, d2Vm

#    Note: C Code disabled as offering only negligible overall speed improvement
#    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#    d = empty(1, 'd');
#
#    code = """
#
#    double f_ = f;
#    double rMax1_ = rMax1;
#    double beta1_ = beta1;
#    double beta2_ = beta2;
#    double dp1_ = dp1;
#    double dp2_ = dp2;
#    double rho_ = rho;
#    double E = exp(1);
#    double nu = pow((rMax2/rMax1_),beta2);
#    double enu = exp(-1.0*nu);
#    double b12 = pow(beta1_,2.0);
#    double b22 = pow(beta2_,2.0);
#    double b13 = pow(beta1_,3.0);
#    double b23 = pow(beta2_,3.0);

#    d = -1.0/(8.0*pow(4.0*beta1_*dp1_/(rho_*E)
#           +(4.0*beta2_*dp2_/rho_)*nu*enu + pow((rMax1_*f_),2.0),1.5))
#           *pow(-1.0*(4.0*b12*dp1_/(rho_*rMax1_*E)) + (4.0*b12*dp1_/(rho_*rMax1_*E))
#           - (4.0*b22*dp2_/rho_)*(nu/rMax1_)*enu
#           + (4.0*b22*dp2_/rho_)*(pow(nu,2.0)/rMax1_)*enu + 2.0*rMax1_*pow(f_,2.0),2.0)
#           + 1.0/(4.0*sqrt((4.0*beta1_*dp1_/(rho_*E)) + (4.0*beta2_*dp2_/rho_)*nu*2.0+enu + pow(rMax1_*f_,2.0))) *
#           ((4.0*b13*dp1_/(rho_*pow(rMax1_,2.0)*E))
#           + (4.0*b12*dp1_/(rho_*pow(rMax1_,2.0)*E))
#           - (12.0*b13*dp1_/(rho_*pow(rMax1_,2.0)*E))
#           - (4.0*b12*dp1_/(rho_*pow(rMax1_,2.0)*E))
#           + (4.0*b13*dp1_/(rho_*pow(rMax1_,2.0)*E))
#           + (4.0*b23*dp2_/rho_)*(nu/pow(rMax1_,2.0))*enu
#           + (4.0*b22*dp2_/rho_)*(nu/pow(rMax1_,2.0))*enu
#           - (12.0*b23*dp2_/rho_)*(pow(nu,2.0)/pow(rMax1_,2.0))*enu
#           - (4.0*b22*dp2_/rho_)*(pow(nu,2.0)/pow(rMax1_,2.0))*enu
#           + (4.0*b23*dp2_/rho_)*(pow(nu,3.0)/pow(rMax1_,2.0))*enu+2.0*pow(f_,2.0));
#
#    """
#
#    err = weave.inline(code,
#                       ['f', 'rMax1', 'rMax2', 'beta1', 'beta2', 'dp1', 'dp2', 'rho', 'd'],
#                       type_converters=converters.blitz,
#                       compiler='gcc')
    
    E = exp(1)
    nu = pow((rMax2/rMax1),beta2);
    
    d2Vm = -1/(8*(4*beta1*dp1/(rho*E) + \
            (4*beta2*dp2/rho)*(nu)*exp(-nu) + (rMax1*f)**2)**(1.5)) * \
            (-(4*(beta1**2)*dp1/(rho*rMax1*E)) + \
            (4*(beta1**2)*dp1/(rho*rMax1*E)) - \
            (4*(beta2**2)*dp2/rho)*(nu/rMax1)*exp(-nu) + \
            (4*(beta2**2)*dp2/rho)*((nu**2)/rMax1)*exp(-nu) + 2*rMax1*f**2)**2 + \
            1/(4*sqrt((4*beta1*dp1/(rho*E)) + \
            (4*beta2*dp2/rho)*nu*2+exp(-nu) + (rMax1*f)**2)) * \
            ((4*(beta1**3)*dp1/(rho*(rMax1**2)*E)) + \
            (4*(beta1**2)*dp1/(rho*(rMax1**2)*E)) - \
            (12*(beta1**3)*dp1/(rho*(rMax1**2)*E)) - \
            (4*(beta1**2)*dp1/(rho*(rMax1**2)*E)) + \
            (4*(beta1**3)*dp1/(rho*(rMax1**2)*E)) + \
            (4*(beta2**3)*dp2/rho)*(nu/(rMax1**2))*exp(-nu) + \
            (4*(beta2**2)*dp2/rho)*(nu/(rMax1**2))*exp(-nu) - \
            (12*(beta2**3)*dp2/rho)*((nu**2)/(rMax1**2))*exp(-nu) - \
            (4*(beta2**2)*dp2/rho)*((nu**2)/(rMax1**2))*exp(-nu) + \
            (4*(beta2**3)*dp2/rho)*((nu**3)/(rMax1**2))*exp(-nu)+2*f**2)

    if d2Vm >= 0.0:
        raise ValueError
    else:
        return d2Vm

