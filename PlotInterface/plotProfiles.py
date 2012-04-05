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


 Title: plotProfiles.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2007-02-13
 Description: Plot the full suite of profiles for a given set of
              conditions for intercomparison.
 Reference:
 SeeAlso:
 Constraints:

 Version: $Rev: 634 $
 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: plotProfiles.py 634 2007-12-16 21:16:58Z carthur $
"""
import os, sys, pdb, logging

import WindfieldInterface.windProfile as windProfile
import Utilities.maputils as maputils
import numpy
from matplotlib import pyplot

class plotProfiles:
    """plotProfiles:
    Description: Plot the radial profile for a given set of cyclone
    parameters.  It generates the profile for radii up to 300 km.
    Results are plotted and the resulting figure can be saved by
    including a filename when calling the plot method.

    Parameters:
    cLat : float
        Cyclone latitude
    cLon : float
        Cyclone longitude
    rMax : float
        Radius of maximum winds
    pEnv : float
        Environmental pressure surrounding cyclone
    pCentre: float
        Central pressure of cyclone
    beta : float
        Holland beta (peakedness) parameter.

    Members:

    Methods:
    plot(profileType, file) :
        Plot the given profile type, if none given then plot all

    Internal Methods:

    """

    def __init__(self, cLat, cLon, rMax, pEnv, pCentre, beta):
        """
        Initialise required fields
        """
        self.R = numpy.array(range(1, 201), 'f')
        self.cLat = cLat
        self.cLon = cLon
        self.rM = rMax
        self.pE = pEnv
        self.pC = pCentre
        self.beta = beta
        self.profile = windProfile.WindProfile(self.R, self.pE, self.pC,
                                               self.rM, self.cLat, self.cLon,
                                               self.beta)

    def __doc__(self):
        """
        Documentation on the function of the class:
        """
        return ' '

    def plot(self, profileType=None, filename=None):
        """plot(filename=None)
        Plot all available radial profiles. Will save the plot to file
        if 'filename' is given.
        """
        pyplot.figure(0)
        pyplot.axes([0.125, 0.2, 0.95-0.125, 0.95-0.3])
        pyplot.hold(True)
        legend = []

        methodList = [method for method in dir(self.profile)
                      if callable(getattr(self.profile, method))]
        if profileType is None:
            for method in methodList:
                if method == '__doc__' or method == '__init__':
                    pass
                else:
                    legend.append(method.capitalize())
                    wind = getattr(self.profile, method)
                    V = wind()
                    pyplot.plot(self.R, abs(V), linewidth=2)
        else:
            legend.append(profileType.capitalize())
            wind = getattr(self.profile, profileType)
            V = wind()
            pylab.plot(self.R,abs(V), linewidth=2)

        pyplot.legend(legend)
        pyplot.grid()
        pyplot.xlabel('Radius (km)', fontsize=14)
        pyplot.ylabel('Wind speed (m/s)', fontsize=14)
        pyplot.title(r'$P_c = %d\hspace{0.5}hPa,\hspace{1} P_e = %d \hspace{0.5} hPa,\hspace{1} R_{max} = %d \hspace{0.5}km$' %
                    (self.pC/100., self.pE/100., self.rM))
        pyplot.savefig(filename)

        if profileType is not None:
            return R, V


if __name__ == "__main__":
    cLat = -12.
    cLon = 130.

    rMax = 30.
    pEnv = 100700.
    pCentre = 95000.
    beta = 1.6
    filename = 'C:/temp/all_profiles.png'
    p = plotProfiles(cLat, cLon, rMax, pEnv, pCentre, beta)
    p.plot(profileType=None, filename=filename)
