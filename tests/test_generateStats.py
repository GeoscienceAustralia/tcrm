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


Title: testWindProfile.py

Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2011-06-10
Description: Unit test for the windProfile class

"""
import os, sys
import unittest
import pickle
from . import NumpyTestCase
import numpy
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from StatInterface import generateStats
from Utilities.files import flStartLog

class TestGenerateStats(NumpyTestCase.NumpyTestCase):

    pkl_file = open(os.path.join(unittest_dir, 'test_data',
                                 'generateStatsTestData.pkl'), 'rb')
    lonLat = pickle.load(pkl_file)
    parameter = pickle.load(pkl_file)  # all_speed
    pkl_file.close()

    configFile = None
    gridLimit = {'xMin': 150.0, 'xMax': 186.0, 'yMin': -31.0, 'yMax': -5.0}
    gridSpace = {'x': 1.0, 'y': 1.0}
    gridInc = {'x': 1.0, 'y': 0.5}
    minSample = 100
    angular = False
    missingValue=min(sys.maxsize, 2147483648) # The second value corresponds to the maximum integer
                                             # encoded in the pickle file.


    def test_generateStats(self):
        """Testing generateStats (warning: will take about 2 mins to run)
        """
        wP = generateStats.GenerateStats(self.parameter, self.lonLat, self.gridLimit,
                                         self.gridSpace, self.gridInc, self.minSample, self.angular,
                                         self.missingValue)

        coeffs_mu_sample = wP.coeffs.mu[0:10]
        coeffs_mu_test = numpy.array([10.61791304, 10.70961832, 9.6531068, 9.99247706, 9.99247706,
                                      10.09581818, 10.42785714, 10.36168067, 10.6747541, 10.55462121])
        self.numpyAssertAlmostEqual(coeffs_mu_sample, coeffs_mu_test)

        coeffs_sig_sample = wP.coeffs.sig[0:10]
        coeffs_sig_test = numpy.array([7.25133853, 7.76308414, 6.53878674, 7.37497593, 7.37497593,
                                       7.42023356, 7.84711405, 7.76943206, 8.02239082, 7.8492764])
        self.numpyAssertAlmostEqual(coeffs_sig_sample, coeffs_sig_test)

        coeffs_alpha_sample = wP.coeffs.alpha[0:10]
        coeffs_alpha_test = numpy.array([0.51818004, 0.55450803, 0.66634809, 0.61266186, 0.61266186,
                                         0.63192755, 0.70984709, 0.64836016, 0.69168147, 0.70101634])
        self.numpyAssertAlmostEqual(coeffs_alpha_sample, coeffs_alpha_test)

        coeffs_phi_sample = wP.coeffs.phi[0:10]
        coeffs_phi_test = numpy.array([0.85527156, 0.83217838, 0.74564081, 0.79034514, 0.79034514,
                                       0.77502747, 0.70435581, 0.76133376, 0.7222027 , 0.71314521])
        self.numpyAssertAlmostEqual(coeffs_phi_sample, coeffs_phi_test)

        coeffs_min_sample = wP.coeffs.min[0:20]
        coeffs_min_test = numpy.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                       0., 1.82, 1.85, 1.85, 1.85])
        self.numpyAssertAlmostEqual(coeffs_min_sample, coeffs_min_test)

        # The land & sea statistics should be identical since reverts to statistics for open ocean
        # when insufficient grid points over land
        self.numpyAssertAlmostEqual(wP.coeffs.mu, wP.coeffs.lmu)
        self.numpyAssertAlmostEqual(wP.coeffs.sig, wP.coeffs.lsig)
        self.numpyAssertAlmostEqual(wP.coeffs.alpha, wP.coeffs.lalpha)
        self.numpyAssertAlmostEqual(wP.coeffs.phi, wP.coeffs.lphi)
        self.numpyAssertAlmostEqual(wP.coeffs.min, wP.coeffs.lmin)

if __name__ == "__main__":
    unittest.main()
