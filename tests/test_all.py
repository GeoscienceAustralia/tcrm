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


Regression testing framework
This module will search for scripts in the same directory named
test_*.py.  Each such script should be a test suite that tests a
module through PyUnit. This script will aggregate all
found test suites into one big test suite and run them all at once.
"""

# Author: Mark Pilgrim
# Modified by Ole Nielsen

import unittest
import os, sys
import pdb

# ---- Uncomment this section when compiling TCRM ----
#from unittests import pathLocate
#from unittests import test_derivatives, test_evd, test_generateStats, \
#                      test_grid, test_gridNC, test_KDEOrigin, test_KDEParameters, \
#                      test_maputils, test_metutils, test_mslp_seasonal_clim, \
#                      test_nctools, test_pressureProfile, test_SamplingOrigin, \
#                      test_SamplingParameters, test_stats, test_vmax, test_windField, \
#                      test_windProfile, test_windVorticity
# ----------------------------------------------------

try:
    import pathLocate
except:
    from unittests import pathLocate

# Add tcrm folder to python path
tcrm_dir = pathLocate.getRootDirectory()
sys.path.append(tcrm_dir)
from Utilities.files import flStartLog

# Switch off minor warning messages
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="pytz")

#List files that should be excluded from the testing process.
#E.g. if they are known to fail and under development

exclude = ['test_GenerateDistributions.py', 'test_lat_long_UTM_conversion.py', 'test_TrackGenerator.py', 'test_nctools.py']
if pathLocate.is_frozen():
    exclude = []

# Test suites still to be produced:
# config.py
# columns.py
# GetType.py
# template.py

def get_test_files(path):

    import sys

    files = os.listdir(path)
    if pathLocate.is_frozen():  #If compiled, then revert to pre-listed test suites
        files = ['test_derivatives', 'test_evd', 'test_generateStats',
                'test_grid', 'test_gridNC', 'test_KDEOrigin', 'test_KDEParameters', 'test_maputils',
                'test_metutils', 'test_mslp_seasonal_clim', 'test_pressureProfile',
                'test_SamplingOrigin', 'test_SamplingParameters', 'test_stats', 'test_vmax', 'test_windField',
                'test_windProfile', 'test_windVorticity']
        files = [k + '.py' for k in files]

    #Check sub directories
    test_files = []
    for file in files:
        if os.path.isdir(file):
            # Sub directory checking switched off
            pass
            ##sys.path.append(file)
            ###print 'Recursing into', file
            ##test_files += get_test_files(path + os.sep + file)
        elif file[:5] == 'test_' and file[-2:] == 'py' and file[-3:]!='pyc':
            #print 'Appending', file
            test_files.append(file)
        elif file[-3:] == 'pyc' and file[:5]=='test_':
            # This stops test_all using the byte-compiled versions of
            # the tests, which causes one of the tests in test_files.py to fail.
            os.remove(os.path.join(path,file))
        else:
            pass
    return test_files



def regressionTest():
    import sys, os, re, unittest
    path = os.path.split(sys.argv[0])[0] or os.getcwd()

    files = get_test_files(path)

    #test = re.compile('^test_[\w]*.py$', re.IGNORECASE)
    #files = filter(test.search, files)

    #try:
    #    files.remove(__file__)  #Remove self from list (Ver 2.3. or later)
    #except:
    #    files.remove('test_all.py')
    print
    print 'Testing:'
    for file in files:
        print '  ' + file
    if globals().has_key('exclude'):
        for file in exclude:
            files.remove(file)
            #print 'WARNING: File '+ file + ' excluded from testing'


    filenameToModuleName = lambda f: os.path.splitext(f)[0]
    #print "files",files
    moduleNames = map(filenameToModuleName, files)
    if pathLocate.is_frozen():
        modules = [sys.modules['unittests.' + k] for k in moduleNames]
    else:
        modules = map(__import__, moduleNames)  
    load = unittest.defaultTestLoader.loadTestsFromModule
    return unittest.TestSuite(map(load, modules))


if __name__ == '__main__':

    from os import sep

    #Attempt to compile all extensions
    #execfile('..' + sep + 'utilities' + sep + 'compile.py')

    #FIXME: Temporary measure
    #os.chdir('..' + sep + 'utilities')
    #execfile('compile.py')
    #os.chdir('..' + sep + 'pyvolution')

    #FIXME: Temporary measure
    #os.chdir('..' + sep + 'triangle')
    #execfile('compile.py')
    #os.chdir('..' + sep + 'pyvolution')

    #os.system('python compile.py')

    #print regressionTest()
    #unittest.main(defaultTest='regressionTest')
    flStartLog('', 'CRITICAL', False)
    suite = regressionTest()
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
