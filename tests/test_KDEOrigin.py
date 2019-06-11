"""
 Title: KDEOriginTestCase.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2006-12-13
 Description: Unit testing module for KDEOrigin.py

 Version: $Rev: 563 $

 $Id: TestKDEOrigin.py 563 2007-10-24 02:52:40Z carthur $
"""
import os, sys, pdb
import pickle
import unittest
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from StatInterface import KDEOrigin
from Utilities.files import flStartLog

class TestKDEOrigin(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        kdeType = 'Epanechnikov'
        gridLimit={'xMin':70, 'xMax':180, 'yMin':-36, 'yMax':0}
        kdeStep = 0.1
        lonLat = pickle.load(open(os.path.join(unittest_dir, 'test_data', 'kde_origin_lonLat.pkl'), 'rb'))
        self.kdeOrigin = KDEOrigin.KDEOrigin(None, gridLimit, kdeStep, lonLat)

    def tearDown(self):
        #forgetAllSingletons()
        pass

    def test_GenerateKDE(self):
        """Testing GenerateKDE for 2D data"""
        pkl_file = open(os.path.join(unittest_dir, 'test_data', 'kdeOrigin_xyz.pkl'), 'rb')
        xp = pickle.load(pkl_file)
        yp = pickle.load(pkl_file)
        zp = pickle.load(pkl_file)

        x, y, z = self.kdeOrigin.generateKDE()
        self.numpyAssertAlmostEqual(xp, x)
        self.numpyAssertAlmostEqual(yp, y)
        self.numpyAssertAlmostEqual(zp, z)

if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(TestKDEOrigin,'test')
    unittest.TextTestRunner(verbosity=2).run(testSuite)
