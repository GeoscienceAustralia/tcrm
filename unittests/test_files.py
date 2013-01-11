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


 Title: test_files.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2009-11-27 9:46:AM
 Description: Unit test for files.py

 Throughout, I convert all directory separators to *NIX style
 separators (i.e. the Windows style '\\' become '/')

 PLEASE READ THE UNITTEST DOCUMENTATION AT
 http://docs.python.org/library/unittest.html
 before adding more tests

 Version :$Rev: 310 $

 $Id: test_files.py 310 2010-07-06 00:31:58Z carthur $
"""
import os, sys, pdb, logging, unittest
import numpy
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from Utilities import files
from Utilities.files import flStartLog


class Test_flModulePath(NumpyTestCase.NumpyTestCase):
    # Set up the test:
    filename = os.path.realpath( __file__ )
    path, fname = os.path.split(filename)
    base, ext = os.path.splitext(fname)
    path = path.replace(os.path.sep,'/')

    def test_flModulePath(self):
        """Test flModulePath returns correct path, base & extension"""
        p,b,e = files.flModulePath()
        self.assertEqual(self.path,p)
        self.assertEqual(self.base,b)
        self.assertEqual(self.ext,e)

class Test_flModuleName(NumpyTestCase.NumpyTestCase):

    def test_flModuleName(self):
        """Test flModuleName returns correct module name"""
        self.assertEqual('test_flModuleName', files.flModuleName())

class Test_flLoadFile(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        self.lat = numpy.arange(0, -21, -1, 'd')
        self.lon = numpy.arange(130, 151, 1, 'd')
        self.testfile = os.path.join(unittest_dir,'test_data/test_files.dat')
        self.testfile2 = os.path.join(unittest_dir,'test_data/test_files2.dat')

    def test_flLoadFile(self):
        """Test flLoadFile loads data correctly"""
        data = files.flLoadFile(self.testfile,comments='%',delimiter=',')
        self.numpyAssertEqual(self.lon, data[:,0])
        self.numpyAssertEqual(self.lat, data[:,1])

    def test_loadNonFile(self):
        """Trying to load a non-file should raise an IOError"""
        self.assertRaises((IOError,), files.flLoadFile, 'string')

    def test_loadBrokenFile(self):
        """Test flLoadFile on a file with differing length columns"""
        self.assertRaises(ValueError, files.flLoadFile,self.testfile2,comments='%',delimiter=',')


class Test_flGetStat(NumpyTestCase.NumpyTestCase):
    """
    Test the flGetStat function
    Note that moving the test_files.dat file will probably change the md5sum, so this may not
    be the best way to test the function...watch this space (CA - 20091127)
    """

    def setUp(self):
        self.testfile = os.path.join(unittest_dir,'test_data/test_files.dat')
        self.testfile3 = os.path.join(unittest_dir,'test_data/test_files3.dat')
        self.testmd5 = 'b8668634af964d4f88b5b83d714a5771'
        self.data = numpy.arange(0,-21,1,'d')

    def test_fileErrors(self):
        """Test flGetStat raises correct exceptions"""
        self.assertRaises(IOError, files.flGetStat, self.data)
        self.assertRaises(IOError, files.flGetStat, self.testfile3)

    def test_flGetStat(self):
        """Test flGetStat function returns correct values"""
        dir,fname,md5sum,moddate = files.flGetStat(self.testfile)
        self.assertEqual(os.path.basename(self.testfile),fname)
        self.assertEqual(self.testmd5,md5sum)



class Test_flConfigFile(NumpyTestCase.NumpyTestCase):
    def test_flConfigFile(self):
        """Test flConfigFile returns correct filename"""
        testconfig = os.path.join(unittest_dir,'test_files.ini').replace(os.path.sep,'/')
        self.assertEqual(testconfig,files.flConfigFile())

class Test_flSize(NumpyTestCase.NumpyTestCase):
    # Set up the test:
    filename = os.path.realpath( __file__ )
    def test_flSize(self):
        """Test flSize returns correct file size value"""
        testsize = os.stat( self.filename ).st_size
        self.assertEqual(testsize,files.flSize(self.filename))

    def test_flSizeError(self):
        """Test flSize raises IOError for non-file"""
        self.assertRaises(IOError, files.flSize, 'string')

class Test_flModDate(NumpyTestCase.NumpyTestCase):
    # Set up the test:
    filename = os.path.realpath( __file__ )
    def test_flModDate(self):
        """Test flModDate function"""
        testdate = time.localtime(os.stat(self.filename).st_mtime)
        testdatestr = time.strftime( '%Y-%m-%dT%H:%M:%S', testdate )
        self.assertEqual( testdatestr, 
                          files.flModDate( self.filename,
                                           '%Y-%m-%dT%H:%M:%S'))




########################################################################
if __name__ == "__main__":
    flStartLog('', 'CRITICAL', False)
    testSuite = unittest.makeSuite(Test_flModulePath,'test')
    testSuite.addTest(Test_flModuleName('test_flModuleName'))
    testSuite.addTest(Test_flLoadFile('test_flLoadFile'))
    testSuite.addTest(Test_flLoadFile('test_loadNonFile'))
    testSuite.addTest(Test_flLoadFile('test_loadBrokenFile'))
    testSuite.addTest(Test_flGetStat('test_fileErrors'))
    testSuite.addTest(Test_flGetStat('test_flGetStat'))
    testSuite.addTest(Test_flConfigFile('test_flConfigFile'))
    testSuite.addTest(Test_flSize('test_flSize'))
    testSuite.addTest(Test_flSize('test_flSizeError'))
    testSuite.addTest(Test_flModDate('test_flModDate'))
    unittest.TextTestRunner(verbosity=2).run(testSuite)
