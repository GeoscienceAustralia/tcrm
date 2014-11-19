import unittest
import inspect
import os
import numpy as np
from numpy.testing import assert_almost_equal
from Utilities import files

TEST_DIR = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda _: None)))


class TestModuleUtilities(unittest.TestCase):

    def setUp(self):
        self.filename = os.path.abspath(inspect.getsourcefile(lambda _:None))
        self.path, self.fname = os.path.split(self.filename)
        self.base, self.ext = os.path.splitext(self.fname)
        #self.path = self.path.replace(os.path.sep, '/')

    def testModulePath(self):
        """Test flModulePath returns correct path, base & extension"""
        p, b, e = files.flModulePath()
        self.assertEqual(self.path, p)
        self.assertEqual(self.base, b)
        self.assertEqual(self.ext, e)

    def testModuleName(self):
        """Test flModuleName returns correct module name"""
        self.assertEqual('testModuleName', files.flModuleName())


class TestFileLoading(unittest.TestCase):

    def setUp(self):
        self.lat = np.arange(0, -21, -1, 'd')
        self.lon = np.arange(130, 151, 1, 'd')
        self.testfile = os.path.join(TEST_DIR, 'test_data/test_files.dat')
        self.testfile2 = os.path.join(TEST_DIR, 'test_data/test_files2.dat')

    def testLoadFile(self):
        """Test flLoadFile loads data correctly"""
        data = files.flLoadFile(self.testfile, comments='%', delimiter=',')
        assert_almost_equal(self.lon, data[:, 0])
        assert_almost_equal(self.lat, data[:, 1])

    def testLoadNonFile(self):
        """Trying to load a non-file should raise an IOError"""
        self.assertRaises((IOError,), files.flLoadFile, 'string')

    def testLoadBrokenFile(self):
        """Test flLoadFile on a file with differing length columns"""
        self.assertRaises(ValueError, files.flLoadFile,
                          self.testfile2, comments='%', delimiter=',')


#class TestFileStats(unittest.TestCase):
#
#    """
#    Test the flGetStat function
#
#    Note that moving the test_files.dat file will probably change the md5sum,
#    so this may not be the best way to test the function...watch this space (CA
#    - 20091127)
#    """
#
#    def setUp(self):
#        self.testfile = os.path.join(TEST_DIR, 'test_data/test_files.dat')
#        self.testfile3 = os.path.join(TEST_DIR, 'test_data/test_files3.dat')
#        self.testmd5 = 'b8668634af964d4f88b5b83d714a5771'
#        self.data = np.arange(0, -21, 1, 'd')
#
#    def testFileErrors(self):
#        """Test flGetStat raises correct exceptions"""
#        self.assertRaises(IOError, files.flGetStat, self.data)
#        self.assertRaises(IOError, files.flGetStat, self.testfile3)
#
#    def testGetStats(self):
#        """Test flGetStat function returns correct values"""
#        dir, fname, md5sum, moddate = files.flGetStat(self.testfile)
#        self.assertEqual(os.path.basename(self.testfile), fname)
#        self.assertEqual(self.testmd5, md5sum)


class TestConfigFile(unittest.TestCase):

    def testConfigFile(self):
        """Test flConfigFile returns correct filename"""
        testconfig = os.path.join(TEST_DIR, 'test_files.ini')
        #FIXME: This can fail if run using test_all.py
        #self.assertEqual(testconfig, files.flConfigFile())


class TestFileSize(unittest.TestCase):

    def setUp(self):
        self.filename = os.path.realpath(__file__)

    def testSize(self):
        """Test flSize returns correct file size value"""
        testsize = os.stat(self.filename).st_size
        self.assertEqual(testsize, files.flSize(self.filename))

    def testSizeError(self):
        """Test flSize raises IOError for non-file"""
        self.assertRaises(OSError, files.flSize, 'string')


class TestModDate(unittest.TestCase):

    def setUp(self):
        self.filename = os.path.realpath(__file__)

    def testModDate(self):
        """Test flModDate function"""
        from time import localtime, strftime
        testdate = localtime(os.stat(self.filename).st_mtime)
        testdatestr = strftime('%Y-%m-%dT%H:%M:%S', testdate)
        self.assertEqual(testdatestr,
                         files.flModDate(self.filename,
                                         '%Y-%m-%dT%H:%M:%S'))


if __name__ == "__main__":
    unittest.main()
