"""

This module tests wind.writer, which is responsible for
saving the windfield dynamical evolution.

"""
# Currently, this module simply invokes the doctests in wind.writer.

# The writer is tightly coupled to the interface exposed elsewhere
# by the wind module and TCRM. Therefore, the test of the writer
# is brittle with regard to plausible refactors of TCRM.
# The use of doctest keeps the test close to the code, where it is
# most likely to be properly maintained.

# TODO: repair broken doctests elsewhere in code base, so that it
# becomes practical to invoke doctests universally (e.g. from nose)
# instead of by this module.

import os
import sys
import doctest
import unittest

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath( __file__ ))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))

import wind.writer

suite = doctest.DocTestSuite(wind.writer)

def load_tests(loader, std, pat): # invoked by unittest discovery process
    return suite

# nosetests discovery process
load_tests.__test__ = False
def test_with_nose():
    assert suite.run(unittest.TestResult()).wasSuccessful()

if __name__ == '__main__':
    unittest.TextTestRunner().run(suite)
