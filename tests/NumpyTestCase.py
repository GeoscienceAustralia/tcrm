"""
 Title: NumpyTestCase.py
 Author: Nariman Habili, nariman.habili@ga.gov.au
 CreationDate: 2006-12-12
 Description: Equality tests for numpy arrays.

 Version: $Rev: 563 $

 $Id: NumpyTestCase.py 563 2007-10-24 02:52:40Z carthur $
"""

from scipy import array
from numpy import allclose, iscomplexobj, alltrue, equal
import unittest

class NumpyTestCase(unittest.TestCase):
    """ExtendsTestCase with equality tests for numpy arrays."""

    def numpyAssertEqual(self, a1, a2):
        """Test for equality of array fields a1 and a2."""

        self.assertEqual(type(a1), type(a2))
        self.assertEqual(a1.shape, a2.shape)
        self.assertEqual(a1.dtype, a2.dtype)
        self.assertTrue(alltrue(equal(a1.ravel(), a2.ravel())))

    def numpyAssertAlmostEqual(self, a1, a2, prec=1.0000000000000001e-005):
        """Test for approximately equality of array fields a1 and a2."""

        self.assertEqual(type(a1), type(a2))
        self.assertEqual(a1.shape, a2.shape)
        self.assertEqual(a1.dtype, a2.dtype)

        if iscomplexobj(a1):
            ar1, ar2 = a1.real.ravel(), a2.real.ravel()
            assert allclose(ar1, ar2, prec)

            ar1, ar2 = a1.imag.ravel(), a2.imag.ravel()
            assert allclose(ar1, ar2, prec)
        else:
            assert allclose(a1, a2, prec)

    def numpyAssertEqualElements(self, a):
        """Test for equality of all elements of array a."""

        if iscomplexobj(a):
            self.assertEqual(a.real.min(), a.real.max())
            self.assertEqual(a.imag.min(), a.imag.max())
        else:
            self.assertEqual(a.min(), a.max())

    def numpyAssertAlmostEqualElements(self, a, prec=1.0000000000000001e-005):
        """Test for equality of all elements of array a."""

        if iscomplexobj(a):
            assert allclose(a.real.min(), a.real.max(), prec)
            assert allclose(a.imag.min(), a.imag.max(), prec)
        else:
            assert allclose(a.min(), a.max(), prec)
