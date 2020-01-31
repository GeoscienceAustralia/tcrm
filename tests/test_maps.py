import unittest

import numpy as np
from . import NumpyTestCase
from matplotlib.colors import LinearSegmentedColormap
from PlotInterface import maps

class TestLevels(NumpyTestCase.NumpyTestCase):

    def test_badInput(self):
        """Non-numeric input raises TypeError"""
        self.assertRaises(TypeError, maps.levels, "str")
        self.assertRaises(TypeError, maps.levels, 10, "str")
        self.assertRaises(TypeError, maps.levels, np.max)

    def test_EmptyRange(self):
        """Minimum & maximum value equal raises ValueError"""
        self.assertRaises(ValueError, maps.levels, 0, 0)
        self.assertRaises(ValueError, maps.levels, 1, 1)

    def test_MinGTMax(self):
        """Minimum greater than maximum automatically swaps values"""
        lvls = np.array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        expo = 0

        rlevs, rexpo = maps.levels(0, minval=1)
        self.numpyAssertAlmostEqual(lvls, rlevs)
        self.assertEqual(expo, rexpo)

    def test_smallLevelValues(self):
        """Test level determination for small input values"""
        lvs = np.array([  0.0000e+00,   1.0000e-06,   2.0000e-06,
                          3.0000e-06,   4.0000e-06,   5.0000e-06,
                          6.0000e-06,   7.0000e-06,   8.0000e-06,
                          9.0000e-06])
        expo = -5
        rlevs, rexpo = maps.levels(0, 10**-5)
        self.numpyAssertAlmostEqual(lvs, rlevs)
        self.assertEqual(expo, rexpo)

    def test_bigLevelValues(self):
        """Test level determination for big input values"""
        lvs = np.array([ 1000000.,  2000000.,  3000000.,  4000000.,
                         5000000.,  6000000.,  7000000.,  8000000.,
                         9000000.])
        expo = 7
        rlevs, rexpo = maps.levels(10**7, 10**6)
        self.numpyAssertAlmostEqual(lvs, rlevs)
        self.assertEqual(expo, rexpo)

class TestSelectColorMap(unittest.TestCase):

    def assertColorMapEqual(self, actual, expected):
        """Test method for equality of LinearSegmentedColormaps"""
        self.assertEqual(actual.N, expected.N)
        self.assertDictEqual(actual._segmentdata,
                             expected._segmentdata)
        for k in list(actual._segmentdata.keys()):
            self.assertListEqual(actual._segmentdata[k],
                                 expected._segmentdata[k])

    def setUp(self):
        import seaborn as sns
        palette = [(1, 1, 1), (0.000, 0.627, 0.235), (0.412, 0.627, 0.235), (0.663, 0.780, 0.282),
        (0.957, 0.812, 0.000), (0.925, 0.643, 0.016), (0.835, 0.314, 0.118),
        (0.780, 0.086, 0.118)]
        div_pal = sns.color_palette("RdBu", 7)
        self.diverging_cmap = sns.blend_palette(div_pal, as_cmap=True)
        self.sequential_cmap = sns.blend_palette(palette, as_cmap=True)

    def test_returnsCmap(self):
        """Test select_colormap returns a cmap object"""
        cmap = maps.selectColormap([-1, 1])
        self.assertIsInstance(cmap, LinearSegmentedColormap)

    def test_badInputType(self):
        """Test select_colormap raises TypeError for wrong input type"""
        self.assertRaises(TypeError, maps.selectColormap, 10)
        self.assertRaises(TypeError, maps.selectColormap, "string")
        self.assertRaises(TypeError, maps.selectColormap, ["str", "str"])

    def test_symmetricDataDiverges(self):
        """Test symmetric data range (around 0) returns divergent color map"""
        cmap = maps.selectColormap([-1, 1])
        self.assertColorMapEqual(cmap, self.diverging_cmap)

    def test_sequentialData(self):
        """Test asymmetric data range returns sequential color map"""
        cmap = maps.selectColormap([0, 1])
        self.assertColorMapEqual(cmap, self.sequential_cmap)

    def test_largePercentThreshold(self):
        """Test select_colormap returns diverging color map for larger threshold"""
        cmap = maps.selectColormap([-7, 10], percent=0.2)
        self.assertColorMapEqual(cmap, self.diverging_cmap)

if __name__ == '__main__':
    unittest.main(verbosity=2)
