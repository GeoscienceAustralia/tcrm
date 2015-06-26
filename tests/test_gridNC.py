import unittest
import inspect
import os
import numpy as np
from numpy.testing import assert_almost_equal
from Utilities import nctools
from Utilities import grid

TEST_DIR = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda _: None)))

class TestSampleGrid(unittest.TestCase):

    """
    Test that the range of methods to load gridded data produce expected
    results. Uses the 0.083 degree land-sea mask dataset as a test dataset.
    """

    def setUp(self):
        self.filename = os.path.join(TEST_DIR, 'test_data', 'landmask.nc')
        # Load the data using grid.grdRead:
        self.lslon, self.lslat, self.lsgrid = grid.grdRead(self.filename)
        # Load the data using nctools.ncLoadFile and nctools.ncGetData:
        ncobj = nctools.ncLoadFile(self.filename)
        self.nclon = nctools.ncGetDims(ncobj, 'lon')
        self.nclat = nctools.ncGetDims(ncobj, 'lat')
        self.ncgrid = nctools.ncGetData(ncobj, 'landmask')
        ncobj.close()
        # Set up an instance of SampleGrid:
        self.sample = grid.SampleGrid(self.filename)
        # Sample the land-sea mask at these points around the globe:
        self.xlon = [100., 130., 180., 250., 300.]
        self.ylat = [-80., -20., 40.]
        # Known point values of the land-sea mask data:
        self.ls = [3., 0., 3., 3., 3., 0., 3., 0., 0., 3., 0., 3., 3., 3., 0.]

    def test_gridNctools(self):
        """Test grid.grdRead produces same result as nctools.ncGetData"""
        assert_almost_equal(self.lslon, self.nclon)
        assert_almost_equal(self.lslat, self.nclat)
        assert_almost_equal(self.lsgrid, self.ncgrid)

    def test_sampleGridArray(self):
        """Test grid.SampleGrid.grid produces flipped result of nctools.ncGetData"""
        assert_almost_equal(self.sample.lon, self.nclon)
        assert_almost_equal(self.sample.lat, self.nclat)
        assert_almost_equal(self.sample.grid, np.flipud(self.ncgrid))

    def test_sampleGridPoints(self):
        """Test grid.SampleGrid returns correct values for known points"""
        i = 0
        for x in self.xlon:
            for y in self.ylat:
                assert_almost_equal(self.sample.sampleGrid(x, y), self.ls[i])
                i += 1

if __name__ == "__main__":
    unittest.main()
