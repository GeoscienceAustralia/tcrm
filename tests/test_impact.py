#def impact(gust, exposure, output='impact.csv')


import os
import sys
import unittest
import tempfile
import shutil
import numpy as np
import xarray
import pandas
from distutils.spawn import find_executable as which # py3: shutil.which

# Add parent folder to python path
unittest_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(unittest_dir, '..')))

import impact

@unittest.skipIf(which('hazimp') is None, "requires hazimp to be already installed")
class test(unittest.TestCase):
    def setUp(self):
        """Build dummy gust.nc and exposure.csv"""

        self.tmpdir = tempfile.mkdtemp()

        lat = np.arange(-10, 11, 1)
        lon = np.arange(100, 150, 1)

        field_equation = lambda lat, lon: lat + 2 * lon # arbitrary

        speed = field_equation(lat[None,:], lon[:, None])
        assert speed.size == lat.size * lon.size # broadcasts

        self.gust = os.path.join(self.tmpdir, 'gust.nc')
        ds = xarray.Dataset({'vmax': (['lon', 'lat'], speed)},
                            coords=dict(lat=lat, lon=lon))
        ds.vmax.attrs['units'] = '0.2s gust at 10m height m/s'
        ds.to_netcdf(self.gust)

        self.exposure = os.path.join(self.tmpdir, 'exposure.csv')
        table = [dict(LID='A', LATITUDE=-2, LONGITUDE=130),
                 dict(LID='B', LATITUDE=-8, LONGITUDE=130),
                 dict(LID='X', LATITUDE=-8, LONGITUDE=105)]
        df = pandas.DataFrame(table).set_index('LID')
        df['WIND_VULNERABILITY_FUNCTION_ID'] = 'dw253'
        df['REPLACEMENT_VALUE'] = 100000
        df.to_csv(self.exposure)

        self.structure_winds = field_equation(df.LATITUDE.values,
                                              df.LONGITUDE.values)
    def tearDown(self):
        """Remove temporary files"""
        shutil.rmtree(self.tmpdir)
    def test_hazimp(self):
        """Cause TCRM module to invoke HazImp, check ran ok"""
        outfile = os.path.join(self.tmpdir, 'test_impact.csv')

        impact.impact(self.gust, self.exposure, outfile)

        result = pandas.read_csv(outfile).set_index('LID')
        winds = result['0.2s gust at 10m height m/s']

        self.assertTrue(len(self.structure_winds) == 3)
        self.assertTrue(len(winds) == 3)
        #self.assertTrue((winds[['A','B','X']].values == self.structure_winds
        #                 ).all())

if __name__ == '__main__':
    unittest.main()
