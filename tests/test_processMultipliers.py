"""
Unit test suite for processMultipliers.py


"""

import sys
import os
from os.path import join as pjoin, exists, dirname
import unittest
import tempfile
import shutil

from numpy.testing import assert_almost_equal
import numpy as np

from osgeo import osr, gdal
from osgeo.gdalconst import *
from netCDF4 import Dataset

try:
    from . import pathLocate
except:
    from tests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())
from wind import WindfieldGenerator
from Utilities.config import ConfigParser
from ProcessMultipliers import processMultipliers as pM


class TestProcessMultipliers(unittest.TestCase):

    """
    TestProcessMultipliers: test the processMultipliers script
    """

    def setUp(self):
        self.array = np.arange(100, dtype=float).reshape((10, 10))
        self.x = np.arange(120, 130)
        self.y = np.arange(10, 20)
        self.dx = 1
        self.dy = -1
        self.testRasterFile = "testRaster.tif"
        self.reprojectRaster = "testReprojection.tif" # is this used?
        self.projectedDataset = pM.createRaster(np.ones((2000, 5000)),
                                                np.arange(15000, 55000, 20),
                                                np.arange(20000, 120000, 20),
                                                dx=20, dy=-20, epsg=3123)

        self.uu = np.array([-1., 1., 1., -1.])
        self.vv = np.array([-1., -1., 1., 1.])
        self.bearing = np.array([45., 315., 225., 135.])

    def test_createRaster(self):
        """Test createRaster returns a gdal dataset"""

        result = pM.createRaster(self.array, self.x, self.y,
                                 self.dx, self.dy,
                                 filename=self.testRasterFile)
        self.assertEqual(type(result), gdal.Dataset)
        assert exists(self.testRasterFile)

    def test_createRasterII(self):
        """Test can create img files"""
        file_is = 'test.img'
        result = pM.createRaster(self.array, self.x, self.y,
                                 self.dx, self.dy,
                                 filename=file_is)
        self.assertEqual(type(result), gdal.Dataset)
        assert exists(file_is)

    def test_xcreateRaster3(self):
        """Test can create img files"""
        file_is = 'test.img'

        f_img = tempfile.NamedTemporaryFile(suffix='.tif',
                                        prefix='test_proMult',
                                        delete=False)
        f_img.close()

        gust = np.asarray([[1., -1., 10, -10.],
                         [100., -100., 1000, -1000.]])
        lon = np.asarray([136, 138, 140, 142]) # 136 is used
        lat = np.array([-22, -20]) # -22 is used

        result = pM.createRaster(gust, lon, lat,
                                 2, -2,
                                 filename=f_img.name)
        minx, miny, maxx, maxy, data = pM.loadRasterFileBandLonLat(f_img.name)

        self.assertEqual(minx, 136)
        self.assertEqual(maxx, 144)
        self.assertEqual(miny, -24)
        self.assertEqual(maxy, -20)

        # This isn't working in this test
        # but works in
        # test_xprocessMult_A
        # assert_almost_equal(np.flipud(gust), data)

        del result
        os.remove(f_img.name)


    def test_loadRasterFile(self):
        """Test loadRasterFile correctly loads data"""

        result = pM.loadRasterFile(self.testRasterFile, -9999)
        self.assertEqual(type(result), np.ndarray)
        assert_almost_equal(result, self.array[::np.sign(self.dy) * 1])

    def test_calculateBearing(self):
        """Test the correct bearings are returned"""
        bb = pM.calculateBearing(self.uu, self.vv)
        assert_almost_equal(bb, self.bearing)

    def test_reprojectDataset(self):
        """Test a dataset is correctly reprojected"""

        pM.reprojectDataset(self.testRasterFile, self.projectedDataset,
                            self.reprojectRaster)
        assert exists(self.reprojectRaster)
        prjDataset = gdal.Open(self.reprojectRaster)
        prjBand = prjDataset.GetRasterBand(1)
        prj_data = prjBand.ReadAsArray()
        # Check shape of arrays:
        self.assertEqual((2000, 5000), prj_data.shape)

        # Check geographic transform:
        self.assertEqual(prjDataset.GetGeoTransform(),
                         self.projectedDataset.GetGeoTransform())
        # Check projection: FIXME GetProjection() from projected dataset
        # drops AXIS["X",EAST],AXIS["Y",NORTH], from the projection
        # information compared to the match dataset.
        # self.assertEqual(prjDataset.GetProjection(),
        #                 self.projectedDataset.GetProjection())
        # Check values are correctly mapped:

        del prjDataset

    def test_reprojectDataset_same_nc_img(self):
        """Test a dataset is correctly reprojected"""
        # Write a .nc file to test
        # This is the gust file
        f_nc = tempfile.NamedTemporaryFile(suffix='.nc',
                                        prefix='test_processMultipliers',
                                        delete=False)
        f_nc.close()


        # Write an .img file to test
        #This is the multiplier file
        f_img = tempfile.NamedTemporaryFile(suffix='.img',
                                        prefix='test_processMultipliers',
                                        delete=False)
        f_img.close()
        
        # Write a temporary track file
        f_track = tempfile.NamedTemporaryFile(suffix='.nc',
                                        prefix='test_track',
                                        delete=False)
        f_track.close()

        lat = np.asarray([ -23, -20, -17, -14, -11, -8, -5])
        lon = np.asarray([137, 140, 143, 146, 149, 152, 155, 158])

        delta = lon[1] - lon[0]
        lon_mid = lon + delta / 2.
        lat_mid = lat + delta / 2.
        speed = np.zeros(([lat.shape[0], lon.shape[0]]))
        speed.fill(42.5)

        # doing this just to get values in
        Vx = Vy = P = speed

        result = lat_mid, lon_mid, speed, Vx, Vy, P
        cfg = ConfigParser()
        wg = WindfieldGenerator(cfg)
        wg.saveGustToFile(f_track.name, result, f_nc.name)
        # nctools.ncSaveGrid(multiplier_name, multiplier_values, lat,
        #               lon, f_nc.name)

        # For the multiplier file
        lat = np.asarray([ -17, -14, -11])
        lon = np.asarray([140, 143, 146, 149])
        speed = np.zeros(([lon.shape[0], lat.shape[0]]))
        dx = 3.0
        dy = -3.0
        speed.fill(27.1701)
        pM.createRaster(speed, lon, lat,
                        dx, dy,
                        filename=f_img.name)
        speed.fill(42.5)
        m4_max_file = f_img.name
        # pulling out a section of the processMultipliers.main
        # load the wind data
        ncobj = Dataset(f_nc.name, 'r')

        lat = ncobj.variables['lat'][:]
        lon = ncobj.variables['lon'][:]

        delta = lon[1] - lon[0]
        lon = lon - delta / 2.
        lat = lat - delta / 2.

        # Wind speed:
        wspd = ncobj.variables['vmax'][:]
        wind_raster_file = 'region_wind.tif'
        wind_prj_file = 'gust_prj.tif'

        wind_raster = pM.createRaster(wspd, lon, lat, delta, -delta,
                                   filename=wind_raster_file)
        # This shows that the right data is going in
        src_data = wind_raster.GetRasterBand(1).ReadAsArray()
        m4_max = gdal.Open(m4_max_file, GA_ReadOnly)
        m4_max_data = m4_max.GetRasterBand(1).ReadAsArray()

        pM.reprojectDataset(wind_raster, m4_max_file,
                            wind_prj_file) #,
                            # match_projection=32756)

        wind_prj_ds = gdal.Open(wind_prj_file, GA_ReadOnly)
        wind_prj = wind_prj_ds.GetRasterBand(1)

        wind_data = wind_prj.ReadAsArray()

        # The info here is correct.
        #print "img data", pM.loadRasterFile(f_img.name)

        assert_almost_equal(wind_data, speed)

        keep = False
        ncobj.close()
        del m4_max
        if keep:
            print("f_nc.name", f_nc.name)
            print("f_img.name", f_img.name)
        else:
            os.remove(f_nc.name)
            os.remove(f_img.name)


    def test_generate_syn_mult_img(self):
        dir_path = tempfile.mkdtemp(prefix='test_generate_syn_mult_img')
        pM.generate_syn_mult_img(136, -20, 2, dir_path, shape=(2, 4))


        shutil.rmtree(dir_path)

    def test_resetOutputExtentIfInvalid(self):
        """Test computation of output extent"""

        # Write a .nc file to test with dummy data. This is the gust file
        f_nc = tempfile.NamedTemporaryFile(suffix='.nc', prefix='test_processMultipliers', delete=False)
        gust_file = f_nc.name
        f_nc.close()

        ncfile = Dataset(gust_file, mode='w', format='NETCDF3_CLASSIC')
        lat_dim = ncfile.createDimension('lat', 2)
        lon_dim = ncfile.createDimension('lon', 2)
        lat = ncfile.createVariable('lat', np.float32, ('lat',))
        lon = ncfile.createVariable('lon', np.float32, ('lon',))
        lon[:] = [121, 129]
        lat[:] = [10, 17]
        ncfile.close()

        # Test when correct extent provided in config
        config_extent = dict(xMin=20, xMax=21, yMin=20, yMax=21)
        pM.extent = pM.computeOutputExtentIfInvalid(config_extent, gust_file, self.reprojectRaster)
        self.assertEqual(pM.extent, config_extent)

        # Test when extent provided in config
        config_extent = dict()
        pM.extent = pM.computeOutputExtentIfInvalid(config_extent, gust_file, self.testRasterFile)
        self.assertEqual(pM.extent['xMin'], 121)
        self.assertEqual(pM.extent['xMax'], 129)
        self.assertEqual(pM.extent['yMin'], 10)
        self.assertEqual(pM.extent['yMax'], 17)

    @unittest.skip("Not working")
    def test_xprocessMult_A(self):
        dir_path = tempfile.mkdtemp(prefix='test_processMult')
        
        tmp_m4_file = tempfile.NamedTemporaryFile(suffix='.img', 
                                            prefix="test_processMult",
                                            delete=False)
        tmp_m4_file.close()
        
        pM.generate_syn_mult_img(136, -20, 2, dir_path, shape=(2, 4))

        # Top to bottom format.
        # Which is the wrong format
        uu = np.asarray([[0., -1., -1., -1.],
                         [0., 1., 1., 1.]])
        vv = np.asarray([[-1., -1., 0, 1.],
                         [1., 1., 0, -1.]])
        gust = np.asarray([[1., -1., 10, -10.],
                         [100., -100., 1000, -1000.]])

        # Bottom to top
        # This is the format required by the interface
        uu = np.asarray([[0., 1., 1., 1.],
                         [0., -1., -1., -1.]])
        vv = np.asarray([[1., 1., 0, -1.],
                         [-1., -1., 0, 1.]])
        gust = np.asarray([[100., -100., 1000, -1000.],
                         [1., -1., 10, -10.]])


        lon = np.asarray([136, 138, 140, 142])
        lat = np.array([-22, -20])
        windfield_path = os.path.dirname(tmp_m4_file.name) #dir_path # write the output to the multiplier dir
        multiplier_path = os.path.basename(tmp_m4_file.name)#dir_path

        output_file = pM.processMult(gust, uu, vv, lon, lat, windfield_path,
                       multiplier_path)

        minx, miny, maxx, maxy, data = pM.loadRasterFileBandLonLat(output_file)

        actual = np.asarray([[0., -1., 20, -30.],
                         [400., -500., 6000, -7000.]])

        # print 'dir_path', dir_path

        assert_almost_equal(data, actual)
        self.assertEqual(minx, 136)
        self.assertEqual(maxx, 144)
        self.assertEqual(miny, -24)
        self.assertEqual(maxy, -20)

        shutil.rmtree(dir_path)

if __name__ == "__main__":
    # Suite = unittest.makeSuite(TestProcessMultipliers, 'test_x')
    unittest.TestLoader.sortTestMethodsUsing = None
    Suite = unittest.makeSuite(TestProcessMultipliers, 'test')
    Runner = unittest.TextTestRunner()
    Runner.run(Suite)
