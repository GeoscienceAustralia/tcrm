"""
Setup file.

To build C extensions in-place:

    python installer/setup.py build_ext -i

"""
import matplotlib
import numpy
import sys
# from setuptools import Extension
from numpy.distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
from os.path import join as pjoin
from glob import glob

opts = {}
py2exe = {}

if 'py2exe' in sys.argv:
    import py2exe

    opts = {
        'py2exe': {
        'includes': ['numpy',
                     'scipy',
                     'netCDF4',
                     'matplotlib.backends',
                     'matplotlib.backends.backend_tkagg',
                     'mpl_toolkits.basemap',
                     'scipy.sparse.csgraph._validation',
                     'scipy.io.matlab.streams',
                     'boto3',
                     'botocore',
                     'netCDF4_utils'],
        'excludes': ['_gtkagg',
                     'wx',
                     'PyQt4',
                     'nose',
                     'curses',
                     'email'],
        'dll_excludes': ['MSVCP90.DLL'],
        'bundle_files': 3,
        'skip_archive': True
        }
    }

    py2exe = {
        'console': ['tcrm.py'],
    }

exts = [
    Extension('Utilities._akima',
              sources=[pjoin('Utilities', 'akima.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=[]),
    Extension('wind.fwind',
              sources=[pjoin('wind', 'fwind.f90')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-g']),
    Extension('Utilities.fmaputils',
                  sources=[pjoin('Utilities', 'maputils.f90')],
                  include_dirs=[pjoin(numpy.get_include(), 'numpy')],
                  extra_compile_args=['-g']),
    Extension('PressureInterface.fpressureProfile',
                  sources=[pjoin('PressureInterface', 'pressureProfile.f90')],
                  include_dirs=[pjoin(numpy.get_include(), 'numpy')],
                  extra_compile_args=['-g'])
]

basemapData = pjoin('mpl_toolkits', 'basemap', 'data')
data = [(basemapData, glob(pjoin(get_python_lib(), basemapData, '*')))] + \
       [('input', glob(pjoin('input', '*')))] + \
       [('MSLP', glob(pjoin('MSLP', '*.nc')))] + \
       [('.', [pjoin('.', 'matplotlibrc')])]

requires = [
    'matplotlib >= 1.1.1',
    'basemap >= 1.0.2',
    'netCDF4 >= 1.0.1',
    'numpy >= 1.7.1',
    'scipy >= 0.12.0']

setup(name='tcrm',
      version='1.0',
      options=opts,
      ext_modules=exts,
      data_files=data,
      install_requires=requires,
      **py2exe)
