"""
Setup file.

To build C extensions in-place:

    python installer/setup.py build_ext -i

"""
import matplotlib
import numpy
import sys
from distutils.core import setup, Extension
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
    Extension('Utilities.Cmap',
              sources=[pjoin('Utilities', 'Cmap.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-std=gnu99']),

    Extension('Utilities.Cstats',
              sources=[pjoin('Utilities', 'Cstats.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-std=gnu99']),

    Extension('Utilities.KPDF',
              sources=[pjoin('Utilities', 'KPDF.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-std=gnu99']),

    Extension('Utilities._akima',
              sources=[pjoin('Utilities', 'akima.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=[])
]

basemapData = pjoin('mpl_toolkits', 'basemap', 'data')
data = matplotlib.get_py2exe_datafiles() + \
       [(basemapData, glob(pjoin(get_python_lib(), basemapData, '*')))] + \
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
