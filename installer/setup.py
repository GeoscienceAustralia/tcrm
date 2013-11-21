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
from os.path import join as pjoin, split as psplit
from glob import glob

opts = {}
py2exe = {}

if 'py2exe' in sys.argv:
    import py2exe

    opts = {
        'py2exe': {
        'includes': ['numpy',
                     'scipy',
                     'matplotlib.backends',
                     'matplotlib.backends.backend_tkagg',
                     'mpl_toolkits.basemap',
                     'scipy.sparse.csgraph._validation',
                     'scipy.io.matlab.streams'],
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
        'console': ['main.py'],
        'windows': ['tcrm.pyw'],
    }

exts = [
    Extension('Utilities.Cmap',
              sources=[pjoin('Utilities', 'Cmap.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-std=c99']),

    Extension('Utilities.Cstats',
              sources=[pjoin('Utilities', 'Cstats.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-std=c99']),

    Extension('Utilities.KPDF',
              sources=[pjoin('Utilities', 'KPDF.c')],
              include_dirs=[pjoin(numpy.get_include(), 'numpy')],
              extra_compile_args=['-std=c99'])
]

basemapData = pjoin('mpl_toolkits', 'basemap', 'data')
data = matplotlib.get_py2exe_datafiles() + \ 
       [(basemapData, glob(pjoin(get_python_lib(), basemapData, '*')))] + \
       [('input', glob('input', '*'))] + \
       [('MSLP', glob('MSLP', '*.nc'))]

setup(name='tcrm',
      version='1.0',
      options=opts,
      ext_modules=exts,
      data_files=data,
      **py2exe)
