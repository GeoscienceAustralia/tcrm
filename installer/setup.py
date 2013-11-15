"""
Setup file.

To build C extensions in-place:

    python installer/setup.py build_ext -i

"""
import matplotlib
import numpy
import glob
import sys

from distutils.core import setup, Extension
from os.path import join as pjoin

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
                     'scipy.sparse.csgraph._validation',
                     'scipy.io.matlab.streams'],
        'excludes': ['_gtkagg',
                     '_ssl',
                     'wx',
                     'PyQt4',
                     'nose',
                     'curses',
                     'email',
                     'doctest',
                     'pdb',
                     'difflib'],
        'dll_excludes': ['MSVCP90.DLL'],
        'ascii': True
        }
    }

    py2exe = {
        'console': ['main.py'],
        'windows': ['tcrm.pyw']
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

setup(name='tcrm',
      version='1.0',
      options=opts,
      ext_modules=exts,
      **py2exe)
