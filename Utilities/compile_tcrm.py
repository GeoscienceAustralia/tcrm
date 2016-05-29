from distutils.core import setup
import matplotlib
import py2exe
import os
import glob

origIsSystemDLL = py2exe.build_exe.isSystemDLL
def isSystemDLL(pathname):
    if os.path.basename(pathname).lower() in ("msvcp71.dll", "msvcrt.dll"):
        return 0
    return origIsSystemDLL(pathname)
py2exe.build_exe.isSystemDLL = isSystemDLL

setup(console=['main.py', 'configeditor.pyw', '.\\unittests\\test_all.py'],
      data_files=[(r'mpl-data\data',glob.glob(r'C:\Python25\Lib\site-packages\mpl_toolkits\basemap\data\*.*'))] + \
                 matplotlib.get_py2exe_datafiles() + \
                 [(r'test_data',glob.glob(r'.\unittests\test_data\*.*'))] + \
                 [(r'MSLP',glob.glob(r'.\MSLP\*.nc'))] + \
                 [(r'output',[])] + \
                 [(r'input',glob.glob(r'.\input\*.*'))],
      options={'py2exe': { "includes" : ["matplotlib.backends",
                                         "matplotlib.backends.backend_qt4agg",
                                         "matplotlib.figure",
                                         "pylab",
                                         "numpy",
                                         "matplotlib.numerix.fft",
                                         "matplotlib.numerix.linear_algebra",
                                         "matplotlib.numerix.random_array",
                                         "matplotlib.backends.backend_tkagg",
                                         "scipy.io.matlab.streams",
                                         "mpl_toolkits.basemap"]}},
      zipfile="lib/shared.zip")
