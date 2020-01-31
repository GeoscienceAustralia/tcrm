"""
Setup file

"""

from setuptools import setup, find_packages

setup(
    name = "TCRM",
    version = '3.0rc3',
    packages=find_packages(), 
    scripts=['tcrm.py', 'tcevent.py'],
    include_package_data=True,
    package_data = {
        '' : ['input/*'],
        'tests/test_data' : ['*'],
        'MSLP/' : ['*.nc']
    },
    
    install_requires = [
    'numpy',
    'scipy',
    'pandas',
    'cartopy'
    'affine',
    'matplotlib',
    'basemap',
    'netcdf4',
    'cftime',
    'configparser',
    'gdal',
    'pycurl',
    'pyproj',
    'seaborn',
    'shapely',
    'simplejson',
    'sqlite',
    'statsmodels',
    'tqdm',
    'xarray',
    'nose',
    'coverage',
    'coveralls'],
    
    # metadata:
    author = "Craig Arthur",
    author_email = "craig.arthur@ga.gov.au",
    description = "Tropical Cyclone Risk Model",
    keywords = "Tropical cyclone risk hazard",
    url = "https://geoscienceastralia.github.io/tcrm",
    
    )
