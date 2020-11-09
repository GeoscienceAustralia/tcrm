"""
Setup file

"""

from setuptools import setup, find_packages

setup(
    name = "TCRM",
    version = '3.1.4',
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
    'matplotlib',
    'basemap',
    'shapely',
    'nose',
    'netcdf4',
    'cftime',
    'coverage',
    'coveralls',
    'pycurl',
    'pyproj',
    'seaborn',
    'simplejson',
    'sqlite',
    'statsmodels',
    'libgdal'
    'gdal',
    'configparser',
    'cartopy',
    'affine',
    'pandas',
    'tqdm',
    'xarray',
    'pthread-stubs',
    'imageio',
    'mpi4py',
    'boto3',
    'botocore'],
    
    # metadata:
    author = "Craig Arthur",
    author_email = "hazards@ga.gov.au",
    description = "Tropical Cyclone Risk Model",
    keywords = "Tropical cyclone risk hazard",
    url = "https://geoscienceaustralia.github.io/tcrm",
    
    )
