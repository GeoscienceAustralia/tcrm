"""
:mod:`nctools` -- NetCDF utility functions
===========================================

.. module:: nctools
    :synopsis: This modlue contains basic utilities for manipulating
               netCDF format files. It relies on the netCDF4 module
               (https://code.google.com/p/netcdf4-python/)

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import logging

from netCDF4 import Dataset
import numpy as np
import time
import getpass

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

ISO_FORMAT = '%Y-%m-%d %H:%M:%S'

def ncLoadFile(filename):
    """
    Load a netCDF file and return a :class:`netCDF4.Dataset` object.

    :param str filename: Path to the netCDF file to open.
    :return: :class:`netCDF4.Dataset` object
    :rtype: :class:`netCDF4.Dataset`

    """

    logger.debug(f"Opening netCDF file {filename} for reading")

    try:
        ncobj = Dataset(filename, mode='r')
    except (IOError, RuntimeError):
        logger.exception(f"Cannot open {filename}")
        raise IOError

    return ncobj

def ncFileInfo(filename, group=None, variable=None, dimension=None):
    """
    Print summary information about a netCDF file.

    Based on ncinfo (https://code.google.com/p/netcdf4-python/source/browse/trunk/utils/ncinfo)

    :param str filename: Path to valid netCDF file.
    :param str group: Name of a `netCDF4.Group` instance to describe.
    :param str variable: Name of a `netCDF4.Variable` instance to describe.
    :param str dimension: Name of a `netCDF4.Dimension` instance to describe.

    """

    def getgrp(g, p):
        grps = p.split("/")
        for gname in grps:
            if gname == "": continue
            g = g.groups[gname]
        return g

    f = Dataset(filename)
    if group is None:
        if variable is None and dimension is None:
            print(f)
        else:
            if variable is not None:
                print((f.variables[variable]))
            if dimension is not None:
                print((f.dimensions[dimension]))
    else:
        if variable is None and dimension is None:
            print((getgrp(f, group)))
        else:
            g = getgrp(f, group)
            if variable is not None:
                print((g.variables[variable]))
            if dimension is not None:
                print((g.dimensions[variable]))
    f.close()


def ncGetDims(ncobj, dim, dtype=float):
    """
    Extract the value of a dimension from a netCDF file. This function
    assumes the file is written following the CF convention, with the
    values of the dimension stored in a 1-d variable with the same name.

    :param ncobj: :class:`netCDF4.Dataset` object
    :type ncobj: :class:`netCDF4.Dataset`
    :param str dim: Name of the desired dimension.
    :param dtype: Data type of the dimension
    :type dtype: :class:`numpy.dtype`

    :return: :class:`numpy.ndarray` of the requested dimension.

    """
    try:
        data = ncobj.variables[dim][:]
    except KeyError:
        logger.exception(f"Dimension {dim} not in file")
        raise

    return np.array(data, copy=True, dtype=dtype)

def ncGetData(ncobj, var):
    """
    Extract data values from a variable in a netCDF file.

    Note that the variable object is a better way to manipulate
    the variables, as the object includes all attributes (e.g. units,
    range, long_name, etc) - use `ncGetVar` for that purpose.

    :param ncobj: :class:`NetCDF4.Dataset` object.
    :type ncobj: :class:`NetCDF4.Dataset`
    :param str var: Name of the variable in the dataset to extract.

    :return: `numpy.masked_array` containing the data of the variable,
              with missing values masked.
    :rtype: `numpy.ndarray`

    :raises KeyError: If variable does not exist in the given file.

    """

    try:
        varobj = ncobj.variables[var]
    except KeyError:
        logger.exception(f"{ncobj.filepath()} does not contain variable {var}")
        raise

    # Automatic conversion of masked values and unpacking
    varobj.set_auto_maskandscale(True)

    # Get the data:
    d = varobj[:]

    data = np.array(d, copy=True, dtype=varobj.dtype)

    return data

def ncGetVar(ncobj, name):
    """
    Return a `netCDF4.variable` object.

    :param ncobj: :class:`netCDF.Group` or :class:`netCDF.Dataset` instance.
    :type ncobj: :class:`netCDF.Group` or :class:`netCDF.Dataset`
    :param str name: Name of the desired variable.

    :return varobj: :class:`netCDF.Variable` instance
    :rtype: :class:`netCDF.Variable`

    """
    try:
        varobj = ncobj.variables[name]
    except KeyError:
        logger.exception(f"{ncobj.filepath()} does not contain variable {name}")
        raise

    return varobj

def ncGetTimes(ncobj, name='time'):
    """
    Get the time data from a netcdf file.

    :param  ncobj: :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance.
    :param str name: Name of the time variable.

    :return times: Array of time dimension values as true Python :class:`datetime` objects.
    :rtype: :class:`numpy.ndarray` of :class:`datetime` objects

    """

    from datetime import datetime
    from cftime import num2pydate

    if hasattr(ncobj, 'Convention'):
        if getattr(ncobj, 'Convention') == "COARDS":
            # COARDS Compliant file - makes examining the data easier.
            times = ncobj.variables['time']
    elif name in ncobj.dimensions:
        times = ncobj.variables[name]
    else:
        logger.debug( f"Unknown time variable name {name}" )

    if hasattr(times, 'units'):
        units = times.units
    if hasattr(times, 'calendar'):
        calendar = times.calendar
    else:
        calendar = 'standard'

    dates = num2pydate(times[:].data, units, calendar)

    return np.array(dates, dtype=datetime)

def ncCreateDim(ncobj, name, values, dtype, atts=None):
    """
    Create a `dimension` instance in a :class:`netcdf4.Dataset` or :class:`netcdf4.Group` instance.

    :param ncobj: :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance.
    :param str name: Name of the dimension.
    :param `numpy.ndarray` values: Dimension values.
    :param `numpy.dtype` dtype: Data type of the dimension.
    :param atts: Attributes to assign to the dimension instance
    :type atts: dict or None

    """

    ncobj.createDimension(name, np.size(values))
    varDim = (name,)
    dimension = ncCreateVar(ncobj, name, varDim, dtype)
    dimension[:] = np.array(values, dtype=dtype)
    if atts:
        dimension.setncatts(atts)

def ncCreateVar(ncobj, name, dimensions, dtype, data=None, atts=None, **kwargs):
    """
    Create a `Variable` instance in a :class:`netCDF4.Dataset` or
    :class:`netCDF4.Group` instance.


    :param ncobj: :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance where the variable will be stored.
    :type ncobj: :class:`netCDF4.Dataset` or :class:`netCDF4.Group`
    :param str name: Name of the variable to be created.
    :param tuple dimensions: dimension names that define the structure of the variable.
    :param dtype: :class:`numpy.dtype` data type.
    :type dtype: :class:`numpy.dtype`
    :param data: :class:`numpy.ndarray` Array holding the data to be stored.
    :type data: :class:`numpy.ndarray` or None.
    :param dict atts: Dict of attributes to assign to the variable.
    :param kwargs: additional keyword args passed directly to the
                   :class:`netCDF4.Variable` constructor

    :return: :class:`netCDF4.Variable` instance
    :rtype: :class:`netCDF4.Variable`

    """
    logger.debug(f"Creating variable {name}")

    var = ncobj.createVariable(name, dtype, dimensions, **kwargs)

    if data:
        var[:] = np.array(data, dtype=dtype)

    if atts:
        var.setncatts(atts)

    return var

def ncSaveGrid(filename, dimensions, variables, nodata=-9999,
                datatitle=None, gatts={}, writedata=True,
                keepfileopen=False, zlib=True, complevel=4, lsd=None):
    """
    Save a gridded dataset to a netCDF file using NetCDF4.

    :param str filename: Full path to the file to write to.
    :param dimensions: :class:`dict`
        The input dict 'dimensions' has a strict structure, to
        permit insertion of multiple dimensions. The dimensions should be keyed
        with the slowest varying dimension as dimension 0.

        ::

            dimesions = {0:{'name':
                            'values':
                            'dtype':
                            'atts':{'long_name':
                                    'units':  ...} },
                         1:{'name':
                            'values':
                            'type':
                            'atts':{'long_name':
                                    'units':  ...} },
                                  ...}

    :param variables: :class:`dict`
        The input dict 'variables' similarly requires a strict structure:

        ::

            variables = {0:{'name':
                            'dims':
                            'values':
                            'dtype':
                            'atts':{'long_name':
                                    'units':
                                    ...} },
                         1:{'name':
                            'dims':
                            'values':
                            'dtype':
                            'atts':{'long_name':
                                    'units':
                                    ...} },
                             ...}

        The value for the 'dims' key must be a tuple that is a subset of
        the dimensions specified above.

    :param float nodata: Value to assign to missing data, default is -9999.
    :param str datatitle: Optional title to give the stored dataset.
    :param gatts: Optional dictionary of global attributes to include in the file.
    :type gatts: `dict` or None
    :param dtype: The data type of the missing value. If not given, infer from other input arguments.
    :type dtype: :class:`numpy.dtype`
    :param bool writedata: If true, then the function will write the provided data
        (passed in via the variables dict) to the file. Otherwise, no data is
        written.

    :param bool keepfileopen:  If True, return a netcdf object and keep the file open, so that data
        can be written by the calling program. Otherwise, flush data to disk and close the file.

    :param bool zlib: If true, compresses data in variables using gzip compression.

    :param integer complevel: Value between 1 and 9, describing level of compression desired.
         Ignored if zlib=False.

    :param integer lsd: Variable data will be truncated to this number of significant digits.

    :return: `netCDF4.Dataset` object (if keepfileopen=True)
    :rtype: :class:`netCDF4.Dataset`

    :raises KeyError: If input dimension or variable dicts do not have required keys.
    :raises IOError: If output file cannot be created.
    :raises ValueError: if there is a mismatch between dimensions and shape of values to write.

    """

    try:
        ncobj = Dataset(filename, 'w', format='NETCDF4', clobber=True)
    except IOError:
        raise IOError(f"Cannot open {filename} for writing")

    # Dict keys required for dimensions and variables
    dimkeys = set(['name', 'values', 'dtype', 'atts'])
    varkeys = set(['name', 'values', 'dtype', 'dims', 'atts'])

    dims = ()
    for d in dimensions.values():
        missingkeys = [x for x in dimkeys if x not in list(d.keys())]
        if len(missingkeys) > 0:
            ncobj.close()
            raise KeyError(f"Dimension dict missing key '{missingkeys}'")

        ncCreateDim(ncobj, d['name'], d['values'], d['dtype'], d['atts'])
        dims = dims + (d['name'],)

    for v in variables.values():
        missingkeys = [x for x in varkeys if x not in list(v.keys())]
        if len(missingkeys) > 0:
            ncobj.close()
            raise KeyError(f"Variable dict missing key '{missingkeys}'")

        if v['values'] is not None:
            if (len(v['dims']) != v['values'].ndim):
                ncobj.close()
                raise ValueError("Mismatch between shape of "
                                 "variable and dimensions")
        if 'least_significant_digit' in v:
            varlsd = v['least_significant_digit']
        else:
            varlsd = lsd

        var = ncobj.createVariable(v['name'], v['dtype'],
                                   v['dims'],
                                   zlib=zlib,
                                   complevel=complevel,
                                   least_significant_digit=varlsd,
                                   fill_value=nodata)

        if (writedata and v['values'] is not None):
            var[:] = np.array(v['values'], dtype=v['dtype'])

        var.setncatts(v['atts'])

    # Additional global attributes:
    gatts['created_on'] = time.strftime(ISO_FORMAT, time.localtime())
    gatts['created_by'] = getpass.getuser()
    gatts['Conventions'] = 'CF-1.6'

    ncobj.setncatts(gatts)

    if datatitle:
        ncobj.setncattr('title', datatitle)

    if keepfileopen:
        return ncobj
    else:
        ncobj.close()
        return
