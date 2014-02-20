#!/usr/bin/env python
"""
:mod: `nctools` -- NetCDF utility functions
===========================================

This modlue contains basic utilities for manipulating
netCDF format files. It relies on the netCDF4 module
(https://code.google.com/p/netcdf4-python/)

"""

import os
import sys 
import logging

from netCDF4 import Dataset
import numpy as np
import time
import getpass

logger = logging.getLogger()

ISO_FORMAT = '%Y-%m-%d %H:%M:%S'

def ncLoadFile(filename):
    """
    Load a netCDF file and return a :class:`netCDF4.Dataset` object.

    Parameters
    ----------

    filename : string
        Path to the netCDF file to open.

    Returns
    -------

    ncobj : :class:`netCDF4.Dataset` object

    """
    logger.debug("Opening netCDF file %s for reading"%filename)
    
    try:
        ncobj = Dataset(filename, mode='r')
    except (IOError, RuntimeError):
        logger.exception("Cannot open %s" % filename)
        raise IOError
    
    return ncobj

def ncFileInfo(filename, group=None, variable=None, dimension=None):
    """Print summary information about a netCDF file.

    Based on ncinfo (https://code.google.com/p/netcdf4-python/source/browse/trunk/utils/ncinfo)

    Parameters
    ----------

    filename : str
        Path to valid netCDF file.

    group : str, optional
        Name of a `netCDF4.Group` instance to describe.

    variable: str, optional
        Name of a `netCDF4.Variable` instance to describe.

    dimension : str, optional
        Name of a `netCDF4.Dimension` instance to describe.
    
    """
    def getgrp(g, p):
        import posixpath
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
            if var is not None:
                print(f.variables[variable])
            if dim is not None:
                print(f.dimensions[dimension])
    else:
        if variable is None and dimension is None:
            print(getgrp(f, group))
        else:
            g = getgrp(f, group)
            if variable is not None:
                print(g.variables[variable])
            if dim is not None:
                print(g.dimensions[variable])
    f.close()
    

def ncGetDims(ncobj, dim, dtype=float):
    """
    Extract the value of a dimension from a netCDF file. This function
    assumes the file is written following the CF convention, with the
    values of the dimension stored in a 1-d variable with the same name.

    Parameters
    ----------
    
    ncobj : :class:`netCDF4.Dataset` object
       
    dim : string
       Name of the desired dimension.

    Returns
    -------

    data : :class:`numpy.ndarray` of the requested dimension.

    """
    try:
        data = ncobj.variables[dim][:]
    except KeyError:
        logger.exception( "Dimension %s not in file" % dim )
        raise
    
    return np.array(data, copy=True, dtype=dtype)

def ncGetData(ncobj, var, missingValue=-9999.):
    """Extract data values from a variable in a netCDF file.
    
    Note that the variable object is a better way to manipulate
    the variables, as the object includes all attributes (e.g. units,
    range, long_name, etc) - use `ncGetVar` for that purpose.

    Parameters
    ----------

    ncobj : :class:`NetCDF4.Dataset` object.

    var : str
        Name of the variable in the dataset to extract.

    missingValue : float or int, optional
        Value to assign to missing data, default is -9999.

    Returns
    -------

    data : :class:`numpy.ndarray`
        Array containing the data of the variable, with missing values
        replaced by `missingValue`.

    Raises
    ------

    KeyError : if variable does not exist in the given file.
    
    """
    
    try:
        varobj = ncobj.variables[var]
    except KeyError:
        logger.exception("File does not contain variable %s"%(var))
        raise

    # Automatic conversion of masked values and unpacking
    varobj.set_auto_maskandscale(True)
    
    # Get the data:
    d = varobj[:]
    
    data = np.array(d, copy=True, dtype=varobj.dtype)

    return data

def ncGetVar(ncobj, name):
    """Return a `netCDF4.variable` object.
    
    Parameters
    ----------
    
    ncobj : :class:`netCDF.Group` or `netCDF.Dataset` instance.
    
    name : string
        name of the desired variable.

    Returns
    -------
    
    varobj : :class:`netCDF.Variable` instance
    """
    try:
        varobj = ncobj.variables[name]
    except KeyError:
        logger.exception("File does not contain variable %s"%(name))
        raise
    
    return varobj

def ncGetTimes(ncobj, name='time'):
    """
    Get the time data from the file

    Parameters
    ----------

    ncobj : :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance
        :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance
        
    descriptor : str, default 'time'
        Name of the time variable.

    Returns
    -------
    
    times : :class:`numpy.ndarray` of :class:`datetime` objects
        Array of time dimension values as :class:`datetime` objects.

    """
    
    from datetime import datetime
    from netCDF4 import num2date
    
    if hasattr(ncobj, 'Convention'):
        if getattr(ncobj, 'Convention') == "COARDS":
            # COARDS Compliant file - makes examining the data easier.
            times = ncobj.variables['time']
    elif name in ncobj.dimensions:
        times = ncobj.variables[name]
    else:
        logger.debug( "Unknown time variable name" )

    if hasattr(times, 'units'):
        units = times.units
    if hasattr(times, 'calendar'):
        calendar = times.calendar
    else:
        calendar ='standard'
        
    dates = num2date(times[:], units, calendar)

    return np.array(dates, dtype=datetime)

def ncCreateDim(ncobj, name, values, dtype, atts=None):
    """
    Create a `dimension` instance in a :class:`netcdf4.Dataset` or
    :class:`netcdf4.Group` instance.

    Parameters
    ----------

    ncobj : :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance

    name : str
        Name of the dimension

    values : :class:`numpy.ndarray`
        Dimension values

    dtype : :class:`numpy.dtype`
        Data type of the dimension

    atts : :dict:, optional
        Attributes to assign to the dimension instance

    Returns
    -------

    None

    """
    ncobj.createDimension(name, np.size(values))
    varDim = (name,)
    dimension = ncCreateVar(ncobj, name, varDim, dtype)
    dimension[:] = np.array(values, dtype=dtype)
    if atts:
        dimension.setncatts(atts)

def ncCreateVar(ncobj, name, dimensions, dtype, data=None, atts=None,
                 zlib=True, complevel=4, lsd=2, nodata=-9999.,):
    """
    Create a `Variable` instance in a :class:`netCDF4.Dataset` or
    :class:`netCDF4.Group` instance.

    Parameters
    ----------
    ncobj : :class:`netCDF4.Dataset` or :class:`netCDF4.Group` instance
        where the variable will be stored

    name : str
        Name of the variable to be created.

    dimensions : :tuple: of dimension names that define the structure
        of the variable

    dtype : `numpy.dtype` data type

    data : :class:`numpy.ndarray`, optional
        array holding the data to be stored. 
        
    atts : :dict:, optional
        dict of attributes to assign to the variable

    zlib : bool, default `True`
         If true, compresses data in variables using gzip compression.

    complevel : integer, default 4
         Value between 1 and 9, describing level of compression desired.
         Ignored if zlib=False
    
    lsd : integer, default 2
        Variable data will be truncated to this number of significant digits.    
    nodata : float or int, default -9999.
         Value representing missing data
    
    Returns
    -------

    var : :class:`netCDF4.Variable` instance

    """
    logger.debug("Creating variable %s" % name)
    
    var = ncobj.createVariable(name, dtype, dimensions,
                               zlib=zlib,
                               complevel=complevel,
                               least_significant_digit=lsd,
                               fill_value=nodata)
    if data:
        var[:] = np.array(data, dtype=dtype)
            
    if atts:
        var.setncatts(atts)
            
    return var

def ncSaveGrid(filename, dimensions, variables, nodata=-9999,
                datatitle=None, gatts={}, dtype='f', writedata=True, 
                keepfileopen=False, zlib=True, complevel=4, lsd=2):
    """
    Save a gridded dataset to a netCDF file using NetCDF4.
    
    Parameters
    ----------

    filename : :class:`str`
        Full path to the file to write to    
    dimensions : :class:`dict`
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

    variables : :class:`dict`
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
    
    nodata : float, optional
        Value to assign to missing data, default is -9999.

    datatitle : str, optional
        Optional title to give the stored dataset
    
    gatts : :class:`dict`, optional
        Optional dictionary of global attributes to include in the file.
    
    dtype : dtype, optional
        The data type of the missing value. If not given, infer from other
        input arguments.
    
    writedata : bool, default `True`
        If true, then the function will write the provided data
        (passed in via the variables dict) to the file. Otherwise, no data is
        written.
    
    keepfileopen : bool, default `False`
        If true, return a netcdf object and keep the file open, so that data
        can be written by the calling program. Otherwise, flush data to disk
        and close the file.
    
    zlib : bool, default `True`
         If true, compresses data in variables using gzip compression.

    complevel : integer, default 4
         Value between 1 and 9, describing level of compression desired.
         Ignored if zlib=False
    
    lsd : integer, default 2
        Variable data will be truncated to this number of significant digits.
    
    Returns
    -------
    
    :class:`netCDF4.Dataset` object (if keepfileopen=True)

    Raises
    ------

    KeyError : if input dimension or variable dicts do not have required keys.

    IOError : if output file cannot be created.

    ValueError : if there is a mismatch between dimensions and shape of values
                 to write.
    """

    try:
        ncobj = Dataset(filename, 'w', format='NETCDF4', clobber=True)
    except IOError:
        raise IOError("Cannot open {0} for writing".format(filename))

    # Dict keys required for dimensions and variables
    dimkeys = set(['name', 'values', 'dtype', 'atts'])
    varkeys = set(['name', 'values', 'dtype', 'dims', 'atts'])
    
    dims = ()
    for i, d in dimensions.iteritems():
        missingkeys = [x for x in dimkeys if x not in d.keys()]
        if len(missingkeys) > 0:
            ncobj.close()
            raise KeyError("Dimension dict missing key '{0}'".format(missingkeys))
        
        name = d['name']
        values = d['values']
        datatype = d['dtype']
        atts = d['atts']
        ncCreateDim(ncobj, name, values, datatype, atts)
        dims = dims + (name,)

    for i, v in variables.iteritems():
        missingkeys = [x for x in varkeys if x not in v.keys()]
        if len(missingkeys) > 0:
            ncobj.close()
            raise KeyError("Dimension dict missing key '{0}'".format(missingkeys))
        name = v['name']
        values = v['values']
        datatype = v['dtype']
        atts = v['atts']
        dimensions = v['dims']

        if len(dimensions) != len(values.shape):
            ncobj.close()
            raise ValueError("Mismatch between shape of variable and dimensions")
        
        var = ncobj.createVariable(name, datatype, dimensions,
                                   zlib=zlib,
                                   complevel=complevel,
                                   least_significant_digit=lsd,
                                   fill_value=nodata)

        if writedata:
            var[:] = np.array(values, dtype=datatype)
    
        var.setncatts(atts)

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
