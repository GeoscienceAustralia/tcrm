#!/usr/bin/env python
"""
:mod: `nctools` -- NetCDF utility functions
===========================================

This modlue contains basic utilities for manipulating
netCDF format files. It relies on the netCDF4 module



Author: Craig Arthur, craig.arthur@ga.gov.au
CreationDate: 2007-05-15
Description: Tools to interrogate a netcdf file. Largely written to read CF-compliant netcdf files
(specifically NCEP/NCAR Reanalysis files).
Reference: CF-conventions - http://cf-pcmdi.llnl.gov/
           NetCDF Interface to Python - http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html

SeeAlso:
Constraints:

Version: 132
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2008-04-08
Modification: Added additional wrapper functions to raise appropriate errors and eliminate
              missing values that were causing errors in subsequent programs

Version: 187
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-11-27 9:35:AM
Modification: Conversion of print statements to logging module calls

Version: 323
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-12-10 8:57:AM
Modification: No longer relies on the old Numeric module - uses numpy for
              all array manipulation

Version: 420
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-09-23 9:42:AM
Modification: Updated ncSaveGrid to permit writing of multiple variables with
               different dimensions to a single file. This increases the
               capabilities of ncSaveGrid to allow 3- and 4-D data to
               be written to file. Retains existing functionality to
               leave a file open without actually writing data to the
               file. Requires modification to the calling program to
               ensure the dimension and variable metadata are structured
               appropriately.

Version: $Rev: 642 $
ModifiedBy: Nicholas Summons, nicholas.summons@ga.gov.au
ModifiedDate: 06/10/11 2:09:PM
Modification: Replaced ScientificPython netcdf library with SciPy equivalent
               (Note: SciPy added netcdf writing support in version 0.9.0).
               This is necessary as ScientificPython does not work with
               latest versions of NumPy. Also allows fall-back to
               ScientificPython if SciPy version lacks netcdf writing
               support.

$Id: nctools.py 642 2012-02-21 07:54:04Z nsummons $

"""

import os
import sys 
import logging

from scipy.io.netcdf import netcdf_file
#from netCDF4 import Dataset as netcdf_file
from netCDF4 import Dataset
import numpy as np
import time
import getpass

logger = logging.getLogger()

if not hasattr(netcdf_file, 'createVariable'):
    try:
        from Scientific.IO.NetCDF import NetCDFFile as netcdf_file
    except ImportError:
        logger.critical("nctools requires either SciPy version >= 0.9.0 or ScientificPython to write NetCDF files.")
        raise

__version__ = '$Id: nctools.py 642 2012-02-21 07:54:04Z nsummons $'
ISO_FORMAT = '%Y-%m-%d %H:%M:%S'

def _ncLoadFile(filename):
    """
    Load a netCDF file and return a :class:`scipy.io.netcdf.netcdf_file`
    object.

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
        ncobj = netcdf_file(filename, 'r')
    except IOError:
        logger.exception("Cannot open %s"%filename)
        raise IOError
    
    return ncobj

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
        ncobj = Dataset(filename, 'r')
    except IOError:
        logger.exception("Cannot open %s" % filename)
        raise IOError
    
    return ncobj


def _ncCreateFile(filename, op='w'):
    """
    Open a netCDF file for writing.

    Parameters
    ----------

    filename : string
        Path to the netCDF file to create.

    op : string, default 'w'
        Type of operation, 'w' for writing, 'a' for append

    Returns
    -------

    ncobj : :class:`netCDF4.Dataset`
        A new :class:`netCDF4.Dataset` object ready for populating
        with dimensions, variables and attributes.

    """
    logger.debug("Opening netCDF file %s for writing" % filename)
    
    try:
       ncobj = netcdf_file(filename, op)
    except IOError:
        logger.exception("Failed to open %s for writing" % filename)
        raise IOError("Failed to open %s for writing" % filename)
    
    return ncobj

def ncCreateFile(filename, op='w'):
    """
    Open a :class:`netCDF4` file for writing.

    Parameters
    ----------

    filename : string
        Path to the netCDF file to create.

    op : string, default 'w'
        Type of operation, 'w' for writing, 'a' for append

    Returns
    -------

    ncobj : :class:`netCDF4.Dataset`
        A new :class:`netCDF4.Dataset` object ready for populating
        with dimensions, variables and attributes.

    """
    logger.debug("Opening netCDF file %s for writing" % filename)
    
    try:
       ncobj = Dataset(filename, op)
    except IOError:
        logger.exception("Failed to open %s for writing" % filename)
        raise IOError("Failed to open %s for writing" % filename)
    
    return ncobj

def _ncFileInfo(ncobj):
    """
    A short helper function to display the contents of a netCDF file
    This should be adapted to allow testing for required attributes.

    """
    print "Dimensions:"
    for dim in ncobj.dimensions.keys():
        try:
            units = getattr(ncobj.variables[dim],'units')
            print "     %s (%s)"%(dim,units)
        except:
            print "     %s"%(dim)

    print "Variables:"
    for var in ncobj.variables.keys():
        if var not in ncobj.dimensions.keys():
            try:
                units = getattr(ncobj.variables[var],'units')
                print "     %s (%s)"%(var,units)
            except:
                print "     %s"%(var)

    print "Global attributes:"
    for att in dir(ncobj):
        # Exclude the attributes associated with the Scientific.IO.NetCDF
        # module
        if att not in ['close', 'createDimension', 'createVariable',
                       'flush','sync']:
            value = getattr(ncobj, att)
            if type(value) is str:
                print att + ": %s"%value
            elif type(value) is int:
                print att + ": %d"%value
            elif type(value) is float:
                print att + ": %d"%value

def ncFileInfo(filename, group=None, variable=None, dimension=None):
    """Print summary information about a netCDF file.

    Based on ncinfo (https://code.google.com/p/netcdf4-python/source/browse/trunk/utils/ncinfo)
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

def _ncGetData(ncobj, var, missingValue=0.0):
    """
    Extract the data values from a variable in a netCDF file.
    Note that the variable object is a better way to manipulate
    the variables, as the object includes all attributes (e.g. units,
    range, long_name, etc) - use ncGetVar() for that purpose.
    Input: ncobj - netCDF file object created by netcdf_file()
           var - string name of the desired variable
           missingValue - float (or int) value to replace missing values
           with in the returned array. It is possible to check if there
           is a _FillValue (or missing_value) attribute for the variable
           (which is required under the CF Conventions).

    Output:
    data - numpy array of the gridded dataset, with null values replaced
           with missingValue

    """
    try:
        varobj = ncobj.variables[var]
    except KeyError:
        logger.exception("File does not contain variable %s"%(var))
        raise
    
    # Get the data:
    d = varobj[:]
    data = np.array(d,copy=True,dtype=float)

    # Get the attributes to determine offsets and scale factors:
    # We make a slightly dangerous assumption that the scale and offset
    # factors will be named as shown (according to CF convention)
    # Data = var*scale_factor + add_offset
    if hasattr(varobj, 'scale_factor'):
        scale = getattr(varobj, 'scale_factor')
        data = data*scale #[0]
    if hasattr(varobj, 'add_offset'):
        offset = getattr(varobj, 'add_offset')
        data = data + offset #[0]

    # Fix to replace NaN values:
    ind = np.where(np.isnan(data))
    data[ind] = missingValue
    ind = np.where(np.isinf(data))
    data[ind] = missingValue

    return data

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

    # FIXME: Automatic conversion of masked values and unpacking
    #varobj.set_auto_maskandscale(True)
    
    # Get the data:
    d = varobj[:]
    
    data = np.array(d, copy=True, dtype=varobj.dtype)

    return data

def ncGetVar(ncobj, var):
    """Return a `netCDF4.variable` object.
    
    Parameters
    ----------
    
    ncobj : :class:`netCDF.Group` or `netCDF.Dataset` instance.
    
    var : string
        name of the desired variable.

    Returns
    -------
    
    varobj : :class:`netCDF.Variable` instance
    """
    try:
        varobj = ncobj.variables[var]
    except KeyError:
        logger.exception("File does not contain variable %s"%(var))
        raise
    
    return varobj

def _ncGetTimes(ncobj, descriptor='time'):
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
    
    times : :class:`numpy.ndarray`
        Array of time dimension values.

    """
    
    if hasattr(ncobj, 'Convention'):
        if getattr(ncobj, 'Convention') == "COARDS":
            # COARDS Compliant file - makes examining the data easier.
            times = ncobj.variables['time'][:]
    elif descriptor in ncobj.dimensions:
        times = ncobj.variables[descriptor][:]
    else:
        logger.debug( "Unknown time descriptor" )

    return np.copy(times)

def ncGetTimes(ncobj, descriptor='time'):
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
    elif descriptor in ncobj.dimensions:
        times = ncobj.variables[descriptor]
    else:
        logger.debug( "Unknown time descriptor" )

    if hasattr(times, 'units'):
        units = times.units
    if hasattr(times, 'calendar'):
        calendar = times.calendar
    else:
        calendar ='standard'
        
    dates = num2date(times[:], units, calendar)

    return np.array(dates, dtype=datetime)

def _ncCreateDim(ncobj, dimName, values, dtype, atts=None):
    """
    Create a dimension in a netcdf object. Also adds the dimensions data
    to the appropriately named variable.
    The atts variable is a dictionary of additional attribute names and
    values to attach to the dimension.
    Input:
    ncobj - netCDF file object created by netcdf_file()
    dimName - string name of the dimension
    values - numpy array of the dimension values
    dtype - string data type (e.g. 'f', 'd', 'i','l','c','b')
    atts - optional dictionary of attributes to assign to the dimension
    """
    ncobj.createDimension(dimName, np.size(values))
    varDim = (dimName,)
    dimension = ncCreateVar(ncobj, dimName, varDim, dtype)
    dimension[:] = np.array(values,dtype=dtype)
    if atts:
        for key,value in atts.items():
            _setattr(dimension, key, value)

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

def _ncCreateVar(ncobj, varName, varDims, dtype, data=None, atts=None):
    """
    Create a variable in a netcdf object
    The atts variable is a dictionary of additional attribute names and
    values to attach to the variable. One could optionally
    assign those attributes to the variable in the program
    Input:
    ncobj - netCDF file object created by netcdf_file()
    varName - string name of the variable to create
    varDims - tuple of dimension names
    dtype - string data type (e.g. 'f', 'd', 'i','l','c','b')
    atts - optional dictionary of attributes to assign to the variable

    Output:
    var - netCDF variable object
    
    """
    
    logger.debug("Creating variable %s"%varName)
    var = ncobj.createVariable(varName, dtype, varDims)
    
    if data != None:
       var[:] = np.array(data,dtype=dtype)
    if atts:
        for key,value in atts.items():
            _setattr(var, key, value)
            
    return var

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

def _ncSaveGrid(filename, dimensions, variables,
                 nodata=-9999, datatitle=None, gatts=None,
                 dtype='f', writedata=True, 
                 keepfileopen=False):

    """
    Save a gridded dataset to a netCDF file.
    Input:
    filename - full path to the file
    dimensions - dictionary - see below for details
    variables - dictionary - see below for details
    nodata - value to assign to missing data
    datatitle - optional title to give the stored dataset
    gatts - optional dictionary of global attributes to include in the file
    dtype - data type of the missing value
    writedata - boolean. If true, then the function will write
                the provided data (passed in via the variables dict)
                to the file. Otherwise, no data is written
    keepfileopen - boolean. If true, return a netcdf object and keep the
                file open, so that data can be written by the calling program
                Otherwise, flush data to disk and close the file.

    Output:
    Netcdf object (if keepfileopen=True)

    The input dict 'dimensions' has a strict structure, to allow us
    to insert multiple dimensions.
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
    The dimensions should be keyed with the slowest varying dimension
    as dimension 0.

    The input dict 'variables' similarly requires a strict structure:
    variables = {0:{'name':
                    'dims':
                    'values':
                    'dtype':
                    'atts':{'long_name':
                            'units':
                            'scale_factor':
                            'add_offset'    ...} },
                 1:{'name':
                    'dims':
                    'values':
                    'dtype':
                    'atts':{'long_name':
                            'units':
                            'scale_factor':
                            'add_offset':   ...} },
                         ...}
    The value for the 'dims' key must be a tuple that is a subset of
    the dimensions specified above.
    """
    nodata = np.array(nodata, dtype=dtype)
    created_time = time.strftime(ISO_FORMAT, time.localtime())
    created_by = getpass.getuser()


    ncobj = ncCreateFile(filename,op='w')
 
   # Create dimensions:
    dims = ()
    for i in dimensions.iterkeys():
        name = dimensions[i]['name']
        values = dimensions[i]['values']
        datatype = dimensions[i]['dtype']
        atts = dimensions[i]['atts']
        ncCreateDim(ncobj, name, values, datatype, atts)
        dims = dims + (name,)

    # Create the variables:
    for i in variables.iterkeys():
        name = variables[i]['name']
        # Note the dimensions *must* be a subset tuple of the dimension
        # names specified above:
        dims = variables[i]['dims']
        values = variables[i]['values']
        datatype = variables[i]['dtype']
        try:
            atts = variables[i]['atts']
        except:
            pass
        var = ncCreateVar(ncobj, name, dims, datatype)

        try:
            datarange = (values.min(),values.max())
        except AttributeError:
            datarange = (-1*sys.maxint,sys.maxint)

        if writedata:
            var[:] = np.array(values,dtype=datatype)
            
        # Add the attributes to the variable:
        _setattr(var,'_FillValue', nodata)
        _setattr(var,'actual_range', datarange)
        # Plus any other attributes:
        if atts:
            for key,value in atts.items():
                _setattr(var, key, value)

    # Add standard global attributes:
    _setattr(ncobj,'created_by', created_by)
    _setattr(ncobj,'created_on', created_time)
    if datatitle:
        _setattr(ncobj,'title',datatitle)

    if gatts is not None:
        for name,value in gatts.iteritems():
            _setattr(ncobj,name,value)

    # Finally, set the Conventions attribute:
    _setattr(ncobj, 'Conventions', 'CF-1.6')
    
    if not keepfileopen:
        ncobj.close()
    else:
        return ncobj

def _setattr(var, key, value):
    """
    Replace empty strings with None to prevent error when scipy writes to netcdf file
    Input: var - netCDF variable object. If setting a global attribute, this should
                 be the netCDF object itself, rather than a variable object.
           key - name of the attribute to set
           value - value to assign to the attribute.
    """
    if value == '':
        value = 'None'
    setattr(var, key, value)


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
    
    Returns:
    --------
    
    NetCDF4 :class:`Dataset` object (if keepfileopen=True)

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
