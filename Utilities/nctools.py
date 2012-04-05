#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Title: nctools.py
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
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
from scipy.io.netcdf import netcdf_file
import numpy
import math
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

def ncLoadFile(filename):
    logger.debug("Opening netCDF file %s for reading"%filename)
    try:
        ncobj = netcdf_file(filename, 'r')
    except IOError:
        logger.exception("Cannot open %s"%filename)
        raise
    return ncobj

def ncCreateFile(filename,op='w'):
    logger.debug("Opening netCDF file %s for writing"%filename)
    try:
       ncobj = netcdf_file(filename,op)
    except IOError:
        logger.exception("Failed to open %s for writing"%filename)
        raise
    return ncobj

def ncFileInfo(ncobj):
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
                print att+": %s"%value
            elif type(value) is int:
                print att+": %d"%value
            elif type(value) is float:
                print att+": %d"%value

def ncGetDims(ncobj, dim):
    """
    Extract the value of a dimension from a netCDF file. This function
    assumes the file is written following the CF convention, with the
    values of the dimension stored in a 1-d variable with the same name.
    Input: ncobj - netCDF file object created by netcdf_file()
           dim - string name of the desired dimenstion
    Output: data - numpy array of the requested dimension.

    """
    try:
        data = ncobj.variables[dim][:]
    except KeyError:
        logger.exception( "Dimension %s not in file"%dim )
        raise
    return numpy.array(data,copy=True,dtype=float)

def ncGetData(ncobj, var, missingValue=0.0):
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
    data = numpy.array(d,copy=True,dtype=float)

    # Get the attributes to determine offsets and scale factors:
    # We make a slightly dangerous assumption that the scale and offset
    # factors will be named as shown (according to CF convention)
    # Data = var*scale_factor + add_offset
    if hasattr(varobj, 'scale_factor'):
        scale = getattr(varobj, 'scale_factor')
        data = data*scale[0]
    if hasattr(varobj, 'add_offset'):
        offset = getattr(varobj, 'add_offset')
        data = data + offset[0]

    # Fix to replace NaN values:
    ind = numpy.where(numpy.isnan(data))
    data[ind] = missingValue
    ind = numpy.where(numpy.isinf(data))
    data[ind] = missingValue

    return data

def ncGetVar(ncobj, var):
    """
    Return a netCDF variable object. The returned object includes
    functions (getValue, assignValue) to either extract the data values
    or assign data values, in addition to the attributes of the variable.
    Input:
    ncobj - netCDF file object created by netcdf_file()
    var - string name of the desired variable

    Output:
    varobj - netCDF variable object
    """
    try:
        varobj = ncobj.variables[var]
    except KeyError:
        logger.exception("File does not contain variable %s"%(var))
        raise
    return varobj

def ncGetTimes(ncobj, descriptor='time'):
    """
    Get the time data from the file
    Input:
    ncobj - netCDF file object created by netcdf_file()
    descriptor - string name of the time variable (default is 'time' -
                 following CF convention)

    Output:
    times - numpy array of time dimension values

    """
    if hasattr(ncobj, 'Convention'):
        if getattr(ncobj, 'Convention') == "COARDS":
        # COARDS Compliant file - makes examining the data easier.
        # All times are stored in hours since the start of the analysis
           times = ncobj.variables['time'][:]
    elif descriptor in ncobj.dimensions:
        times = ncobj.variables[descriptor][:]
    else:
        logger.debug( "Unknown time descriptor" )

    return numpy.copy(times)

def ncCreateDim(ncobj, dimName, values, dtype, atts=None):
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
    ncobj.createDimension(dimName, numpy.size(values))
    varDim = (dimName,)
    dimension = ncCreateVar(ncobj, dimName, varDim, dtype)
    dimension[:] = numpy.array(values,dtype=dtype)
    if atts:
        for key,value in atts.items():
            _setattr(dimension, key, value)

def ncCreateVar(ncobj, varName, varDims, dtype, data=None, atts=None):
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
       var[:] = numpy.array(data,dtype=dtype)
    if atts:
        for key,value in atts.items():
            _setattr(var, key, value)
    return var



def ncSaveGrid(filename, lon, lat, data, varname, units,
                 lonunits='degrees_east', latunits='degrees_north',
                 nodata=-9999,longname=None,datatitle=None,dtype='f',
                 writedata=True, keepfileopen=False):

    """
    Save a gridded dataset to a netCDF file.
    Input:
    filename - full path to the file
    lon - array of the longitude points in the grid
    lat - array of the latitude points in the grid
    data - 2D array of the data to be stored. Initially, only 1 record
           is stored in the file
    varname - string which will be the short name of the variable
    units - string describing the units the data is stored in

    Output:
    None
    """
    nodata = numpy.array(nodata, dtype=dtype)
    created_time = time.strftime(ISO_FORMAT, time.localtime())
    created_by = getpass.getuser()
    ncobj = ncCreateFile(filename,op='w')
    latrange = [min(lat), max(lat)]
    lonrange = [min(lon), max(lon)]
    datarange = [data.min(),data.max()]

    # Define attributes for the dimensions:
    lat_atts = {'units':latunits,
                'long_name':'latitude',
                'actual_range': latrange}
    lon_atts = {'units':lonunits,
                'long_name':'longitude',
                'actual_range': lonrange}

    # Create dimensions:
    ncCreateDim(ncobj, 'lat', lat,'f',lat_atts)
    ncCreateDim(ncobj, 'lon', lon,'f',lon_atts)

    # The order of these is significant under the CF conventions
    dims = ('lat','lon',)

    # Create the variable:
    var = ncCreateVar(ncobj, varname, dims, 'f')
    if writedata:
        var[:] = numpy.array(data,dtype=dtype)

    # Add the attributes to the variable:
    _setattr(var,'units',units)
    _setattr(var,'_FillValue', nodata)
    _setattr(var,'actual_range', datarange)
    if longname:
        _setattr(var,'long_name',longname)

    # Add global attributes:
    _setattr(ncobj,'created_by', created_by)
    _setattr(ncobj,'created_on', created_time)
    if datatitle:
        _setattr(ncobj,'title',datatitle)

    # Finally, set the convention attribute:
    _setattr(ncobj,'Conventions','CF-1.4')
    if not keepfileopen:
        ncobj.close()
    else:
        return ncobj

def _ncSaveGrid(filename, dimensions, variables,
                 nodata=-9999,datatitle=None,gatts=None,
                 dtype='f',writedata=True, 
                 keepfileopen=False,packed=False):

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
    packed - boolean. If true, attempts to pack the data according to
             recommended best practices
            (http://www.unidata.ucar.edu/software/netcdf/docs/BestPractices.html)
            At this time, not implemented in working form.
    Output:
    Netcdf object (if keepfileopen=True)

    The input dict 'dimensions' has a strict structure, to allow us
    to insert multiple dimensions.
    dimesions = {0:{'name':
                    'values':
                    'type':
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
    nodata = numpy.array(nodata, dtype=dtype)
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
        """
        if packed:
            var = ncCreateVar(ncobj, name, dims, 'i')
        else:
            var = ncCreateVar(ncobj, name, dims, datatype)
        """
        try:
            datarange = [values.min(),values.max()]
        except AttributeError:
            datarange = [-1*sys.maxint,sys.maxint]

        if writedata:
            var[:] = numpy.array(values,dtype=datatype)
            """
            if packed:
                # Pack the data using a standard packing algorithm:
                n = 8
                scale = (values.max() - values.min())/((2**n)-1)
                offset = values.min() #+ scale*(2**(n-1))
                packdata =((values - offset)/scale).astype(int)
                #pdb.set_trace()
                var[:] = numpy.array(packdata,dtype=int)
                _setattr(var, 'scale_factor', scale)
                _setattr(var, 'add_offset', offset)
            else:
                var[:] = numpy.array(values,dtype=datatype)
            """
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

    # Finally, set the convention attribute:
    _setattr(ncobj,'Conventions','CF-1.4')
    if not keepfileopen:
        ncobj.close()
    else:
        return ncobj

def _setattr(var, key, value):
    # Replace empty strings with None to prevent error when scipy writes to netcdf file
    if value == '':
        value = 'None'
    setattr(var, key, value)
