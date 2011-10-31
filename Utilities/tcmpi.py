#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

 Title: tcgpi.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2007-05-15
 Description: Determine genesis potential index.

 Reference:
 Emanuel, K.A. and D.S. Nolan (2004) Tropical cyclone activity and the global climate system
    Proc. 26th Conf. on Hurricanes & Trop. Meteor.

 SeeAlso:
 Constraints:

 Version: $Rev: 512 $
 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: tcmpi.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, sys, math, pdb
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
sys.path.append('/home/carthur/sandpit/climate/lib/python/')
sys.path.append(os.environ.get('CDATLIB'))
import cdms, cdtime, cdutil, regrid, MA, MV
import Numeric
import my_tool as myutils
import nctools
import met
import pcmin
import vorticity
from datetime import * # NWS: Module appears to be unused?  Also name conflict between datetime.time & module 'time' when imported this way.
import time
from numpy import *
import numpy, pylab
__version__ = '$Id: tcmpi.py 512 2011-10-31 07:20:38Z nsummons $'
def tcmpi(infile, verbose=False):
    sLat = -50.
    nLat = 25
    dLat = 2.0
    sLon = 90.0
    nLon = 45
    dLon = 2.0

    # This is a CDML file for the corresponding dataset
    fin = cdms.open(infile)
    vars = ['ts', 'ta', 'psl', 'hus']
    for v in vars:
        if v not in fin.variables.keys():
            if verbose:
                print "%s is not in dataset - bailing out"%var
            return 1
    path,base = os.path.split(infile)
    base,ext = os.path.splitext(base)
    model,scen,run,region = base.split('.')

    startdate = cdtime.reltime(fin['time'][0], fin['time'].units)
    enddate = cdtime.reltime(fin['time'][-1], fin['time'].units)
    startyr = str(startdate.tocomp().year)
    endyr = str(enddate.tocomp().year)

    outputfile = os.path.join(path,'mpi.'+model+'.'+scen+'.nc')
    # Delete existing output file to avoid overwrite errors
    if outputfile and os.path.isfile(outputfile):
        os.unlink(outputfile)

    lat = MA.arange(sLat, sLat+nLat*dLat, dLat,typecode='f')
    lon = MA.arange(sLon, sLon+nLon*dLon, dLon,typecode='f')
    try:
        plev = met.convert(fin['plev'][:], fin['plev'].units, 'hPa')
    except TypeError:
        if verbose:
            raise
        return 1
    times = fin['time'][:]
    nlevs = len(plev)
    ntimes = len(times)

    # Create the uniform grid to regrid the data to:
    outgrid = cdms.createUniformGrid(sLat, nLat, dLat, sLon, nLon, dLon)

    # We assume that all variables have their WCRP CMIP3 ids:
    if verbose:
        print "Regridding all variables..."
    try:
        tos = fin['ts'].regrid(outgrid)
    except AttributeError:
        if verbose:
            raise
        return 1
    try:
        ta = fin['ta'].regrid(outgrid)
    except AttributeError:
        if verbose:
            raise
        return 1
    try:
        psl = fin['psl'].regrid(outgrid)
    except AttributeError:
        if verbose:
            raise
        return 1
    try:
        hus = fin['hus'].regrid(outgrid)
    except AttributeError:
        if verbose:
            raise
        return 1

    # Convert data to the required units:
    tosdata = met.convert(numpy.array(tos.getValue()), tos.units, 'C')
    tadata = met.convert(numpy.array(ta.getValue()), ta.units, 'C')
    psldata = met.convert(numpy.array(psl.getValue()), psl.units, 'hPa')
    shdata = numpy.array(hus.getValue())

    if verbose:
        print "Calculating intermediate fields..."


    mrdata = numpy.zeros(shape(tadata))
    pmin = MV.zeros(shape(psldata), typecode='f')
    vmax = MV.zeros(shape(psldata), typecode='f')


    if verbose:
        print "Calculating indices..."
    for n in range(ntimes):
        for i in range(nLat):
            for j in range(nLon):
                try:
                    mrdata[n,:,i,j] = met.spHumToMixRat(shdata[n,:,i,j], 'kgkg')
                except IndexError:
                    fin.close()
                    if verbose:
                        raise
                    return 1
                try:
                    pmin[n,i,j], vmax[n,i,j], status = pcmin.pcmin(tosdata[n,i,j], psldata[n,i,j], plev,
                                                                   tadata[n,:,i,j], mrdata[n,:,i,j],
                                                                   idiss=1, sig=1.0, b=2.0)
                except IndexError:
                    fin.close()
                    if verbose:
                        raise
                    return 1
                if (status !=1) or (vmax[n,i,j] < 10e-4):
                    pmin[n,i,j] = vmax[n,i,j] = 1.0e36

    if verbose:
        print "Output file: %s"%outputfile
    fout = cdms.open(outputfile, 'w')
    # Create the axes for the output dataset:
    fout.createAxis('lat', lat)
    fout.createAxis('lon', lon)
    fout.createAxis('time', times,cdms.Unlimited)

    if verbose:
        print "Creating output variables..."
    vardims = [fout['time'], fout['lat'], fout['lon']]
    fout.createVariable('vmax', 'd', vardims, fill_value=1.0e36)
    fout.createVariable('pmin', 'd', vardims, fill_value=1.0e36)

    fout['vmax'].units = 'm/s'
    fout['vmax'].long_name = 'Maximum theoretical wind speed'

    fout['pmin'].units = 'hPa'
    fout['pmin'].long_name = 'Minimum theoretical central pressure'

    # Set the units for the axes:
    fout['time'].units = fin['time'].units
    fout['lon'].units = fin['lon'].units
    fout['lat'].units = fin['lat'].units

    # Calendar:
    fout['time'].calendar = fin['time'].calendar

    # Standard names:
    fout['lon'].standard_name = fin['lon'].standard_name
    fout['lat'].standard_name = fin['lat'].standard_name

    # Copy across the global attibutes from the input dataset:
    fout.history = fin.history
    fout.title = fin.title
    #fout.source = fin.source
    #fout.experiment = fin.experiment
    #fout.comment = fin.comment

    # A couple of additional attributes:
    #fout.source = cfg['Output.Source']
    fout.CreatedDate = time.ctime()
    fout.CreatedBy = myutils.getUser()
    #fout.sync()

    if verbose:
        print "Writing variables to %s"%outputfile
    fout['vmax'].assignValue(vmax)
    fout['pmin'].assignValue(pmin)

    fin.close()
    fout.close()
    return 0

if __name__ == "__main__":
    infile = sys.argv[1]
    rc = tcmpi(infile, verbose=True)
