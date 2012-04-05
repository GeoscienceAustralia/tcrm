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

 Title: tcgpi.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2007-05-15
 Description: Determine genesis potential index.

 Reference:

 Emanuel, K.A. and D.S. Nolan (2004) Tropical cyclone activity and the
  global climate system, Proc. 26th Conf. on Hurricanes & Trop. Meteor.

 SeeAlso:
 Constraints:

 Version: $Rev: 686 $
 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: tcgpi.py 686 2012-03-29 04:24:59Z carthur $
"""
import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

sys.path.append(os.environ.get('CDATLIB'))

import cdms, cdtime, cdutil, regrid, MA, MV
import Numeric
import my_tool as myutils
import nctools
import metutils
import pcmin
import vorticity
from datetime import *  # NWS: Module appears to be unused?  Also name conflict between datetime.time & module 'time' when imported this way.
import time
import numpy, pylab

__version__ = '$Id: tcgpi.py 686 2012-03-29 04:24:59Z carthur $'
logger = logging.getLogger()

def getVar(obj, outgrid, varname, verbose=False):
    try:
        var = obj[varname].regrid(outgrid)
    except AttributeError:
        if verbose:
            raise
        return None
    else:
        return var

def tcgpi(infile, verbose=False, sLat=-50.0, nLat=25, dLat=2.0,
          sLon=90.0, nLon=45, dLon=1.0):

    # This is a CDML file for the corresponding dataset
    fin = cdms.open(infile)
    vars = ['ts', 'ta', 'psl', 'hus', 'ua', 'va']
    for v in vars:
        if v not in fin.variables.keys():
            logger.critical("%s is not in dataset - bailing out"%v)
            return 1
    path,base = os.path.split(infile)
    base,ext = os.path.splitext(base)
    model,scen,run,region = base.split('.')

    startdate = cdtime.reltime(fin['time'][0], fin['time'].units)
    enddate = cdtime.reltime(fin['time'][-1], fin['time'].units)
    startyr = str(startdate.tocomp().year)
    endyr = str(enddate.tocomp().year)

    outputfile = os.path.join(path,'gpi.'+model+'.'+scen+'.nc')
    # Delete existing output file to avoid overwrite errors
    if outputfile and os.path.isfile(outputfile):
        os.unlink(outputfile)

    lat = MA.arange(sLat, sLat+nLat*dLat, dLat, typecode='f')
    lon = MA.arange(sLon, sLon+nLon*dLon, dLon, typecode='f')
    try:
        plev = metutils.convert(fin['plev'][:], fin['plev'].units,'hPa')
    except TypeError:
        logger.critical("Input pressure data is of incorrect type")
        return 1
    times = fin['time'][:]
    nlevs = len(plev)
    ntimes = len(times)

    # Create the uniform grid to regrid the data to:
    outgrid = cdms.createUniformGrid(sLat, nLat, dLat, sLon, nLon, dLon)

    # We assume that all variables have their WCRP CMIP3 ids:
    logger.debug("Regridding all variables...")
    #tos = getVar(fin, outgrid, 'ts' )

    try:
        tos = fin['ts'].regrid(outgrid)
    except AttributeError:
        return 1
    try:
        ta = fin['ta'].regrid(outgrid)
    except AttributeError:
        return 1
    try:
        psl = fin['psl'].regrid(outgrid)
    except AttributeError:
        return 1
    try:
        hus = fin['hus'].regrid(outgrid)
    except AttributeError:
        return 1
    try:
        ua = fin['ua'].regrid(outgrid)
    except AttributeError:
        return 1
    try:
        va = fin['va'].regrid(outgrid)
    except AttributeError:
        return 1
    # Convert data to the required units:
    tosdata = metutils.convert(numpy.array(tos.getValue()), tos.units,'C')
    tadata = metutils.convert(numpy.array(ta.getValue()), ta.units,'C')
    psldata = metutils.convert(numpy.array(psl.getValue()), psl.units,'hPa')
    shdata = hus.getValue()
    uadata = ua.getValue()
    vadata = va.getValue()

    mrdata = numpy.zeros(shape(tadata), dtype='d')
    rhdata = numpy.zeros(shape(psldata), dtype='d')
    shear = numpy.zeros(shape(psldata[0,:,:]), dtype='d')
    uu = numpy.zeros(shape(psldata[0,:,:]), dtype='d')
    vv = numpy.zeros(shape(psldata[0,:,:]), dtype='d')
    vort = numpy.zeros(shape(psldata[0,:,:]), dtype='d')
    pmin = MV.zeros(shape(psldata), typecode='d')
    vmax = MV.zeros(shape(psldata), typecode='d')
    gpi = MV.zeros(shape(psldata), typecode='d')
    k = numpy.where(plev == 700)[0][0]
    l = numpy.where(plev == 850)[0][0]
    m = numpy.where(plev == 200)[0][0]

    logger.debug("Calculating indices...")
    for n in xrange(ntimes):
        uu = uadata[n,l,:,:]
        vv = vadata[n,l,:,:]
        zeta = vorticity.absolute(numpy.transpose(uu), numpy.transpose(vv),
                                  lon, lat)
        vort[1:-1,1:-1] = zeta[:,:]
        du = uadata[n,l,:,:]-uadata[n,m,:,:]
        dv = vadata[n,l,:,:]-vadata[n,m,:,:]
        shear = Numeric.sqrt(du**2+dv**2)

        for i in xrange(nLat):
            for j in xrange(nLon):
                try:
                    mrdata[n,:,i,j] = metutils.spHumToMixRat(shdata[n,:,i,j], 'kgkg')
                except IndexError:
                    fin.close()
                    return 1
                try:
                    rhdata[n,i,j] = metutils.spHumToRH(shdata[n,k,i,j], tadata[n,k,i,j], 700.0)
                except IndexError:
                    fin.close()
                    return 1
                try:
                    pmin[n,i,j], vmax[n,i,j], status = pcmin.pcmin(tosdata[n,i,j], psldata[n,i,j], plev, tadata[n,:,i,j],
                                                                   mrdata[n,:,i,j], idiss=1, sig=1.0, b=2.0)
                except IndexError:
                    fin.close()
                    return 1
                if (status !=1) or (vmax[n,i,j] < 10e-4):
                    pmin[n,i,j] = vmax[n,i,j] = gpi[n,i,j] = 1.0e36
                else:
                    gpi[n,i,j] = float(((abs(100000.0*vort[i,j]))**1.5)*((rhdata[n,i,j]/50.0)**3.)* \
                                       ((vmax[n,i,j]/70.0)**3.)/((1.0+0.1*shear[i,j])**2.0))

    logger.debug("Output file: %s"%outputfile)
    fout = cdms.open(outputfile, 'w')
    # Create the axes for the output dataset:
    fout.createAxis('lat', lat)
    fout.createAxis('lon', lon)
    fout.createAxis('time', times, cdms.Unlimited)

    logger.debug("Creating output variables...")
    vardims = [fout['time'], fout['lat'], fout['lon']]
    fout.createVariable('vmax', 'd', vardims, fill_value=1.0e36)
    fout.createVariable('pmin', 'd', vardims, fill_value=1.0e36)
    fout.createVariable('gpi', 'd', vardims, fill_value=1.0e36)

    fout['vmax'].units = 'm/s'
    fout['vmax'].long_name = 'Maximum theoretical wind speed'

    fout['pmin'].units = 'hPa'
    fout['pmin'].long_name = 'Minimum theoretical central pressure'
    fout['gpi'].units = ' '
    fout['gpi'].long_name = 'Genesis Potential Index'

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

    logger.debug("Writing variables to %s"%outputfile)

    fout['vmax'].assignValue(vmax)
    fout['pmin'].assignValue(pmin)
    fout['gpi'].assignValue(gpi)

    fin.close()
    fout.close()
    return 0

if __name__ == "__main__":
    infile = sys.argv[1]
    rc = tcgpi(infile, verbose=True)
