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
"""
import os, sys, math, pdb
sys.path.append('/home/carthur/sandpit/climate/lib/python/')
sys.path.append(os.environ.get('CDATLIB'))
import cdms, MV, cdutil
import Numeric
import pyclimate.pydcdflib
import time


def ttest(mu,sigma,n,value,significance=0.05):
    """
Compute the Student t-test for unequal sample sizes, given a mean (mu),
variance (sigma), number of samples (n) and a value
    """
    if value < 10e10:   # Some arbitrary threshold?
        s = Numeric.sqrt(sigma*((1.0/n)+1)) # Pooled variance
        theta = (mu - value) / s
        t = pyclimate.pydcdflib.CDFT()
        t.which = 2
        t.p = significance/2.0
        t.df = n-1
        pyclimate.pydcdflib.pycdft(t)
        tlow = t.t
        t.p = 1.0 - significance/2.0
        pyclimate.pydcdflib.pycdft(t)
        thigh = t.t
        if (theta > thigh) or (theta < tlow):
            return 1
        else:
            return 0
    else:
        return 0

if __name__ == "__main__":
    variables = ['pmin','vmax']
    scenarios = ['a1b','a2','b1']
    years = ['2010','2030','2050','2070','2090']
    for v in variables:
        for s in scenarios:
            for y in years:
                counts = '/nas/NHIP/garnaut/AUST/SRES'+s.upper()+'/'+v+ \
                         '.all.model.ensemble_count.'+y+'.nc'
                means = '/nas/NHIP/garnaut/AUST/SRES'+s.upper()+'/'+v+ \
                        '.all.model.mean.'+s+'.'+y+'.nc'
                sigmas = '/nas/NHIP/garnaut/AUST/SRES'+s.upper()+'/'+v+ \
                         '.all.model.variance.'+s+'.'+y+'.nc'
                values = '/nas/NHIP/garnaut/AUST/SRES'+s.upper()+'/'+v+ \
                         '.all.model.ensemble.'+s+'.'+y+'.nc'
                outputfile = '/nas/NHIP/garnaut/AUST/SRES'+s.upper()+'/' \
                             +v+'.all.model.sig.'+s+'.'+y+'.nc'

                countfid = cdms.open(counts)
                meanfid = cdms.open(means)
                sigmafid = cdms.open(sigmas)
                valuefid = cdms.open(values)

                nrec = len(valuefid['record'][:])
                nlon = len(valuefid['lon'][:])
                nlat = len(valuefid['lat'][:])
                sig = MV.zeros(Numeric.shape(valuefid[v]),typecode='i')

                lon = valuefid['lon'][:]
                lat = valuefid['lat'][:]
                records = valuefid['record'][:]

                for rec in xrange(nrec):
                    for i in xrange(nlat):
                        for j in xrange(nlon):
                            value = valuefid[v][rec,i,j]
                            if abs(value) < 10e10:
                                count = countfid['count'][i,j]
                                mu = meanfid[v][i,j]
                                sigma = sigmafid[v][0,i,j]
                                sig[rec,i,j] = ttest(mu, sigma, count, value,
                                                     significance=0.05)
                            else:
                                sig[rec,i,j] = 0

                fout = cdms.open(outputfile,'w')
                # Create the axes for the output dataset:
                fout.createAxis('lat',lat)
                fout.createAxis('lon',lon)
                fout.createAxis('record',records,cdms.Unlimited)

                vardims = [fout['record'],fout['lat'],fout['lon']]
                fout.createVariable('significance', 'i', vardims,
                                    fill_value=-99999)
                fout['significance'].long_name = 'Significance'
                fout['lon'].units = valuefid['lon'].units
                fout['lat'].units = valuefid['lat'].units
                # Standard names:
                fout['lon'].standard_name = valuefid['lon'].standard_name
                fout['lat'].standard_name = valuefid['lat'].standard_name
                fout.CreatedDate = time.ctime()
                fout['significance'].assignValue(sig)

                fout.close()
                valuefid.close()
                meanfid.close()
                sigmafid.close()
                countfid.close()
