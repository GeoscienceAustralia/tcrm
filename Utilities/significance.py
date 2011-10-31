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

 Title: significance.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 05/23/08 5:08:PM
 Description:

 Version :$Rev: 512 $

 $Id: significance.py 512 2011-10-31 07:20:38Z nsummons $
"""
import os, sys, pdb
sys.path.append(os.environ.get('CDATLIB'))
#sys.path.append(os.environ.get('PYTHONSTARTUP'))
import cdms, cdutil, genutil, cdtime
from pyclimate.tools import *
import Numeric
from numpy import *
from log_utils import *
__version__ = '$Id: significance.py 512 2011-10-31 07:20:38Z nsummons $'
def significance(data,testfun=ttest, sig=0.05):
    """significance(data, testfun=ttest, sig=0.05):
    Calculate the statistical significance of the difference between the
    mean of a set of datasets and the individual datasets. This will
    loop through the directory sttructure looking for folders of the
    given year, and select the files containing the given variable in
    its name (and hopefully contained within the file!).
    """
    result = empty(shape(data))
    testdata = Numeric.empty((1, shape(data)[1], shape(data)[2]))
    for n in range(shape(data)[0]):
        testdata[0] = data[n]
        result[n] = test(data, testdata, sig)
    return result

def ttest(mu,sigma, n, value, significance=0.05):
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

if __name__=="__main__":
    inputDir = sys.argv[1]
    year = sys.argv[2]
    variable = sys.argv[3]
    filelist = []
    for root, dirs, files in os.walk(inputDir):
        if year in root:
            for f in files:
                if f.startswith(variable):
                    filelist.append(os.path.join(root, f))
    nfiles = len(filelist)

    # Open a test file to get the dimensions first!
    # We also assume that the files all have the same dimensions
    fid1 = cdms.open(filelist[0])
    nlon = len(fid1['lon'])
    nlat = len(fid1['lat'])
    fid1.close()

    data = Numeric.empty((nfiles, nlat, nlon))
    result = Numeric.empty((nfiles, nlat, nlon))
    i = 0
    for f in filelist:
        fid = cdms.open(f)
        data[i] = fid[variable].getValue()
        fid.close()
        i += 1
    result = significance(data, testfun=ttest, sig=0.05)
