"""
:mod:`tsmultipliers` -- apply site-exposure multipliers to time series output
=============================================================================

.. module:: tsmultipliers
    :synopsis: Multiply the wind speed in a timeseries file by the
               appropriate multiplier values. Still a very rudimentary
               process.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import sys
import logging


import numpy as np
from files import flLoadFile, flSaveFile

def tsmultiply(inputFile):
    """
    Apply multipliers to a single file. Values are combined then written
    back to the source file. 

    :param str inputFile: Path to the input timeseries file. This will need
                          to contain the values of the three multipliers
                          (topography, terrain and shielding) for each of
                          eight directions.
    """
    tsdata = flLoadFile(inputFile, delimiter=',')
    tstep = tsdata[:,0]
    lon = tsdata[:,1]
    lat = tsdata[:,2]
    gust = tsdata[:,3]
    uu = tsdata[:,4]
    vv = tsdata[:,5]
    bear = tsdata[:,6]
    pressure = tsdata[:,7]
    bear = np.mod((180./np.pi)*np.arctan2(-uu,-vv),360.)

    # Multipliers are stored in the data file:
    mse = tsdata[0,8]
    msne = tsdata[0,9]
    msn = tsdata[0,10]
    msnw = tsdata[0,11]
    msse = tsdata[0,12]
    mss = tsdata[0,13]
    mssw = tsdata[0,14]
    msw = tsdata[0,15]
    mze = tsdata[0,16]
    mzne = tsdata[0,17]
    mzn = tsdata[0,18]
    mznw = tsdata[0,19]
    mzse = tsdata[0,20]
    mzs = tsdata[0,21]
    mzsw = tsdata[0,22]
    mzw = tsdata[0,23]
    mhe = tsdata[0,24]
    mhne = tsdata[0,25]
    mhn = tsdata[0,26]
    mhnw = tsdata[0,27]
    mhse = tsdata[0,28]
    mhs = tsdata[0,29]
    mhsw = tsdata[0,30]
    mhw = tsdata[0,31]

    # Combine multipliers into a single value:
    mce = mse * mze * mhe
    mcne = msne * mzne * mhne
    mcn = msn * mzn * mhn
    mcnw = msnw * mznw * mhnw
    mcse = msse * mzse * mhse
    mcs = mss * mzs * mhs
    mcsw = mssw * mzsw * mhsw
    mcw = msw * mzw * mhw

    # Apply multipliers:
    """ii = np.where((bear < 22.5) | (bear >= 337.5))
    gust[ii] *= mcn
    ii = np.where((bear >= 22.5) & (bear < 67.5))
    gust[ii] *= mcne
    ii = np.where((bear >= 67.5) & (bear < 112.5))
    gust[ii] *= mce
    ii = np.where((bear >= 112.5) & (bear < 157.5))
    gust[ii] *= mcse
    ii = np.where((bear >= 157.5) & (bear < 202.5))
    gust[ii] *= mcs
    ii = np.where((bear >= 202.5) & (bear < 247.5))
    gust[ii] *= mcsw
    ii = np.where((bear >= 247.5) & (bear < 292.5))
    gust[ii] *= mcw
    ii = np.where((bear >= 292.5) & (bear < 337.5))
    gust[ii] *= mcnw

    """
    ii = np.where((bear >= 0.0) & (bear < 45.))
    gust[ii] *= (1./45.)*(mcn*(bear[ii]-0.0) +
                          mcne*(45. - bear[ii]))
    ii = np.where((bear >= 45.0) & (bear < 90.))
    gust[ii] *= (1./45.)*(mcne*(bear[ii] - 45.0) +
                          mce*(90. - bear[ii]))
    ii = np.where((bear >= 90.0) & (bear < 135.))
    gust[ii] *= (1./45.)*(mce*(bear[ii] - 90.0) +
                          mcse*(135. - bear[ii]))
    ii = np.where((bear >= 135.0) & (bear < 180.))
    gust[ii] *= (1./45.)*(mcse*(bear[ii] - 135.0) +
                          mcs*(180. - bear[ii]))
    ii = np.where((bear >= 180.0) & (bear < 225.))
    gust[ii] *= (1./45.)*(mcs*(bear[ii] - 180.0) +
                          mcsw*(225. - bear[ii]))
    ii = np.where((bear >= 225.0) & (bear < 270.))
    gust[ii] *= (1./45.)*(mcsw*(bear[ii] - 225.0) +
                          mcw*(270. - bear[ii]))
    ii = np.where((bear >= 270.0) & (bear < 315.))
    gust[ii] *= (1./45.)*(mcw*(bear[ii] - 270.0) +
                          mcnw*(315. - bear[ii]))
    ii = np.where((bear >= 315.0) & (bear <= 360.))
    gust[ii] *= (1./45.)*(mcnw*(bear[ii] - 315.0) +
                          mcn*(360. - bear[ii]))


    ii = np.where(gust==0)
    bear[ii]=0

    data = np.transpose([tstep, lon, lat, gust, uu, vv, bear])
    header = 'Time,Longitude,Latitude,Speed,UU,VV,Bearing'
    flSaveFile(inputFile, data, header, ',', '%f')

if __name__=='__main__':
    try:
        inputPath = sys.argv[1]
    except:
        inputPath = eval(raw_input("Enter input path: "))
    ls = os.listdir(inputPath)
    for f in ls:
        inputFile = os.path.join(inputPath,f)
        tsmultiply(inputFile)

