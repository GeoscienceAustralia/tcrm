"""
:mod:`tsmultipliers` -- apply site-exposure multipliers to time series output
=============================================================================

.. module:: tsmultipliers
    :synopsis: Multiply the wind speed in a timeseries file by the
               appropriate multiplier values. Still a very rudimentary process. 

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

import numpy
from files import flLoadFile, flSaveFile

__version__ = '$Id: tsmultipliers.py 642 2012-02-21 07:54:04Z nsummons $'

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
    bearing = tsdata[:,6]
    bearing = numpy.mod((180./numpy.pi)*numpy.arctan2(-uu,-vv),360.)

    # Multipliers are stored in the data file:
    mse = tsdata[0,7]
    msne = tsdata[0,8]
    msn = tsdata[0,9]
    msnw = tsdata[0,10]
    msse = tsdata[0,11]
    mss = tsdata[0,12]
    mssw = tsdata[0,13]
    msw = tsdata[0,14]
    mze = tsdata[0,15]
    mzne = tsdata[0,16]
    mzn = tsdata[0,17]
    mznw = tsdata[0,18]
    mzse = tsdata[0,19]
    mzs = tsdata[0,20]
    mzsw = tsdata[0,21]
    mzw = tsdata[0,22]
    mhe = tsdata[0,23]
    mhne = tsdata[0,24]
    mhn = tsdata[0,25]
    mhnw = tsdata[0,26]
    mhse = tsdata[0,27]
    mhs = tsdata[0,28]
    mhsw = tsdata[0,29]
    mhw = tsdata[0,30]

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
    """ii = numpy.where((bearing < 22.5) | (bearing >= 337.5))
    gust[ii] *= mcn
    ii = numpy.where((bearing >= 22.5) & (bearing < 67.5))
    gust[ii] *= mcne
    ii = numpy.where((bearing >= 67.5) & (bearing < 112.5))
    gust[ii] *= mce
    ii = numpy.where((bearing >= 112.5) & (bearing < 157.5))
    gust[ii] *= mcse
    ii = numpy.where((bearing >= 157.5) & (bearing < 202.5))
    gust[ii] *= mcs
    ii = numpy.where((bearing >= 202.5) & (bearing < 247.5))
    gust[ii] *= mcsw
    ii = numpy.where((bearing >= 247.5) & (bearing < 292.5))
    gust[ii] *= mcw
    ii = numpy.where((bearing >= 292.5) & (bearing < 337.5))
    gust[ii] *= mcnw

    """
    ii = numpy.where((bearing >= 0.0) & (bearing<45.))
    gust[ii] *= (1./45.)*(mcn*(bearing[ii]-0.0) + mcne*(45. - bearing[ii]))
    ii = numpy.where((bearing >= 45.0) & (bearing < 90.))
    gust[ii] *= (1./45.)*(mcne*(bearing[ii]-45.0) + mce*(90. - bearing[ii]))
    ii = numpy.where((bearing >= 90.0) & (bearing < 135.))
    gust[ii] *= (1./45.)*(mce*(bearing[ii]-90.0) + mcse*(135. - bearing[ii]))
    ii = numpy.where((bearing >= 135.0) & (bearing < 180.))
    gust[ii] *= (1./45.)*(mcse*(bearing[ii]-135.0) + mcs*(180. - bearing[ii]))
    ii = numpy.where((bearing >= 180.0) & (bearing < 225.))
    gust[ii] *= (1./45.)*(mcs*(bearing[ii]-180.0) + mcsw*(225. - bearing[ii]))
    ii = numpy.where((bearing >= 225.0) & (bearing < 270.))
    gust[ii] *= (1./45.)*(mcsw*(bearing[ii]-225.0) + mcw*(270. - bearing[ii]))
    ii = numpy.where((bearing >= 270.0) & (bearing < 315.))
    gust[ii] *= (1./45.)*(mcw*(bearing[ii]-270.0) + mcnw*(315. - bearing[ii]))
    ii = numpy.where((bearing >= 315.0) & (bearing <= 360.))
    gust[ii] *= (1./45.)*(mcnw*(bearing[ii]-315.0) + mcn*(360. - bearing[ii]))


    ii = numpy.where(gust==0)
    bearing[ii]=0

    data = numpy.transpose([tstep, lon, lat, gust, uu, vv, bearing])
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

