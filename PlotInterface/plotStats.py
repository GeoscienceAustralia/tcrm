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


 Title: plotStats.py - plot basic statistical info based on processed data
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2010-05-10
 Description:

 Version: $Rev: 762 $

 $Id: plotStats.py 762 2011-11-25 05:02:07Z nsummons $
"""

import os, sys, pdb, logging

from matplotlib import pyplot
from scipy.stats import linregress, probplot
from scipy.stats import scoreatpercentile as percentile
import numpy

import Utilities.files as files
import Utilities.config as config
import Utilities.nctools as nctools
import Utilities.stats as stats

__version__ = '$Id: plotStats.py 762 2011-11-25 05:02:07Z nsummons $'



def _linreg(data):
    """
    Calculate the linear regression of the data against itself (lag-1)
    Returns the slope, intercept, correlation, two-tailed
    probability and standard error of the estimate.
    """
    tData = numpy.array([data[1:],data[:-1]])
    i = numpy.where((tData[0,:]<sys.maxint) & (tData[1,:]<sys.maxint))[0]
    m,c,r,pr,err = linregress(tData[:,i])
    return m, c, r, pr, err

def _quantiles(data,dist='normal',mean=0.0,sigma=1.0):
    """
    Calculate normalised anomalies and associated quantile estimates
    """
    d = data.compress(data<sys.maxint)
    m = numpy.average(d)
    sd = numpy.std(d)
    nd = (d-m)/sd
    qd = [percentile(nd,i) for i in range(1,100)]
    f = getattr(numpy.random,dist)
    nf = f(mean,sigma,1000)
    qn = [percentile(nf,i) for i in range(1,100)]
    return qd, qn

def _plotHistogram(data,xlims,binInt,outputImg,xlabel):
    """
    Plot a simple histogram of the data values
    """
    pyplot.rcdefaults()
    fig = pyplot.figure()
    fig.clf()
    ax1 = fig.add_subplot(111)
    d = stats.statRemoveNum(numpy.array(data), sys.maxint)
    bins = numpy.arange(xlims[0],xlims[1]+1,binInt)
    #pdf,bins = numpy.histogram(d,bins,normed=True)
    #cdf = numpy.cumsum(pdf)
    n,b,p = pyplot.hist(d,bins,align='left',normed=True)
    ax1.set_ylabel('Probability')
    #ax2 = pyplot.twinx()
    #pyplot.plot(b[:-1],cdf,'k-')
    #ax2.set_ylabel('Cumulative probability')

    ax1.set_xlabel(xlabel)
    pyplot.savefig(outputImg)
    pyplot.rcdefaults()



def plotPressure(pAllData,pRateData,outputPath):
    """
    Plot the input pressure values lagged against themselves,
    and the same for the changes in pressure.
    """
    pyplot.rcdefaults()

    pyplot.rcParams['figure.figsize'] = (7,12)
    pyplot.figure(1)
    pyplot.subplot(211)
    pyplot.plot(pAllData[1:],pAllData[:-1],'k.',markersize=1)
    m,c,r,p,e = _linreg(pAllData)
    
    x = numpy.arange(880.,1021.,1.)
    y = m*x+c
    #pyplot.plot(x,y,'r-')
    pyplot.plot(x,x,'k-')
    pyplot.text(900,1010, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.xlim(880.,1020.)
    pyplot.ylim(880.,1020.)
    pyplot.xticks(numpy.arange(880.,1021.,20.))
    pyplot.yticks(numpy.arange(880.,1021.,20.))
    pyplot.ylabel(r"$p (t)$",fontsize=16)
    pyplot.xlabel(r"$p (t-1)$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Pressure")

    pyplot.subplot(212)

    pyplot.plot(pRateData[1:],pRateData[:-1],'k.',markersize=1)
    m,c,r,p,e = _linreg(pRateData)
    pyplot.text(-8.,8., "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.xlim(-10.,10.)
    pyplot.ylim(-10.,10.)
    pyplot.xticks(numpy.arange(-10.,11.,2.5))
    pyplot.yticks(numpy.arange(-10.,11.,2.5))
    pyplot.ylabel(r"$\partial p/\partial t (t)$",fontsize=16)
    pyplot.xlabel(r"$\partial p/\partial t (t-1)$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Pressure rate of change")

    pyplot.savefig(os.path.join(outputPath,'prs_corr.png'))
    pyplot.savefig(os.path.join(outputPath,'prs_corr.eps'))

    #_plotHistogram(pAllData,(880.,1020.),5,os.path.join(outputPath,'pressure_hist.png'),'Pressure')
    #_plotHistogram(pAllData,(880.,1020.),5,os.path.join(outputPath,'pressure_hist.eps'),'Pressure')
    #_plotHistogram(pRateData,(-10.,10.),0.5,os.path.join(outputPath,'pressureRate_hist.png'),'Pressure rate of change (hPa/hr)')
    #_plotHistogram(pRateData,(-10.,10.),0.5,os.path.join(outputPath,'pressureRate_hist.eps'),'Pressure rate of change (hPa/hr)')

def plotBearing(bAllData,bRateData,outputPath):
    pyplot.rcParams['figure.figsize'] = (7,12)
    pyplot.figure(4)
    pyplot.subplot(211)
#    pyplot.plot(numpy.cos(bAllData[1:]),numpy.cos(bAllData[:-1]),'kx')
#    pyplot.xlim(-1.,1.)
#    pyplot.ylim(-1.,1.)
#    pyplot.xticks(numpy.arange(-1.,1.1,1.))
#    pyplot.yticks(numpy.arange(-1.,1.1,1.))
    bAllData_t0 = bAllData[1:]
    bAllData_tm1 = bAllData[:-1]
    bAllData_skip = (bAllData_t0>=sys.maxint) | (bAllData_tm1>=sys.maxint)
    bAllData_t0 = bAllData_t0.compress(bAllData_skip == False)
    bAllData_tm1 = bAllData_tm1.compress(bAllData_skip == False)
    pyplot.plot(numpy.cos(numpy.radians(bAllData_t0)),numpy.cos(numpy.radians(bAllData_tm1)),'k.',markersize=1)
    pyplot.xlim(-1.,1.)
    pyplot.ylim(-1.,1.)
    pyplot.xticks(numpy.arange(-1.,1.1,1.))
    pyplot.yticks(numpy.arange(-1.,1.1,1.))
    m,c,r,p,e = _linreg(numpy.cos(numpy.radians(bAllData)).compress(bAllData < sys.maxint))

    x = numpy.arange(-1.0,1.1,0.1)
    y = m*x+c
    #pyplot.plot(x,y,'r-')
    pyplot.plot(x,x,'k-')

    pyplot.text(-0.8,0.8, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.ylabel(r"$cos(\theta (t))$",fontsize=16)
    pyplot.xlabel(r"$cos(\theta (t-1))$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Bearing")

    pyplot.subplot(212)
    pyplot.plot(bRateData[1:],bRateData[:-1],'k.',markersize=1)
    m,c,r,p,e = _linreg(numpy.cos(numpy.radians(bRateData)))
    pyplot.text(-25,25, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)

    pyplot.xlim(-30.,30.)
    pyplot.ylim(-30.,30.)
    pyplot.xticks(numpy.arange(-30.,31.,10.))
    pyplot.yticks(numpy.arange(-30.,31.,10.))
    pyplot.ylabel(r"$\partial \theta /\partial t (t)$",fontsize=16)
    pyplot.xlabel(r"$\partial \theta /\partial t (t-1)$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Bearing rate of change")
    pyplot.savefig(os.path.join(outputPath,'bear_corr.png'))
    #pyplot.savefig(os.path.join(outputPath,'bear_corr.eps'))


def plotSpeed(sAllData,sRateData,outputPath):
    pyplot.rcParams['figure.figsize'] = (7,12)
    pyplot.figure(5)
    pyplot.subplot(211)
    pyplot.plot(sAllData[1:],sAllData[:-1],'k.',markersize=1)

    m,c,r,p,e = _linreg(sAllData)

    x = numpy.arange(0.,101.,1.)
    y = m*x+c
    #pyplot.plot(x,y,'r-')
    pyplot.plot(x,x,'k-')
    pyplot.text(10,90, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.xlim(0.,100.)
    pyplot.ylim(0.,100.)
    pyplot.xticks(numpy.arange(0,101.,20))
    pyplot.yticks(numpy.arange(0,101.,20))
    #pyplot.grid(True)
    pyplot.ylabel(r"$v (t)$",fontsize=16)
    pyplot.xlabel(r"$v (t-1)$",fontsize=16)

    pyplot.title("Speed")

    pyplot.subplot(212)
    pyplot.plot(sRateData[1:],sRateData[:-1],'k.',markersize=1)
    m,c,r,p,e = _linreg(sRateData)
    pyplot.text(-25,25, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.xlim(-30.,30.)
    pyplot.ylim(-30.,30.)
    pyplot.xticks(numpy.arange(-30.,31.,10.))
    pyplot.yticks(numpy.arange(-30.,31.,10.))

    pyplot.ylabel(r"$\partial v /\partial t (t)$",fontsize=16)
    pyplot.xlabel(r"$\partial v /\partial t (t-1)$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Speed rate of change")
    pyplot.savefig(os.path.join(outputPath,'spd_corr.png'))
    pyplot.savefig(os.path.join(outputPath,'spd_corr.eps'))
    pyplot.rcdefaults()

def plotSpeedBear(sAllData,bAllData,outputPath):
    pyplot.rcParams['figure.figsize'] = (7,7)
    pyplot.figure(6)
    pyplot.subplot(111)
    ii = numpy.where((sAllData<sys.maxint) & (bAllData<sys.maxint))
    pyplot.polar((numpy.pi/2. - numpy.radians(bAllData[ii])),sAllData[ii],'k.',markersize=2)
    thetalabels=(90-numpy.arange(0,360,45))
    ii = numpy.where(thetalabels<0)
    thetalabels[ii]+=360
    lines, labels = pyplot.rgrids(numpy.arange(20.,101.,20.), labels=None, angle=67.5)
    lines, labels = pyplot.thetagrids(numpy.arange(0.,360.,45.),thetalabels)
    pyplot.ylim(0,100.)
    pyplot.grid(True)
    r = numpy.corrcoef(bAllData[ii],sAllData[ii])[1,0]
    pyplot.text(45,125, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.title("Speed vs bearing")
    pyplot.savefig(os.path.join(outputPath,'spd_bear_corr.png'))
    pyplot.savefig(os.path.join(outputPath,'spd_bear_corr.eps'))
    pyplot.rcdefaults()


def plotFrequency(years, frequency, outputPath):
    pyplot.rcdefaults()
    pyplot.rcParams['figure.figsize'] = (14,5)
    pyplot.figure(7)
    pyplot.plot(years,frequency,'k-',linewidth=3)
    xmax = 5*int((1+years.max()/5))
    xmin = 5*int((years.min()/5))
    ymax = 5*int(1+frequency.max()/5)
    pyplot.xlim(xmin,xmax)
    pyplot.ylim(0.0,ymax)
    pyplot.xticks(numpy.arange(xmin,xmax,5))
    pyplot.yticks(numpy.arange(5,ymax,5))

    pyplot.xlabel("Year")
    pyplot.ylabel("Frequency")

    m,c,r,pr,err = linregress(numpy.array([years,frequency]))
    x = numpy.arange(xmin,xmax)
    y = m*x+c

    pyplot.plot(x,y,'r--')
    #pyplot.text(xmin+5,ymax-5, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)

    pyplot.xlim(xmin,xmax)
    pyplot.ylim(0.0,ymax)
    pyplot.xticks(numpy.arange(xmin,xmax,5))
    pyplot.yticks(numpy.arange(0,ymax,5))

    pyplot.grid(True)
    pyplot.title("Annual frequency (%d - %d)"%(years.min(),years.max()))
    pyplot.savefig(os.path.join(outputPath,'frequency.png'))
    pyplot.savefig(os.path.join(outputPath,'frequency.eps'))
    pyplot.rcdefaults()

def plotLonLat(lonData,latData,indicator,outputPath):
    """
    Short description:
    """
    pyplot.rcParams['figure.figsize'] = (7,12)
    pyplot.figure(8)

    dlon = lonData[1:]-lonData[:-1]
    dlat = latData[1:]-latData[:-1]
    j = numpy.where(indicator[1:]==0)
    dlon = dlon[j]
    dlat = dlat[j]

    pyplot.subplot(211)
    pyplot.plot(dlon[1:],dlon[:-1],'k.',markersize=1)
    m,c,r,p,e = _linreg(dlon)
    pyplot.text(-3,3, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.xlim(-4.,4.)
    pyplot.ylim(-4.,4.)
    pyplot.xticks(numpy.arange(-4.,4.1,1.))
    pyplot.yticks(numpy.arange(-4.,4.1,1.))
    pyplot.ylabel(r"$\Delta lon (t)$",fontsize=16)
    pyplot.xlabel(r"$\Delta lon (t-1)$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Longitude rate of change")

    pyplot.subplot(212)
    pyplot.plot(dlat[1:],dlat[:-1],'k.',markersize=1)
    m,c,r,p,e = _linreg(dlat)
    pyplot.text(-3,3, "r = %5.3f"%r,ha='center',va='center',color='r',size=14)
    pyplot.xlim(-4.,4.)
    pyplot.ylim(-4.,4.)
    pyplot.xticks(numpy.arange(-4.,4.1,1.))
    pyplot.yticks(numpy.arange(-4.,4.1,1.))
    pyplot.ylabel(r"$\Delta lat (t)$",fontsize=16)
    pyplot.xlabel(r"$\Delta lat (t-1)$",fontsize=16)
    #pyplot.grid(True)
    pyplot.title("Latitude rate of change")

    pyplot.savefig(os.path.join(outputPath,'lonlat_corr.png'))
    pyplot.savefig(os.path.join(outputPath,'lonlat_corr.eps'))
    pyplot.rcdefaults()



def quantile(data,outputPath,parameterName,dist='normal',mean=0.0,sigma=1.0):
    """
    Generate a probability plot of the given data; data should be an array of anomalies

    """
    pyplot.rcParams['figure.figsize'] = (8,7)
    pyplot.figure(9)
    pyplot.clf()
    d = data.compress(data<sys.maxint)
    m = numpy.average(d)
    sd = numpy.std(d)
    nd = (d-m)/sd
    (osm, osr), (slope, intercept, r) = probplot(nd,plot=pyplot)
    #pyplot.xticks(numpy.arange(-3.,3.1,1.))
    #pyplot.yticks(numpy.arange(-3.,3.1,1.))
    #pyplot.xlabel("Normal")
    pyplot.ylabel(parameterName)
    pyplot.title("Q-Q plot - %s"%parameterName)
    pyplot.xlim((-5,5))
    pyplot.ylim((-5,5))
    pos = 2, -4.8
    pyplot.text(2,-4.9,r"$r^2=%1.4f$" % r, fontsize=12)

    pyplot.savefig(os.path.join(outputPath,'qqplot_%s.png'%parameterName))
    pyplot.savefig(os.path.join(outputPath,'qqplot_%s.eps'%parameterName))
    pyplot.rcdefaults()

def minPressureHist(index,pAllData,outputPath):
    """ Plot a histogram of the minimum central pressures from the input
        dataset.
    """
    pyplot.rcParams['figure.figsize'] = (8,7)
    pyplot.figure(10)
    pyplot.clf()
    pcarray = []
    index = index.astype(int)
    for i in range(len(index)-1):
        if index[i] == 1:
            pcarray.append(pAllData[i])
        else:
            if pAllData[i] is not None:
                if pAllData[i]<pcarray[-1]:
                    pcarray[-1] = pAllData[i]

    pbins = numpy.arange(850.,1020.,5)
    n,b,p = pyplot.hist(numpy.array(pcarray),pbins,normed=False,lw=2,ec='k',fc='w')
    pyplot.xlabel("Minimum central pressure (hPa)")
    pyplot.ylabel("Count")
    pyplot.title("Distribution of minimum central pressure")
    pyplot.savefig(os.path.join(outputPath,"min_pressure_hist.png"))
    pyplot.savefig(os.path.join(outputPath,"min_pressure_hist.eps"))
    x = numpy.zeros((len(b),2))
    x[:,0] = b
    x[1:,1] = n
    files.flSaveFile(os.path.join(outputPath,'min_pressure_hist.csv'),x,delimiter=',',fmt='%6.2f')
    pyplot.rcdefaults()

def minPressureLat(pAllData,latData,outputPath,latMin=-40.,latMax=0.):
    """
    Plot the minimum central pressures as a function of latitude
    """
    rLat = numpy.round(latData,0)
    lats = numpy.arange(latMin,latMax+0.1,1)
    minP = numpy.zeros(len(lats))
    n = 0
    for l in lats:
        i = numpy.where(rLat==l)[0]
        if len(i>0):
            pvals = pAllData[i]
            pvals = stats.statRemoveNum(pvals,0)
            if len(pvals)>0:
                minP[n] = pvals.min()
            else:
                minP[n] = 1020.
        else:
            minP[n] = 1020.
        n += 1
    pyplot.figure(11)
    pyplot.plot(lats,minP,'r-',linewidth=2,label=r'Min $P_{centre}$')
    pyplot.xlim(latMin,latMax)
    pyplot.ylim(800,1020)

    pyplot.xlabel('Latitude',fontsize=10)
    pyplot.ylabel('Minimum central pressure (hPa)',fontsize=10)
    pyplot.legend(loc=3)
    pyplot.grid(True)

    pyplot.savefig(os.path.join(outputPath,"min_pressure_lat.png"))
    pyplot.savefig(os.path.join(outputPath,"min_pressure_lat.eps"))
    x = numpy.zeros((len(lats),2))
    x[:,0] = lats
    x[:,1] = minP
    files.flSaveFile(os.path.join(outputPath,'min_pressure_lat.csv'),x,delimiter=',',fmt='%6.2f')

def julianDay(julianDayObs,julianDayGenesis,outputPath):
    pyplot.rcParams['figure.figsize'] = (10,6)
    pyplot.figure(12)
    pyplot.clf()
    pyplot.bar(julianDayObs[:,0],julianDayObs[:,1])
    pyplot.xlim((1,365))
    pyplot.xlabel('Day of year',fontsize=10)
    pyplot.ylabel('Number of observations',fontsize=10)
    pyplot.savefig(os.path.join(outputPath,"julian_day_observations.png"))
    pyplot.clf()
    pyplot.bar(julianDayGenesis[:,0],julianDayGenesis[:,1])
    pyplot.xlim((1,365))
    pyplot.xlabel('Day of year',fontsize=10)

    pyplot.ylabel('Number of genesis events',fontsize=10)
    pyplot.savefig(os.path.join(outputPath,"julian_day_genesis.png"))

if __name__=="__main__":
    configFile = sys.argv[1]
    dataPath = config.cnfGetIniValue(configFile, 'Output','Path',os.getcwd())
    outputPath = os.path.join(dataPath,'plots')
    pyplot.rcParams['figure.figsize'] = (7,12)

    pRateData = files.flLoadFile(os.path.join(dataPath,'pressure_rate'))
    pAllData = files.flLoadFile(os.path.join(dataPath,'all_pressure'))
    bRateData = files.flLoadFile(os.path.join(dataPath,'bearing_rate'))
    bAllData = files.flLoadFile(os.path.join(dataPath,'all_bearing'))
    sRateData = files.flLoadFile(os.path.join(dataPath,'speed_rate'))
    sAllData = files.flLoadFile(os.path.join(dataPath,'all_speed'))
    freq = files.flLoadFile(os.path.join(dataPath,'frequency'))
    years = freq[:,0]
    frequency = freq[:,1]

    plotPressure(pAllData,pRateData,outputPath)
    plotBearing(bAllData,bRateData,outputPath)
    plotSpeed(sAllData,sRateData,outputPath)
    plotFrequency(years,frequency,outputPath)
    plotSpeedBear(sAllData,bAllData,outputPath)
