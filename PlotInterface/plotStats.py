import os 
import logging
import sys

from os.path import join as pjoin
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import linregress, probplot
from scipy.stats import scoreatpercentile as percentile
import numpy as np

import Utilities.files as files
import Utilities.config as config
import Utilities.nctools as nctools
import Utilities.stats as stats

import seaborn as sns
sns.set(style="ticks")

log = logging.getLogger(__name__)
pyplot.ioff()

def get_bins(data, allpos=False):
    """
    Establish nice bin widths for plotting histograms - as a default it sets
    the width to 0.25 times the standard deviation of the data
    """
    binwidth = 0.25*np.std( data )
    dmax = np.max(np.fabs( data ))

    lims = (int(dmax/binwidth) + 1) * binwidth
    if allpos:

        bins = np.arange(0, lims+binwidth, binwidth)
    else:
        bins = np.arange(-lims, lims+binwidth, binwidth)

    return bins


def _linreg(data):
    """
    Calculate the linear regression of the data against itself (lag-1)
    Returns the slope, intercept, correlation, two-tailed
    probability and standard error of the estimate.
    """
    tData = np.array([data[1:], data[:-1]])
    i = np.where((tData[0, :] < sys.maxint) & (tData[1, :] < sys.maxint))[0]
    m, c, r, pr, err = linregress(tData[:, i])
    return m, c, r, pr, err


class PlotData(object):
    """
    Base class for plotting summaries of input data.

    """
    def __init__(self, output_path, output_format, context="poster"):
        """
        Initialise statistical plotting routines, fixing the output path and
        the format of all images generated
        """
        self.outpath = output_path
        self.fmt = output_format
        self.fignum = 0
        sns.set_context(context)

    def figurenum(self):
        """
        Increment the figure number and return the new value
        """
        self.fignum += 1
        return self.fignum

    def savefig(self, filename, **kwargs):
        """
        Provide a shortcut to pyplot.savefig() that adds the output path to the
        filename.
        All keyword args are passed without alteration
        """
        outputfile = ".".join([filename, self.fmt])
        pyplot.savefig(pjoin(self.outpath, outputfile),
                       format=self.fmt, **kwargs)
        pyplot.close()

    def scatterHistogram(self, x, y, labels, img_name, allpos=False):
        """
        Create a scatter plot with marginal distributions
        
        :param x: `numpy.ndarray` of x values.
        :param y: `numpy.ndarray` of y values.
        :param list labels: A length-2 list with the x and y labels as strings.
        :param str img_name: Name to use in the file name for saving the 
                             figure.
        
        """

        i = np.where((x < sys.maxint) & (y < sys.maxint))[0]
        x = x[i]
        y = y[i]
        jp = sns.jointplot(x, y, kind='reg',
                           joint_kws={'scatter_kws':
                                      {'color':'slategray', 
                                       'alpha':0.5}})
        xlabel, ylabel = labels
        jp.set_axis_labels(xlabel, ylabel)
        pyplot.tight_layout()
        self.savefig(img_name)

    def plotPressure(self, pAllData, pRateData):
        """
        Plot the input pressure values lagged against themselves,
        and the same for the changes in pressure.
        """

        fig, (ax1, ax2) = pyplot.subplots(2,1,figsize=(7,12))

        i = np.where((pAllData < sys.maxint) & (pAllData != 0.))[0]
        j = np.where((np.abs(pRateData) < 100.))[0]

        sns.regplot(pAllData[i][1:], pAllData[i][:-1],
                    fit_reg=True, ax=ax1, dropna=True)
        sns.regplot(pRateData[j][1:], pRateData[j][:-1], 
                    fit_reg=True, ax=ax2, dropna=True)

        m, c, r, p, e = _linreg(pAllData)

        ax1.text(900, 1010, "r = %5.3f"%r, ha='center', va='center')
        ax1.set_xlim(880., 1020.)
        ax1.set_ylim(880., 1020.)
        ax1.set_ylabel(r"$p (t)$")
        ax1.set_xlabel(r"$p (t-1)$")
        ax1.set_title("Pressure")

        m, c, r, p, e = _linreg(pRateData)
        ax2.text(-8., 8., "r = %5.3f"%r, ha='center', va='center')
        ax2.set_xlim(-10., 10.)
        ax2.set_ylim(-10., 10.)
        ax2.set_ylabel(r"$\Delta p/\Delta t (t)$") 
        ax2.set_xlabel(r"$\Delta p/\Delta t (t-1)$")
        ax2.set_title("Pressure rate of change")
        fig.tight_layout()

        self.savefig('prs_corr')


    def plotSpeed(self, v, dv):
        """
        Plot the input speed values lagged against themselves,
        and the same for the changes in speed.
        """

        fig, (ax1, ax2) = pyplot.subplots(2,1,figsize=(7,12))

        i = np.where((v < sys.maxint))[0]
        j = np.where((dv < sys.maxint))[0]
        
        sns.regplot(v[i][:-1], v[i][1:], fit_reg=True, ax=ax1, dropna=True)
        sns.regplot(dv[j][:-1], dv[j][1:], fit_reg=True, ax=ax2, dropna=True)

        m, c, r, p, e = _linreg(v)

        ax1.text(10, 90, "r = %5.3f"%r, ha='center', va='center')

        ax1.set_xlim(0., 100.)
        ax1.set_ylim(0., 100.)
        ax1.set_ylabel(r"$v (t)$")
        ax1.set_xlabel(r"$v (t-1)$")
        ax1.set_title("Speed")


        m, c, r, p, e = _linreg(dv)
        ax2.text(-25, 25, "r = %5.3f"%r, ha='center', va='center')
        ax2.set_xlim(-30.,30.)
        ax2.set_ylim(-30.,30.)

        ax2.set_ylabel(r"$\partial v /\partial t (t)$")
        ax2.set_xlabel(r"$\partial v /\partial t (t-1)$")

        ax2.set_title("Speed rate of change")
        fig.tight_layout()
        self.savefig("spd_corr")

        labels = [r'$v(t)$', r'$v(t-1)$']
        self.scatterHistogram(v[1:], v[:-1], labels,
                              'speed_scatterHist', allpos=True)
        labels = [r'$\frac{\Delta v}{\Delta t}(t)$', 
                  r'$\frac{\Delta v}{\Delta t}(t-1)$']
        self.scatterHistogram(dv[1:], dv[:-1], labels, 
                              'speedRate_scatterHist', allpos=True)

    def plotBearing(self, bAllData, bRateData):
        """
        Plot the (cosine of) input bearing values lagged against themselves,
        and the same for the changes in bearing.
        """

        pyplot.figure(self.figurenum(), figsize=(7, 12))
        pyplot.subplot(211)
    #    pyplot.plot(np.cos(bAllData[1:]),np.cos(bAllData[:-1]),'kx')
    #    pyplot.xlim(-1.,1.)
    #    pyplot.ylim(-1.,1.)
    #    pyplot.xticks(np.arange(-1.,1.1,1.))
    #    pyplot.yticks(np.arange(-1.,1.1,1.))
        bAllData_t0 = bAllData[1:]
        bAllData_tm1 = bAllData[:-1]
        bAllData_skip = (bAllData_t0>=sys.maxint) | (bAllData_tm1>=sys.maxint)
        bAllData_t0 = bAllData_t0.compress(bAllData_skip == False)
        bAllData_tm1 = bAllData_tm1.compress(bAllData_skip == False)
        pyplot.plot(np.cos(np.radians(bAllData_t0)),
                    np.cos(np.radians(bAllData_tm1)),
                    'k.', markersize=1)
        pyplot.xlim(-1., 1.)
        pyplot.ylim(-1., 1.)
        pyplot.xticks(np.arange(-1.,1.1,1.))
        pyplot.yticks(np.arange(-1.,1.1,1.))
        m, c, r, p, e = _linreg(np.cos(np.radians(bAllData)). \
                                compress(bAllData < sys.maxint))

        x = np.arange(-1.0,1.1,0.1)
        y = m*x+c
        pyplot.plot(x,y,'r-')
        pyplot.plot(x,x,'k-')

        pyplot.text(-0.8, 0.8, "r = %5.3f"%r,
                    ha='center', va='center', color='r', size=14)
        pyplot.ylabel(r"$cos(\theta (t))$", fontsize=16)
        pyplot.xlabel(r"$cos(\theta (t-1))$", fontsize=16)
        #pyplot.grid(True)
        pyplot.title("Bearing")

        pyplot.subplot(212)
        pyplot.plot(bRateData[1:],bRateData[:-1],'k.',markersize=1)
        m,c,r,p,e = _linreg(np.cos(np.radians(bRateData)))
        pyplot.text(-25, 25, "r = %5.3f"%r,
                    ha='center', va='center', color='r', size=14)

        pyplot.xlim(-30., 30.)
        pyplot.ylim(-30., 30.)
        pyplot.xticks(np.arange(-30., 31., 10.))
        pyplot.yticks(np.arange(-30., 31., 10.))
        pyplot.ylabel(r"$\partial \theta /\partial t (t)$", fontsize=16)
        pyplot.xlabel(r"$\partial \theta /\partial t (t-1)$", fontsize=16)
        #pyplot.grid(True)
        pyplot.title("Bearing rate of change")
        self.savefig('bear_corr')
        #pyplot.savefig(os.path.join(outputPath,'bear_corr.eps'))
        labels = [r'$\theta(t)$', 
                  r'$\theta(t-1)$']
        self.scatterHistogram(bAllData[1:], bAllData[:-1], labels, 
                              'bear_scatterHist', allpos=True)
        labels = [r'$\frac{\Delta \theta}{\Delta t}(t)$', 
                  r'$\frac{\Delta \theta}{\Delta t}(t-1)$']
        self.scatterHistogram(bRateData[1:], bRateData[:-1], labels,
                              'bearRate_scatterHist', allpos=True)

    def julianDay(self, julianDayObs, julianDayGenesis):
        """
        Plot bar graphs of the number of TC observations and genesis events
        """

        pyplot.figure(self.figurenum(), figsize=(10, 6))
        pyplot.clf()
        pyplot.bar(julianDayObs[:, 0], julianDayObs[:, 1])
        pyplot.xlim((1, 365))
        pyplot.xlabel('Day of year', fontsize=10)
        pyplot.ylabel('Number of observations', fontsize=10)
        self.savefig("julian_day_obs")

        pyplot.figure(self.figurenum())
        pyplot.clf()
        pyplot.bar(julianDayGenesis[:, 0], julianDayGenesis[:, 1])
        pyplot.xlim((1, 365))
        pyplot.xlabel('Day of year', fontsize=10)

        pyplot.ylabel('Number of genesis events', fontsize=10)
        self.savefig("julian_day_genesis")

    def plotLonLat(self, lonData, latData, indicator):
        """
        Plot the input lat/lon values lagged against themselves,
        and the same for the changes in lat/lon.
        """

        pyplot.figure(self.figurenum(), figsize=(7, 12))

        dlon = lonData[1:] - lonData[:-1]
        dlat = latData[1:] - latData[:-1]
        j = np.where(indicator[1:] == 0)
        dlon = dlon[j]
        dlat = dlat[j]

        # Correct change in longitude where the value jumps across the 180E
        # meridian
        k = np.where(dlon<-180.)
        dlon[k] += 360.

        pyplot.subplot(211)
        pyplot.plot(dlon[1:], dlon[:-1], 'k.', markersize=1)
        m, c, r, p, e = _linreg(dlon)
        pyplot.text(-3, 3, "r = %5.3f"%r, ha='center',
                    va='center', color='r', size=14)
        pyplot.xlim(-4., 4.)
        pyplot.ylim(-4., 4.)
        pyplot.xticks(np.arange(-4., 4.1, 1.))
        pyplot.yticks(np.arange(-4., 4.1, 1.))
        pyplot.ylabel(r"$\Delta lon (t)$", fontsize=16)
        pyplot.xlabel(r"$\Delta lon (t-1)$", fontsize=16)
        #pyplot.grid(True)
        pyplot.title("Longitude rate of change")

        pyplot.subplot(212)
        pyplot.plot(dlat[1:], dlat[:-1], 'k.', markersize=1)
        m, c, r, p, e = _linreg(dlat)
        pyplot.text(-3, 3, "r = %5.3f"%r, ha='center',
                    va='center', color='r', size=14)
        pyplot.xlim(-4., 4.)
        pyplot.ylim(-4., 4.)
        pyplot.xticks(np.arange(-4., 4.1, 1.))
        pyplot.yticks(np.arange(-4., 4.1, 1.))
        pyplot.ylabel(r"$\Delta lat (t)$", fontsize=16)
        pyplot.xlabel(r"$\Delta lat (t-1)$", fontsize=16)
        pyplot.title("Latitude rate of change")

        self.savefig('lonlat_corr')

        labels = [r'$\Delta \phi(t)$', r'$\Delta \phi(t-1)$']
        self.scatterHistogram(dlon[1:], dlon[:-1], labels, 'dlon_scatterHist')
        labels = [r'$\Delta \lambda(t)$', r'$\Delta \lambda(t-1)$']
        self.scatterHistogram(dlat[1:], dlat[:-1], labels, 'dlat_scatterHist')

    def plotFrequency(self, years, frequency):
        """Plot annual frequency of TCs, plus linear trend line"""

        pyplot.figure(self.figurenum(), figsize=(14, 5))
        pyplot.plot(years, frequency) 
        xmax = 5*int((1 + years.max()/5))
        xmin = 5*int((years.min()/5))
        ymax = 5*int(1 + frequency.max()/5)
        pyplot.xlim(xmin, xmax)
        pyplot.ylim(0.0, ymax)
        pyplot.xticks(np.arange(xmin, xmax, 5))
        pyplot.yticks(np.arange(5, ymax, 5))

        pyplot.xlabel("Year")
        pyplot.ylabel("Frequency")

        m, c, r, pr, err = linregress(np.array([years, frequency]))
        x = np.arange(xmin,xmax)
        y = m*x + c

        pyplot.plot(x, y, 'r--', linewidth=0.5)

        pyplot.xlim(xmin, xmax)
        pyplot.ylim(0.0, ymax)
        #pyplot.xticks(np.arange(xmin, xmax, 5))
        pyplot.yticks(np.arange(0, ymax, 5))

        pyplot.title("Annual frequency (%d - %d)"%(years.min(), years.max()))
        self.savefig('frequency')

        return

    def minPressureLat(self, pAllData, latData, latMin=-40., latMax=0.):
        """
        Plot the minimum central pressures as a function of latitude
        """
        rLat = np.round(latData, 0)
        lats = np.arange(latMin, latMax + 0.1, 1)
        minP = np.zeros(len(lats))
        n = 0
        for l in lats:
            i = np.where(rLat == l)[0]
            if len(i > 0):
                pvals = pAllData[i]
                pvals = stats.statRemoveNum(pvals, 0)
                if len(pvals)>0:
                    minP[n] = pvals.min()
                else:
                    minP[n] = 1020.
            else:
                minP[n] = 1020.
            n += 1
        pyplot.figure(self.figurenum())
        pyplot.plot(lats, minP, label=r'Min $P_{centre}$')
        pyplot.xlim(latMin, latMax)
        pyplot.ylim(800, 1020)

        pyplot.xlabel('Latitude')
        pyplot.ylabel('Minimum central pressure (hPa)')
        pyplot.legend(loc=3)
        pyplot.grid(True)

        self.savefig("min_pressure_lat")

        x = np.zeros((len(lats), 2))
        x[:, 0] = lats
        x[:, 1] = minP
        files.flSaveFile(os.path.join(self.outpath, 'min_pressure_lat.csv'), x,
                         delimiter=',', fmt='%6.2f')

    def minPressureHist(self, index, pAllData):
        """
        Plot a histogram of the minimum central pressures from the input
        dataset.
        """

        pyplot.figure(self.figurenum(), figsize=(8, 7))
        pyplot.clf()
        pcarray = []
        index = index.astype(int)
        for i in range(len(index) - 1):
            if index[i] == 1:
                pcarray.append(pAllData[i])
            else:
                if pAllData[i] is not None:
                    if pAllData[i] < pcarray[-1]:
                        pcarray[-1] = pAllData[i]

        pbins = np.arange(850., 1020., 5)
        n, b, p = pyplot.hist(np.array(pcarray), pbins, normed=False,
                              lw=2, ec='k', fc='w')
        pyplot.xlabel("Minimum central pressure (hPa)")
        pyplot.ylabel("Count")
        pyplot.title("Distribution of minimum central pressure")
        self.savefig("min_pressure_hist")

        x = np.zeros((len(b), 2))
        x[:, 0] = b
        x[1:, 1] = n
        files.flSaveFile(os.path.join(self.outpath, 'min_pressure_hist.csv'), x,
                         delimiter=',', fmt='%6.2f')

    def plotSpeedBear(self, sAllData, bAllData):
        """
        Plot speed and bearing against each other
        """
        pyplot.figure(self.figurenum(), figsize=(7, 7))
        pyplot.subplot(111)
        ii = np.where((sAllData < sys.maxint) & (bAllData < sys.maxint))
        pyplot.polar((np.pi/2. - np.radians(bAllData[ii])), sAllData[ii],
                     'k.', markersize=2)
        thetalabels=(90 - np.arange(0, 360, 45))
        ii = np.where(thetalabels < 0)
        thetalabels[ii] += 360
        lines, labels = pyplot.rgrids(np.arange(20., 101., 20.),
                                      labels=None, angle=67.5)
        lines, labels = pyplot.thetagrids(np.arange(0., 360., 45.),
                                          thetalabels)
        pyplot.ylim(0, 100.)
        pyplot.grid(True)
        r = np.corrcoef(bAllData[ii], sAllData[ii])[1, 0]
        pyplot.text(45, 125, "r = %5.3f"%r, ha='center',
                    va='center', color='r', size=14)
        pyplot.title("Speed vs bearing")
        self.savefig('spd_bear_corr')

    def quantile(self, data, parameterName, dist='normal', mean=0.0, sigma=1.0):
        """
        Generate a probability plot of the given data; data should be an
        array of anomalies

        """
        pyplot.figure(self.figurenum(), figsize=(8, 7))
        pyplot.clf()
        d = data.compress(data < sys.maxint)
        m = np.average(d)
        sd = np.std(d)
        nd = (d-m)/sd
        (osm, osr), (slope, intercept, r) = probplot(nd, dist=dist, plot=pyplot)
        #pyplot.xticks(np.arange(-3.,3.1,1.))
        #pyplot.yticks(np.arange(-3.,3.1,1.))
        #pyplot.xlabel("Normal")
        pyplot.ylabel(parameterName)
        pyplot.title("Q-Q plot - %s" % parameterName)
        pyplot.xlim((-5, 5))
        pyplot.ylim((-5, 5))
        pos = (2, -4.8)
        pyplot.text(2, -4.9, r"$r^2=%1.4f$" % r, fontsize=12)

        self.savefig('qqplot_%s' % parameterName)

if __name__=="__main__":
    from os.path import join as pjoin

    configFile = sys.argv[1]
    dataPath = config.cnfGetIniValue(configFile, 'Output', 'Path', os.getcwd())
    inputPath = pjoin(dataPath, 'process')
    outputPath = pjoin(dataPath, 'plots')
    pyplot.rcParams['figure.figsize'] = (7, 12)

    pRateData = files.flLoadFile(pjoin(inputPath, 'pressure_rate'))
    pAllData = files.flLoadFile(pjoin(inputPath, 'all_pressure'))
    bRateData = files.flLoadFile(pjoin(inputPath, 'bearing_rate'))
    bAllData = files.flLoadFile(pjoin(inputPath, 'all_bearing'))
    sRateData = files.flLoadFile(pjoin(inputPath, 'speed_rate'))
    sAllData = files.flLoadFile(pjoin(inputPath, 'all_speed'))
    freq = files.flLoadFile(pjoin(inputPath, 'frequency'))
    years = freq[:, 0]
    frequency = freq[:, 1]

    plotting = PlotData(outputPath, "png")
    plotting.plotPressure(pAllData, pRateData)
    plotting.plotBearing(bAllData, bRateData)
    plotting.plotSpeed(sAllData, sRateData)
    plotting.plotFrequency(years, frequency)
    plotting.plotSpeedBear(sAllData, bAllData)
