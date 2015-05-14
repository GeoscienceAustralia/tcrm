import os 
import logging
import sys

from os.path import join as pjoin
from matplotlib import pyplot

from scipy.stats import linregress, probplot, frechet_l
import numpy as np

import Utilities.files as files
import Utilities.config as config
import Utilities.stats as stats

import seaborn as sns
sns.set_style("ticks")

log = logging.getLogger(__name__)
pyplot.ioff()

def linreg(data):
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

    def scatterHistogram(self, x, y, labels, name):
        """
        Create a scatter plot with marginal distributions
        
        :param x: `numpy.ndarray` of x values.
        :param y: `numpy.ndarray` of y values.
        :param list labels: A length-2 list with the x and y labels as strings.
        :param str img_name: Name to use in the file name for saving the 
                             figure.
        
        """

        i = np.where((x < sys.maxint) & (y < sys.maxint))[0]
        xx = x[i]
        yy = y[i]
        jp = sns.jointplot(xx, yy, kind='reg',
                           joint_kws={'scatter_kws':
                                      {'color':'slategray', 
                                       'alpha':0.5}})

        jp.set_axis_labels(*labels)
        pyplot.tight_layout()

        self.savefig(name)

    def plotRegression(self, x, name, labels, transform=lambda x: x):
        """
        A generic function to plot a lag-1 autoregression of the 
        variable, with a joint probability plot including 
        marginal distribution plots.

        :param x: `numpy.ndarray` of the data variable to plot.
        :param str name: Name of the variable. Used for saving the 
                         resulting image.
        :param list labels: A list of the x- and y-labels to apply.
        :param func transform: A transform function to apply to 
                               the data. Default is to leave the 
                               data unchanged.
        
        """
        
        fig, ax = pyplot.subplots(1, 1)
        x_t = x[1:]
        x_tm = x[:-1]
        skip = (x_t >= sys.maxint) | (x_tm >= sys.maxint)
        x_t = x_t.compress(skip==False)
        x_tm = x_tm.compress(skip==False)
        
        x_t = transform(x_t)
        x_tm = transform(x_tm)
        
        sns.regplot(x_t, x_tm, fit_reg=True, ax=ax, dropna=True)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        fig.tight_layout()

        self.savefig(name)

        self.scatterHistogram(x_t, x_tm, labels, name+'_scatter')

    def plotPressure(self, pAllData, pRateData):
        """
        Plot the input pressure values lagged against themselves,
        and the same for the changes in pressure.
        """

        fig, (ax1, ax2) = pyplot.subplots(2, 1, figsize=(7, 12))

        i = np.where((pAllData < sys.maxint) & (pAllData != 0.))[0]
        j = np.where((np.abs(pRateData) < 100.))[0]

        sns.regplot(pAllData[i][1:], pAllData[i][:-1],
                    fit_reg=True, ax=ax1, dropna=True)
        sns.regplot(pRateData[j][1:], pRateData[j][:-1], 
                    fit_reg=True, ax=ax2, dropna=True)

        params = linreg(pAllData)

        ax1.text(900, 1010, "r = %5.3f"%params[2], ha='center', va='center')
        ax1.set_xlim(880., 1020.)
        ax1.set_ylim(880., 1020.)
        ax1.set_ylabel(r"$p (t)$")
        ax1.set_xlabel(r"$p (t-1)$")
        ax1.set_title("Pressure")

        params = linreg(pRateData)
        ax2.text(-8., 8., "r = %5.3f"%params[2], ha='center', va='center')
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

        fig, (ax1, ax2) = pyplot.subplots(2, 1, figsize=(7, 12))

        i = np.where((v < sys.maxint))[0]
        j = np.where((dv < sys.maxint))[0]
        
        sns.regplot(v[i][:-1], v[i][1:], 
                    fit_reg=True, ax=ax1, dropna=True)
        sns.regplot(dv[j][:-1], dv[j][1:], 
                    fit_reg=True, ax=ax2, dropna=True)

        params = linreg(v)

        ax1.text(10, 90, "r = %5.3f"%params[2], ha='center', va='center')

        ax1.set_xlim(0., 100.)
        ax1.set_ylim(0., 100.)
        ax1.set_ylabel(r"$v (t)$")
        ax1.set_xlabel(r"$v (t-1)$")
        ax1.set_title("Speed")


        params = linreg(dv)
        ax2.text(-25, 25, "r = %5.3f"%params[2], ha='center', va='center')
        ax2.set_xlim(-30., 30.)
        ax2.set_ylim(-30., 30.)

        ax2.set_ylabel(r"$\partial v /\partial t (t)$")
        ax2.set_xlabel(r"$\partial v /\partial t (t-1)$")

        ax2.set_title("Speed rate of change")
        fig.tight_layout()
        self.savefig("spd_corr")

        labels = [r'$v(t)$', r'$v(t-1)$']
        self.scatterHistogram(v[1:], v[:-1], labels,
                              'speed_scatterHist')
        labels = [r'$\frac{\Delta v}{\Delta t}(t)$', 
                  r'$\frac{\Delta v}{\Delta t}(t-1)$']
        self.scatterHistogram(dv[1:], dv[:-1], labels, 
                              'speedRate_scatterHist')

        
    def plotBearing(self, bAllData, bRateData):
        """
        Plot the (cosine of) input bearing values lagged against themselves,
        and the same for the changes in bearing.
        """

        fig, (ax1, ax2) = pyplot.subplots(2, 1, figsize=(7, 12))

        bAllData_t0 = bAllData[1:]
        bAllData_tm1 = bAllData[:-1]
        bAllData_skip = (bAllData_t0>=sys.maxint) | (bAllData_tm1>=sys.maxint)
        bAllData_t0 = bAllData_t0.compress(bAllData_skip == False)
        bAllData_tm1 = bAllData_tm1.compress(bAllData_skip == False)

        sns.regplot(np.cos(np.radians(bAllData_t0)),
                    np.cos(np.radians(bAllData_tm1)), 
                    fit_reg=True, ax=ax1, dropna=True)
                    
        ax1.set_xlim(-1., 1.)
        ax1.set_ylim(-1., 1.)
        ax1.set_xticks(np.arange(-1., 1.1, 1.))
        ax1.set_yticks(np.arange(-1., 1.1, 1.))
        params = linreg(np.cos(np.radians(bAllData)). \
                               compress(bAllData < sys.maxint))

        x = np.arange(-1.0, 1.1, 0.1)
        y = params[0]*x + params[1]
        pyplot.plot(x, y, 'r-')
        pyplot.plot(x, x, 'k-')

        ax1.text(-0.8, 0.8, "r = %5.3f"%params[2],
                    ha='center', va='center', color='r')
        ax1.set_ylabel(r"$cos(\theta (t))$")
        ax1.set_xlabel(r"$cos(\theta (t-1))$")
        #pyplot.grid(True)
        ax1.set_title("Bearing")

        sns.regplot(bRateData[1:], bRateData[:-1],
                    fit_reg=True, ax=ax2, dropna=True)
        params = linreg(np.cos(np.radians(bRateData)))
        ax2.text(-25, 25, "r = %5.3f"%params[2],
                    ha='center', va='center', color='r')

        ax2.set_xlim(-30., 30.)
        ax2.set_ylim(-30., 30.)
        ax2.set_xticks(np.arange(-30., 31., 10.))
        ax2.set_yticks(np.arange(-30., 31., 10.))
        ax2.set_ylabel(r"$\partial \theta /\partial t (t)$")
        ax2.set_xlabel(r"$\partial \theta /\partial t (t-1)$")
        ax2.set_title("Bearing rate of change")
        self.savefig('bear_corr')

        labels = [r'$cos(\theta(t))$', 
                  r'$cos(\theta(t-1))$']
        self.scatterHistogram(np.cos(np.radians(bAllData[1:])), 
                              np.cos(np.radians(bAllData[:-1])), 
                              labels, 'bear_scatterHist')
        labels = [r'$\frac{\Delta \theta}{\Delta t}(t)$', 
                  r'$\frac{\Delta \theta}{\Delta t}(t-1)$']
        self.scatterHistogram(bRateData[1:], bRateData[:-1], labels,
                              'bearRate_scatterHist')

    def julianDay(self, julianDayObs, julianDayGenesis):
        """
        Plot bar graphs of the number of TC observations and genesis events
        """

        pyplot.figure(self.figurenum())
        pyplot.clf()
        ax = sns.barplot(julianDayObs[:, 0], julianDayObs[:, 1])
        ax.set_xlim((1, 365))
        ax.set_xlabel('Day of year')
        ax.set_ylabel('Number of observations')
        self.savefig("julian_day_obs")

        pyplot.figure(self.figurenum())
        pyplot.clf()
        ax = sns.barplot(julianDayGenesis[:, 0], julianDayGenesis[:, 1])
        ax.set_xlim((1, 365))
        ax.set_xlabel('Day of year')

        ax.set_ylabel('Number of genesis events')
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
        k = np.where(dlon < -180.)
        dlon[k] += 360.

        pyplot.subplot(211)
        pyplot.plot(dlon[1:], dlon[:-1], 'k.', markersize=1)
        params = linreg(dlon)
        pyplot.text(-3, 3, "r = %5.3f"%params[2], ha='center',
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
        params = linreg(dlat)
        pyplot.text(-3, 3, "r = %5.3f"%params[2], ha='center',
                    va='center', color='r', size=14)
        pyplot.xlim(-4., 4.)
        pyplot.ylim(-4., 4.)
        pyplot.xticks(np.arange(-4., 4.1, 1.))
        pyplot.yticks(np.arange(-4., 4.1, 1.))
        pyplot.ylabel(r"$\Delta lat (t)$", fontsize=16)
        pyplot.xlabel(r"$\Delta lat (t-1)$", fontsize=16)
        pyplot.title("Latitude rate of change")
        pyplot.tight_layout()
        self.savefig('lonlat_corr')

        labels = [r'$\Delta \phi(t)$', r'$\Delta \phi(t-1)$']
        self.scatterHistogram(dlon[1:], dlon[:-1], labels, 'dlon_scatterHist')
        labels = [r'$\Delta \lambda(t)$', r'$\Delta \lambda(t-1)$']
        self.scatterHistogram(dlat[1:], dlat[:-1], labels, 'dlat_scatterHist')

    def plotFrequency(self, years, frequency):
        """Plot annual frequency of TCs, plus linear trend line"""

        pyplot.figure(self.figurenum(), figsize=(14, 5))
        pyplot.plot(years[1:-1], frequency[1:-1]) 
        xmax = 5*int((1 + years.max()/5))
        xmin = 5*int((years.min()/5))
        ymax = 5*int(1 + frequency.max()/5)
        pyplot.xlim(xmin, xmax)
        pyplot.ylim(0.0, ymax)
        
        pyplot.xlabel("Year")
        pyplot.ylabel("Frequency")

        params = linregress(np.array([years, frequency]))
        x = np.arange(xmin,xmax)
        y = params[0]*x + params[1]

        pyplot.plot(x, y, 'r--')
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
                if len(pvals) > 0:
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
        files.flSaveFile(os.path.join(self.outpath, 'min_pressure_lat.csv'), 
                         x,  delimiter=',', fmt='%6.2f')

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
        pcarray = np.array(pcarray)
        pc = np.take(pcarray, np.where(pcarray<sys.maxint))
        ax = sns.distplot(pc, bins=pbins, fit=frechet_l, 
                          kde_kws={'label':'KDE'},
                          fit_kws={'color':'r',
                                   'label':'Fitted distribution'})

        sns.despine(ax=ax, offset=10, trim=True)
        ax.set_xlabel("Minimum central pressure (hPa)")
        ax.set_ylabel("Probability")
        ax.set_title("Distribution of minimum central pressure")
        ax.legend(loc=0)
        pyplot.tight_layout()
        self.savefig("min_pressure_hist")

    def plotSpeedBear(self, sAllData, bAllData):
        """
        Plot speed and bearing against each other
        """
        pyplot.figure(self.figurenum(), figsize=(7, 7))
        pyplot.subplot(111)
        ii = np.where((sAllData < sys.maxint) & (bAllData < sys.maxint))
        pyplot.polar((np.pi/2. - np.radians(bAllData[ii])), sAllData[ii],
                     'k.', markersize=2)
        thetalabels = (90 - np.arange(0, 360, 45))
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

    def quantile(self, data, parameterName, dist='normal'):
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

if __name__ == "__main__":

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
