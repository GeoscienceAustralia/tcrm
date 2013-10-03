import os, logging, sys

from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import linregress, probplot
from scipy.stats import scoreatpercentile as percentile
import numpy

import Utilities.files as files
import Utilities.config as config
import Utilities.nctools as nctools
import Utilities.stats as stats

__version__ = '$Id: plotStats.py 810 2012-02-21 07:52:50Z nsummons $'

log = logging.getLogger(__name__)

def get_bins(data, allpos=False):
    """
    Establish nice bin widths for plotting histograms - as a default it sets
    the width to 0.25 times the standard deviation of the data
    """
    binwidth = 0.25*numpy.std( data )
    dmax = numpy.max(numpy.fabs( data ))

    lims = (int(dmax/binwidth) + 1) * binwidth
    if allpos:

        bins = numpy.arange(0, lims+binwidth, binwidth)
    else:
        bins = numpy.arange(-lims, lims+binwidth, binwidth)

    return bins


def _linreg(data):
    """
    Calculate the linear regression of the data against itself (lag-1)
    Returns the slope, intercept, correlation, two-tailed
    probability and standard error of the estimate.
    """
    tData = numpy.array([data[1:], data[:-1]])
    i = numpy.where((tData[0, :] < sys.maxint) & (tData[1, :] < sys.maxint))[0]
    m, c, r, pr, err = linregress(tData[:, i])
    return m, c, r, pr, err


class PlotData:
    """
    Base class for plotting summaries of input data.

    """
    def __init__(self, output_path, output_format):
        """
        Initialise statistical plotting routines, fixing the output path and
        the format of all images generated
        """
        self.outpath = output_path
        self.fmt = output_format
        self.fignum = 0

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
        pyplot.savefig(os.path.join(self.outpath, outputfile),
                       format=self.fmt, **kwargs)
        pyplot.close()

    def scatterHistogram(self, x, y, img_name, allpos=False):
        """
        Generate a combined scatter/histogram plot of data
        """
        pyplot.figure(self.figurenum(), figsize=(8, 8))
        i = numpy.where((x < sys.maxint) & (y < sys.maxint))[0]
        axScatter = pyplot.subplot(111)
        axScatter.scatter(x[i], y[i], color='k', s=3)
        axScatter.set_aspect(1.)
        for yticks in axScatter.get_yticklabels():
            yticks.set_fontsize('small')
        for xticks in axScatter.get_xticklabels():
            xticks.set_fontsize('small')

        # Create new axes instance on the right of the current axes
        divider = make_axes_locatable(axScatter)
        axHistX = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter)
        axHistY = divider.append_axes("right", 1.2, pad=0.2, sharey=axScatter)

        # Make labels invisible:
        pyplot.setp(axHistX.get_xticklabels() + axHistY.get_yticklabels(), visible=False)

        xbins = get_bins(x[i], allpos)
        ybins = get_bins(y[i], allpos)

        axHistX.hist(x[i], bins=xbins, fc='k', ec=None)
        axHistY.hist(y[i], bins=ybins, fc='k', ec=None, orientation='horizontal')

        for xticks in axHistX.get_xticklabels():
            xticks.set_visible(False)
        for yticks in axHistX.get_yticklabels():
            yticks.set_fontsize('x-small')

        for yticks in axHistY.get_yticklabels():
            yticks.set_visible(False)
        for xticks in axHistY.get_xticklabels():
            xticks.set_fontsize('x-small')
            xticks.set_rotation(90)

        outputfile = '%s.%s' % (img_name, self.fmt)
        pyplot.savefig(os.path.join(self.outpath, outputfile), format=self.fmt)
        pyplot.close()

    def plotPressure(self, pAllData, pRateData):
        """
        Plot the input pressure values lagged against themselves,
        and the same for the changes in pressure.
        """
        pyplot.rcdefaults()

        pyplot.rcParams['figure.figsize'] = (7, 12)
        pyplot.figure(self.figurenum())
        pyplot.subplot(211)
        pyplot.plot(pAllData[1:], pAllData[:-1], 'k.', markersize=1)

        m, c, r, p, e = _linreg(pAllData)

        x = numpy.arange(880., 1021., 1.)
        y = m*x + c

        pyplot.plot(x, y, 'r-')
        pyplot.plot(x, x, 'k-')
        pyplot.text(900, 1010, "r = %5.3f"%r, ha='center', va='center',
                    color='r', size=14)
        pyplot.xlim(880., 1020.)
        pyplot.ylim(880., 1020.)
        pyplot.xticks(numpy.arange(880., 1021., 20.))
        pyplot.yticks(numpy.arange(880., 1021., 20.))
        pyplot.ylabel(r"$p (t)$", fontsize=16)
        pyplot.xlabel(r"$p (t-1)$", fontsize=16)
        #pyplot.grid(True, linewidth=0.5)
        pyplot.title("Pressure")

        pyplot.subplot(212)

        pyplot.plot(pRateData[1:], pRateData[:-1], 'k.', markersize=1)
        m, c, r, p, e = _linreg(pRateData)
        pyplot.text(-8., 8., "r = %5.3f"%r, ha='center', va='center',
                    color='r', size=14)
        pyplot.xlim(-10., 10.)
        pyplot.ylim(-10., 10.)
        pyplot.xticks(numpy.arange(-10., 11., 2.5))
        pyplot.yticks(numpy.arange(-10., 11., 2.5))
        pyplot.ylabel(r"$\partial p/\partial t (t)$", fontsize=16)
        pyplot.xlabel(r"$\partial p/\partial t (t-1)$", fontsize=16)
        #pyplot.grid(True)
        pyplot.title("Pressure rate of change")

        outputfile = 'prs_corr.%s' % self.fmt
        pyplot.savefig(os.path.join(self.outpath, outputfile), format=self.fmt)
        pyplot.close()

    def plotSpeed(self, speedData, speedRate):
        """
        Plot the input speed values lagged against themselves,
        and the same for the changes in speed.
        """
        pyplot.rcParams['figure.figsize'] = (7, 12)
        pyplot.figure(self.figurenum())
        pyplot.subplot(211)
        pyplot.plot(speedData[1:], speedData[:-1], 'k.', markersize=1)

        m, c, r, p, e = _linreg(speedData)

        x = numpy.arange(0., 101., 1.)
        y = m*x+c
        pyplot.plot(x, y, 'r-')
        pyplot.plot(x, x, 'k-')
        pyplot.text(10, 90, "r = %5.3f"%r, ha='center', va='center',
                    color='r', size=14)
        pyplot.xlim(0., 100.)
        pyplot.ylim(0., 100.)
        pyplot.xticks(numpy.arange(0, 101., 20))
        pyplot.yticks(numpy.arange(0, 101., 20))
        #pyplot.grid(True)
        pyplot.ylabel(r"$v (t)$", fontsize=14)
        pyplot.xlabel(r"$v (t-1)$", fontsize=14)

        pyplot.title("Speed")

        pyplot.subplot(212)
        pyplot.plot(speedRate[1:], speedRate[:-1], 'k.', markersize=1)
        m, c, r, p, e = _linreg(speedRate)
        pyplot.text(-25, 25, "r = %5.3f"%r,
                    ha='center', va='center', color='r', size=14)
        pyplot.xlim(-30.,30.)
        pyplot.ylim(-30.,30.)
        pyplot.xticks(numpy.arange(-30.,31.,10.))
        pyplot.yticks(numpy.arange(-30.,31.,10.))

        pyplot.ylabel(r"$\partial v /\partial t (t)$", fontsize=14)
        pyplot.xlabel(r"$\partial v /\partial t (t-1)$", fontsize=14)
        #pyplot.grid(True)
        pyplot.title("Speed rate of change")
        self.savefig("spd_corr")

        pyplot.rcdefaults()

        self.scatterHistogram(speedData[1:], speedData[:-1],
                              'speed_scatterHist', allpos=True)

        self.scatterHistogram(speedRate[1:], speedRate[:-1],
                              'speedRate_scatterHist', allpos=True)

    def plotBearing(self, bAllData, bRateData):
        """
        Plot the (cosine of) input bearing values lagged against themselves,
        and the same for the changes in bearing.
        """
        pyplot.rcParams['figure.figsize'] = (7, 12)
        pyplot.figure(self.figurenum())
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
        pyplot.plot(numpy.cos(numpy.radians(bAllData_t0)),
                    numpy.cos(numpy.radians(bAllData_tm1)),
                    'k.', markersize=1)
        pyplot.xlim(-1., 1.)
        pyplot.ylim(-1., 1.)
        pyplot.xticks(numpy.arange(-1.,1.1,1.))
        pyplot.yticks(numpy.arange(-1.,1.1,1.))
        m, c, r, p, e = _linreg(numpy.cos(numpy.radians(bAllData)). \
                                compress(bAllData < sys.maxint))

        x = numpy.arange(-1.0,1.1,0.1)
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
        m,c,r,p,e = _linreg(numpy.cos(numpy.radians(bRateData)))
        pyplot.text(-25, 25, "r = %5.3f"%r,
                    ha='center', va='center', color='r', size=14)

        pyplot.xlim(-30., 30.)
        pyplot.ylim(-30., 30.)
        pyplot.xticks(numpy.arange(-30., 31., 10.))
        pyplot.yticks(numpy.arange(-30., 31., 10.))
        pyplot.ylabel(r"$\partial \theta /\partial t (t)$", fontsize=16)
        pyplot.xlabel(r"$\partial \theta /\partial t (t-1)$", fontsize=16)
        #pyplot.grid(True)
        pyplot.title("Bearing rate of change")
        self.savefig('bear_corr')
        #pyplot.savefig(os.path.join(outputPath,'bear_corr.eps'))
        self.scatterHistogram(bAllData[1:], bAllData[:-1],
                              'bear_scatterHist', allpos=True)

        self.scatterHistogram(bRateData[1:], bRateData[:-1],
                              'bearRate_scatterHist', allpos=True)

    def julianDay(self, julianDayObs, julianDayGenesis):
        """
        Plot bar graphs of the number of TC observations and genesis events
        """
        pyplot.rcParams['figure.figsize'] = (10, 6)
        pyplot.figure(self.figurenum())
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
        pyplot.rcParams['figure.figsize'] = (7, 12)
        pyplot.figure(self.figurenum())

        dlon = lonData[1:] - lonData[:-1]
        dlat = latData[1:] - latData[:-1]
        j = numpy.where(indicator[1:] == 0)
        dlon = dlon[j]
        dlat = dlat[j]

        # Correct change in longitude where the value jumps across the 180E
        # meridian
        k = numpy.where(dlon<-180.)
        dlon[k] += 360.

        pyplot.subplot(211)
        pyplot.plot(dlon[1:], dlon[:-1], 'k.', markersize=1)
        m, c, r, p, e = _linreg(dlon)
        pyplot.text(-3, 3, "r = %5.3f"%r, ha='center',
                    va='center', color='r', size=14)
        pyplot.xlim(-4., 4.)
        pyplot.ylim(-4., 4.)
        pyplot.xticks(numpy.arange(-4., 4.1, 1.))
        pyplot.yticks(numpy.arange(-4., 4.1, 1.))
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
        pyplot.xticks(numpy.arange(-4., 4.1, 1.))
        pyplot.yticks(numpy.arange(-4., 4.1, 1.))
        pyplot.ylabel(r"$\Delta lat (t)$", fontsize=16)
        pyplot.xlabel(r"$\Delta lat (t-1)$", fontsize=16)
        pyplot.title("Latitude rate of change")

        self.savefig('lonlat_corr')
        pyplot.rcdefaults()

        self.scatterHistogram(dlon[1:], dlon[:-1], 'dlon_scatterHist')

        self.scatterHistogram(dlat[1:], dlat[:-1], 'dlat_scatterHist')

    def plotFrequency(self, years, frequency):
        """Plot annual frequency of TCs, plus linear trend line"""
        pyplot.rcdefaults()
        pyplot.rcParams['figure.figsize'] = (14, 5)
        pyplot.figure(self.figurenum())
        pyplot.plot(years, frequency, 'k-', linewidth=3)
        xmax = 5*int((1 + years.max()/5))
        xmin = 5*int((years.min()/5))
        ymax = 5*int(1 + frequency.max()/5)
        pyplot.xlim(xmin, xmax)
        pyplot.ylim(0.0, ymax)
        pyplot.xticks(numpy.arange(xmin, xmax, 5))
        pyplot.yticks(numpy.arange(5, ymax, 5))

        pyplot.xlabel("Year")
        pyplot.ylabel("Frequency")

        m, c, r, pr, err = linregress(numpy.array([years, frequency]))
        x = numpy.arange(xmin,xmax)
        y = m*x + c

        pyplot.plot(x, y, 'r--', linewidth=0.5)

        pyplot.xlim(xmin, xmax)
        pyplot.ylim(0.0, ymax)
        pyplot.xticks(numpy.arange(xmin, xmax, 5))
        pyplot.yticks(numpy.arange(0, ymax, 5))

        pyplot.grid(True)
        pyplot.title("Annual frequency (%d - %d)"%(years.min(), years.max()))
        self.savefig('frequency')
        pyplot.rcdefaults()
        return

    def minPressureLat(self, pAllData, latData, latMin=-40., latMax=0.):
        """
        Plot the minimum central pressures as a function of latitude
        """
        rLat = numpy.round(latData, 0)
        lats = numpy.arange(latMin, latMax + 0.1, 1)
        minP = numpy.zeros(len(lats))
        n = 0
        for l in lats:
            i = numpy.where(rLat == l)[0]
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
        pyplot.plot(lats, minP, 'r-', linewidth=2, label=r'Min $P_{centre}$')
        pyplot.xlim(latMin, latMax)
        pyplot.ylim(800, 1020)

        pyplot.xlabel('Latitude', fontsize=10)
        pyplot.ylabel('Minimum central pressure (hPa)', fontsize=10)
        pyplot.legend(loc=3)
        pyplot.grid(True)

        self.savefig("min_pressure_lat")

        x = numpy.zeros((len(lats), 2))
        x[:, 0] = lats
        x[:, 1] = minP
        files.flSaveFile(os.path.join(self.outpath, 'min_pressure_lat.csv'), x,
                         delimiter=',', fmt='%6.2f')

    def minPressureHist(self, index, pAllData):
        """
        Plot a histogram of the minimum central pressures from the input
        dataset.
        """
        pyplot.rcParams['figure.figsize'] = (8, 7)
        pyplot.figure(self.figurenum())
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

        pbins = numpy.arange(850., 1020., 5)
        n, b, p = pyplot.hist(numpy.array(pcarray), pbins, normed=False,
                              lw=2, ec='k', fc='w')
        pyplot.xlabel("Minimum central pressure (hPa)")
        pyplot.ylabel("Count")
        pyplot.title("Distribution of minimum central pressure")
        self.savefig("min_pressure_hist")

        x = numpy.zeros((len(b), 2))
        x[:, 0] = b
        x[1:, 1] = n
        files.flSaveFile(os.path.join(self.outpath, 'min_pressure_hist.csv'), x,
                         delimiter=',', fmt='%6.2f')
        pyplot.rcdefaults()

    def plotSpeedBear(self, sAllData, bAllData):
        """
        Plot speed and bearing against each other
        """
        pyplot.rcParams['figure.figsize'] = (7, 7)
        pyplot.figure(self.figurenum())
        pyplot.subplot(111)
        ii = numpy.where((sAllData < sys.maxint) & (bAllData < sys.maxint))
        pyplot.polar((numpy.pi/2. - numpy.radians(bAllData[ii])), sAllData[ii],
                     'k.', markersize=2)
        thetalabels=(90 - numpy.arange(0, 360, 45))
        ii = numpy.where(thetalabels < 0)
        thetalabels[ii] += 360
        lines, labels = pyplot.rgrids(numpy.arange(20., 101., 20.),
                                      labels=None, angle=67.5)
        lines, labels = pyplot.thetagrids(numpy.arange(0., 360., 45.),
                                          thetalabels)
        pyplot.ylim(0, 100.)
        pyplot.grid(True)
        r = numpy.corrcoef(bAllData[ii], sAllData[ii])[1, 0]
        pyplot.text(45, 125, "r = %5.3f"%r, ha='center',
                    va='center', color='r', size=14)
        pyplot.title("Speed vs bearing")
        self.savefig('spd_bear_corr')
        pyplot.rcdefaults()

    def quantile(self, data, parameterName, dist='normal', mean=0.0, sigma=1.0):
        """
        Generate a probability plot of the given data; data should be an
        array of anomalies

        """
        pyplot.rcParams['figure.figsize'] = (8, 7)
        pyplot.figure(self.figurenum())
        pyplot.clf()
        d = data.compress(data < sys.maxint)
        m = numpy.average(d)
        sd = numpy.std(d)
        nd = (d-m)/sd
        (osm, osr), (slope, intercept, r) = probplot(nd, plot=pyplot)
        #pyplot.xticks(numpy.arange(-3.,3.1,1.))
        #pyplot.yticks(numpy.arange(-3.,3.1,1.))
        #pyplot.xlabel("Normal")
        pyplot.ylabel(parameterName)
        pyplot.title("Q-Q plot - %s" % parameterName)
        pyplot.xlim((-5, 5))
        pyplot.ylim((-5, 5))
        pos = (2, -4.8)
        pyplot.text(2, -4.9, r"$r^2=%1.4f$" % r, fontsize=12)

        self.savefig('qqplot_%s' % parameterName)
        pyplot.rcdefaults()

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
