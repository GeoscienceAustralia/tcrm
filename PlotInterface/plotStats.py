"""
:mod:`plotStats` -- statistical plotting routines
=================================================

.. module:: plotStats
    :synopsis: A collection of statistical plotting routines.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import sys

from os.path import join as pjoin
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

from PlotInterface.figures import LaggedRegressionFigure
from scipy.stats import linregress, probplot, frechet_l
import numpy as np

import seaborn as sns
sns.set_style("ticks")

plt.ioff()


MONTH_LOCATOR = mdates.MonthLocator(bymonthday=15, interval=1)
YEAR_LOCATOR = mdates.YearLocator()
TICK_LOCATOR = ticker.MultipleLocator(base=5)
MONTH_FORMATTER = mdates.DateFormatter("%b")
YEAR_FORMATTER = mdates.DateFormatter("%Y")

def linreg(data):
    """
    Calculate the linear regression of the data against itself (lag-1)
    Returns the slope, intercept, correlation, two-tailed
    probability and standard error of the estimate.
    """
    tData = np.array([data[1:], data[:-1]])
    i = np.where((tData[0, :] < sys.maxsize) & (tData[1, :] < sys.maxsize))[0]
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
        Provide a shortcut to plt.savefig() that adds the output path to the
        filename.
        All keyword args are passed without alteration
        """
        outputfile = ".".join([filename, self.fmt])
        plt.savefig(pjoin(self.outpath, outputfile),
                       format=self.fmt, **kwargs)
        plt.close()

    def scatterHistogram(self, xdata, ydata, labels, name,
                         transform=lambda x: x, **kwargs):
        """
        Create a scatter plot with marginal distributions

        :param x: `numpy.ndarray` of x values.
        :param y: `numpy.ndarray` of y values.
        :param list labels: A length-2 list with the x and y labels as strings.
        :param str name: Name to use in the file name for saving the
                         figure.
        :param func transform: A function to transform the data.
        :param kwargs: Additional keyword arguments passed to
                       `seaborn.jointplot`.


        """

        i = np.where((xdata < sys.maxsize) & (ydata < sys.maxsize))[0]
        xx = transform(xdata[i])
        yy = transform(ydata[i])
        jp = sns.jointplot(xx, yy, kind='reg',
                           joint_kws={'scatter_kws':
                                      {'color':'slategray',
                                       'alpha':0.5}},
                           **kwargs)

        jp.set_axis_labels(*labels)
        plt.tight_layout()

        self.savefig(name)

    def barPlot(self, xdata, ydata, name, labels):
        """
        Bar plot, with added trend line or mean line

        :param x: `numpy.ndarray` of x-values.
        :param y: `numpy.ndarray` of y-values.
        :param str name: Name of the parameter, which will be used
                         in the filename.
        :param list labels: A list of the x- and y-axis labels to apply.
        """

        fig, ax = plt.subplots(1, 1)
        sns.barplot(xdata, ydata, ax=ax)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.axhline(np.mean(ydata))
        sns.despine()
        fig.tight_layout()
        self.savefig(name)

    def plotRegression(self, xdata, name, labels, transform=lambda x: x):
        """
        A generic function to plot a lag-1 autoregression of the
        variable, with a joint probability plot including
        marginal distribution plots.

        :param x: `numpy.ndarray` of the data variable to plot.
        :param str name: Name of the variable. Used for saving the
                         resulting image.
        :param list labels: A list of the x- and y-axis labels to apply.
        :param func transform: A transform function to apply to
                               the data. Default is to leave the
                               data unchanged.

        """
        figure = LaggedRegressionFigure()
        figure.add(xdata, labels[0], name, transform)
        self.savefigure(figure, name)

        x_t = xdata[1:]
        x_tm = xdata[:-1]
        skip = (x_t >= sys.maxsize) | (x_tm >= sys.maxsize)
        x_t = x_t.compress(skip==False)
        x_tm = x_tm.compress(skip==False)

        x_t = transform(x_t)
        x_tm = transform(x_tm)

        self.scatterHistogram(x_t, x_tm, labels, name+'_scatter')


    def minPressureHist(self, index, pAllData):
        """
        Plot a histogram of the minimum central pressures from the input
        dataset.
        """

        plt.figure(self.figurenum(), figsize=(8, 7))
        plt.clf()
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
        pc = np.take(pcarray, np.where(pcarray<sys.maxsize))
        ax = sns.distplot(pc, bins=pbins, fit=frechet_l,
                          kde_kws={'label':'KDE'},
                          fit_kws={'color':'r',
                                   'label':'Fitted distribution'})

        sns.despine(ax=ax, offset=10, trim=True)
        ax.set_xlabel("Minimum central pressure (hPa)")
        ax.set_ylabel("Probability")
        ax.set_title("Distribution of minimum central pressure")
        ax.legend(loc=0)
        plt.tight_layout()
        self.savefig("min_pressure_hist")

    def plotSpeedBear(self, sAllData, bAllData):
        """
        Plot speed and bearing against each other
        """
        plt.figure(self.figurenum(), figsize=(7, 7))
        plt.subplot(111)
        ii = np.where((sAllData < sys.maxsize) & (bAllData < sys.maxsize))
        plt.polar((np.pi/2. - np.radians(bAllData[ii])), sAllData[ii],
                     'k.', markersize=2)
        thetalabels = (90 - np.arange(0, 360, 45))
        ii = np.where(thetalabels < 0)
        thetalabels[ii] += 360
        lines, labels = plt.rgrids(np.arange(20., 101., 20.),
                                      labels=None, angle=67.5)
        lines, labels = plt.thetagrids(np.arange(0., 360., 45.),
                                          thetalabels)
        plt.ylim(0, 100.)
        plt.grid(True)
        r = np.corrcoef(bAllData[ii], sAllData[ii])[1, 0]
        plt.text(45, 125, "r = %5.3f"%r, ha='center',
                    va='center', color='r', size=14)
        plt.title("Speed vs bearing")
        self.savefig('spd_bear_corr')

    def quantile(self, data, parameterName, dist='norm'):
        """
        Generate a probability plot of the given data; data should be an
        array of anomalies

        """
        #plt.figure(self.figurenum(), figsize=(8, 7))
        #plt.clf()
        d = data.compress(data < sys.maxsize)
        m = np.average(d)
        sd = np.std(d)
        nd = (d-m)/sd
        (osm, osr), (slope, intercept, r) = probplot(nd, dist=dist, plot=plt)

        plt.ylabel(parameterName)
        plt.title("Q-Q plot - %s" % parameterName)
        plt.xlim((-5, 5))
        plt.ylim((-5, 5))
        pos = (2, -4.8)
        plt.text(2, -4.9, r"$r^2=%1.4f$" % r, fontsize=12)

        self.savefig('qqplot_%s' % parameterName)


class PlotPressure(PlotData):

    def plotPressure(self, data):
        labels = [r'$p_c(t)$', r'$p_c(t-1)$']
        self.scatterHistogram(data[1:], data[:-1],
                              labels, 'pressure',
                              xlim=(850, 1020), ylim=(850, 1020))

    def plotPressureRate(self, data):
        labels = [r'$\frac{\delta p_c}{\delta t}(t)$',
                  r'$\frac{\delta p_c}{\delta t}(t-1)$']
        self.scatterHistogram(data[1:], data[:-1],
                              labels, 'pressure_rate')
        self.quantile(data, 'pressure_rate', 'logistic')

    def plotMinPressure(self, index, pAllData):
        """
        Plot the distribution of minimum central pressure, and
        include a fitted distribution (presently uses the
        `scipy.stats.frechet_l` distribution).

        :param index: `numpy.ndarray` of 1/0 that indicates the start of
                      separate TC tracks.
        :param pAllData: `numpy.ndarray` of pressure observations from
                         TC events.

        """
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
        pc = np.take(pcarray, np.where(pcarray<sys.maxsize))
        ax = sns.distplot(pc, bins=pbins, fit=frechet_l,
                          kde_kws={'label':'KDE'},
                          fit_kws={'color':'r',
                                   'label':'Fitted distribution'})

        sns.despine(ax=ax, offset=10, trim=True)
        ax.set_xlabel("Minimum central pressure (hPa)")
        ax.set_ylabel("Probability")
        ax.set_title("Distribution of minimum central pressure")
        ax.legend(loc=0)
        plt.tight_layout()
        self.savefig("min_pressure_hist")

class PlotBearing(PlotData):

    def plotBearing(self, data):
        labels = [r'$\theta(t)$', r'$\theta(t-1)$']
        def transform(x):
            return np.cos(x/2.)

        self.scatterHistogram(data[1:], data[:-1], labels, 'bearing',
                              transform=transform, xlim=(-1, 1),
                              ylim=(-1, 1))

    def plotBearingRate(self, data):
        labels = [r'$\frac{\delta \theta(t)}{\delta t}$',
                  r'$\frac{\delta \theta(t-1)}{\delta t}$']
        self.scatterHistogram(data[1:], data[:-1],
                              labels, 'bearing_rate')
        self.quantile(data, 'bearing_rate', 'logistic')

class PlotSpeed(PlotData):

    def plotSpeed(self, data):
        labels = [r'$v (t)$', r'$v (t-1)$']
        self.scatterHistogram(data[1:], data[:-1], labels, 'speed',
                              xlim=(0, 100), ylim=(0, 100))

    def plotSpeedRate(self, data):
        labels = [r'$\frac{\delta v(t)}{\delta t}$',
                  r'$\frac{\delta v(t-1)}{\delta t}$']
        self.scatterHistogram(data[1:], data[:-1],
                              labels, 'speed_rate')
        self.quantile(data, 'speed_rate', 'logistic')

class PlotFrequency(PlotData):

    def plotFrequency(self, years, frequency):
        """
        Plot annual count of events within the domain.

        TODO: Automatically adjust the x-tickmarks to be spaced
              nicely (i.e. no overlap of labels).
              Offer option of drawing the mean value (using
              `axes.axhline`) or a linear trend line.
        """
        labels = ["Year", "Number"]
        self.barPlot(years.astype(int), frequency, "frequency", labels)

class PlotDays(PlotData):

    def plotJulianDays(self, julianDayObs, julianDayGenesis):
        """
        Plot bar graphs of the number of TC observations per day of year
        and genesis events per day of year.

        TODO: Format the tick marks to represent monthly intervals.
              Add a KDE over both bar plots, remembering this is cyclic
              data (so KDE[-1] ~ KDE[0]).
        """

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        dateLocator = mdates.MonthLocator(bymonthday=15, interval=1)
        dateFormat = mdates.DateFormatter("%b")
        dates = mdates.num2date(julianDayObs[:, 0])

        sns.barplot(julianDayObs[:, 0], julianDayObs[:, 1], ax=ax1)
        ax1.set_xlim((1, 365))
        ax1.set_ylabel('Observations')

        sns.barplot(julianDayGenesis[:, 0], julianDayGenesis[:, 1], ax=ax2)
        ax2.set_xlim((1, 365))
        ax2.set_xlabel('Month')
        ax2.set_ylabel('Genesis events')

        ax2.xaxis.set_major_locator(MONTH_LOCATOR)
        ax2.xaxis.set_major_formatter(MONTH_FORMATTER)

        sns.despine()
        self.savefig("julian_day")

class PlotLonLat(PlotData):

    def plotLonLat(self, lonData, latData, indicator):
        dlon = lonData[1:] - lonData[:-1]
        dlat = latData[1:] - latData[:-1]
        j = np.where(indicator[1:] == 0)
        dlon = dlon[j]
        dlat = dlat[j]
        k = np.where(dlon < -180.)
        dlon[k] += 360.
        labels = [r"$\Delta \phi (t-1)$", r"$\Delta \phi (t)$"]
        self.scatterHistogram(dlon[1:], dlon[:-1], labels, 'delta_lon')

        labels = [r"$\Delta \lambda (t-1)$", r"$\Delta \lambda (t)$"]
        self.scatterHistogram(dlat[1:], dlat[:-1], labels, 'delta_lat')

        labels = [r"$\Delta \phi (t)$", r"$\Delta \lambda (t)$"]
        self.scatterHistogram(dlon, dlat, labels, 'lonlat_corr')

        self.quantile(dlon, 'delta_lon', 'logistic')
        self.quantile(dlat, 'delta_lat', 'logistic')

