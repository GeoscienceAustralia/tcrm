from __future__ import division

import sys
import numpy as np
import scipy.stats as stats
import wind.windmodels as windmodels

from matplotlib import pyplot
from matplotlib.figure import Figure
from matplotlib.artist import setp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

import seaborn
seaborn.set(style="ticks")

class WindProfileFigure(Figure):

    def __init__(self, lat, lon, eP, cP, rMax, beta, beta1=1.5, beta2=1.4):

        Figure.__init__(self)
        self.R = np.array(range(1, 201), 'f')
        self.lat = lat
        self.lon = lon
        self.rMax = rMax
        self.eP = eP
        self.cP = cP
        self.beta = beta
        self.beta1 = beta1
        self.beta2 = beta2

    def plot(self, profileType=None):
        profiles = []

        if profileType:
            profiles.append(profileType)
        else:
            for p in windmodels.PROFILES:
                profiles.append(p)

        ax = self.add_subplot(1, 1, 1)
        ax.hold(True)
        legend = []

        for name in profiles:
            try:
                cls = windmodels.profile(name)
                params = windmodels.profileParams(name)
                values = [getattr(self, p) for p in params if hasattr(self, p)]
                profile = cls(self.lat, self.lon, self.eP, self.cP,
                              self.rMax, *values)
                V = profile.velocity(self.R)
                ax.plot(self.R, abs(V), linewidth=2)
                legend.append(name.capitalize())
            except TypeError:
                pass

        ax.legend(legend)
        ax.grid()
        ax.set_xlabel('Radius (km)', fontsize=14)
        ax.set_ylabel('Wind speed (m/s)', fontsize=14)
        ax.set_title((r'$P_c = %d\hspace{0.5}hPa,\hspace{1} P_e' +
                      r'= %d \hspace{0.5} hPa,\hspace{1} R_{max}' +
                      r'= %d \hspace{0.5}km$') %
                    (self.cP / 100., self.eP / 100., self.rMax))

class ScatterHistogramFigure(Figure):

    def __init__(self):
        Figure.__init__(self)

    def render(self, axes, x, y):
        color = axes._get_lines.color_cycle

        axes.scatter(x, y, color=color.next(), s=3)
        axes.set_aspect('equal')

        divider = make_axes_locatable(axes)
        axesTop = divider.append_axes('top', 1.2, pad=0.1, sharex=axes)
        axesRight = divider.append_axes('right', 1.2, pad=0.2, sharey=axes)

        # Make labels invisible:
        setp(axesTop.get_xticklabels() + axesRight.get_yticklabels(),
                    visible=False)

        axesTop.hist(x, bins=20, fc='k', ec=None, normed=True)
        axesRight.hist(y, bins=20, fc='k', ec=None, orientation='horizontal', normed=True)
        [t.set_rotation(-90) for t in axesRight.get_xticklabels()]


    def plot(self, xdata, ydata):
        i = np.where((xdata < sys.maxint) & (ydata < sys.maxint))[0]
        x = xdata[i]
        y = ydata[i]

        axes = self.add_subplot(1, 1, 1)
        self.render(axes, x, y)

class QuantileFigure(Figure):
    
    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []

    def percentiles(self, data):
        """Calculate percentile values from 1 to 100"""
        p = np.array([stats.scoreatpercentile(data, q) for q in range(1, 101)])
        return p

    def percentilerange(self, data):
        """Calculate a range of percentile values"""
        b = np.zeros((1000,100))
        for n in xrange(1000):
            b[n, :] = np.array([np.random.choice(data) for _ in xrange(100)])
        
        d = stats.scoreatpercentile(b, range(1,101), axis=1)
        u = stats.scoreatpercentile(d, 95, axis=1)
        l = stats.scoreatpercentile(d, 5, axis=1)
        return u, l

    def add(self, xdata, ydata, axisrange, xlabel='x', ylabel='y',
            title='Q-Q plot'):

        self.subfigures.append((xdata, ydata, axisrange, 
                                xlabel, ylabel, title))

    def formatAxes(self, axes, limits):
        axes.set_ylim(limits)
        axes.set_xlim(limits)

        #axes.set_xticks(np.linspace(*limits, num=11))
        #axes.set_yticks(np.linspace(*limits, num=11))

    def addGrid(self, axes):
        axes.grid(True, which='major', color='k', linestyle='-', linewidth=0.1)


    def subplot(self, axes, subfigure):
        xdata, ydata, axisrange, xlabel, ylabel, title = subfigure

        xper = self.percentiles(xdata)
        yper = self.percentiles(ydata)

        xupper, xlower = self.percentilerange(xdata)
        yupper, ylower = self.percentilerange(ydata)
        
        axes.scatter(xper, yper, s=5, c='r', marker='o', 
                     edgecolor='none', alpha=0.5)

        axes.plot(xupper, ylower, color='0.5', linestyle='solid', linewidth=0.5)
        axes.plot(xlower, yupper, color='0.5', linestyle='solid', linewidth=0.5)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        
        xx = np.arange(*axisrange)
        axes.plot(xx, xx, 'b--', linewidth=0.5)

        self.formatAxes(axes, axisrange)
        self.addGrid(axes)

    def plot(self):
        n = len(self.subfigures)
        r = int(np.ceil(np.sqrt(n)))
        c = int(np.ceil(n / r))
        w, h = self.get_size_inches()
        self.set_size_inches(w * c, r * h)
        for i, subfigure in enumerate(self.subfigures):
            axes = self.add_subplot(r, c, i)
            self.subplot(axes, subfigure)

class RegressionFigure(Figure):

    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []

    def add(self, xdata, ydata, xlabel='x', ylabel='y',
        title='x', transform=lambda x: x):
        self.subfigures.append((xdata, ydata, xlabel, ylabel, title, transform))

    def prepareData(self, xdata, ydata, transform):
        i = np.where((xdata < sys.maxint) & (ydata < sys.maxint))[0]
        xdata, ydata = transform(xdata), transform(ydata)
        return xdata[i], ydata[i]

    def formatAxes(self, axes, xdata, ydata):
        def filterOutliers(data, m=12.):
            d = np.abs(data - np.median(data))
            mdev = np.median(d)
            s  = d/mdev if mdev else 0
            return data[s<m]

        xmin = filterOutliers(xdata).min()
        xmax = filterOutliers(xdata).max()
        ymin = filterOutliers(ydata).min()
        ymax = filterOutliers(ydata).max()

        axes.set_xlim(xmin, xmax)
        axes.set_ylim(ymin, ymax)
        axes.set_xticks(np.linspace(xmin, xmax, 7))
        axes.set_yticks(np.linspace(ymin, ymax, 7))

        return xmin, xmax

    def subplot(self, axes, subfigure):
        xdata, ydata, xlabel, ylabel, title, transform = subfigure
        color = axes._get_lines.color_cycle

        xdata, ydata = self.prepareData(xdata, ydata, transform)

        k = color.next()
        axes.plot(xdata, ydata, '.', alpha=0.8, color=k)

        xmin, xmax = self.formatAxes(axes, xdata, ydata)

        m, c, r, p, e = stats.linregress(np.array([xdata, ydata]))
        x = np.array([xmin, xmax])
        y = m * x + c

        k = color.next()
        axes.plot(x, x, '-', color=k)

        k = color.next()
        axes.plot(x, y, '-', color=k, label='r = %5.3f' % r)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)

        legend = axes.legend(loc=2)
        legend.get_frame().set_alpha(0.5)

    def plot(self):
        n = len(self.subfigures)
        r = int(np.ceil(np.sqrt(n)))
        c = int(np.ceil(n / r))
        w, h = self.get_size_inches()
        self.set_size_inches(w * c, r * h)
        for i, subfigure in enumerate(self.subfigures):
            axes = self.add_subplot(r, c, i + 1)
            self.subplot(axes, subfigure)
        


class LaggedRegressionFigure(RegressionFigure):

    def add(self, data, label='x', title='x', transform=lambda x: x):
        super(LaggedRegressionFigure, self).add(data[1:],
                                                data[:-1],
                                                r'$' + label + r'(t-1)$',
                                                r'$' + label + r'(t)$',
                                                title,
                                                transform)

    def subplot(self, axes, subfigure):
        super(LaggedRegressionFigure, self).subplot(axes, subfigure)
        axes.set_aspect('equal')

class LineRegressionFigure(RegressionFigure):

    def add(self, x, y, xlabel='x', ylabel='y', title='x'):
        super(LineRegressionFigure, self).add(x, y,
                                              xlabel,
                                              ylabel,
                                              title)

    def subplot(self, axes, subfigure):
        xdata, ydata, xlabel, ylabel, title, transform = subfigure
        color = axes._get_lines.color_cycle

        xdata, ydata = self.prepareData(xdata, ydata, transform)

        k = color.next()
        axes.plot(xdata, ydata, '-', alpha=0.8, color=k)

        xmin, xmax = self.formatAxes(axes, xdata, ydata)

        m, c, r, p, e = stats.linregress(np.array([xdata, ydata]))
        x = np.array([xmin, xmax])
        y = m * x + c

        k = color.next()
        axes.plot(x, x, '-', color=k)

        k = color.next()
        axes.plot(x, y, '-', color=k, label='r = %5.3f' % r)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)

        legend = axes.legend(loc=2)
        legend.get_frame().set_alpha(0.5)

class PressureFigure(LaggedRegressionFigure):

    def plot(self, pressures, pressureRates):
        self.add(pressures, 'p', 'Pressure')
        self.add(pressureRates, r'\Delta p', 'Pressure rate of change')
        super(PressureFigure, self).plot()


class SpeedFigure(LaggedRegressionFigure):

    def plot(self, speeds, speedRates):
        self.add(speeds, 'v', 'Speed')
        self.add(speedRates, r'\Delta v', 'Speed rate of change (Acceleration)')
        super(SpeedFigure, self).plot()


class BearingFigure(LaggedRegressionFigure):

    def plot(self, bearings, bearingRates):
        def transform(z):
            return np.cos(np.radians(z))
        self.add(bearings, r'\cos(\theta)', 'Bearing', transform)
        self.add(bearingRates, r'\Delta\theta',
                 'Bearing rate of change')
        super(BearingFigure, self).plot()


class FrequencyFigure(LineRegressionFigure):

    def plot(self, years, frequency):
        self.add(np.array(years[1:-1], int), frequency[1:-1], 'Year', 'Frequency', 'Annual Frequency')
        super(FrequencyFigure, self).plot()


def saveFigure(figure, filename):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename)


def saveWindProfilesFigure(lat, lon, eP, cP, rMax, beta,
                           filename='docs/windprofiles.png'):
    fig = WindProfileFigure(lat, lon, eP, cP, rMax, beta)
    fig.plot()
    saveFigure(fig, filename)


def savePressureFigure(pressures, pressureRates, filename='docs/prs_corr.png'):
    fig = PressureFigure()
    fig.plot(pressures, pressureRates)
    saveFigure(fig, filename)
    fig = ScatterHistogramFigure()
    fig.plot(pressures[1:], pressures[:-1])
    saveFigure(fig, 'docs/pressuresh.png')
    fig = ScatterHistogramFigure()
    fig.plot(pressureRates[1:], pressureRates[:-1])
    saveFigure(fig, 'docs/pressureratesh.png')


def saveSpeedFigures(speeds, speedRates):
    fig = SpeedFigure()
    fig.plot(speeds, speedRates)
    saveFigure(fig, 'docs/spd_corr.png')
    fig = ScatterHistogramFigure()
    fig.plot(speeds[1:], speeds[:-1])
    saveFigure(fig, 'docs/speedsh.png')
    fig = ScatterHistogramFigure()
    fig.plot(speedRates[1:], speedRates[:-1])
    saveFigure(fig, 'docs/speedratesh.png')

def saveBearingFigure(bearings, bearingRates):
    fig = BearingFigure()
    fig.plot(bearings, bearingRates)
    saveFigure(fig, 'docs/bear_corr.png')
    fig = ScatterHistogramFigure()
    fig.plot(bearings[1:], bearings[:-1])
    saveFigure(fig, 'docs/bearingsh.png')
    fig = ScatterHistogramFigure()
    fig.plot(bearingRates[1:], bearingRates[:-1])
    saveFigure(fig, 'docs/bearingratesh.png')
    


def saveFrequencyFigure(years, frequency, filename='docs/freq_corr.png'):
    fig = FrequencyFigure()
    fig.plot(years, frequency)
    saveFigure(fig, filename)


def main():
    import Utilities.files as files
    from os.path import join as pjoin, dirname, normpath
    baseDir = normpath(pjoin(dirname(__file__), '..'))
    inputPath = pjoin(baseDir, 'output', 'process')

    saveWindProfilesFigure(-12., 130., 100700., 95000., 30., 1.6)

    pressureRates = files.flLoadFile(pjoin(inputPath, 'pressure_rate'))
    pressures = files.flLoadFile(pjoin(inputPath, 'all_pressure'))
    savePressureFigure(pressures, pressureRates)

    speedRates = files.flLoadFile(pjoin(inputPath, 'speed_rate'))
    speeds = files.flLoadFile(pjoin(inputPath, 'all_speed'))
    saveSpeedFigures(speeds, speedRates)

    bearingRates = files.flLoadFile(pjoin(inputPath, 'bearing_rate'))
    bearings = files.flLoadFile(pjoin(inputPath, 'all_bearing'))
    saveBearingFigure(bearings, bearingRates)

    freq = files.flLoadFile(pjoin(inputPath, 'frequency'))
    saveFrequencyFigure(np.array(freq[:,0],int), freq[:,1])

if __name__ == "__main__":
    main()
