import sys
import numpy as np
import scipy.stats as stats

from matplotlib.figure import Figure
from matplotlib.artist import setp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap

class CurveFigure(Figure):

    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []
        self.subplots_adjust(hspace=0.5)

    def add(self, xdata, ydata, xlabel='x', ylabel='y', title='x'):
        self.subfigures.append((xdata, ydata, xlabel, ylabel, title))
    
    def subplot(self, axes, subfigure):
        xdata, ydata, xlabel, ylabel, title = subfigure

        axes.plot(xdata, ydata, '-', color='k')
        axes.set_xlabel(xlabel)
        axes.set_ticklabels(xdata.astype(int))
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

    def addGrid(self, axes):
        axes.grid(True, which='both', color='k', linestyle=':', linewidth=0.2)

    def addRange(self, axes, xdata, ymin, ymax):
        axes.fill_between(xdata, ymax, ymin, facecolor='0.75',
                          edgecolor='0.99', alpha=0.7)


    def plot(self):
        """
        Plot the subfigures. 

        For a number of subplots, they are organised into an approximately
        square arrangement.
        """

        n = len(self.subfigures)
        r = int(np.ceil(np.sqrt(n)))
        c = int(np.ceil(n / r))

        w, h = self.get_size_inches()
        self.set_size_inches(w * c, r * h)
        for i, subfigure in enumerate(reversed(self.subfigures)):
            x = i % c + 1
            y = i // c + 1
            axes = self.add_subplot(n, x, y)
            self.subplot(axes, subfigure)

class SemilogxCurve(CurveFigure):

    def subplot(self, axes, subfigure):
        xdata, ydata, xlabel, ylabel, title = subfigure

        axes.semilogx(xdata, ydata, '-', color='k', subsx=xdata)
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

class RangeCurve(CurveFigure):

    def add(self, xdata, ymean, ymax, ymin, xlabel='x', ylabel='y', title='x'):
        self.subfigures.append((xdata, ymean, ymax, ymin,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        xdata, ymean, ymax, ymin, xlabel, ylabel, title = subfigure

        axes.plot(xdata, ymean, color='k', lw=2)
        self.addRange(axes, xdata, ymin, ymax)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

class RangeCompareCurve(CurveFigure):

    def add(self, xdata, y1, y2, y2max, y2min, xlabel='x', ylabel='y', title='x'):
        self.subfigures.append((xdata, y1, y2, y2max, y2min,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        xdata, y1, y2, y2max, y2min, xlabel, ylabel, title = subfigure

        axes.plot(xdata, y1, color='r', lw=2, label="")
        axes.plot(xdata, y2, color='k', lw=2, label="")
        self.addRange(axes, xdata, y2min, y2max)
        
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

class SemilogRangeCurve(CurveFigure):

    def add(self, xdata, ymean, ymax, ymin, xlabel, ylabel, title):
        self.subfigures.append((xdata, ymean, ymax, ymin,
                                xlabel, ylabel, title))
    def addGrid(self, axes):
        xlims = axes.get_xlim()
        xmax = np.max(xlims)
        xticks = np.array([])
        x = np.arange(1,10, dtype=int)
        for i in range(int(np.log10(xmax)+1)):
            xticks = np.append(xticks, x*np.power(10, i).astype(int))

        xticks = np.append(xticks, xmax)
        axes.set_xticks(xticks)
        axes.autoscale(True, axis='x', tight=True)
        axes.grid(True, which='both', color='k', linestyle=':', linewidth=0.2)
        
    def subplot(self, axes, subfigure):
        xdata, ymean, ymax, ymin, xlabel, ylabel, title = subfigure

        axes.semilogx(xdata, ymean, color='k', lw=2, subsx=xdata)
        if (ymin[0] > 0) and (ymax[0] > 0):
            self.addRange(axes, xdata, ymin, ymax)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title, fontsize='small')
        ylim = (0., np.ceil(ymean.max()/10.)*10.)
        axes.set_ylim(ylim)
        self.addGrid(axes)

class SemilogRangeCompareCurve(CurveFigure):

    def add(self, xdata, y1, y2, y2max, y2min, xlabel='x', ylabel='y', title='x'):
        self.subfigures.append((xdata, y1, y2max, y2min,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        xdata, y1, y2, y2max, y2min, xlabel, ylabel, title = subfigure

        axes.semilogx(xdata, y1, color='r', lw=2, label="", subsx=xdata)
        axes.semilogx(xdata, y2, color='k', lw=2, label="", subsx=xdata)
        self.addRange(axes, xdata, y2min, y2max)
        
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)


class HazardCurve(SemilogRangeCurve):

    def plot(self, years, wspd, wspdupper, wspdlower, xlabel, ylabel, title):
        self.add(years, wspd, wspdupper, wspdlower, xlabel, ylabel, title)
        super(HazardCurve, self).plot()

class DistributionCurve(RangeCompareCurve):

    def plot(self, x, y1, y2, y2max, y2min, xlabel, ylabel, title):
        self.add(x, y1, y2, y2max, y2min, xlabel, ylabel, title)
        super(DistributionCurve, self).plot()


def saveFigure(figure, filename):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename)

def saveHazardCurve(years, wspd, wspdupper, wspdlower, xlabel, ylabel, title, filename):
    fig = HazardCurve()
    fig.plot(years, wspd, wspdupper, wspdlower, xlabel, ylabel, title)
    saveFigure(fig, filename)

def saveDistributionCurve(x, y1, y2, y2max, y2min, xlabel, ylabel, title, filename):
    fig = DistributionCurve()
    fig.plot(x, y1, y2, y2max, y2min, xlabel, ylabel, title)
    saveFigure(fig, filename)
    
