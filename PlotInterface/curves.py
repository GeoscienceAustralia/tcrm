"""
:mod:`curves` -- plot curve figures
===================================

.. module:: curves
    :synopsis: provide methods to plot curves (hazard curves,
               distributions, logarithmic axes, etc).

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

The routines here make use of the themes from
`seaborn <https://seaborn.pydata.org/api.html>`_ to
define the line styles and annotations (font sizes, etc.).

"""


import logging
import numpy as np

import seaborn
seaborn.set_style("ticks")

from matplotlib.figure import Figure
from matplotlib.ticker import LogLocator, FormatStrFormatter

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

class CurveFigure(Figure):
    """
    Base class for plotting line figures.

    """

    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []
        self.subplots_adjust(hspace=0.5)

    def add(self, xdata, ydata, xlabel='x', ylabel='y', title='x'):
        """
        Add a new subplot to the list of subplots to be created.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param ydata: y values of the data points.
        :type ydata: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """

        self.subfigures.append((xdata, ydata, xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        """
        Draw a line plot on an :class:`matplotlib.axes` instance. The
        data and labels are contained in a tuple. The basic plot is a
        black line plot, with the x-ticks presented as integer values.
        A grid is added with a call to :meth:`RangeFigure.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.

        """

        xdata, ydata, xlabel, ylabel, title = subfigure

        axes.plot(xdata, ydata, '-')
        axes.set_xlabel(xlabel)
        axes.set_ticklabels(xdata.astype(int))
        axes.autoscale(True, axis='x', tight=True)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

    def addGrid(self, axes):
        """
        Add a graticule to the subplot axes. The x-axis is autoscaled,
        the grid is a dotted black line, linewidth 0.2.

        :param axes: :class:`axes` instance.

        """

        axes.autoscale(True, axis='x', tight=True)
        axes.grid(True, which='both', color='k', linestyle=':', linewidth=0.2)

    def addRange(self, axes, xdata, ymin, ymax):
        """
        Add a shaded range of values to the plot.

        :param axes: :class:`matplotlib.axes` instance.
        :param xdata: x-data values.
        :type xdata: `numpy.ndarray`
        :param ymin: Lower range of y values.
        :type ymin: `numpy.ndarray`
        :param ymax: Upper range of y values.
        :type ymax: `numpy.ndarray`

        """

        axes.fill_between(xdata, ymax, ymin, alpha=0.5)


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
        for i, subfigure in enumerate(self.subfigures):
            axes = self.add_subplot(r, c, i + 1)
            self.subplot(axes, subfigure)
        self.tight_layout()

class SemilogCurve(CurveFigure):
    """
    Extend the basic :class:`CurveFigure` to use a logarithmic scale
    on the x-axis.
    """

    def subplot(self, axes, subfigure):
        """
        Draw a line plot on an :class:`matplotlib.axes` instance, with a
        logarithmic scale on the x-axis. Data and labels are
        contained in a tuple. x-ticks presented as integer values.
        A grid is added with a call to :meth:`RangeFigure.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.

        """

        xdata, ydata, xlabel, ylabel, title = subfigure

        axes.semilogx(xdata, ydata, '-', subsx=xdata)
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

    def addGrid(self, axes):
        """
        Add a logarithmic graticule to the subplot axes.

        :param axes: :class:`axes` instance.

        """

        axes.xaxis.set_major_locator(LogLocator())
        axes.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        axes.xaxis.set_minor_locator(LogLocator(subs=[.1, .2, .3, .4, .5, .6, .7, .8, .9]))
        axes.autoscale(True, axis='x', tight=True)
        axes.grid(True, which='major', linestyle='-', linewidth=0.5)
        axes.grid(True, which='minor', linestyle='-', linewidth=0.5)

class RangeCurve(CurveFigure):
    """
    A line plot, with additional range (e.g. confidence interval)
    """

    def add(self, xdata, ymean, ymax, ymin, xlabel='x', ylabel='y', title='x'):
        """
        Add a new subplot to the list of subplots to be created.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param ymean: y values of the data points.
        :type ymean: `numpy.ndarray`

        :param ymax: y values of the upper range data points.
        :type  ymax: `numpy.ndarray`

        :param ymin: y values of the lower range data points.
        :type  ymin: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """
        self.subfigures.append((xdata, ymean, ymax, ymin,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        """
        Draw a line plot on an :class:`matplotlib.axes` instance. The
        data and labels are contained in a tuple. This will plot a
        line plot and a grey range (for displaying uncertainty). A
        grid is added with a call to :meth:`RangeFigure.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.

        """

        xdata, ymean, ymax, ymin, xlabel, ylabel, title = subfigure

        axes.plot(xdata, ymean, lw=2)
        self.addRange(axes, xdata, ymin, ymax)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        axes.autoscale(True, axis='x', tight=True)
        self.addGrid(axes)

class RangeCompareCurve(CurveFigure):
    """
    A line plot, with additional range (e.g. confidence interval) and a second
    line for comparison. e.g. if you have model output (with an upper and lower
    estimate) and an observed value to compare against.

    """

    def add(self, xdata, y1, y2, y2max, y2min, xlabel='x',
            ylabel='y', title='x'):
        """
        Add a new subplot to the list of subplots to be created.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param y1: y values of the primary data points.
        :type  y1: `numpy.ndarray`
        :param y2: y values of the secondary data points.
        :type  y2: `numpy.ndarray`

        :param ymax: y values of the upper range data points.
        :type  ymax: `numpy.ndarray`

        :param ymin: y values of the lower range data points.
        :type  ymin: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """
        self.subfigures.append((xdata, y1, y2, y2max, y2min,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        """
        Draw a line plot on an :class:`matplotlib.axes` instance. The
        data and labels are contained in a tuple. This will plot two
        lines and a grey range (for displaying uncertainty). A grid is
        added with a call to :meth:`RangeFigure.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.
        """
        xdata, y1, y2, y2max, y2min, xlabel, ylabel, title = subfigure

        axes.plot(xdata, y1, color='r', lw=2, label="")
        axes.plot(xdata, y2, color='k', lw=2, label="")
        self.addRange(axes, xdata, y2min, y2max)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        axes.autoscale(True, axis='x', tight=True)
        self.addGrid(axes)

class SemilogRangeCurve(SemilogCurve):
    """
    A line plot on a semilog-x plot with additional range (e.g. confidence
    interval).

    """

    def add(self, xdata, ymean, ymax, ymin, xlabel, ylabel, title):
        """
        Add a new subplot to the list of subplots to be created.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param ymean: y values of the data points.
        :type ymean: `numpy.ndarray`

        :param ymax: y values of the upper range data points.
        :type  ymax: `numpy.ndarray`

        :param ymin: y values of the lower range data points.
        :type  ymin: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """

        self.subfigures.append((xdata, ymean, ymax, ymin,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        """
        Draw a line and range plot on an :class:`matplotlib.axes`
        instance, with a logarithmic scale on the x-axis. Data and
        labels are contained in a tuple. x-ticks presented as integer
        values. A grid is added with a call to
        :meth:`SemilogCurve.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.

        """
        xdata, ymean, ymax, ymin, xlabel, ylabel, title = subfigure

        axes.semilogx(xdata, ymean, lw=2, subsx=xdata)
        if (ymin[0] > 0) and (ymax[0] > 0):
            self.addRange(axes, xdata, ymin, ymax)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        ylim = (0., np.max([100, np.ceil(ymean.max()/10.)*10.]))
        axes.set_ylim(ylim)
        self.addGrid(axes)

class SemilogRangeScatterCurve(SemilogCurve):
    """
    A line plot on a semilog-x plot with additional range (e.g. confidence
    interval).

    """

    def add(self, xdata, events, ymean, ymax, ymin, xlabel, ylabel, title, fit):
        """
        Add a new subplot to the list of subplots to be created.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param ymean: y values of the data points.
        :type ymean: `numpy.ndarray`

        :param ymax: y values of the upper range data points.
        :type  ymax: `numpy.ndarray`

        :param ymin: y values of the lower range data points.
        :type  ymin: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """

        self.subfigures.append((xdata, events, ymean, ymax, ymin,
                                xlabel, ylabel, title, fit))

    def subplot(self, axes, subfigure):
        """
        Draw a line and range plot on an :class:`matplotlib.axes`
        instance, with a logarithmic scale on the x-axis. Data and
        labels are contained in a tuple. x-ticks presented as integer
        values. A grid is added with a call to
        :meth:`SemilogCurve.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.

        """
        xdata, events, ymean, ymax, ymin, xlabel, ylabel, title, fit = subfigure

        emprp = self.empReturnPeriod(events)
        log.debug("xvalues = {0} length".format(len(emprp)))
        log.debug("xvalues = {0}".format(emprp))

        axes.semilogx(xdata, ymean, lw=2, subsx=xdata, 
                      label = 'Fitted hazard curve ({0})'.format(fit))
        axes.scatter(emprp[emprp > 1], events[emprp > 1], s=100,
                color='r', label = 'Empirical ARI')
        axes.legend(loc = 2)

        if (ymin[0] > 0) and (ymax[0] > 0):
            self.addRange(axes, xdata, ymin, ymax)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        ylim = (0., np.max([100, np.ceil(ymean.max()/10.)*10.]))
        axes.set_ylim(ylim)
        self.addGrid(axes)

    def empReturnPeriod(self, data, npyr=365.25):
        """
        Returns the empirically-based recurrence interval (in years) for a set
        of observations.
        It is assumed the data are daily observations. If the observations are not
        daily, there are two options: set the ``npyr`` variable, or backfill the
        ``data`` variable with zero values to match the assumed length of the
        record.
        The highest return period should be (approximately) len(``data``)/``npyr``.
        :param data: :class:`numpy.ndarray` containing the observed values (with
                      missing values removed).
        :param float npy: Number of observations per year (default=365.25)
        :returns: Recurrence intervals for the observed data.
        :rtype: :class:`numpy.ndarray`
        """

        log.debug("Calculating xvalues for the scatter plots")
        nobs = len(data)
        # Empirical return periods:
        emprp = 1. / (1. - np.arange(1, nobs + 1, 1) / (nobs + 1)) / npyr
        return emprp

class SemilogRangeCompareCurve(SemilogCurve):
    """
    A line plot on a semilog-x plot with additional range (e.g. confidence
    interval) and a second line for comparison. e.g. if you have model
    output (with an upper and lower estimate) and an observed value to
    compare against.

    """
    def add(self, xdata, y1, y2, y2max, y2min, xlabel='x',
            ylabel='y', title='x'):
        """
        Add a new subplot to the list of subplots to be created.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param y1: y values of the primary data points.
        :type  y1: `numpy.ndarray`
        :param y2: y values of the secondary data points.
        :type  y2: `numpy.ndarray`

        :param ymax: y values of the upper range data points.
        :type  ymax: `numpy.ndarray`

        :param ymin: y values of the lower range data points.
        :type  ymin: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """
        self.subfigures.append((xdata, y1, y2, y2max, y2min,
                                xlabel, ylabel, title))

    def subplot(self, axes, subfigure):
        """
        Draw a range-compare plot on an :class:`matplotlib.axes`
        instance, with a logarithmic scale on the x-axis. Data and
        labels are contained in a tuple. x-ticks presented as integer
        values. A grid is added with a call to
        :meth:`RangeFigure.addGrid`.

        :param axes: :class:`matplotlib.axes` instance.
        :param tuple subfigure: Holds the data and labels to be added
                                to the subplot.

        """
        xdata, y1, y2, y2max, y2min, xlabel, ylabel, title = subfigure

        axes.semilogx(xdata, y1, color='r', lw=2, label="", subsx=xdata)
        axes.semilogx(xdata, y2, color='k', lw=2, label="", subsx=xdata)
        self.addRange(axes, xdata, y2min, y2max)
        ylim = (0., np.max([100, np.ceil(y2.max()/10.)*10.]))
        axes.set_ylim(ylim)

        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        self.addGrid(axes)

class HazardCurve(SemilogRangeCurve):

    def plot(self, years, wspd, wspdupper, wspdlower, xlabel, ylabel, title):
        """
        Plot a hazard curve, including upper and lower conffidence
        range estimates. This uses the :meth:`SemilogRangeCurve`
        method for plotting the curve on a semilog x-axis.
        :param years: Return period years to plot.
        :param wspd: Return period wind speed values.
        :param wspdupper: Upper confidence range wind speed values.
        :param wspdlower: Lower confidence range wind speed values.
        :type years: `numpy.ndarray`
        :type wspd: `numpy.ndarray`
        :type wspdupper: `numpy.ndarray`
        :type wspdlower: `numpy.ndarray`
        :param str xlabel: x-axis label.
        :param str ylabel: y-axis label.
        :param str title: Plot title.
        """

        self.add(years, wspd, wspdupper, wspdlower, xlabel, ylabel, title)
        super(HazardCurve, self).plot()

class HazardScatterCurve(SemilogRangeScatterCurve):

    def plot(self, years, events, wspd, wspdupper, wspdlower, xlabel, ylabel, title, fit):
        """
        Plot a hazard curve, including upper and lower conffidence
        range estimates. This uses the :meth:`SemilogRangeCurve`
        method for plotting the curve on a semilog x-axis.

        :param years: Return period years to plot.
        :param events: `numpy.ndarray`
        :param wspd: Return period wind speed values.
        :param wspdupper: Upper confidence range wind speed values.
        :param wspdlower: Lower confidence range wind speed values.
        :type years: `numpy.ndarray`
        :type wspd: `numpy.ndarray`
        :type wspdupper: `numpy.ndarray`
        :type wspdlower: `numpy.ndarray`

        :param str xlabel: x-axis label.
        :param str ylabel: y-axis label.
        :param str title: Plot title.
        :param str fit: Label for the legend containing the EVD fit type

        """

        self.add(years, events, wspd, wspdupper, wspdlower, xlabel, ylabel, title, fit)
        super(HazardScatterCurve, self).plot()

class DistributionCurve(RangeCompareCurve):

    def plot(self, x, y1, y2, y2max, y2min, xlabel, ylabel, title):
        """
        Plot a distribution curve, which compares two curves, plus an
        upper/lower confidence range.

        :param xdata: x values of the data points.
        :type xdata: `numpy.ndarray`

        :param y1: y values of the primary data points.
        :type  y1: `numpy.ndarray`
        :param y2: y values of the secondary data points.
        :type  y2: `numpy.ndarray`

        :param ymax: y values of the upper range data points.
        :type  ymax: `numpy.ndarray`

        :param ymin: y values of the lower range data points.
        :type  ymin: `numpy.ndarray`

        :param str xlabel: Label for the x-axis
        :param str ylabel: Label for the y-axis
        :param str title: Plot title

        """
        self.add(x, y1, y2, y2max, y2min, xlabel, ylabel, title)
        super(DistributionCurve, self).plot()


def saveFigure(figure, filename):
    """
    Save a figure canvas to file.

    :param figure: :class:`Figure` instance.
    :param str filename: Path to the location to store the image.
    """

    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename, bbox_inches='tight', dpi=300)

def saveHazardCurve(years, events, wspd, wspdupper, wspdlower,
                    xlabel, ylabel, title, filename, fit):
    """
    Plot a hazard curve, including upper and lower conffidence
    range estimates. This uses the :meth:`SemilogRangeCurve`
    method for plotting the curve on a semilog x-axis.

    :param years: Return period years to plot.
    :param events: `numpy.ndarray`
    :param wspd: Return period wind speed values.
    :param wspdupper: Upper confidence range wind speed values.
    :param wspdlower: Lower confidence range wind speed values.
    :type years: `numpy.ndarray`
    :type wspd: `numpy.ndarray`
    :type wspdupper: `numpy.ndarray`
    :type wspdlower: `numpy.ndarray`

    :param str xlabel: x-axis label.
    :param str ylabel: y-axis label.
    :param str title: Plot title.
    :param str filename: Path to save teh figure to.
    :param str fit: Label for the legend containing the EVD fit type

    """

    fig = HazardScatterCurve()
    fig.plot(years, events, wspd, wspdupper, wspdlower, xlabel, ylabel, title, fit)
    saveFigure(fig, filename)

def saveDistributionCurve(x, y1, y2, y2max, y2min,
                          xlabel, ylabel, title, filename):
    """
    Plot and save a distribution curve, which compares two curves,
    plus an upper/lower confidence range.

    :param xdata: x values of the data points.
    :type xdata: `numpy.ndarray`

    :param y1: y values of the primary data points.
    :type  y1: `numpy.ndarray`
    :param y2: y values of the secondary data points.
    :type  y2: `numpy.ndarray`

    :param ymax: y values of the upper range data points.
    :type  ymax: `numpy.ndarray`

    :param ymin: y values of the lower range data points.
    :type  ymin: `numpy.ndarray`

    :param str xlabel: Label for the x-axis
    :param str ylabel: Label for the y-axis
    :param str title: Plot title
    :param str filename: Path to save the figure to

    """
    fig = DistributionCurve()
    fig.plot(x, y1, y2, y2max, y2min, xlabel, ylabel, title)
    saveFigure(fig, filename)

