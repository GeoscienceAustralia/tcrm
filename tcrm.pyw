#!/usr/bin/env pythonw
"""
TCRM User Interface
"""
import Tkinter as tk
import numpy as np
import subprocess
import threading
import logging
import time
import json
import ttk
import sys
import os
import io

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure as MatplotlibFigure
from matplotlib.widgets import RectangleSelector
from Utilities.config import ConfigParser
from matplotlib.patches import Rectangle
from mpl_toolkits.basemap import Basemap
from os.path import join as pjoin, isdir
from collections import OrderedDict
from itertools import tee, izip
from Queue import Queue, Empty
from contextlib import closing

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

# tk or ttk?
Frame = ttk.Frame
Labelframe = ttk.Labelframe
Entry = ttk.Entry
Scale = ttk.Scale
Label = ttk.Label
Notebook = ttk.Notebook
Combobox = ttk.Combobox
Checkbutton = ttk.Checkbutton
Text = tk.Text
Scrollbar = tk.Scrollbar
Button = ttk.Button
PanedWindow = ttk.PanedWindow
Canvas = tk.Canvas

ON_POSIX = 'posix' in sys.builtin_module_names
POLL_INTERVAL = 1000  # ms
TEXT_FONT = ('Helvetica', '10')

ICON = """
R0lGODlhIAAgAPcAAAAAAAAAMwAAZgAAmQAAzAAA/wArAAArMwArZgArmQArzAAr/wBVAA
BVMwBVZgBVmQBVzABV/wCAAACAMwCAZgCAmQCAzACA/wCqAACqMwCqZgCqmQCqzACq/wDV
AADVMwDVZgDVmQDVzADV/wD/AAD/MwD/ZgD/mQD/zAD//zMAADMAMzMAZjMAmTMAzDMA/z
MrADMrMzMrZjMrmTMrzDMr/zNVADNVMzNVZjNVmTNVzDNV/zOAADOAMzOAZjOAmTOAzDOA
/zOqADOqMzOqZjOqmTOqzDOq/zPVADPVMzPVZjPVmTPVzDPV/zP/ADP/MzP/ZjP/mTP/zD
P//2YAAGYAM2YAZmYAmWYAzGYA/2YrAGYrM2YrZmYrmWYrzGYr/2ZVAGZVM2ZVZmZVmWZV
zGZV/2aAAGaAM2aAZmaAmWaAzGaA/2aqAGaqM2aqZmaqmWaqzGaq/2bVAGbVM2bVZmbVmW
bVzGbV/2b/AGb/M2b/Zmb/mWb/zGb//5kAAJkAM5kAZpkAmZkAzJkA/5krAJkrM5krZpkr
mZkrzJkr/5lVAJlVM5lVZplVmZlVzJlV/5mAAJmAM5mAZpmAmZmAzJmA/5mqAJmqM5mqZp
mqmZmqzJmq/5nVAJnVM5nVZpnVmZnVzJnV/5n/AJn/M5n/Zpn/mZn/zJn//8wAAMwAM8wA
ZswAmcwAzMwA/8wrAMwrM8wrZswrmcwrzMwr/8xVAMxVM8xVZsxVmcxVzMxV/8yAAMyAM8
yAZsyAmcyAzMyA/8yqAMyqM8yqZsyqmcyqzMyq/8zVAMzVM8zVZszVmczVzMzV/8z/AMz/
M8z/Zsz/mcz/zMz///8AAP8AM/8AZv8Amf8AzP8A//8rAP8rM/8rZv8rmf8rzP8r//9VAP
9VM/9VZv9Vmf9VzP9V//+AAP+AM/+AZv+Amf+AzP+A//+qAP+qM/+qZv+qmf+qzP+q///V
AP/VM//VZv/Vmf/VzP/V////AP//M///Zv//mf//zP///wAAAAAAAAAAAAAAACH5BAEAAP
wALAAAAAAgACAAAAj/AB9UqCCQ4MCCCA8qNMiw4MAKDg5SMEjB4cOLGDPu28ix48aMIEF6
HBktpMmHI0eevGhhIIcKHFJ6hAnzJQebNjfQLAKkyA8LRSp0rBCUww8jR5PyXNqzpxGBRo
wUidrRCBCpWKdezer0B8KKP37A6fjGSFk1cIyoMWsWCJy0DRcK7Pj27Ru7dX8QjPigYgUK
fiNC7CiQL8EiHoP6PehQ4IOOhitCS6lMYMWJgSEK5ZhwrEw4juMe7Hh5sMyNi0vzpUD4IG
LODh5zHGj46+aNhj3v0/RQ075o+/QaLChY9kbiuu2sflB7omaDDkgX1M0bY9/hyytChvjg
DUc4FWPHz77c9wHz8cwhl/e+UVPE99zfO5gfPzrs+cnPz6c/X//++ZDtx94+mfQn3n8I/h
fgfAPCkeCDCXaEwHwIDFjgfDiE9QOCDUyIgIeQHTDfGBwtM+EBHeHggIcINCDiAS52dECL
LXakDA7AcfRihwjM2OGMMs54AIqnRSNkjzAO+WKQPxKZUgBJujgkkkMGqeQBAYxEz5VUXo
llJpqACaWXAByQSSZpeHkllGWOeUCbBgwJ5ZhzYvmmnW2WCUAABuxpAJ+AAtDnn4Pyqeeb
gBoQEAA7
"""

DEFAULTS = """
[Actions]
DownloadData=True
DataProcess=True
ExecuteStat=True
ExecuteTrackGenerator=True
ExecuteWindfield=True
PlotData=True

[DataProcess]
InputFile=input/Allstorms.ibtracs_wmo.v03r05.csv
Source=IBTRACS
StartSeason=1981

[Input]
MSLPGrid=MSLP/mslp_annual_climatology.nc

[Region]
gridLimit={'xMin':112.0,'xMax':160.0,'yMin':-47.0,'yMax':16.0}

[StatInterface]
kdeType=Biweight
kde2DType=Gaussian
kdeStep=0.2
gridSpace={'x':1.0,'y':1.0}
gridInc={'x':1.0,'y':0.5}
minSamplesCell=40

[TrackGenerator]
NumSimulations=1024
YearsPerSimulation=1
gridSpace={'x':1.0,'y':1.0}
gridInc={'x':1.0,'y':0.5}
SeasonSeed=2
TrackSeed=1111

[WindfieldInterface]
NumberofFiles=1024
TrackPath=output/tracks
Margin=2.0
Resolution=0.05
Source=TCRM
profileType=powell
windFieldType=kepert
beta=1.5
beta1=1.5
beta2=1.4
thetaMax=70.0

[Output]
Path=output

[Logging]
LogFile=output/log/default.log
LogLevel=INFO
Verbose=True
ProgressBar=False

[Process]
DatFile=output/process/dat/default.dat
ExcludePastProcessed=False

[RMW]
GetRMWDistFromInputData=False
mean=50.0
sigma=0.6

[TCRM]
Columns=index,age,lon,lat,speed,bearing,pressure,penv,rmax
FieldDelimiter=,
NumberOfHeadingLines=1
SpeedUnits=kph
PressureUnits=hPa
"""

CHOICES = {
    'dataSources': ['IBTRACS'],
    'kdeKernel': ['Biweight', 'Epanechnikov', 'Triangular', 'Gaussian'],
    'kde2DKernel': ['Epanechnikov', 'Gaussian'],
    'profileType': ['Jelesnianski', 'Holland', 'Willoughby',
                    'Rankine', 'Schloemer', 'DoubleHolland', 'Powell',
                    'NewHolland'],
    'windfieldType': ['Hubbert', 'McConochie', 'Kepert']
}

json.encoder.FLOAT_REPR = lambda f: ('%.2f' % f)


def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def colorGenerator(saturation=0.5):
    def hsvToRgb(h, s, v):
        hi = int(h * 6)
        f = h * 6 - hi
        p = v * (1 - s)
        q = v * (1 - f * s)
        t = v * (1 - (1 - f) * s)
        if hi == 0:
            r, g, b = v, t, p
        if hi == 1:
            r, g, b = q, v, p
        if hi == 2:
            r, g, b = p, v, t
        if hi == 3:
            r, g, b = p, q, v
        if hi == 4:
            r, g, b = t, p, v
        if hi == 5:
            r, g, b = v, p, q
        return (int(r * 256), int(g * 256), int(b * 256))

    golden = 0.618033988749895
    h = 0.1

    while True:
        h += golden
        h %= 1
        rgb = hsvToRgb(h, saturation, 0.90)
        txt = '#%02x%02x%02x' % rgb
        yield txt


class Observable(object):

    """
    Base class for observable objects.
    """

    def __init__(self):
        """
        Initialise the class.
        """
        self.callbacks = []

    def addCallback(self, callback):
        """
        Add a function to be called when the 'notify' event occurs.
        """
        self.callbacks.append(callback)

    def notify(self, value):
        """
        Notify all subscribed callback functions. This function is
        typically called by subclasses of this object when they wish
        to notify that their state has changed. All registered
        callbacks are called with 'value' as their argument.
        """
        for callback in self.callbacks:
            callback(value)


class ObservableDict(Observable):

    """
    An observable dictionary object.
    """

    def __init__(self, theDict=None):
        """
        Construct an observable dictionary from `theDict`. If `theDict`
        is not provided then a new `dict` is created.
        """
        super(ObservableDict, self).__init__()
        self.theDict = theDict or {}

    def __setitem__(self, key, value):
        """
        Set an item in the dict with `dict[key] = value`.
        """
        prev = self.theDict.get(key)
        if value != prev:
            self.theDict[key] = value
            self.notify(self)

    def __getitem__(self, key):
        """
        Get an item in the dict with `dict[key]`.
        """
        return self.theDict.get(key)

    def __delitem__(self, key):
        del self.theDict[key]

    def __repr__(self):
        """
        Return a string representation of the dictionary.
        """
        return self.theDict.__repr__()

    def items(self):
        return self.theDict.items()

    def keys(self):
        return self.theDict.keys()


class ObservableVariable(Observable):

    """
    An observable variable. The `variable` is internally stored as a
    tk.StringVar object so that it can be used by Tk controls. Basic
    python objects are serialised to JSON format so that they can be
    easily edited through Tk controls as well.
    """

    def __init__(self, value):
        """
        Construct a variable with initial value `value`.
        """
        super(ObservableVariable, self).__init__()
        self.variable = tk.StringVar()
        self.set(value)
        self.variable.trace('w', self.changed)
        self.delay = 1000
        self.noPendingNotification = True

    def changed(self, *args):
        """
        Delays the 'notify' event by `self.delay` milliseconds. This is
        automatically called by the superclass when the variable is
        changed.
        """
        if self.noPendingNotification:
            self.after_idle(self.doNotify)  # FIXME
            self.noPendingNotification = False

    def doNotify(self, *args):
        """
        Perform a notify event. Called by the `changed` method to
        trigger the actual `notify` event.
        """
        try:
            self.notify(self.get())
        except ValueError:
            pass
        self.noPendingNotification = True

    def set(self, value):
        """
        Set the value of the variable.
        """
        if isinstance(value, str):
            text = value
        else:
            text = json.dumps(value)
        prev = self.variable.get()
        if text != prev:
            self.variable.set(text)

    def get(self):
        """
        Get the value of the variable.
        """
        return json.loads(self.variable.get())


class View(object, Frame):

    """
    An abstract view.
    """

    def __init__(self, parent, **kwargs):
        Frame.__init__(self, parent)


class GriddedView(object):

    def __init__(self, parent, **kwargs):
        self.parent = parent

    def grid(self, *args, **kwargs):
        column, row = kwargs['column'], kwargs['row']
        self.label.grid(column=column, row=row, sticky='W', padx=2, pady=2)
        self.control.grid(column=column + 1, row=row, sticky='E', padx=2,
                          pady=2)

    def rowconfigure(self, *args, **kwargs):
        row = kwargs.pop('row')
        self.parent.rowconfigure(row=row, **kwargs)

    def after_idle(self, func):
        self.parent.after_idle(func)


class LabeledView(object, Labelframe):

    """
    A labeled view.
    """

    def __init__(self, parent, name='', **kwargs):
        Labelframe.__init__(self, parent, text=name)


class ObservableEntry(Frame, ObservableVariable):

    """
    An observable entry box control with an associated label
    to the left of the control.

    :type  parent: object
    :param parent: the parent gui control.

    :type  name: str
    :param name: the name of the variable to put in the text label.

    :type  value: str
    :param value: the initial value to put in the entry box.
    """

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')
        value = kwargs.pop('value', '')

        Frame.__init__(self, parent)
        ObservableVariable.__init__(self, value=value)

        self.label = Label(self, text=name + ':')
        self.label.grid(column=0, row=0, sticky='W', padx=2)

        self.entry = Entry(self, **kwargs)
        self.entry.grid(column=1, row=0, sticky='E', padx=2)
        self.entry.config(textvariable=self.variable)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)


class ObservableSpinbox(Frame, ObservableVariable):

    """
    An observable spin box control with an associated label
    to the left of the control.
    """

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')
        width = kwargs.pop('width', 10)
        lower = kwargs.pop('lower', 0)
        upper = kwargs.pop('upper', 100)
        value = kwargs.pop('value', 0)

        Frame.__init__(self, parent)
        ObservableVariable.__init__(self, value)

        self.label = Label(self, width=width, text=name + ':')
        self.label.grid(column=0, row=0, sticky='NSEW')

        self.spinbox = tk.Spinbox(self, from_=lower, to=upper, increment=10)
        self.spinbox.grid(column=1, row=0, sticky='NSEW')
        self.spinbox.config(textvariable=self.variable)

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)


class ObservableScale(Frame, ObservableVariable):

    """
    An observable entry box and scale control with an associated label
    to the left.

    :type  parent: object
    :param parent: the parent gui control.

    :type  name: str
    :param name: the name of the variable to put in the text label.

    :type  value: int
    :param value: the initial value to put in the entry box.

    :type  lower: int
    :param lower: the lower value limit for the scale control.

    :type  upper: int
    :param upper: the upper value limit for the scale control.

    :type  delay: int
    :param delay: the amount of time to delay the 'notify' event.
    """

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')
        lower = kwargs.pop('lower', 0)
        upper = kwargs.pop('upper', 100)
        value = kwargs.pop('value', 0)
        delay = kwargs.pop('delay', 1000)

        Frame.__init__(self, parent)
        ObservableVariable.__init__(self, value)

        self.label = Label(self, text=name + ':')
        self.label.grid(column=0, row=0, sticky='EW', padx=2)

        self.value = Entry(self, textvariable=self.variable, **kwargs)
        self.value.grid(column=1, row=0, sticky='W', padx=2)

        self.scale = Scale(self, from_=lower, to=upper, value=value,
                           length=50)
        self.scale.grid(column=2, row=0, sticky='EW', padx=2)
        self.scale.config(variable=self.variable, command=self.format)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.rowconfigure(0, weight=1)

        self.noPendingNotification = True
        self.delay = delay

    def format(self, *args):
        """
        Format the variable by rounding it to 2 decimal places. This is
        called automatically when the scale control changes the
        variable. Subclasses can overwrite this to provide custom
        formatting.
        """
        self.set(round(self.get(), 2))

    def changed(self, *args):
        """
        Delays the 'notify' event by `self.delay` milliseconds. This is
        automatically called by the superclass when the variable is
        changed.
        """
        if self.noPendingNotification:
            self.after(self.delay, self.doNotify)
            self.noPendingNotification = False

    def doNotify(self, *args):
        """
        Perform a notify event. Called by the `changed` method to
        trigger the actual `notify` event.
        """
        try:
            self.notify(self.get())
        except ValueError:
            pass
        self.noPendingNotification = True


class ObservableCombobox(Frame, ObservableVariable):

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')
        width = kwargs.pop('width', 15)
        values = kwargs.pop('values', ['a', 'b', 'c'])

        Frame.__init__(self, parent)
        ObservableVariable.__init__(self, value='')

        self.label = Label(self, width=width, text=name + ':')
        self.label.grid(column=0, row=0, sticky='W')

        self.combo = Combobox(self, width=width, values=values)
        self.combo.grid(column=1, row=0, sticky='E')
        self.combo.config(textvariable=self.variable)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)

    def grid(self, *args, **kwargs):
        Frame.grid(self, *args, **kwargs)


class GriddedCombobox(GriddedView, ObservableVariable):

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')

        GriddedView.__init__(self, parent)
        ObservableVariable.__init__(self, value='')

        self.label = Label(parent, text=name + ':')
        self.control = Combobox(parent, textvariable=self.variable, width=9,
                                state='readonly', **kwargs)

    def setValues(self, values):
        self.control.config(values=values)


class GriddedEntry(GriddedView, ObservableVariable):

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')
        value = kwargs.pop('value', '')

        GriddedView.__init__(self, parent)
        ObservableVariable.__init__(self, value=value)

        self.label = Label(parent, text=name + ':')
        self.control = Entry(parent, textvariable=self.variable, width=11,
                             **kwargs)


class ObservableCheckbutton(Frame, ObservableVariable):

    """
    An observable check button control.

    :type  parent: object
    :param parent: the parent gui control.

    :type  name: str
    :param name: the name of the variable to put in the text label.

    :type  value: bool
    :param value: the initial state.
    """

    def __init__(self, parent, **kwargs):
        value = kwargs.pop('value', True)

        Frame.__init__(self, parent)
        ObservableVariable.__init__(self, value=value)

        self.button = Checkbutton(self, **kwargs)
        self.button.grid(column=1, row=0, sticky='E', padx=2)
        if not ('offvalue' in kwargs and 'onvalue' in kwargs):
            self.button.config(offvalue=False, onvalue=True)
        self.button.config(variable=self.variable)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)


class MapView(Frame):

    """
    A Map view control.
    """

    def __init__(self, parent, **kwargs):
        drawCoastlines = kwargs.pop('drawCoastlines', True)
        fillContinents = kwargs.pop('fillContinents', True)
        continentColor = kwargs.pop('continentColor', '0.8')
        coastlineWidth = kwargs.pop('coastlineWidth', 0.3)
        projection = kwargs.pop('projection', 'mill')
        resolution = kwargs.pop('resolution', 'c')
        figSize = kwargs.pop('figSize', (1.3, 1))

        Frame.__init__(self, parent)
        figure = MatplotlibFigure(figsize=figSize)
        self.canvas = FigureCanvasTkAgg(figure, master=self)

        self.axes = figure.add_axes([0, 0, 1, 1], aspect='auto')

        self.basemap = Basemap(llcrnrlon=0, llcrnrlat=-80, urcrnrlon=360,
                               urcrnrlat=80, projection=projection,
                               resolution=resolution,
                               ax=self.axes, fix_aspect=False)
        if drawCoastlines:
            self.basemap.drawcoastlines(linewidth=coastlineWidth)

        if fillContinents:
            self.basemap.fillcontinents(color=continentColor)

        widget = self.canvas.get_tk_widget()
        widget.config(highlightthickness=0, relief=tk.GROOVE, borderwidth=1)
        widget.pack(fill='both', expand=False)


class MapRegionSelector(MapView, ObservableVariable):

    """
    A map region selection widget.
    """

    def __init__(self, parent, **kwargs):
        value = kwargs.pop('value', 0)
        self.regionFaceColor = kwargs.pop('regionFaceColor', 'cyan')
        self.regionEdgeColor = kwargs.pop('regionEdgeColor', 'black')
        self.regionAlpha = kwargs.pop('regionAlpha', 0.5)

        MapView.__init__(self, parent, **kwargs)
        ObservableVariable.__init__(self, value=value)

        self.selector = RectangleSelector(self.axes, self.onSelected,
                                          drawtype='box', useblit=True,
                                          minspanx=5, minspany=5,
                                          spancoords='pixels')
        self.addCallback(self.onValueChanged)
        self.region = None

    def drawRegion(self, xMin, xMax, yMin, yMax):
        x0, y0 = self.basemap(xMin, yMin)
        x1, y1 = self.basemap(xMax, yMax)
        w, h = (x1 - x0), (y1 - y0)
        if self.region is not None:
            self.region.remove()
        self.region = Rectangle((x0, y0), w, h, ec=self.regionEdgeColor,
                                fc=self.regionFaceColor,
                                alpha=self.regionAlpha)
        self.axes.add_patch(self.region)
        self.canvas.draw()

    def onValueChanged(self, value):
        limits = self.get()
        xMin = limits['xMin']
        xMax = limits['xMax']
        yMin = limits['yMin']
        yMax = limits['yMax']
        self.drawRegion(xMin, xMax, yMin, yMax)

    def onSelected(self, click, release):
        x1, y1 = click.xdata, click.ydata
        x2, y2 = release.xdata, release.ydata
        xMin, yMin = self.basemap(min(x1, x2), min(y1, y2), inverse=True)
        xMax, yMax = self.basemap(max(x1, x2), max(y1, y2), inverse=True)
        self.set({'xMin': xMin, 'xMax': xMax, 'yMin': yMin, 'yMax': yMax})


class DictEntry(Labelframe, ObservableVariable):

    """
    A control for editing dictionaries.
    """

    def __init__(self, parent, **kwargs):
        name = kwargs.pop('name', '')
        columns = kwargs.pop('columns', 2)
        keys = kwargs.pop('keys', None)

        Labelframe.__init__(self, parent, text=name)
        ObservableVariable.__init__(self, value='')

        self.entries = {key: ObservableEntry(self, name=key, **kwargs)
                        for key in keys}
        for i, key in enumerate(keys):
            c, r = i % columns, i / columns
            entry = self.entries[key]
            entry.addCallback(self.changed)
            entry.grid(row=r, column=c, sticky='N', padx=2)
            self.rowconfigure(r, weight=0)
            self.columnconfigure(c, weight=1)

        self.keys = keys
        self.allowNotify = False

    def changed(self, *args):
        if self.allowNotify:
            self.notify(self.get())

    def set(self, value):
        if isinstance(value, dict):
            self.allowNotify = False
            for key in self.keys:
                entry = self.entries[key]
                entry.set(value[key])
            self.allowNotify = True
            self.changed()

    def get(self):
        value = {}
        for key in self.keys:
            entry = self.entries[key]
            value[key] = entry.get()
        return value


class MapRegionGrid(MapView):

    """
    A widget that shows the current region and the current grid.
    """

    def __init__(self, parent, **kwargs):
        self.gridStyle = kwargs.pop('gridStyle', '-')
        self.gridColor = kwargs.pop('gridColor', '0.7')
        MapView.__init__(self, parent, **kwargs)
        self.region = None
        self.step = None

    def set(self, value):
        if isinstance(value, dict):
            region = value
            step = self.step
        else:
            region = self.region
            step = value
        if (self.region != region) or (self.step != step):
            if region:
                self.setRegionLimits(region)
                self.axes.grid(False)
                self.canvas.draw()
                self.region = region
            if step:
                self.setGridStep(step)
                self.step = step
                self.after_idle(self.canvas.draw)

    def setRegionLimits(self, limits):
        xMin = limits['xMin']
        xMax = limits['xMax']
        yMin = limits['yMin']
        yMax = limits['yMax']
        x0, y0 = self.basemap(xMin, yMin)
        x1, y1 = self.basemap(xMax, yMax)
        self.axes.set_xlim((x0, x1))
        self.axes.set_ylim((y0, y1))

    def setGridStep(self, step):
        if step < 0.5 or self.region is None:
            return
        xMin = self.region['xMin']
        xMax = self.region['xMax']
        yMin = self.region['yMin']
        yMax = self.region['yMax']
        xs = np.arange(min(xMin, xMax), max(xMin, xMax), step)
        ys = np.arange(min(yMin, yMax), max(yMin, yMax), step)
        xs, ys = self.basemap(*np.meshgrid(xs, ys))
        self.axes.xaxis.set_ticks(xs.flatten())
        self.axes.yaxis.set_ticks(ys.flatten())
        self.axes.xaxis.set_ticklabels([])
        self.axes.yaxis.set_ticklabels([])
        self.axes.grid(True, which='both', linestyle=self.gridStyle,
                       color=self.gridColor)


class LogView(logging.Handler, Frame):

    """
    Log view.
    """

    def __init__(self, parent, **kwargs):
        logging.Handler.__init__(self)
        Frame.__init__(self, parent)

        self.console = Text(self, **kwargs)
        self.console.config(state=tk.DISABLED, wrap='none')
        self.console.config(font=TEXT_FONT)

        self.scroll = Scrollbar(self, orient=tk.VERTICAL)
        self.scroll.config(command=self.console.yview)

        self.console.config(yscrollcommand=self.scroll.set)

        self.console.grid(column=0, row=0, sticky='NSEW', padx=2, pady=2)
        self.scroll.grid(column=0, row=0, sticky='NSEW')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.setFormatter(logging.Formatter('%(asctime)s: %(message)s\n'))
        log.addHandler(self)

    def emit(self, message):
        self.console.config(state=tk.NORMAL)
        self.console.insert(tk.END, self.format(message))
        self.console.config(state=tk.DISABLED)
        self.console.see(tk.END)


class RegionSettingsView(Frame):

    """
    Region settings.
    """

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.region = DictEntry(self, keys=('xMin', 'xMax', 'yMin', 'yMax'),
                                name='Region', width=5, justify='center')
        self.region.grid(column=0, row=0, padx=2, pady=2, sticky='N')

        self.step = ObservableScale(self, name='Grid step', width=5,
                                    justify='center')
        self.step.grid(column=0, row=1, padx=2, pady=2, sticky='N')

        # Locality

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)


class CalibrationSettingsView(Frame):

    """
    Model calibration settings.
    """

    def __init__(self, parent):
        Frame.__init__(self, parent)

        # Data

        self.source = GriddedCombobox(self, name='Data Source',
                                      justify='center')
        self.source.grid(column=0, row=0, sticky='N')

        self.start = GriddedEntry(self, name='Season Start Year',
                                  justify='center')
        self.start.grid(column=0, row=1, sticky='N')

        # Stats

        self.kdeKernel = GriddedCombobox(self, name='KDE Kernel',
                                         justify='center')
        self.kdeKernel.grid(column=0, row=2, sticky='N')

        self.kde2DKernel = GriddedCombobox(self, name='KDE 2D Kernel',
                                           justify='center')
        self.kde2DKernel.grid(column=0, row=3, sticky='N')

        self.kdeStep = GriddedEntry(self, name='KDE Step', justify='center')
        self.kdeStep.grid(column=0, row=4, sticky='N')

        self.minSamples = GriddedEntry(self, name='Min Cell Samples',
                                       justify='center')
        self.minSamples.grid(column=0, row=5, sticky='N')

        self.gridSpace = DictEntry(self, keys=('x', 'y'),
                                   name='Grid Space', width=5,
                                   justify='center')
        self.gridSpace.grid(column=0, columnspan=2, row=6, sticky='N')

        self.gridInc = DictEntry(self, keys=('x', 'y'), name='Grid Increment',
                                 width=5, justify='center')
        self.gridInc.grid(column=0, columnspan=2, row=7, sticky='N')

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        for i in range(7):
            self.rowconfigure(i, weight=0)


class TrackSettingsView(Frame):

    """
    Track settings.
    """

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.gridSpace = DictEntry(self, name='Grid Space', keys=('x', 'y'),
                                   justify='center', width=5)
        self.gridInc = DictEntry(self, name='Grid Increment', keys=('x', 'y'),
                                 justify='center', width=5)
        self.nSims = GriddedEntry(self, name='Simulations', justify='center')
        self.nYears = GriddedEntry(self, name='Years Per Sim',
                                   justify='center')
        self.seasonSeed = GriddedEntry(self, name='Season Seed',
                                       justify='center')
        self.trackSeed = GriddedEntry(self, name='Track Seed',
                                      justify='center')

        for i, control in enumerate([self.nSims, self.nYears, self.seasonSeed,
                                     self.trackSeed]):
            control.grid(column=0, row=i)
            self.rowconfigure(i, weight=0)

        self.gridSpace.grid(column=0, columnspan=2, row=5, sticky='N')
        self.gridInc.grid(column=0, columnspan=2, row=6, sticky='N')

        self.rowconfigure(5, weight=0)
        self.rowconfigure(6, weight=0)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)


class WindfieldSettingsView(Frame):

    """
    Wind field settings.
    """

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.margin = GriddedEntry(self, name='Margin', justify='center')
        self.margin.grid(column=0, row=1, sticky='N')

        self.resolution = GriddedEntry(self, name='Resolution',
                                       justify='center')
        self.resolution.grid(column=0, row=1, sticky='N')

        self.profileType = GriddedCombobox(self, name='Profile Type',
                                           justify='center')
        self.profileType.grid(column=0, row=2, sticky='N')

        self.windfieldType = GriddedCombobox(self, name='Windfield Type',
                                             justify='center')
        self.windfieldType.grid(column=0, row=3, sticky='N')

        self.modelParams = GriddedEntry(self, name='Model Parameters',
                                        justify='center')
        self.modelParams.grid(column=0, row=4, sticky='N')

        self.thetaMax = GriddedEntry(self, name='Max Theta',
                                     justify='center')
        self.thetaMax.grid(column=0, row=5, sticky='N')

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        for i in range(5):
            self.rowconfigure(i, weight=0)


class StageProgressView(Canvas, ObservableVariable):

    """
    Stage progress view.
    """

    def __init__(self, parent, stages):
        bgColor = '#%02x%02x%02x' % \
            tuple([c / 255 for c in
                   parent.winfo_rgb('SystemButtonFace')])

        Canvas.__init__(self, parent, bg=bgColor, height=30, width=250,
                        borderwidth=0, highlightthickness=0,
                        selectborderwidth=0)
        ObservableVariable.__init__(self, value='')

        self.stages = OrderedDict()
        self.colors = []
        for stage, color in zip(stages, colorGenerator()):
            self.stages[stage] = ObservableVariable(False)
            self.colors.append(color)

        self.bind("<Configure>", self.resize)

    def resize(self, event):
        self.render(event.width, event.height)

    def render(self, w, h):
        self.delete('all')
        n = len(self.stages)
        xs = np.linspace(5, w - 1 - 10, n + 1)
        stages = self.stages.keys()
        for stage, pts, color in zip(stages, pairwise(xs), self.colors):
            x0, x1 = pts
            poly = [x0 + 5, h / 2, x0, 0, x1, 0, x1 + 5, h / 2, x1, h - 1, x0,
                    h - 1]
            self.create_polygon(poly, fill=color)
            self.create_text((x0 + x1) / 2 + 2.5, h / 2, text=stage,
                             font=TEXT_FONT)
        self.addtag_all('all')


class SubprocessOutputView(Frame, Observable):

    def __init__(self, parent, **kwargs):
        Frame.__init__(self, parent)
        Observable.__init__(self)

        self.console = Text(self, **kwargs)
        self.console.config(state=tk.DISABLED, wrap='none')
        self.console.config(font=TEXT_FONT)
        self.console.config(yscrollcommand=self.onScroll)

        self.console.tag_config('important', background='yellow')

        self.scroll = Scrollbar(self, orient=tk.VERTICAL)
        self.scroll.config(command=self.console.yview)

        self.console.grid(column=0, row=0, sticky='NSEW', padx=2, pady=2)
        self.scroll.grid(column=0, row=0, sticky='NSE')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.lockScroll = True
        self.alarm = None

        self.poll()

    def onScroll(self, *args):
        if args[1] == '1.0':
            self.lockScroll = True
        else:
            self.lockScroll = False
        self.scroll.set(*args)

    def emit(self, message, *tags):
        self.console.config(state=tk.NORMAL)
        self.console.insert(tk.END, message, *tags)
        self.console.config(state=tk.DISABLED)
        if self.lockScroll:
            self.console.see(tk.END)

    def clear(self):
        self.console.config(state=tk.NORMAL)
        self.console.delete(1.0, tk.END)
        self.console.config(state=tk.DISABLED)
        self.console.see(tk.END)

    def poll(self):
        self.notify(self)
        self.alarm = self.after(POLL_INTERVAL, lambda: self.poll())


class MainView(Frame):

    def __init__(self, parent):
        Frame.__init__(self)

        paned = PanedWindow(self, orient=tk.HORIZONTAL)
        paned.grid(column=0, row=0, sticky='NSEW')

        # Left

        frame = Frame(paned)
        paned.add(frame)

        self.region = MapRegionSelector(frame, figSize=(1.3, 1.2))
        self.region.grid(column=0, row=0, padx=2, pady=2, sticky='NEW')

        notebook = Notebook(frame)
        notebook.grid(column=0, row=1, sticky='N', padx=2, pady=2)

        # self.regionSettings = RegionSettingsView(notebook)
        self.calibration = CalibrationSettingsView(notebook)
        self.track = TrackSettingsView(notebook)
        self.wind = WindfieldSettingsView(notebook)

        # notebook.add(self.regionSettings, text='Region')
        notebook.add(self.calibration, text='Calibration')
        notebook.add(self.track, text='Tracks')
        notebook.add(self.wind, text='Windfields')

        self.stage = StageProgressView(frame,
                                       ['Calibration', 'Tracks', 'Windfields'])
        self.stage.grid(column=0, row=2, sticky='NWE', padx=2, pady=2)

        self.run = Button(frame, text='Start')
        self.run.grid(column=0, row=3, sticky='N', padx=2, pady=2)

        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=0)
        frame.rowconfigure(1, weight=0)
        frame.rowconfigure(2, weight=0)
        frame.rowconfigure(3, weight=0)

        # Right

        rightFrame = Frame(paned)
        paned.add(rightFrame)

        self.view = MapRegionGrid(rightFrame, figSize=(7, 7),
                                  continentColor='#cdcbc1',
                                  coastlineWidth=0.8,
                                  resolution='i')
        self.view.grid(column=0, row=0, sticky='NSEW')

        self.output = SubprocessOutputView(rightFrame, width=80, height=5)
        self.output.grid(column=0, row=1, sticky='NSEW')

        rightFrame.columnconfigure(0, weight=1)
        rightFrame.rowconfigure(0, weight=2)
        rightFrame.rowconfigure(1, weight=1)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)



class TropicalCycloneRiskModel(object):

    def __init__(self):
        self.cmd = [sys.executable, '-u', 'tcrm.py', '-v', 'default.ini']
        self.running = False
        self.output = Queue()
        self.newFiles = Queue()
        self.process = None
        self.monitorThread = None
        self.outputDirectories = []

    def run(self):
        """
        Start TCRM.
        """
        if (not self.running) and (self.process is None):
            # create the subprocess that runs TCRM and redirect stdout
            self.process = subprocess.Popen(' '.join(self.cmd),
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            shell=True,
                                            bufsize=1)
                                            
            # monitor the process by a thread
            self.monitorThread = threading.Thread(target=self.monitor)
            self.monitorThread.daemon = True
            self.monitorThread.start()

            self.running = True

    def monitor(self):
        """
        Monitor function. Run in a thread and continuously steps
        coroutines to check for output, new files, etc.

        The logic here is that if there is a new line of output from
        TCRM then there is a chance that there are some new output files
        created. This seems better than polling every x number of
        seconds.
        """
        lineOfOutput = self.enqueueOutput()
        additionOfFiles = self.enqueueFileAdditions()
        try:
            while True:
                next(lineOfOutput)
                next(additionOfFiles)
        except StopIteration:
            # the thread was probably killed so exit nicely
            pass

    def isAlive(self):
        return self.running

    def enqueueOutput(self):
        """
        Coroutine to enqueue new output.
        """
        with closing(self.process.stdout) as out:
            for line in iter(out.readline, b''):
                self.output.put(line)
                yield

        self.running = False

    def enqueueFileAdditions(self):
        """
        Coroutine to enqueue new files.
        """
        def files():
            lst = []
            for directory in self.outputDirectories:
                if isdir(directory):
                    fns = os.listdir(directory)
                    lst.extend([pjoin(directory, fn) for fn in fns])
            return lst

        prev = files()

        while True:
            now = files()
            self.newFiles.put([f for f in now if not f in prev])
            prev = now
            yield

    def getOutput(self):
        """
        Get all new process output since last call.
        """
        lines = []
        while True:
            try:
                line = self.output.get_nowait()
                lines.append(line)
            except Empty:
                break
        return ''.join(lines)

    def getFileAdditions(self):
        """
        Get all new file additions since last call.
        """
        additions = []
        while True:
            try:
                files = self.newFiles.get_nowait()
                additions.extend(files)
            except Empty:
                break
        return additions

    def saveFlatConfig(self, flatConfig, destFile='default.ini'):
        with open(destFile, 'w') as dest:
            config = ConfigParser()

            for key, value in flatConfig.items():
                section, option = key.split('_')
                config.set(section, option, value)

            config.write(dest)

    def quit(self):
        """
        Nicely kill the TCRM process.
        """
        if self.process is None:
            return

        self.process.terminate()
        
        for i in range(5):
            if self.process.poll() is None:
                self.process.kill()
                time.sleep(0.5)
            else:
                break

        self.monitorThread.shutdown = True
        self.monitorThread.join()

        self.process = None
        self.running = False

        
class Controller(tk.Tk):

    """
    Controller.
    """

    def __init__(self):
        tk.Tk.__init__(self)
        self.withdraw()
        self.title('Tropical Cyclone Risk Model')
        self.tk.call('wm', 'iconphoto', self._w, tk.PhotoImage(data=ICON))

        self.protocol('WM_DELETE_WINDOW', self.onQuit)

        self.running = False
        self.view = MainView(self)
        self.tcrm = TropicalCycloneRiskModel()

        view = self.view

        view.run.config(command=self.toggleRun)
        view.output.addCallback(self.onWantOutput)
        view.output.addCallback(self.onCheckAlive)

        self.loadSettings()

        mappings = [
            (view.region, 'Region_gridlimit'),
            # (view.regionSettings.region, 'Region_gridlimit'),
            (view.view, 'Region_gridlimit'),
            (view.calibration.source, 'DataProcess_source'),
            (view.calibration.start, 'DataProcess_startseason'),
            (view.calibration.kdeKernel, 'StatInterface_kdetype'),
            (view.calibration.kde2DKernel, 'StatInterface_kde2dtype'),
            (view.calibration.kdeStep, 'StatInterface_kdestep'),
            (view.calibration.gridSpace, 'StatInterface_gridspace'),
            (view.calibration.gridInc, 'StatInterface_gridinc'),
            (view.calibration.minSamples, 'StatInterface_minsamplescell'),
            (view.track.gridSpace, 'TrackGenerator_gridspace'),
            (view.track.gridInc, 'TrackGenerator_gridinc'),
            (view.track.nSims, 'TrackGenerator_numsimulations'),
            (view.track.nYears, 'TrackGenerator_yearspersimulation'),
            (view.track.seasonSeed, 'TrackGenerator_seasonseed'),
            (view.track.trackSeed, 'TrackGenerator_trackseed'),
            (view.wind.margin, 'WindfieldInterface_margin'),
            (view.wind.resolution, 'WindfieldInterface_resolution'),
            (view.wind.profileType, 'WindfieldInterface_profiletype'),
            (view.wind.windfieldType, 'WindfieldInterface_windfieldtype'),
            (view.wind.thetaMax, 'WindfieldInterface_thetamax'),
        ]

        def makeCallback(key):
            return lambda x: self.onControlChanged(x, key)

        self.notifyWhom = {}
        for control, key in mappings:
            controls = self.notifyWhom.setdefault(key, [])
            controls.append(control)
            if hasattr(control, 'addCallback'):
                control.addCallback(makeCallback(key))

        comboboxes = [
            (view.calibration.source, 'dataSources'),
            (view.calibration.kdeKernel, 'kdeKernel'),
            (view.calibration.kde2DKernel, 'kde2DKernel'),
            (view.wind.profileType, 'profileType'),
            (view.wind.windfieldType, 'windfieldType')
        ]
        for control, key in comboboxes:
            values = CHOICES[key]
            control.setValues([v.lower() for v in values])

        self.onSettingsChanged(self.settings)
        self.settings.addCallback(self.onSettingsChanged)

        view.pack(fill='both', expand=True)
        self.after_idle(self.show)

    def loadSettings(self, configFile=None):
        config = ConfigParser()
        if configFile:
            config.read(configFile)
        else:
            config.readfp(io.BytesIO(DEFAULTS))

        outputPath = config.get('Output', 'Path')
        self.tcrm.outputDirectories.append(outputPath)

        flatConfig = {}
        for section in config.sections():
            for name, value in config.items(section):
                flatConfig['%s_%s' % (section, name)] = value

        self.settings = ObservableDict(flatConfig)

    def onSettingsChanged(self, settings):
        for key, value in settings.items():
            try:
                controls = self.notifyWhom[key]
                for control in controls:
                    control.set(value)
            except KeyError:
                pass

    def onControlChanged(self, value, key):
        self.settings[key] = value

    def onWantOutput(self, control):
        if self.running:
            output = self.tcrm.getOutput()
            control.emit(output)
            added = self.tcrm.getFileAdditions()
            for f in added:
                control.emit('new file: %s' % f)

    def onCheckAlive(self, control):
        if self.running:
            if not self.tcrm.isAlive():
                self.doFlushAndReset()

    def doFlushAndReset(self):
        if self.running:
            self.running = False
            self.view.run.config(text='Start')
            self.view.output.emit('TCRM STOPPED', ('important'))

    def toggleRun(self):
        if self.running:
            self.tcrm.quit()
            self.after(3, self.doFlushAndReset)
        else:
            self.view.output.clear()
            self.tcrm.saveFlatConfig(self.settings)
            self.tcrm.run()
            self.view.run.config(text='Stop')
            self.running = True

    def onQuit(self):
        self.tcrm.quit()
        self.destroy()

    def show(self):
        self.update()
        self.deiconify()

Controller().mainloop()
