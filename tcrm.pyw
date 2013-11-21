#!/usr/bin/env pythonw
"""
TCRM User Interface
"""
import matplotlib as mp
import Tkinter as tk
import numpy as np
import subprocess
import threading
import logging
import time
import json
import ttk
import sys

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure as MatplotlibFigure
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import Rectangle
from mpl_toolkits.basemap import Basemap
from Queue import Queue, Empty
from contextlib import closing

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

json.encoder.FLOAT_REPR = lambda f: ('%.2f' % f)

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

ON_POSIX = 'posix' in sys.builtin_module_names
POLL_INTERVAL = 1000

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
            self.after_idle(self.doNotify)
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
        self.label.grid(column=0, row=0, sticky='EW',  padx=2)

        self.value = Entry(self, textvariable=self.variable, **kwargs)
        self.value.grid(column=1, row=0, sticky='W', padx=2)

        self.scale = Scale(self, from_=lower, to=upper, value=value, length=50)
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
        width = kwargs.pop('width', 10)
        values = kwargs.pop('values', ['a', 'b', 'c'])

        Frame.__init__(self, parent)
        ObservableVariable.__init__(self, value='')

        self.label = Label(self, width=width, text=name + ':')
        self.label.grid(column=0, row=0, sticky='NSEW')

        self.combo = Combobox(self, values=values)
        self.combo.grid(column=1, row=0, sticky='NSEW')
        self.combo.config(textvariable=self.variable)

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=2)
        self.rowconfigure(0, weight=1)


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
        if not (kwargs.has_key('offvalue') and kwargs.has_key('onvalue')):
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
        figSize = kwargs.pop('figSize', (3, 1.3))
        bgColor = '#%02x%02x%02x' % \
                  tuple([c / 255 for c in
                         parent.winfo_rgb('SystemButtonFace')])

        Frame.__init__(self, parent)

        figure = MatplotlibFigure(figsize=figSize, facecolor=bgColor)

        self.axes = figure.add_subplot(111, aspect='auto')

        self.basemap = Basemap(llcrnrlon=0, llcrnrlat=-80, urcrnrlon=360,
                               urcrnrlat=80, projection=projection,
                               ax=self.axes)
        if drawCoastlines:
            self.basemap.drawcoastlines(linewidth=coastlineWidth)

        if fillContinents:
            self.basemap.fillcontinents(color=continentColor)

        self.axes.set_aspect('auto')
        figure.set_tight_layout(True)

        self.canvas = FigureCanvasTkAgg(figure, master=self)
        widget = self.canvas.get_tk_widget()
        widget.config(highlightthickness=0)
        widget.config(background=bgColor, relief=tk.GROOVE, borderwidth=1)
        widget.grid(column=0, row=0, sticky='NSEW')

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)


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
            entry.grid(row=r, column=c, sticky='EW', padx=2, pady=2)
            self.rowconfigure(r, weight=1)
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

    def __init__(self, parent, **kwargs):
        logging.Handler.__init__(self)
        Frame.__init__(self, parent)

        self.console = Text(self, **kwargs)
        self.console.config(state=tk.DISABLED, wrap='none')
        self.console.config(font=('Helvetica', 8))

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


class GridSettingsView(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.region = DictEntry(self, keys=('xMin', 'xMax', 'yMin', 'yMax'),
                                name='Region', width=5, justify='center')
        self.region.grid(column=0, row=0, padx=2, pady=2, sticky='N')

        self.step = ObservableScale(self, name='Grid step', width=5,
                                    justify='center')
        self.step.grid(column=0, row=1, padx=2, pady=2, sticky='N')

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)


class TrackSettingsView(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.gridSpace = DictEntry(self, name='Grid space', keys=('x', 'y'),
                                   width=7, justify='center')
        self.gridInc = DictEntry(self, name='Grid increment', keys=('x', 'y'),
                                 width=7, justify='center')
        self.nSims = ObservableEntry(self, name='Simulations',
                                     width=7, justify='center')
        self.nYears = ObservableEntry(self, name='Years per sim',
                                      width=7, justify='center')
        self.seasonSeed = ObservableEntry(self, name='Season seed', width=7,
                                          justify='center')
        self.trackSeed = ObservableEntry(self, name='Track seed', width=7,
                                         justify='center')
        for i, control in enumerate([self.gridSpace, self.gridInc, self.nSims,
                                     self.nYears, self.seasonSeed,
                                     self.trackSeed]):
            control.grid(column=0, row=i, sticky='WE', padx=4, pady=4)
            self.rowconfigure(i, weight=1)
        self.columnconfigure(0, weight=1)


class StageProgressViewOld(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.calibrate = ObservableCheckbutton(self, text='Calibrate',
                                               takefocus=False)
        self.trackSim = ObservableCheckbutton(self, text='Track simulation',
                                              takefocus=False)
        self.windFields = ObservableCheckbutton(self, text='Windfields',
                                                takefocus=False)

        self.calibrate.grid(column=0, row=0, sticky='NW', padx=4, pady=4)
        self.trackSim.grid(column=0, row=1, sticky='NW', padx=4, pady=4)
        self.windFields.grid(column=0, row=2, sticky='NW', padx=4, pady=4)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)


class StageProgressView(Canvas):

    def __init__(self, parent):
        Canvas.__init__(self, parent, width=200, height=100, borderwidth=0)
        self.bind("<Configure>", self.resize)
        self.config(offset='3,3')

        #self.render(200, 100)

    def resize(self, event):
        w, h = event.width, event.height
        print(w,h)
        self.render(w, h)

    def render(self, w, h):
        x = 3
        self.create_rectangle(x, x, w-x, h-x)


class SubprocessOutputView(Frame,Observable):

    def __init__(self, parent, **kwargs):
        Frame.__init__(self, parent)
        Observable.__init__(self)

        self.console = Text(self, **kwargs)
        self.console.config(state=tk.DISABLED, wrap='none')
        self.console.config(font=('Helvetica', 8))
        self.console.config(yscrollcommand=self.onScroll)

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

    def emit(self, message):
        self.console.config(state=tk.NORMAL)
        self.console.insert(tk.END, message)
        self.console.config(state=tk.DISABLED)
        if self.lockScroll:
            self.console.see(tk.END)

    def poll(self):
        self.notify(self)
        self.alarm = self.after(POLL_INTERVAL, lambda: self.poll())

class MainView(Frame):

    def __init__(self, parent):
        Frame.__init__(self)

        paned = PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True)

        # Left

        frame = Frame(paned)
        paned.add(frame)

        self.region = MapRegionSelector(frame, figSize=(2, 1.3))
        self.region.grid(column=0, row=0, padx=2, pady=2, sticky='NEW')

        notebook = Notebook(frame)
        notebook.grid(column=0, row=1, sticky='N', padx=2, pady=2)

        self.gridSettings = GridSettingsView(notebook)
        self.track = TrackSettingsView(notebook)

        notebook.add(self.gridSettings, text='Region')
        notebook.add(self.track, text='Simulation')

        self.stage = StageProgressView(frame)
        self.stage.grid(column=0, row=2, sticky='N', padx=2, pady=2)

        self.run = Button(frame, text='Run')
        self.run.grid(column=0, row=3, sticky='WN', padx=2, pady=2)

        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=0)
        frame.rowconfigure(1, weight=0)
        frame.rowconfigure(2, weight=0)
        frame.rowconfigure(3, weight=0)

        # Right

        right = Frame(paned)
        paned.add(right)

        rightPane = PanedWindow(right, orient=tk.VERTICAL)
        rightPane.pack(fill=tk.BOTH, expand=True)

        self.view = MapRegionGrid(rightPane, figSize=(7, 5),
                                  continentColor='#cdcbc1',
                                  coastlineWidth=0.8)
        rightPane.add(self.view)

        self.output = SubprocessOutputView(right, width=80, height=3)
        rightPane.add(self.output)


class TropicalCycloneRiskModel(object):

    def __init__(self):
        self.output = Queue()
        self.process = None
        self.monitorThread = None

    def run(self):
        if self.process is None:
            self.process = subprocess.Popen(['python', '-u', 'test.py'],
                                            stdout=subprocess.PIPE,
                                            shell=True,
                                            bufsize=1)
            self.monitorThread = threading.Thread(target=self.enqueueOutput)
            self.monitorThread.daemon = True
            self.monitorThread.start()

    def enqueueOutput(self):
        with closing(self.process.stdout) as out:
            for line in iter(out.readline, b''):
                self.output.put(line)

    def getOutput(self):
        lines = []
        while True:
            try:
                line = self.output.get_nowait()
                lines.append(line)
            except Empty:
                break
        return ''.join(lines)

    def quit(self):
        if self.process is None:
            return True

        self.process.terminate()

        for i in range(5):
            if self.process.poll() is None:
                self.process.kill()
                self.process.wait()
                time.sleep(1)
            else:
                break


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

        view = MainView(self)
        tcrm = TropicalCycloneRiskModel()

        view.run.config(command=self.onRun)
        view.output.addCallback(self.onWantOutput)

        self.settings = ObservableDict({
            'GRID_LIMITS': {
                'xMin': 113,
                'xMax': 116,
                'yMin': 10.5,
                'yMax': -25
            },
            'GRID_STEP': 20.
        }) # Model

        mappings = [(view.region, 'GRID_LIMITS'),
                    (view.gridSettings.region, 'GRID_LIMITS'),
                    (view.gridSettings.step, 'GRID_STEP'),
                    (view.view, 'GRID_LIMITS'),
                    (view.view, 'GRID_STEP')]

        def makeCallback(key):
            return lambda x: self.onControlChanged(x, key)

        self.notifyWhom = {}
        for control, key in mappings:
            controls = self.notifyWhom.setdefault(key, [])
            controls.append(control)
            if hasattr(control, 'addCallback'):
                control.addCallback(makeCallback(key))

        self.onSettingsChanged(self.settings)
        self.settings.addCallback(self.onSettingsChanged)

        view.pack(fill='both', expand=True)
        self.after_idle(self.show)

        self.view = view
        self.tcrm = tcrm

    def onSettingsChanged(self, settings):
        log.info('Settings changed')
        for key, value in settings.items():
            controls = self.notifyWhom[key]
            for control in controls:
                control.set(value)

    def onControlChanged(self, value, key):
        self.settings[key] = value

    def onWantOutput(self, control):
        output = self.tcrm.getOutput()
        control.emit(output)

    def onRun(self):
        print('onRun')
        self.tcrm.run()

    def onQuit(self):
        self.tcrm.quit()
        self.destroy()

    def show(self):
        self.update()
        self.deiconify()

Controller().mainloop()
