"""
:mod:`maps` -- base class for generating maps
=============================================

.. module:: maps
    :synopsis: Generate map images using `mpl_toolkits.basemap`

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

Note: Many of the defaults (e.g. coastline colours, font sizes for 
      the graticule) are set up for GA's style. These can be overridden
      in this base class if users want a different style for their output.

"""

from __future__ import division

import numpy as np

from matplotlib.figure import Figure
from mpl_toolkits.basemap import Basemap

from Utilities.smooth import smooth

import seaborn

def levels(maxval, minval=0):
    """
    Calculate a nice number of levels between `minval` and `maxval`
    for plotting contour intervals.

    :param float maxval: Maximum value in the field to be plotted.
    :param float minval: Minimum value in the field to be plotted.
    
    :returns: Array of levels and the exponential by which the 
              values are scaled
    :rtype: `tuple`

    """

    min_levels = 7.0
    level_opts = np.array([5.0, 1.0, 0.5, 0.25, 0.2, 0.1])
    exponent = int(np.floor(np.log10(maxval)))
    significand = (maxval - minval) * 10**-exponent
    level_step = level_opts[np.where((significand/level_opts) > \
                                     min_levels)[0][0]]
    levels = np.arange(minval, maxval, level_step*(10.0**exponent))

    return levels, exponent

class MapFigure(Figure):

    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []
        self.cmap = seaborn.light_palette('#0A437A', as_cmap=True)
        
    def add(self, data, xgrid, ygrid, title, levels, cbarlab, map_kwargs):
        """
        Add a new subfigure to the collection of subfigures to be 
        plotted.
        
        :param data: `numpy.ndarray` of the data to be plotted.
        :param xgrid: `numpy.ndarray` representing the points on the 
                      x-axis. For gridded arrays of data, `xgrid` must
                      have the same shape as `data`.
        :param ygrid: `numpy.ndarray` representing the points on the 
                      y-axis. For gridded arrays of data, `ygrid` must
                      have the same shape as `data`.
        :param str title: Subfigure title
        :param levels: `numpy.ndarray` of levels used for plotting contours.
        :param str cbarlabel: Label for the colorbar (if one is to be added).
        :param dict map_kwargs: `dict` of keyword arguments used for setting
                                up the `Basemap` instance.

        """
        
        self.subfigures.append((data, xgrid, ygrid, title, 
                                levels, cbarlab, map_kwargs))

    def labelAxes(self, axes, xlabel='Longitude', ylabel='Latitude'):
        """
        Add labels to the x- and y-axis on the current `matplotlib.Axes`
        instance.
        
        :param axes: Current `matplotlib.Axes` instance to annotate.
        :param str xlabel: Label for the x-axis. Defaults to 'Longitude'.
        :param str ylabel: Label for the y-axis. Defaults to 'Latitude'.
        
        """

        axes.set_xlabel(xlabel, labelpad=20,) 
        axes.set_ylabel(ylabel, labelpad=35,) 

    def addGraticule(self, axes, mapobj, dl=10.):
        """
        Add a graticule to the map instance (tick marks and 
        lines for parallels and meridians).

        :param axes: Current `matplotlib.Axes` instance being worked on.
        :param mapobj: Current `Basemap` instance to annotate.
        :param float dl: Optional increment between parallels/meridians. 
                         Default value is 10 degrees.

        """

        xmin = mapobj.llcrnrlon
        xmax = mapobj.urcrnrlon
        ymin = mapobj.llcrnrlat
        ymax = mapobj.urcrnrlat

        meridians = np.arange(dl*np.floor(xmin / dl),
                                dl*np.ceil(xmax / dl) + dl, dl)
        parallels = np.arange(dl*np.floor(ymin / dl),
                                dl*np.ceil(ymax / dl) + dl, dl)

        mapobj.drawparallels(parallels, linewidth=0.25,
                             labels=[1, 0, 0, 1], style="italic") 
        mapobj.drawmeridians(meridians, linewidth=0.25,
                             labels=[1, 0, 0, 1], style='italic') 

        axes.tick_params(direction='in', length=4, width=1)

    def addCoastline(self, mapobj):
        """
        Draw coastlines and a map background to the current
        `Basemap` instance.

        :param mapobj: Current `Basemap` instance to add coastlines to.

        """

        mapobj.drawcoastlines(linewidth=.5, color="k")
        mapobj.drawmapboundary(fill_color="#BEE8FF")

    def fillContinents(self, mapobj):
        """
        Fill continents with a base color in the current
        `Basemap` instance.

        :param mapobj: Current `Basemap` instance to color fill the 
        continents on.

        """
        mapobj.fillcontinents(color="#FFDAB5",
                              lake_color="#BEE8FF",
                              zorder=0)

    def addMapScale(self, mapobj):
        """
        Add a map scale to the curent `Basemap` instance. This
        automatically determines a 'nice' length forthe scale bar - 
        chosen to be approximately 20% of the map width at the centre
        of the map.

        :param mapobj: Current `Basemap` instance to add the scale bar to.

        
        """
        lonmin = mapobj.lonmin
        lonmax = mapobj.lonmax
        latmin = mapobj.latmin
        latmax = mapobj.latmax
        midlon = (lonmax - lonmin) / 2.
        midlat = (latmax - latmin) / 2.

        xmin = mapobj.llcrnrx
        xmax = mapobj.urcrnrx
        ymin = mapobj.llcrnry
        ymax = mapobj.urcrnry

        xloc = xmin + 0.15 * abs(xmax - xmin)
        yloc = ymin + 0.1 * abs(ymax - ymin)

        lonloc, latloc = mapobj(xloc, yloc, inverse=True)

        # Set scale length to nearest 100-km for 20% of map width
        scale_length = 100*int((0.2 * (xmax - xmin) / 1000.)/100)
        mapobj.drawmapscale(lonloc, latloc, midlon, midlat, scale_length,
                            barstyle='fancy', zorder=10)

    def createMap(self, axes, xgrid, ygrid, map_kwargs):
        """
        Create a :class:`basemap` object, and convert the x/y grid
        to map coordinates
        """

        mapobj = Basemap(ax=axes, **map_kwargs)
        mx, my = mapobj(xgrid, ygrid)
        return mapobj, mx, my

    def subplot(self, axes, subfigure):
        """
        Add a new subplot to the current figure.

        :param axes: `matplotlib.Axes` instance that will be used for the
                     subplot. Usually created by a call to `self.add_subplot()`
        :param tuple subfigure: A tuple containing the data, x- and y-grid,
                                title, levels, colorbar label and 
                                a map_kwargs `dict`.

        """

        data, xgrid, ygrid, title, levels, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        CS = mapobj.contour(mx, my, data, levels=levels, 
                            colors='k', linewidth=1.5)
        axes.clabel(CS, inline=1)
        axes.set_title(title)
        self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.addMapScale(mapobj)

    def plot(self):
        """
        Plot all subfigures onto a figure in a (roughly) square
        arrangement. 
        
        """

        n = len(self.subfigures)
        r = int(np.ceil(np.sqrt(n)))
        c = int(np.ceil(n / r))

        w, h = self.get_size_inches()
        self.set_size_inches(w * c , r * h)
        for i, subfigure in enumerate(self.subfigures):
            axes = self.add_subplot(r, c, i+1)
            self.subplot(axes, subfigure)
        
        
class FilledContourMapFigure(MapFigure):

    def add(self, data, xgrid, ygrid, title, levels, cbarlab, map_kwargs):
        super(FilledContourMapFigure, self).add(data, xgrid, ygrid,
                                          title, levels, cbarlab, map_kwargs)

    def subplot(self, axes, subfigure):
        data, xgrid, ygrid, title, levels, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)
        CS = mapobj.contourf(mx, my, data, levels=levels, extend='both', cmap=self.cmap)
        CB = mapobj.colorbar(CS, location='right', pad='5%', ticks=levels[::2],
                             fig=self, ax=axes, extend='both')
        CB.set_label(cbarlab)
        axes.set_title(title)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.addMapScale(mapobj)

class MaskedContourMapFigure(FilledContourMapFigure):
    """Filled contour plot with ocean areas masked"""
    def subplot(self, axes, subfigure):
        from mpl_toolkits.basemap import maskoceans
        data, xgrid, ygrid, title, levels, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        masked_data = maskoceans(xgrid, ygrid, data, inlands=False)
        CS = mapobj.contourf(mx, my, masked_data, levels=levels, extend='both', cmap=self.cmap)
        CB = mapobj.colorbar(CS, location='right', pad='5%', ticks=levels[::2],
                             fig=self, ax=axes, extend='both')
        CB.set_label(cbarlab)
        axes.set_title(title)
        self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.addMapScale(mapobj)

class ArrayMapFigure(MapFigure):
    """Array plot"""
    def add(self, data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs):
        super(ArrayMapFigure, self).add(data, xgrid, ygrid, title, datarange,
                                        cbarlab, map_kwargs)

    def subplot(self, axes, subfigure):
        data, xgrid, ygrid, title, \
            datarange, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        vmin = datarange[0]
        vmax = datarange[1]
        CS = mapobj.pcolormesh(mx, my, data, vmin=vmin, vmax=vmax, cmap=self.cmap)
        CB = mapobj.colorbar(CS, location='right', pad='5%',
                             fig=self, ax=axes)
        CB.set_label(cbarlab)
        axes.set_title(title)
        self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.addMapScale(mapobj)

class MaskedArrayMapFigure(ArrayMapFigure):
    """Array plot with ocean areas masked"""
    def subplot(self, axes, subfigure):
        from mpl_toolkits.basemap import maskoceans
        data, xgrid, ygrid, title, \
            datarange, cbarlabel, map_kwargs = subfigure

        masked_data = maskoceans(xgrid, ygrid, data, inlands=False)
        subfigure = (masked_data, xgrid, ygrid, title, 
                     datarange, cbarlab, map_kwargs)
        super(MaskedArrayMapFigure, self).subplot(axes, subfigure)

class BarbMapFigure(MapFigure):

    def add(self, xdata, ydata, xgrid, ygrid, title, 
            levels, cbarlab, map_kwargs):
        self.subfigures.append((mapobj, xdata, ydata, xgrid, ygrid,
                                title, levels, cbarlab, map_kwargs))

    def subplot(self, axes, subfigure):
        xdata, ydata, xgrid, ygrid, title, \
            levels, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        mag = np.sqrt(xdata*xdata + ydata*ydata)
        CS = mapobj.contourf(mx, my, mag, levels, cmap=self.cmap)
        CB = mapobj.colorbar(CS, location='right', pad='5%',
                             fig=self, ax=axes, ticks=levels[::2])
        CB.set_label(cbarlab)
        mapobj.barbs(xgrid, ygrid, xdata, ydata, length=5, linewidth=0.5)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)

class ScatterMapFigure(MapFigure):

    def add(self, data, xgrid, ygrid, title, levels, cbarlab, map_kwargs):
        self.subfigures.append((data, xgrid, ygrid, title, levels, 
                                cbarlab, map_kwargs))

    def subplot(self, axes, subfigure):
        data, xgrid, ygrid, title, levels, cbarlab, map_kwargs = subfigure
        xp, yp = data
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        self.addCoastline(mapobj)
        self.fillContinents(mapobj)
        
        mxp, myp = mapobj(xp, yp)
        mapobj.scatter(mxp, myp)
        axes.set_title(title)
        self.addGraticule(axes, mapobj)
        self.addMapScale(mapobj)
        
        

class HazardMap(MaskedContourMapFigure):

    def plot(self, data, xgrid, ygrid, title, levels, cbarlab, map_kwargs):
        # Smooth the data to reduce 'lines-on-a-map' inferences:
        dx = np.mean(np.diff(xgrid))
        data = smooth(data, int(1/dx))
        self.add(data, xgrid, ygrid, title, levels, cbarlab, map_kwargs)
        super(HazardMap, self).plot()



class WindfieldMap(FilledContourMapFigure):

    def plot(self, data, xgrid, ygrid, title, levels, cbarlab, map_kwargs):
        self.add(data, xgrid, ygrid, title, levels, cbarlab, map_kwargs)
        super(WindfieldMap, self).plot()

class ArrayMap(ArrayMapFigure):

    def plot(self, data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs):
        self.add(data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs)
        super(ArrayMap, self).plot()

def saveFigure(figure, filename):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename,dpi=300)

def saveHazardMap(data, xgrid, ygrid, title, levels, cbarlab, 
                  map_kwargs, filename):
    fig = HazardMap()

    fig.plot(data, xgrid, ygrid, title, levels, cbarlab, map_kwargs)
    saveFigure(fig, filename)

def saveWindfieldMap(data, xgrid, ygrid, title, levels, 
                     cbarlab, map_kwargs, filename):
    fig = WindfieldMap()
    fig.plot(data, xgrid, ygrid, title, levels, cbarlab, map_kwargs)
    saveFigure(fig, filename)

def saveArrayMap(data, xgrid, ygrid, title, datarange, cbarlab, 
                 map_kwargs, filename):
    fig = ArrayMap()
    fig.plot(data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs)
    saveFigure(fig, filename)

def main():
    from netCDF4 import Dataset
    from os.path import join as pjoin, dirname, normpath
    baseDir = normpath(pjoin(dirname(__file__), '..'))
    inputPath = pjoin(baseDir, 'output', 'hazard')

    hazardFile = pjoin(inputPath, 'hazard.nc')
    f = Dataset(hazardFile, 'r')

    xdata = f.variables['lon'][:]
    ydata = f.variables['lat'][:]
    tdata = f.variables['years'][:]

    vdata = f.variables['wspd'][:]

    [xgrid, ygrid] = np.meshgrid(xdata, ydata)
    map_kwargs = dict(llcrnrlon=xdata.min(),
                      llcrnrlat=ydata.min(),
                      urcrnrlon=xdata.max(),
                      urcrnrlat=ydata.max(),
                      projection='merc',
                      resolution='i')
    title = "Return period wind speed"
    cbarlab = "Wind speed (%s)"%f.variables['wspd'].units
    levels = np.arange(30, 101., 5.)
    datarange = (np.floor(vdata[2,:,:].min()), np.ceil(vdata[4,:,:].max()))

    figure = ArrayMapFigure()
    figure.add(vdata[2,:,:], xgrid, ygrid, title, datarange, cbarlab, map_kwargs)
    figure.add(vdata[4,:,:], xgrid, ygrid, title, datarange, cbarlab, map_kwargs)
    figure.plot()
    saveFigure(figure, 'docs/hazard_map.png')

if __name__ == "__main__":
    main()
