"""
:mod:`maps` -- base class for generating maps
=============================================

.. module:: maps
    :synopsis: Generate map images using `cartopy`

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

Note: Many of the defaults (e.g. coastline colours, font sizes for
      the graticule) are set up for GA's style. These can be overridden
      in this base class if users want a different style for their output.

"""



import numpy as np
import numpy.ma as ma

from matplotlib.figure import Figure
import cartopy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import seaborn as sns

def levels(maxval, minval=0):
    """
    Calculate a nice number of levels between `minval` and `maxval`
    for plotting contour intervals.

    :param float maxval: Maximum value in the field to be plotted.
    :param float minval: Minimum value in the field to be plotted.

    :returns: Array of levels and the exponential by which the
              values are scaled
    :rtype: `tuple`

    :raises: ValueError if the minimum and maximum values are equal.

    """

    if not isinstance(maxval, (int, float)):
        raise TypeError("maxval must be numeric value")

    if not isinstance(minval, (int, float)):
        raise TypeError("minval must be numeric value")

    if maxval == minval:
        raise ValueError("Minimum and maximum value are equal")

    if maxval < minval:
        mm = maxval
        maxval = minval
        minval = mm

    min_levels = 7.0
    level_opts = np.array([5.0, 1.0, 0.5, 0.25, 0.2, 0.1])
    exponent = int(np.floor(np.log10(maxval)))
    significand = (maxval - minval) * 10**-exponent
    level_step = level_opts[np.where((significand/level_opts) > \
                                     min_levels)[0][0]]
    lvls = np.arange(minval, maxval, level_step*(10.0**exponent))

    return lvls, exponent

def selectColormap(data_range, percent=0.1):
    """
    Determine whether to use a sequential or diverging color map, based on the
    defined data range to be used for the plot. Note the diverging color map
    only works when the data range spans zero (and it won't automatically put
    the neutral colour at zero).
    
    Red to green colour palette using recommendations from  ISO22324 (2015).
    Unsaturated colour palette:
    [(0.486, 0.722, 0.573), (0.447, 0.843, 0.714), (0.875, 0.882, 0.443), 
    (0.969, 0.906, 0.514), (0.933, 0.729, 0.416), (0.906, 0.522, 0.373), 
    (0.937, 0.522, 0.616)]

    :param data_range: array-like containing either the minimum and maximum
                       levels of the data range, or the array of levels (e.g.
                       for contour maps).
    :param float percent: Threshold for switching from diverging to sequential
                          colormap

    :returns: A `matplotlib.colors.LinearSegmentedColormap` instance used for
              setting the `Figure.cmap` attribute.

    :raises: TypeError for non-numeric input, or if the `data_range` is not
             array-like.
    """
    if not isinstance(data_range, (list, np.ndarray, tuple)):
        raise TypeError("Data range must be a list or array of numeric values")

    if not isinstance(max(data_range), (int, float)):
        raise TypeError("Data range must be a list or array of numeric values")

    x = (abs(max(data_range)) - abs(min(data_range)))/ \
        (max(data_range) - min(data_range))
    if abs(x) < percent:
        palette = sns.color_palette("RdBu", 7)
        cmap = sns.blend_palette(palette, as_cmap=True)
    else:
        palette = [(1, 1, 1),
                   (0.000, 0.627, 0.235),
                   (0.412, 0.627, 0.235), 
                   (0.663, 0.780, 0.282),
                   (0.957, 0.812, 0.000),
                   (0.925, 0.643, 0.016),
                   (0.835, 0.314, 0.118),
                   (0.780, 0.086, 0.118)]
        cmap = sns.blend_palette(palette, as_cmap=True)

    return cmap


class MapFigure(Figure):
    """
    A base class for all map figures. Implements methods to annotate
    maps in a consistent fashion.

    """

    def __init__(self):
        Figure.__init__(self)
        self.subfigures = []
        palette = sns.color_palette("YlOrRd", 7)
        self.cmap = sns.blend_palette(palette, as_cmap=True)

        #self.canvas = FigureCanvas

    def add(self, data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs):
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
        :param lvls: `numpy.ndarray` of levels used for plotting contours.
        :param str cbarlabel: Label for the colorbar (if one is to be added).
        :param dict map_kwargs: `dict` of keyword arguments used for setting
                                up the `GeoAxes` instance.

        """

        self.subfigures.append((data, xgrid, ygrid, title,
                                lvls, cbarlab, map_kwargs))

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

    def addGraticule(self, axes, mapobj):
        """
        Add a graticule to the map instance (tick marks and
        lines for parallels and meridians). Grid spacing is automatically
        determined from the maximum span of either longitude or latitude.

        :param axes: Current `matplotlib.Axes` instance being worked on.
        :param mapobj: Current `GeoAxes` instance to annotate.

        """
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        xmin, xmax, ymin, ymax = mapobj.get_extent()

        dx = abs(xmin - xmax)
        dy = abs(ymin - ymax)
        dd = max(dx, dy)
        gr_opts = np.array([30., 10., 5., 4., 2.])
        min_gr = 5
        try:
            dl = gr_opts[np.where((dd/gr_opts) >= min_gr)[0][0]]
        except IndexError:
            dl = 2.

        meridians = np.arange(xmin // dl * dl, xmax + dl, dl)
        parallels = np.arange(ymin // dl * dl, ymax + dl, dl)

        gl = mapobj.gridlines(xlocs=meridians, ylocs=parallels,
                              draw_labels=True)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.top_labels = False
        gl.right_labels = False

    def addCoastline(self, mapobj):
        """
        Draw coastlines and a map background to the current
        `GeoAxes` instance.

        :param mapobj: Current `GeoAxes` instance to add coastlines to.

        """
        mapobj.coastlines(resolution='10m')

    def fillContinents(self, mapobj, fillcolor="#FFDAB5"):
        """
        Fill continents with a base color in the current
        `GeoAxes` instance.

        :param mapobj: Current `GeoAxes` instance to color fill the
        continents on.

        """
        mapobj.add_feature(cartopy.feature.LAND, color=fillcolor)

    def maskOceans(self, mapobj, fillcolor="#66ccff"):
        """
        Mask oceans with a blue background in the current `GeoAxes` instance.

        :param mapobj: Current `GeoAxes` instance to colour fill the oceans
        :param fillcolor: Optional colour to use (default #66CCFF)
        """

        mapobj.add_feature(cartopy.feature.OCEAN, color=fillcolor)

    def addMapScale(self, mapobj):
        """
        Add a map scale to the curent `Basemap` instance. This
        automatically determines a 'nice' length forthe scale bar -
        chosen to be approximately 20% of the map width at the centre
        of the map.

        :param mapobj: Current `Basemap` instance to add the scale bar to.

        """
        return # TODO: migrate to cartopy - see https://stackoverflow.com/questions/32333870/

        midlon = (mapobj.lonmax - mapobj.lonmin) / 2.
        midlat = (mapobj.latmax - mapobj.latmin) / 2.

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
        Configure the map extent
        """
        x0 = map_kwargs['llcrnrlon']
        y0 = map_kwargs['llcrnrlat']
        x1 = map_kwargs['urcrnrlon']
        y1 = map_kwargs['urcrnrlat']

        axes.set_extent([x0, x1, y0, y1])

        return axes, xgrid, ygrid

    def subplot(self, axes, subfigure):
        """
        Add a new subplot to the current figure.

        :param axes: `matplotlib.Axes` instance that will be used for the
                     subplot. Usually created by a call to `self.add_subplot()`
        :param tuple subfigure: A tuple containing the data, x- and y-grid,
                                title, levels, colorbar label and
                                a map_kwargs `dict`.

        """

        data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)

        CS = mapobj.contour(mx, my, data, levels=lvls,
                            colors='k', linewidth=1.5)
        axes.clabel(CS, inline=1)
        axes.set_title(title)
        self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.fillContinents(mapobj)
        self.addMapScale(mapobj)

    def plot(self):
        """
        Plot all subfigures onto a figure in a (roughly) square
        arrangement.

        """

        nfig = len(self.subfigures)
        rows = int(np.ceil(np.sqrt(nfig)))
        cols = int(np.ceil(nfig / rows))

        width, height = self.get_size_inches()
        self.set_size_inches(width * cols, rows * height)
        for i, subfigure in enumerate(self.subfigures):
            if subfigure[-1]['projection'] == 'merc':
                axes = self.add_subplot(rows, cols, i+1,
                                        projection=cartopy.crs.PlateCarree())
            else:
                raise NotImplementedError(
                          "Only Mercator projection supported")
            self.subplot(axes, subfigure)

    def save(self, filename):
        """
        Save the figure to file.

        :param str filename: Full path (incl. extension) to save the figure to.

        """
        canvas = FigureCanvas(self)
        self.tight_layout()
        canvas.print_figure(filename, dpi=300, bbox_inches='tight')


class FilledContourMapFigure(MapFigure):
    """
    A filled contour map figure.

    Useful for plotting continuous values across a spatial domain.

    """

    def add(self, data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs):
        """
        Add a filled contour subfigure to the list of subfigures to be
        plotted.

        :param data: :class:`numpy.ndarray` of (continuous)  values to be
                     plotted.
        :param xgrid: :class:`numpy.ndarray` of x-coordinates of the data
                      points.
        :param ygrid: :class:`numpy.ndarray` of y-coordinates of the data
                      points.
        :param str title: Title to put on the subfigure.
        :param lvls: :class:`numpy.ndarra` or `list` of values that
                       define the contour levels to use.
        :param str cbarlab: Label for the color bar.
        :param dict map_kwargs: A dict containing keyword arguments for
                                setting up the :class:`GeoAxes` instance.

        """
        super(FilledContourMapFigure, self).add(data, xgrid, ygrid,
                                                title, lvls, cbarlab,
                                                map_kwargs)

    def subplot(self, axes, subfigure):
        """
        Generate the subplot on the figure.

        :param axes: :class:`matplotlib.axes` instance where the plot will
                     be added.
        :param tuple subfigure: A tuple that holds all the required elements
                                of the plot, including the data, x- and y-grids,
                                title, contour levels, colorbar label and
                                map keyword arguments.

        """
        data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)
        cmap = selectColormap(lvls)
        CS = mapobj.contourf(mx, my, data, levels=lvls,
                             extend='both', cmap=cmap)
        CB = self.colorbar(CS, ticks=lvls[::2], ax=axes, extend='both')
        CB.set_label(cbarlab)
        axes.set_title(title)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.addMapScale(mapobj)

class MaskedContourMapFigure(FilledContourMapFigure):
    """
    Create a filled contour plot with ocean areas masked.

    """
    def subplot(self, axes, subfigure):
        """
        Generate the subplot on the figure.

        :param axes: :class:`matplotlib.axes` instance where the plot will
                     be added.
        :param tuple subfigure: A tuple that holds all the required elements
                                of the plot, including the data, x- and y-grids,
                                title, contour levels, colorbar label and
                                map keyword arguments.

        """

        data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)
        cmap = selectColormap(lvls)
        CS = mapobj.contourf(mx, my, data, levels=lvls,
                             extend='both', cmap=cmap)
        CB = self.colorbar(CS, ticks=lvls[::2], ax=axes, extend='both')
        CB.set_label(cbarlab)
        axes.set_title(title)
        self.labelAxes(axes)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.maskOceans(mapobj)
        self.addMapScale(mapobj)

class ArrayMapFigure(MapFigure):
    """Array plot"""
    def add(self, data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs):

        super(ArrayMapFigure, self).add(data, xgrid, ygrid, title, datarange,
                                        cbarlab, map_kwargs)

    def subplot(self, axes, subfigure):
        """
        Generate the subplot on the figure.

        :param axes: :class:`matplotlib.axes` instance where the plot will
                     be added.
        :param tuple subfigure: A tuple that holds all the required elements
                                of the plot, including the data, x- and y-grids,
                                title, data range, colorbar label and
                                map keyword arguments.

        """
        data, xgrid, ygrid, title, \
            datarange, cbarlab, map_kwargs = subfigure
        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)
        cmap = selectColormap(datarange)
        vmin = datarange[0]
        vmax = datarange[1]
        CS = mapobj.pcolormesh(mx, my, data, vmin=vmin,
                               vmax=vmax, cmap=cmap)
        CB = self.colorbar(CS, orientation='horizontal', pad=0.1,
                             ax=axes)
        CB.set_label(cbarlab)
        axes.set_title(title)
        self.addGraticule(axes, mapobj)
        self.addCoastline(mapobj)
        self.addMapScale(mapobj)

class MaskedArrayMapFigure(ArrayMapFigure):
    """Array plot with ocean areas masked"""
    def subplot(self, axes, subfigure):
        """
        Generate the subplot on the figure.

        :param axes: :class:`matplotlib.axes` instance where the plot will
                     be added.
        :param tuple subfigure: A tuple that holds all the required elements
                                of the plot, including the data, x- and y-grids,
                                title, data range, colorbar label and
                                map keyword arguments.

        """
        data, xgrid, ygrid, title, \
            datarange, cbarlab, map_kwargs = subfigure

        subfigure = (data, xgrid, ygrid, title,
                     datarange, cbarlab, map_kwargs)
        super(MaskedArrayMapFigure, self).subplot(axes, subfigure)

#class BarbMapFigure(MapFigure):
#    """
#    Map figure with velocity indicated by wind barbs
#
#    """
#    def add(self, xdata, ydata, xgrid, ygrid, title,
#            lvls, cbarlab, map_kwargs):
#        self.subfigures.append((xdata, ydata, xgrid, ygrid,
#                                title, lvls, cbarlab, map_kwargs))
#
#    def subplot(self, axes, subfigure):
#        """
#        Generate the subplot on the figure.
#
#        :param axes: :class:`matplotlib.axes` instance where the plot will
#                     be added.
#        :param tuple subfigure: A tuple that holds all the required elements
#                                of the plot, including the data, x- and y-grids,
#                                title, contour levels, colorbar label and
#                                map keyword arguments.
#
#        """
#        xdata, ydata, xgrid, ygrid, title, \
#            lvls, cbarlab, map_kwargs = subfigure
#        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)
#
#        mag = np.sqrt(xdata*xdata + ydata*ydata)
#        CS = mapobj.contourf(mx, my, mag, lvls, cmap=self.cmap)
#        CB = mapobj.colorbar(CS, location='right', pad='5%',
#                             fig=self, ax=axes, ticks=lvls[::2])
#        CB.set_label(cbarlab)
#        mapobj.barbs(xgrid, ygrid, xdata, ydata, length=5, linewidth=0.5, 
#                     latlon=True)
#        axes.set_title(title)
#        self.addGraticule(axes, mapobj)
#        self.addCoastline(mapobj)
#
#class ScatterMapFigure(MapFigure):
#    """
#    Scatter plot over a map figure.
#
#    """
#    def add(self, data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs):
#        self.subfigures.append((data, xgrid, ygrid, title, lvls,
#                                cbarlab, map_kwargs))
#
#    def subplot(self, axes, subfigure):
#        """
#        Generate the subplot on the figure.
#
#        :param axes: :class:`matplotlib.axes` instance where the plot will
#                     be added.
#        :param tuple subfigure: A tuple that holds all the required elements
#                                of the plot, including the data, x- and y-grids,
#                                title, contour levels, colorbar label and
#                                map keyword arguments. ``data`` is an n-by-2
#                                :class:`numpy.ndarray` containing the x- and y-
#                                coordinates of the scatter points.
#
#        """
#        data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs = subfigure
#        xp, yp = data
#        mapobj, mx, my = self.createMap(axes, xgrid, ygrid, map_kwargs)
#
#        self.addCoastline(mapobj)
#        self.fillContinents(mapobj)
#
#        mxp, myp = mapobj(xp, yp)
#        mapobj.scatter(mxp, myp)
#        axes.set_title(title)
#        self.addGraticule(axes, mapobj)
#        self.addMapScale(mapobj)

class HazardMap(FilledContourMapFigure):
    """
    A map for presenting return level data. 

    """
    def plot(self, data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs):
        #dx = np.mean(np.diff(xgrid))
        #dmask = data.mask
        #data = ma.array(data, mask=dmask)
        self.add(data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs)
        self.cmap = sns.light_palette("orange", as_cmap=True)
        super(HazardMap, self).plot()

class WindfieldMap(FilledContourMapFigure):
    """
    Plot a wind field using filled contours. Only presents the magnitude of
    the wind field, not the direction.

    """
    def plot(self, data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs):
        self.add(data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs)
        self.cmap = sns.light_palette("orange", as_cmap=True)
        super(WindfieldMap, self).plot()

class ArrayMap(ArrayMapFigure):
    """
    Plot a gridded array of values over a map
    """
    def plot(self, data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs):
        self.add(data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs)
        super(ArrayMap, self).plot()

def saveFigure(figure, filename):
    """
    Save a :class:`MapFigure` instance to file.

    :param figure: :class:`MapFigure` instance.
    :param str filename: Path (including extension) to save the figure to.

    """
    canvas = FigureCanvas(figure)
    canvas.print_figure(filename, dpi=300, bbox_inches='tight')

def saveHazardMap(data, xgrid, ygrid, title, lvls, cbarlab,
                  map_kwargs, filename):
    """
    Plot and save a hazard map.

    """
    fig = HazardMap()

    fig.plot(data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs)
    fig.save(filename)

def saveWindfieldMap(data, xgrid, ygrid, title, lvls,
                     cbarlab, map_kwargs, filename):
    """
    Plot and save a map of wind speeds as a filled contour map.
    """
    fig = WindfieldMap()
    fig.plot(data, xgrid, ygrid, title, lvls, cbarlab, map_kwargs)
    fig.save(filename)

#def saveArrayMap(data, xgrid, ygrid, title, datarange, cbarlab,
#                 map_kwargs, filename):
#    """
#    Plot and save a map of an array of data values.
#
#    """
#    fig = ArrayMap()
#    fig.plot(data, xgrid, ygrid, title, datarange, cbarlab, map_kwargs)
#    fig.save(filename)

def demo():
    """
    Perform a basic demonstration of the plotting routines.

    """
    from netCDF4 import Dataset
    from os.path import join as pjoin, dirname, normpath
    baseDir = normpath(pjoin(dirname(__file__), '..'))
    inputPath = pjoin(baseDir, 'output', 'hazard')

    hazardFile = pjoin(inputPath, 'hazard.nc')
    try:
        ncobj = Dataset(hazardFile, 'r')
    except:
        raise "{0} doesn't exist".format(hazardFile)

    xdata = ncobj.variables['lon'][:]
    ydata = ncobj.variables['lat'][:]
    tdata = ncobj.variables['ari'][:]
    vdata = ncobj.variables['wspd'][:]
    ncobj.close()
    [xgrid, ygrid] = np.meshgrid(xdata, ydata)
    map_kwargs = dict(llcrnrlon=xdata.min(),
                      llcrnrlat=ydata.min(),
                      urcrnrlon=xdata.max(),
                      urcrnrlat=ydata.max(),
                      projection='merc',
                      resolution='i')
    cbarlab = "Wind speed (%s)" % ncobj.variables['wspd'].units
    lvls, exp = levels(np.max(vdata[8, :, :]))

    figure = MaskedContourMapFigure()
    figure.add(vdata[2, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[2]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[4, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[4]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[6, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[6]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[8, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[8]),
               lvls, cbarlab, map_kwargs)
    figure.plot()
    figure.save('docs/hazard_map.png')
    figure = FilledContourMapFigure()
    figure.add(vdata[2, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[2]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[4, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[4]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[6, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[6]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[8, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[8]),
               lvls, cbarlab, map_kwargs)
    figure.plot()
    figure.save('docs/hazard_map.png')
    figure = MapFigure()
    figure.add(vdata[2, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[2]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[4, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[4]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[6, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[6]),
               lvls, cbarlab, map_kwargs)
    figure.add(vdata[8, :, :], xgrid, ygrid,
               "{0:.0f}-year return period wind speed".format(tdata[8]),
               lvls, cbarlab, map_kwargs)
    figure.plot()
    figure.save('docs/hazard_map_contour.png')

    saveHazardMap(vdata[2, :, :], xgrid, ygrid,
                  "{0:.0f}-year return period wind speed".format(tdata[2]),
                  lvls, cbarlab, map_kwargs, 'docs/hazard_map_ex.png')
if __name__ == "__main__":
    demo()
