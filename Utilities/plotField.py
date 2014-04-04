"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Title: plotField.py
Author: Craig Arthur, craig.arthur@ga.gov.au
Created: 2009-10-14
Description: Plot the contents of a gridded ascii file over a given
             domain.  One can optionally apply a smoothing function.

Version: 164
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-10-14
Modification: First version

Version: 351
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-10-22 4:35:PM
Modification: Provides a generic interface to plotting data on a base map,
            with wrappers to permit plotting wind fields, pressure fields,
            or to use a configuration file to control the program.

Version: 411
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2010-10-19 1:07:PM
Modification: Apply masked arrays to data being plotted (if data is loaded from
              a netCDF file).

Version: $Rev: 642 $
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2011-06-08 1:00:PM
Modification: Added plotArray and plotBarb functions.

Id: $Id: plotField.py 642 2012-02-21 07:54:04Z nsummons $
"""

import os, sys, logging

logger = logging.getLogger()

import getopt, traceback

from numpy import *
from matplotlib import pyplot, cm, rcParams

try:
    from mpl_toolkits.basemap import Basemap
    NO_BASEMAP = False
except ImportError:
    NO_BASEMAP = True
    logger.warn('Basemap package not installed. Disabling some plots')

# GA modules:
from files import flStartLog, flConfigFile, flLogFatalError
from config import cnfGetIniValue
import nctools
import grid
import colours
from smooth import smooth

filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)

__version__ = '$Id: plotField.py 642 2012-02-21 07:54:04Z nsummons $'


def _usage():
    """
    Short function to describe how to use this function using the
    command line interface.
    """
    print "plotField.py: "
    print __version__
    print "Plot the contents of a gridded ascii or netcdf file over a basemap."
    print "Options:"
    print "-c, --config: configuration file                             "
    print "-h, --help: Prints this message.                             "
    print "-l, --loglevel: Log level ('DEBUG','INFO','WARN','CRITICAL') "
    print "-v, --verbose: print log messages to STDOUT.                 "
    print ""
    """
    print "Configuration file settings:"
    print "[Input]"
    print "File = full path to input file"
    print "Format= format of the input file (either .txt for ascii grd format files,"
    print " or .nc for netcdf files"
    print "Variable=name of the variable to plot (for netcdf files only)"
    print "Record=integer number of the record dimension (e.g. time) to "
    print "       select for plotting"
    print ""
    print "[Output]"
    print "File=Full path to output image file to create"
    print "Format=image format to generate (eps, png, gif)"
    print "ColourMap=one of the pre-defined colour scales from Matplotlib,"
    print "          or alternatively a custom colour scale"
    print "Smoothing=the number of grid-points to include in a gaussian"
    print "          filter applied to the data. Set to False to turn "
    print "          the smoothing off"
    print "Label=Title for the plot"
    print "ScaleMin=minimum value for the colour scale"
    print "ScaleMax=maximum value for the colour scale"
    print "ScaleInt=interval for the colour scale"
    print "MaskLand=Optionally mask land areas. If a shapefile name is
    print "         provided, that will act as the mask. Defaluts to False"
    print "MaskOcean=Optionally mask ocean areas. Defaults to False"

    """

def _process(argv):
    """
    A wrapper function to provide an interface between the command line
    args and the actual plotField function. This function reads settings
    from a configuration file and then passes those arguments to the
    plotField function.
    """
    if len(argv)==0:
        _usage()
        sys.exit(2)

    logLevel = 'INFO'
    verbose = False
    configFile = flConfigFile()

    try:
        opts, args = getopt.getopt(argv, "c:hl:v", ["config=", "help",
                                                   "loglevel=", "verbose"])
    except getopt.GetoptError:
        _usage()
        sys.exit(2)
    for opt,arg in opts:
        if opt in ("-h", "--help"):
            _usage()
            sys.exit(2)
        elif opt in ("-c", "--config"):
            configFile = arg
        elif opt in ("-l", "--loglevel"):
            logLevel = arg
        elif opt in ("-v", "--verbose"):
            verbose = True

    flStartLog(cnfGetIniValue(configFile, 'Logging', 'LogFile', flConfigFile('.log')),
               cnfGetIniValue(configFile, 'Logging', 'LogLevel', logLevel),
               cnfGetIniValue(configFile, 'Logging', 'Verbose', verbose))

    # Input data:
    inputFile = cnfGetIniValue(configFile, 'Input', 'File')
    inputFormat = cnfGetIniValue(configFile, 'Input', 'Format', os.path.splitext(inputFile)[-1])
    varname = cnfGetIniValue(configFile,'Input','Variable','')
    record = cnfGetIniValue(configFile,'Input','Record',0)
    lvl = cnfGetIniValue(configFile,'Input','Level',0)

    # Output settings - the default is to use the input filename, with
    # the extension replaced by the image format:
    # The smoothing is optional. Set it to the number of grid points to
    # smooth over (recommend the reciprocal of the data resolution in degrees).
    imgfmt = cnfGetIniValue(configFile, 'Output', 'Format','png')
    outputFile = cnfGetIniValue(configFile, 'Output', 'File',
                                "%s.%s" % (os.path.splitext(inputFile)[0], imgfmt))
    smoothing = cnfGetIniValue(configFile, 'Output', 'Smoothing', False)
    cmapName = cnfGetIniValue(configFile, 'Output', 'ColourMap', 'gist_ncar')
    label = cnfGetIniValue(configFile, 'Output', 'Label', '')
    mask = cnfGetIniValue(configFile, 'Output', 'MaskLand', False)
    maskocean = cnfGetIniValue(configFile, 'Output', 'MaskOcean', False)
    fill = cnfGetIniValue(configFile, 'Output', 'FillContours', True)
    title = cnfGetIniValue(configFile,'Plot','Title',None)
    # Load data:
    if inputFormat == '.txt':
        # Attempt to load the dataset:
        try:
            lon,lat,data = grid.grdRead(inputFile)
        except:
            logger.critical("Cannot load input file: %s"%inputFile)
            raise
    elif inputFormat == '.nc':
        try:
            ncobj = nctools.ncLoadFile(inputFile)
            lon = nctools.ncGetDims(ncobj,'lon')
            lat = nctools.ncGetDims(ncobj,'lat')
            data = nctools.ncGetData(ncobj,varname)
            mv = getattr(ncobj.variables[varname],'_FillValue')
            ncobj.close()
        except:
            logger.critical("Cannot load input file: %s"%inputFile)
            raise
        if len(shape(data))==3:
            data = data[record,:,:]
        elif len(shape(data))==4:
            data = data[record,lvl,:,:]

        # Create a masked array:
        datamask = (data==mv)
        data = ma.array(data,mask=datamask)

    else:
        logger.critical("Unknown data format")
        raise IOError

    # Set defaults for the extent of the map to match the data in the
    # input file:
    llLon = min(lon)
    urLon = max(lon)
    llLat = min(lat)
    urLat = max(lat)
    res = 'l'
    dl = 10.

    # Domain settings - can override the default settings:
    domain = cnfGetIniValue(configFile, 'Domain', 'Name', None)
    if domain is not None:
        llLon = cnfGetIniValue(configFile, domain, 'LowerLeftLon', min(lon))
        llLat = cnfGetIniValue(configFile, domain, 'LowerLeftLat', min(lat))
        urLon = cnfGetIniValue(configFile, domain, 'UpperRightLon', max(lon))
        urLat = cnfGetIniValue(configFile, domain, 'UpperRightLat', max(lat))
        res = cnfGetIniValue(configFile, domain, 'Resolution', res)
        dl = cnfGetIniValue(configFile, domain, 'GridInterval', dl)

    [x,y] = meshgrid(lon, lat)

    # Set the scale:
    scaleMin = cnfGetIniValue(configFile, 'Output', 'ScaleMin', 0)
    scaleMax = cnfGetIniValue(configFile, 'Output', 'ScaleMax', 101)
    scaleInt = cnfGetIniValue(configFile, 'Output', 'ScaleInt', 10)
    levels = arange(scaleMin, scaleMax, scaleInt)
    plotField(x,y,data, llLon, llLat, urLon, urLat, res, dl, levels,
              cmapName, smoothing, title=title, xlab='Longitude',
              ylab='Latitude', clab=label, maskland=mask,
              maskocean=maskocean,outputFile=outputFile,fill=fill)

    logger.info("Completed %s"%sys.argv[0])

def plotField(x, y, data, llLon=None, llLat=None, urLon=None, urLat=None,
              res='i', dl=10., levels=10, cmap='jet', smoothing=False,
              title=None, xlab='Longitude', ylab='Latitude', clab=None,
              maskland=False,maskocean=False,outputFile=None,fill=True):
    """
    plotField:
    The main function. This function sets up the map instance (using
    Basemap), contours the data over the map, applies labels, titles
    and a colourbar.
    """
    if NO_BASEMAP:
        return

    pyplot.rcdefaults()
    pyplot.figure()
    if (len(x.shape)<2)and(len(y.shape)<2):
        [xx,yy] = meshgrid(x,y)
    else:
        xx = x
        yy = y
    if llLon is not None:
        llcrnrlon = llLon
    else:
        llcrnrlon = x.min()

    if llLat is not None:
        llcrnrlat = llLat
    else:
        llcrnrlat = y.min()

    if urLon is not None:
        urcrnrlon = urLon
    else:
        urcrnrlon = x.max()

    if urLat is not None:
        urcrnrlat = urLat
    else:
        urcrnrlat = y.max()

    # Apply smoothing:
    if smoothing:
        logger.debug("Applying smoothing to data")
        if hasattr(data,'mask'):
            # Masked array
            data.data[:] = smooth(data.data, smoothing)
        else:
            data = smooth(data, smoothing)

    # Set the colour map:
    if hasattr(cm, cmap):
        cmap = getattr(cm, cmap)
    else:
        cmap = colours.colourMap(cmap, 'stretched')

    meridians = arange(dl*floor(llcrnrlon/dl), dl*ceil(urcrnrlon/dl), dl)
    parallels = arange(dl*floor(llcrnrlat/dl), dl*ceil(urcrnrlat/dl), dl)
    logger.debug("Generating map object")
    # Create the map object:
    try:
        m = Basemap(projection='cyl',
                resolution=res,
                llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat)
    except IOError:
        flLogFatalError(traceback.format_exc().splitlines())

    # Plot the field - by default, we specify the range to extend both
    # above and below the max & min values specified.
    if fill:
        if maskocean:
            try:
                from mpl_toolkits.basemap import maskoceans
            except ImportError:
                logger.debug("Maskoceans module unavailable, skipping this command")
                pass
            else:
                datam = maskoceans(x,y,data,inlands=False)
                m.contourf(x, y, datam, levels, extend='both', cmap=cmap)
        else:
            m.contourf(xx, yy, data, levels, extend='both', cmap=cmap)
        cb = pyplot.colorbar(shrink=0.5, orientation='horizontal', extend='both',
                             ticks=levels[::2], pad=0.05)
        cb.set_label(clab, fontsize=10)
        if cb.orientation=='horizontal':
            for t in cb.ax.get_xticklabels():
                t.set_fontsize(8)
    else:
        CS = m.contour(xx,yy,data,levels,colors='k',linewidth=0.5)
        pyplot.clabel(CS, CS.levels,inline=True,fontsize=7)

    # Add meridians, parallels, coastlines, etc:
    m.drawcoastlines(linewidth=0.5)
    if maskland:
        m.fillcontinents()
    m.drawparallels(parallels, labels=[1,0,0,0], fontsize=9, linewidth=0.2)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=9, linewidth=0.2)
    pyplot.grid(True)

    if title is not None:
        pyplot.title(title)

    # Save the image:
    if outputFile:
        pyplot.savefig(outputFile)
    else:
        pyplot.show()

def plotWindfield(x, y, data, llLon=None, llLat=None, urLon=None, urLat=None,
                  res='i', dl=10., cmapName='jet', smoothing=None,
                  title='Cyclone wind field', xlab='Longitude',
                  ylab='Latitude', clab='Wind speed (m/s)', outputFile=None):
    """
    Plot the wind field for a given grid.
    x, y and data are arrays of the same size.
    (x and y are typically generated by [x,y] = meshgrid(lon,lat))
    This is essentially a wrapper that allows us to set the levels and titles
    without having to code them into the main TCRM functions.
    """

    if llLon is None:
        # Assume in this case that none of the remaining corner
        # co-ords are defined either.
        llLon = min(x[0,:])
        urLon = max(x[0,:])
        llLat = min(y[:,0])
        urLat = max(y[:,0])

    levels = arange(30, 101, 5)

    plotField(x, y, data, llLon, llLat, urLon, urLat, res, dl, levels,
              cmapName, smoothing, title=title, xlab='Longitude',
              ylab='Latitude', clab=clab, outputFile=outputFile)

def plotPressurefield(x, y, data, llLon=None, llLat=None, urLon=None,
                      urLat=None, res='i', dl=10., cmapName='jet_r',
                      smoothing=None, title='Cyclone pressure field',
                      xlab='Longitude', ylab='Latitude',
                      clab='Minimum pressure (hPa)', outputFile=None):
    """
    Plot the pressure field for a given grid.
    x, y and data are arrays of the same size.
    (x and y are typically generated by [x,y] = meshgrid(lon,lat))
    """

    # Convert to hectopascals if necessary (crudely, we check the minimum
    # pressure value to determine whether the data is in hPa or Pa):
    pMin = pressure.min()
    if pMin > 5000.:
        pressure = metutils.convert(pressure, 'Pa', 'hPa')

    if llLon is None:
        llLon = min(xGrid[0,:])
        urLon = max(xGrid[0,:])
        llLat = min(yGrid[:,0])
        urLat = max(yGrid[:,0])

    levels = arange(900, 1020, 5)

    plotField(x, y, data, llLon, llLat, urLon, urLat, res, dl, levels,
              cmapName, smoothing, title=title, xlab='Lonigtude',
              ylab='Latitude', clab=clab, outputFile=outputFile)


def plotArray(x,y,data,llLon=None, llLat=None, urLon=None,
              urLat=None, res='i', dl=10., datarange=(-1.,1.),
              cmap='jet_r',title=None, xlab='Longitude',
              ylab='Latitude', clab=None, maskland=False,
              maskocean=False, outputFile=None):
    """
    Plot a grid of values using pyplot.pcolormesh() to create a pseudocolour
    plot of a 2-dimensional array, with a basemap included.
    """
    #pyplot.rcdefaults()

    if NO_BASEMAP:
        return

    pyplot.figure()
    if (len(x.shape)<2)and(len(y.shape)<2):
        [xx,yy] = meshgrid(x,y)
    else:
        xx = x
        yy = y
    if llLon is not None:
        llcrnrlon = llLon
    else:
        llcrnrlon = x.min()

    if llLat is not None:
        llcrnrlat = llLat
    else:
        llcrnrlat = y.min()

    if urLon is not None:
        urcrnrlon = urLon
    else:
        urcrnrlon = x.max()

    if urLat is not None:
        urcrnrlat = urLat
    else:
        urcrnrlat = y.max()
    meridians = arange(dl*floor(llcrnrlon/dl), dl*ceil(urcrnrlon/dl), dl)
    parallels = arange(dl*floor(llcrnrlat/dl), dl*ceil(urcrnrlat/dl), dl)
    m = Basemap(projection='cyl',
                resolution=res,
                llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat)

    # Set the colour map:
    if hasattr(cm, cmap):
        cmap = getattr(cm, cmap)
    else:
        cmap = colours.colourMap(cmap, 'stretched')

    if maskocean:
        try:
            from mpl_toolkits.basemap import maskoceans
        except ImportError:
            logger.debug("Maskoceans module unavailable, skipping this command")
            pass
        else:
            datam = maskoceans(xx,yy,data,inlands=False)
            m.pcolormesh(xx,yy,datam,edgecolors='None',vmin=datarange[0],vmax=datarange[1],cmap=cmap)
    else:
        m.pcolormesh(xx,yy,data,edgecolors='None',vmin=datarange[0],vmax=datarange[1],cmap=cmap)
    cb = pyplot.colorbar(shrink=0.5, orientation='horizontal', extend='both',pad=0.05)
    if cb.orientation=='horizontal':
        for t in cb.ax.get_xticklabels():
            t.set_fontsize(10)
    if clab is not None:
        cb.set_label(clab)

    m.drawcoastlines(linewidth=0.5)
    if maskland:
        m.fillcontinents(color='white')
    m.drawparallels(parallels, labels=[1,0,0,0], fontsize=9, linewidth=0.2)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=9, linewidth=0.2)
    pyplot.grid(True)

    # Save the image:
    if outputFile:
        pyplot.savefig(outputFile)
    else:
        pyplot.show()

def plotBarb(x,y,u,v,llLon=None,llLat=None,urLon=None,
             urLat=None, res='i', dl=10.,cmapName='jet_r',
             title=None, xlab='Longitude',ylab='Latitude',
             clab=None, plotMagnitude=False, outputFile=None):
    """
    Plot vector data on a grid using wind barbs. Options exist (plotMagnitude)
    to add a filled contour map indicating the magnitude of the data,
    though wind barbs already perform that action through the use of
    barbs and flags.
    """

    pyplot.rcdefaults()
    pyplot.figure()
    if (len(x.shape)<2)and(len(y.shape)<2):
        [xx,yy] = meshgrid(x,y)
    else:
        xx = x
        yy = y

    if llLon is not None:
        llcrnrlon = llLon
    else:
        llcrnrlon = x.min()

    if llLat is not None:
        llcrnrlat = llLat
    else:
        llcrnrlat = y.min()

    if urLon is not None:
        urcrnrlon = urLon
    else:
        urcrnrlon = x.max()

    if urLat is not None:
        urcrnrlat = urLat
    else:
        urcrnrlat = y.max()

    meridians = arange(dl*floor(llcrnrlon/dl), dl*ceil(urcrnrlon/dl), dl)
    parallels = arange(dl*floor(llcrnrlat/dl), dl*ceil(urcrnrlat/dl), dl)
    m = Basemap(projection='cyl',
                resolution=res,
                llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon,
                llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat)

    if plotMagnitude:
        mag = sqrt(u*u+v*v)
        m.contourf(xx, yy, ws, levels, extend='both', cmap=cmap)
        cb = pyplot.colorbar(shrink=0.5, orientation='horizontal', extend='both',
                     ticks=levels[::2], pad=0.05)

        cb.set_label(clab, fontsize=10)
        if cb.orientation=='horizontal':
            for t in cb.ax.get_xticklabels():
                t.set_fontsize(8)

    dx = numpy.mean(numpy.diff(lon))
    xdim = urcrnrlon - llcrnrlon
    nx = xdim/dx
    skip = int(numpy.round(nx/25))
    m.barbs(xx[0::skip,0::skip],yy[0::skip,0::skip],
            u[0::skip,0::skip],v[0::skip,0::skip],
            length=5,linewidth=0.5)

    m.drawcoastlines()

    # Save the image:
    if outputFile:
        pyplot.savefig(outputFile)
    else:
        pyplot.show()

#To run this program from the command line, call the main program:
if __name__ == "__main__":
    _process(sys.argv[1:])
