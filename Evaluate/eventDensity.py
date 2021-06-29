#!/usr/bin/env python
# coding: utf-8
"""
:mod:`eventFrequency` -- calculate spatial mean frequency
=========================================================

.. module:: eventFrequency
   :synopsis: Calculate frequency of TC tracks on a grid (TCs/year)

.. moduleauthor: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import time
import logging
import getpass

from os.path import join as pjoin
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import xarray as xr

from cartopy import crs as ccrs
from Utilities.track import ncReadTrackData
from Utilities.config import ConfigParser
from PlotInterface.maps import selectColormap

from shapely.geometry import LineString, Polygon
from shapely.geometry import box as sbox

import seaborn as sns

ISO_FORMAT = '%Y-%m-%d %H:%M:%S'

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

sns.set_style('whitegrid')

palette = [(1.000, 1.000, 1.000),
           (0.000, 0.627, 0.235),
           (0.412, 0.627, 0.235),
           (0.663, 0.780, 0.282),
           (0.957, 0.812, 0.000),
           (0.925, 0.643, 0.016),
           (0.835, 0.314, 0.118),
           (0.780, 0.086, 0.118)]
cmap = sns.blend_palette(palette, as_cmap=True)


def filter_tracks_domain(df, minlon=90, maxlon=180, minlat=-40, maxlat=0):
    """
    Takes a `DataFrame` and filters on the basis of whether the track interscts
    the given domain, which is specified by the minimum and maximum longitude and 
    latitude.
    
    NOTE: This assumes the tracks and bounding box are in the same geographic 
    coordinate system (i.e. generally a latitude-longitude coordinate system). 
    It will NOT support different projections (e.g. UTM data for the bounds and
    geographic for the tracks).
    
    NOTE: This doesn't work if there is only one point for the track. 
    
    :param df: :class:`pandas.DataFrame` that holds the TCLV data
    :param float minlon: minimum longitude of the bounding box
    :param float minlat: minimum latitude of the bounding box
    :param float maxlon: maximum longitude of the bounding box
    :param float maxlat: maximum latitude of the bounding box
    """
    domain = sbox(minlon, minlat, maxlon, maxlat, ccw=False)
    tracks = df.groupby('num')
    tempfilter = tracks.filter(lambda x: len(x) > 1)
    filterdf = tempfilter.groupby('num').filter(lambda x: LineString(zip(x['lon'], x['lat'])).intersects(domain))
    return filterdf

def readTracks(trackFile):
    """
    Read all the tracks from a given track file (in TCRM-netcdf format), add a
    geometry and return as a :class:`geopandas.GeoDataFrame`. Includes
    calculation of the pressure deficit, a TC intensity category, and normalised
    intensity.

    :param str trackFile: path to a TCRM-netcdf format track file

    :returns: :class:`geopandas.GeoDataFrame` containing all the tracks from the
    file.
    """
    tracks = ncReadTrackData(trackFile)
    trackgdf = []
    for t in tracks:
        segments = []
        for n in range(len(t.data) - 1):
            segment = LineString([[t.Longitude[n], t.Latitude[n]],
                                  [t.Longitude[n+1], t.Latitude[n+1]]])
            
            segments.append(segment)
        gdf = gpd.GeoDataFrame.from_records(t.data[:-1])
        gdf['geometry'] = segments
        gdf['category'] = pd.cut(gdf['CentralPressure'], 
                                bins=[0, 930, 955, 970, 985, 990, 1020], 
                                labels=[5,4,3,2,1,0])
        # Calculate pressure difference and normalised intensity
        gdf['pdiff'] = gdf.EnvPressure - gdf.CentralPressure
        gdf['ni'] = gdf.pdiff / gdf.pdiff.max()
        trackgdf.append(gdf)
    trackgdf = pd.concat(trackgdf)
    return trackgdf


def createGrid(xmin, xmax, ymin, ymax, width, length):
    """
    Create a grid over which to perform the analysis.
    """
    cols = list(np.arange(xmin, xmax + width, width))
    rows = list(np.arange(ymin, ymax + length, length))
    gridid = 0
    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(Polygon([(x, y),
                                    (x + width, y),
                                    (x + width, y + length),
                                    (x, y + length)]))
    gridid = np.arange(len(polygons))
    grid = gpd.GeoDataFrame({'gridid': gridid,
                             'geometry': polygons})
    return grid


def gridCount(grid, tracks, dims):
    """
    Given a grid and a collection of tracks, calculate the number of
    unique events that enter each box in the grid. Return the data as
    an array with shape defined by `dims`

    :param grid: `geopandas.GeoDataFrame` containing the grid
    boxes. Assumes this is a rectangular grid of boxes 

    :param tracks: `geopandas.GeoDataFrame` of the tracks, with the
    geometry stored as a LineString

    """

    df = gpd.sjoin(grid, tracks)
    df2 = df.groupby('gridid')['CycloneNumber'].nunique()
    dfcount = grid.merge(df2, how='left', left_on='gridid', right_index=True)
    dfcount.rename(columns={'CycloneNumber':'count'}, inplace=True)
    dfcount['count'] = dfcount['count'].fillna(0)
    return dfcount['count'].values.reshape(dims)


def calculateGridFrequency(configFile, plot=False):
    """
    Set up a grid and process all available track files for the track
    frequency.
    
    :param str configFile: Path to TCRM configuration file
    :param bool plot: If `True` plot the resulting data and save image

    :returns: None

    """

    config = ConfigParser()
    config.read(configFile)
    outputPath = config.get('Output', 'Path')
    trackPath = pjoin(outputPath, 'tracks')
    plotPath = pjoin(outputPath, 'plots', 'stats')
    dataPath = pjoin(outputPath, 'process')
    
    # Implicitly assume we are working with TCRM-format track files
    filelist = [f for f in os.listdir(trackPath) if f.endswith('nc')]
    nfiles = len(filelist)
    print(f"There are {nfiles} track files")

    gridLimit = config.geteval('Region', 'gridLimit')
    gridSpace = config.geteval('Region', 'gridSpace')

    # The gridSpace is refined for finer resolution than the
    # statistics used in the track generation.
    lon = np.arange(gridLimit['xMin'], gridLimit['xMax'], gridSpace['x']/5.)
    lat = np.arange(gridLimit['yMin'], gridLimit['yMax'], gridSpace['y']/5.)
    X, Y = np.meshgrid(lon, lat)
    dims = X.shape

    grid = createGrid(gridLimit['xMin'], gridLimit['xMax'], gridLimit['yMin'],
                      gridLimit['yMax'], gridSpace['x']/5., gridSpace['y']/5.)
    grarray = np.empty((nfiles, *dims))

    for sim, f in enumerate(filelist):
        q, r = np.divmod(sim*10, nfiles)
        if r==0: print(f"{q*10}% complete")
        tracks = readTracks(pjoin(trackPath, f))
        grarray[sim, :, :] = gridCount(grid.copy(), tracks, dims)

    outputFile = pjoin(dataPath, "mean_track_density.nc")
    log.info(f"Saving track density data to {outputFile}")
    # TODO: Add data range, standard_name (?) attributes to the var
    da = xr.DataArray(np.nanmean(grarray, axis=0).T,
                      coords=[lat, lon],
                      dims=['lat', 'lon'],
                      attrs=dict(long_name="Mean annual TC frequency",
                                 units="year-1"))
    # TODO: More attributes required here
    ds = xr.Dataset({"frequency": da},
                    attrs=dict(description="Mean annual TC frequency",
                               created_by=getpass.getuser(),
                               created_on=time.strftime(ISO_FORMAT, time.localtime())))
    ds.to_netcdf(outputFile)

    if plot:
        plotData(plotPath, ds)
    return

def plotData(plotPath, ds):
    """
    Plot the data using in-built plotting routines from xarray
    
    :param str plotPath: destination folder for the output image
    :param ds: :class:`xarray.Dataset` containing the frequency data

    :returns: None

    """
    prj = ccrs.PlateCarree()
    fig, ax = plt.subplots(1, 1, figsize=(8, 6),
                           subplot_kw={'projection':prj})
    ds.frequency.plot(ax=ax, transform=prj)
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False

    plt.savefig(pjoin(plotPath, "mean_TC_frequency.png"), bbox_inches="tight")
    
if __name__ == "__main__":
    configFile = "/g/data/w85/QFES_SWHA/configuration/tcrm/hazard/QLD.GROUP2.RCP85_2081-2100.ini"
    calculateGridFrequency(configFile, plot=True)




