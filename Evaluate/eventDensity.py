#!/usr/bin/env python
# coding: utf-8

import os
import logging

from os import walk
from os.path import join as pjoin
from pathlib import Path

import cftime
from datetime import datetime
datefmt = "%Y-%m-%d %H:%M"

import matplotlib.pyplot as plt
import geopandas as gpd
import xarray as xr

from cartopy import crs as ccrs
from Utilities import track
from shapely.geometry import LineString, Point, Polygon
import shapely.geometry as sg
from shapely.geometry import box as sbox
import numpy as np
import pandas as pd

import scipy.stats as stats
import seaborn as sns

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

sns.set_style('whitegrid')
sns.set_context('talk')

palette = [(1.000, 1.000, 1.000), (0.000, 0.627, 0.235), (0.412, 0.627, 0.235), 
           (0.663, 0.780, 0.282), (0.957, 0.812, 0.000), (0.925, 0.643, 0.016), 
           (0.835, 0.314, 0.118), (0.780, 0.086, 0.118)]
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
    tempfilter.head()
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
    tracks = track.ncReadTrackData(trackFile)
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

def createGrid(xmin, xmax, ymin, ymax, wide, length):
    """
    Create a grid over which to perform the analysis.
    """
    cols = list(np.arange(xmin, xmax + wide, wide))
    rows = list(np.arange(ymin, ymax+length, length))
    gridid = 0
    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(Polygon([(x, y),
                                    (x + wide, y),
                                    (x + wide, y + length),
                                    (x, y + length)]))
    gridid = np.arange(len(polygons))
    grid = gpd.GeoDataFrame({'gridid': gridid,
                             'geometry': polygons})
    return grid


datapath = "C:/WorkSpace/data/tcha/tracks"

filelist = [f for f in os.listdir(datapath) if f.endswith('nc')]
nfiles = len(filelist)
log.info(f"There are {nfiles} track files")

minlon = 130
maxlon = 160
minlat = -30
maxlat = -5
dx = .2
dy = .2

lon = np.arange(minlon, maxlon, dx)
lat = np.arange(minlat, maxlat, dy)
xx, yy = np.meshgrid(lon, lat)
grid = createGrid(minlon, maxlon, minlat, maxlat, dx, dy)
dims = (int((maxlon - minlon)/dx), int((maxlat - minlat)/dy))

grarray = np.empty((nfiles, *dims))
for sim, f in enumerate(filelist):
    #if sim >= 500: break
    q, r = np.divmod(sim*10, len(filelist))
    if r==0:
        log.info(f"{q*10}% complete")
    griddf = grid.copy()
    tracks = readTracks(pjoin(datapath,"tracks",f))
    dfjoin = gpd.sjoin(griddf, tracks)
    df2 = dfjoin.groupby('gridid')['CycloneNumber'].nunique()
    dfcount = griddf.merge(df2, how='left', left_on='gridid', right_index=True)
    dfcount.rename(columns={'CycloneNumber':'count'}, inplace=True)
    dfcount['count'] = dfcount['count'].fillna(0)
    grarray[sim, :, :] = dfcount['count'].values.reshape(dims)

ax = plt.axes(projection=ccrs.PlateCarree())
cb = plt.contourf(xx, yy, np.nanmean(grarray, axis=0).T, cmap=cmap, transform=ccrs.PlateCarree())
plt.colorbar(cb)
ax.coastlines()
ax.gridlines()
plt.savefig("mean_TC_frequency.png", bbox_inches="tight")

da = xr.DataArray(np.nanmean(grarray, axis=0), coords=[lon, lat], dims=['lon', 'lat'],
                  attrs=dict(long_name="Mean annual TC frequency",
                             units="1/year"))
ds = xr.Dataset({"frequency": da},
                attrs=dict(description="Mean annual TC frequency"))

ds.to_netcdf(os.path.join(os.path.dirname(datapath),"mean_track_density.nc"))



