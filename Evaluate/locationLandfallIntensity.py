#!/usr/bin/env python
# coding: utf-8

# # Landfall rates for simulated tropical cyclones 
# 
# This notebook looks at the rate of landfalling tropical cyclones simulated in
# the TCHA. "Landfall" is a loose definition, where we have set out a series of
# gates around the coastline, which are set 50 km off the coast, and each gate
# is 200 km wide. 
# 
# We look at landfall as this is where the impacts of TCs are felt. We will
# compare to historical rates at a later point. 
# 
# The gates are defined in a vector shapefile as line segments. We'll look at the 
# counts, and the mean intensity at landfall to see how well it corresponds to 
# historical landfall rates and intensity. 
# 
# Note however, we'd expect the historical record to be somewhat different to the 
# mean values determined here. Arguably, the historical record should be (statistically) 
# indistiguishable from the range of scenarios. It should appear to be like any other 
# member of the distribution.
# 
# The simulations used here are 35-year simulations of TC activity in the Australian region.
# This corresponds to the length of historical record used as the input record
# to the TCHA (1981-2016). The TCs are simulated with a 3-hour timestep. 

import os
import sys
from os import walk
from os.path import join as pjoin
import matplotlib.pyplot as plt
from Utilities import track

from shapely.geometry import LineString, Point

import numpy as np
import pandas as pd
import geopandas as gpd

import logging
import seaborn as sns
sns.set_style('whitegrid')

LOGGER = logging.getLogger(__name__)

# Start with reading in the gates into a `GeoDataFrame`, and adding some additional 
# attributes. This `GeoDataFrame` will be duplicated for each simulation, then 
# aggregated for the summary statistics.

def loadGateFile(gateFile):
    try:
        gates = gpd.read_file(gateFile)
    except:
        LOGGER.error(f"Cannot open {gateFile}")
        LOGGER.error("Check the file exists and is a valid shapefile")
        sys.exit()

    # Add some extra attributes
    gates['sim'] = 0
    gates['count'] = 0
    gates['meanlfintensity'] = np.nan
    gates['minlfintensity'] = np.nan
    gates['cat1'] = 0
    gates['cat2'] = 0
    gates['cat3'] = 0
    gates['cat4'] = 0
    gates['cat5'] = 0

    return gates

# Define a function to read in the track files. The track files are the
# TCRM format files, which are netCDF files, with a heirarchical structure. 
# The `tracks.ncReadTrackData` function does this, but we still need to convert 
# each track to a `GeoDataFrame` and add a `geometry` attribute that 
# represents the line of each time step in the track history. We also assign
# a category attribute to each segment, which at this point in time is simply 
# based on categorising the central pressure according to the Bureau of 
# Meteorology's TC intensity scale.

def readTracks(trackFile):
    """
    Read a track file and create a `GeoPandas.GeoDataFrame` from the data, with
    separate polyline features for each unique track in the file. Also adds a
    nominal intensity category based on the central pressure value.
    
    :param str trackFile: Path to a TCRM-format track file
    
    :returns: `GeoPandas.GeoDataFrame` of tracks
    """

    LOGGER.debug(f"Loading track data from {trackFile}")
    tracks = track.ncReadTrackData(trackFile)
    trackgdf = []
    for t in tracks:
        segments = []
        for n in range(len(t.data) - 1):
            segment = LineString([[t.Longitude[n], t.Latitude[n]],[t.Longitude[n+1], t.Latitude[n+1]]])
            segments.append(segment)
        gdf = gpd.GeoDataFrame.from_records(t.data[:-1])
        gdf['geometry'] = segments
        gdf['category'] = pd.cut(gdf['CentralPressure'], 
                                bins=[0, 930, 955, 970, 985, 990, 1020], 
                                labels=[5,4,3,2,1,0])
        trackgdf.append(gdf)
    trackgdf = pd.concat(trackgdf)
    return trackgdf

def isLeft(line, point):
    """
    Test whether a point is to the left of a (directed) line segment. 
    
    :param line: :class:`Shapely.geometry.LineString` of the line feature being tested
    :param point: :class:`Shapely.geometry.Point` being tested
    
    :returns: `True` if the point is to the left of the line segment, `False` otherwise 
    """
    start = Point(line.coords[0])
    end = Point(line.coords[1])

    det = (end.x - start.x) * (point.y - start.y) - (end.y - start.y) * (point.x - start.x)
    if det > 0: return True
    if det <= 0: return False

def enters(feature, line):
    start = Point(line.coords[0])
    end = Point(line.coords[1])
    if not start.within(feature):
        return True
    else:
        return False

def isLandfall(feat, tracks):
    """
    Determine the count of tracks crossing a gate segment. 
    
    :param feat: `shapely.Geometry.Polygon` 
    :param tracks: `GeoPandas.GeoDataFrame` of tracks, where the `geometry` field is a
                   `LineSegment`
                   
    """
    LOGGER.debug(f"Determining landfalls for track collection")

    crossings = tracks.crosses(feat.geometry)
    landfall = []
    for t in tracks[crossings].itertuples():
        if enters(feat.geometry, Point(t.geometry.coords[0])):
            landfall.append(True)
        else:
            landfall.append(False)

    return tracks[crossings][landfall]


# This function counts the number of track segments that cross the coastal gates. 

def countCrossings(gdf, tracks, sim):
    """
    Count the crossing rate of all gates for all tracks in a given simulation.
    
    :param gdf: `GeoDataFrame` containing the landfall gates
    :param tracks: `GeoDataFrame` containing `Track` objects
    :param int sim: Ordinal simulation number

    """
    gdf['sim'] = sim
    for i, feat in enumerate(gdf.itertuples(index=False)):
        ncrossings = 0
        l = isLandfall(feat, tracks)
        ncrossings = len(l)
        if ncrossings > 0:
            gdf['count'].iloc[i] = ncrossings
            gdf['meanlfintensity'].iloc[i] = l['CentralPressure'].mean()
            gdf['minlfintensity'].iloc[i] = l['CentralPressure'].min()
            cathist, bins = np.histogram(l['category'].values, bins=[0,1,2,3,4,5, 6])
            gdf['cat1'].iloc[i] = cathist[0]
            gdf['cat2'].iloc[i] = cathist[1]
            gdf['cat3'].iloc[i] = cathist[2]
            gdf['cat4'].iloc[i] = cathist[3]
            gdf['cat5'].iloc[i] = cathist[4]
        else:
            gdf['count'].iloc[i] = 0
            gdf['meanlfintensity'].iloc[i] = np.nan
            gdf['minlfintensity'].iloc[i] = np.nan
            gdf['cat1'].iloc[i] = 0
            gdf['cat2'].iloc[i] = 0
            gdf['cat3'].iloc[i] = 0
            gdf['cat4'].iloc[i] = 0
            gdf['cat5'].iloc[i] = 0
            
    return gdf

# A bunch of helper functions to return quantiles
def q10(x): return x.quantile(0.1)
def q90(x): return x.quantile(0.9)
def q25(x): return x.quantile(0.25)
def q75(x): return x.quantile(0.75)



# Get the list of track files from the input data path (this'll become a 
# config option in a scripted version of theis notebook). Just find all 
# files in the directory that end with "nc" - assume that you're pointing
# to a directory that only has track files.

def loadLandfallRates(datapath, gdf):

    if not os.path.isdir(datapath):
        LOGGER.error(f"{datapath} is not a valid directory")
        raise IOError(f"{datapath} is not a valid directory")
    filelist = []
    for (dirpath, dirnames, filenames) in walk(datapath):
        filelist.extend([fn for fn in filenames if fn.endswith('nc')])
        break
    nfiles = len(filelist)
    if nfiles==0:
        LOGGER.error(f"There are no valid trackfiles in the input path {datapath}")
        raise IOError(f"There are no valid trackfiles in the input path {datapath}")
    LOGGER.info(f"There are {nfiles} track files in {datapath}")


    # Now we loop through all the gates and determine the landfall rates for each simulation. 

    gatedflist = []
    for sim, f in enumerate(filelist):
        q, r = np.divmod(sim*10, len(filelist))
        if r==0:
            LOGGER.info(f"Loaded {sim} events ({q*10}%)")
        gatedf = gdf.copy()
        tracks = readTracks(pjoin(datapath,f))
        gatedf = countCrossings(gatedf, tracks, sim)
        gatedflist.append(gatedf)

    gatesummary = pd.concat(gatedflist)
    LOGGER.info("Grouping data")


    gs = gatesummary.groupby('LGA_CODE19').agg({'count':['sum',np.nanmean, np.nanstd, 'min', 'max', q10, q90],
                                        'cat1': ['sum',np.nanmean, np.nanstd],
                                        'cat2': ['sum',np.nanmean, np.nanstd],
                                        'cat3': ['sum',np.nanmean, np.nanstd],
                                        'cat4': ['sum',np.nanmean, np.nanstd],
                                        'cat5': ['sum',np.nanmean, np.nanstd],
                                        'meanlfintensity':[np.nanmean, q10, q25, q75, q90],
                                        'minlfintensity':[np.nanmean,'min', np.nanstd]}, as_index=False)

    gs.columns = ['_'.join(col).strip() for col in gs.columns.values]
    gs.reset_index(col_level=1)
    gs.columns = gs.columns.get_level_values(0)


    gatedata = gates[['LGA_CODE19', 'LGA_NAME19', 'geometry']].join(gs, on='LGA_CODE19')
    
    gatedata.rename(columns={
        'count_sum':'sum',
        'count_nanmean': 'count_mean',
        'count_nanstd': 'count_std',
        'cat1_nanmean': 'cat1_mean',
        'cat2_nanmean': 'cat2_mean',
        'cat3_nanmean': 'cat3_mean',
        'cat4_nanmean': 'cat4_mean',
        'cat5_nanmean': 'cat5_mean',
        'meanlfintensity_nanmean': 'mlfimean',
        'meanlfintensity_q10': 'mlfiq10',
        'meanlfintensity_q25': 'mlfiq25',
        'meanlfintensity_q75': 'mlfiq75',
        'meanlfintensity_q90': 'mlfiq90',
        'minlfintensity_nanmean': 'minlfi',
        'minlfintensity_min': 'minlfim'
    }, inplace=True)

    return gatedata

def plotLandfallIntensity(df, datapath):

    LOGGER.info("Plotting landfall intensity")
    width=0.4
    fig, ax = plt.subplots(3,1,figsize=(12,16),sharex=True)
    cat12 = np.add(df['cat1_mean'], df['cat2_mean']).tolist()
    cat123 = np.add(cat12, df['cat3_mean']).tolist()
    cat1234 = np.add(cat123, df['cat4_mean']).tolist()
    ax[0].bar(df['LGA_CODE19'], df['cat1_mean'], color='b', label="Cat 1")
    ax[0].bar(df['LGA_CODE19'], df['cat2_mean'], bottom=df['cat1_mean'], color='g', label='Cat 2')
    ax[0].bar(df['LGA_CODE19'], df['cat3_mean'], bottom=cat12, color='y', label='Cat 3')
    ax[0].bar(df['LGA_CODE19'], df['cat4_mean'], bottom=cat123, color='orange', label='Cat 4')
    ax[0].bar(df['LGA_CODE19'], df['cat5_mean'], bottom=cat1234, color='r', label='Cat 5')

    ax[0].legend()
    ax[0].set_ylabel("Mean number of TCs")
    ax[1].plot(df['LGA_CODE19'], df['minlfi'], 
               label='Minimum landfall intensity')
    ax[1].plot(df['LGA_CODE19'], df['mlfimean'], color='r', 
               label='Mean landfall intensity')
    ax[1].fill_between(df['LGA_CODE19'], df['mlfiq10'],
                       df['mlfiq90'], color='r', alpha=0.25)

    ax[1].legend(loc=2)
    ax[1].set_ylim((900, 1020))
    ax[1].set_ylabel("Pressure (hPa)")
    ax[2].plot(df['LGA_CODE19'], df['sum']/10000)
    ax[2].fill_between(df['LGA_CODE19'], df['count_q90']/10000,
                       df['count_q10']/10000, alpha=0.25)
    ax[2].set_xlim((0,81))
    ax[2].set_xticks(np.arange(0,80))
    ax[2].set_yticks(np.arange(0,5.1,.5))
    ax[2].set_xticklabels(df['LGA_NAME19'][::2], rotation='vertical')
    ax[2].set_ylabel("Mean proportion of landfall")
    plt.savefig(os.path.join(datapath, "landfall_intensity_lga.png"), bbox_inches='tight')

def plotIntensityDistribution(df, plotpath):
    """
    Plot a distribution of intenstity at landfall

    :param df: `DataFrame` containing the gate data
    :param str plotpath: Path to where the plots will be saved/
    """
    LOGGER.info("Plotting landfall intensity distribution")
    width=0.4
    fig, ax = plt.subplots(1,1, figsize=(12,6), sharex=True)
    cat12 = np.add(df['cat1_mean'], df['cat2_mean']).tolist()
    cat123 = np.add(cat12, df['cat3_mean']).tolist()
    cat1234 = np.add(cat123, df['cat4_mean']).tolist()
    ax.bar(df['LGA_CODE19'], df['cat1_mean'], color='b', label="Cat 1")
    ax.bar(df['LGA_CODE19'], df['cat2_mean'], bottom=df['cat1_mean'], color='g', label='Cat 2')
    ax.bar(df['LGA_CODE19'], df['cat3_mean'], bottom=cat12, color='y', label='Cat 3')
    ax.bar(df['LGA_CODE19'], df['cat4_mean'], bottom=cat123, color='orange', label='Cat 4')
    ax.bar(df['LGA_CODE19'], df['cat5_mean'], bottom=cat1234, color='r', label='Cat 5')

    ax.legend()
    ax.set_ylabel("Number of TCs")
    ax.set_xlim((0,81))
    ax.set_xticks(np.arange(0,81))
    ax.set_yticks(np.arange(0,2.5,.2))
    ax.set_xticklabels(gatedata['LGA_NAME19'], rotation='vertical')
    ax.set_ylabel("Mean rate of landfall")
    plt.savefig(os.path.join(plotpath, "mean_landfall_rate_intensity_lga.png"), bbox_inches='tight')

def plotLandfallMap(df, plotpath):

    bins = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0]
    ax = df.plot(column='count_mean', 
                 figsize=(9,13.5), 
                 cmap='OrRd', 
                 legend=True,
                 edgecolor='0.5',
                 scheme='UserDefined', 
                 classification_kwds={'bins':bins})
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Mean TC occurence (TC/year")
    plt.savefig(os.path.join(plotpath, "mean_landfall_rate_intensity_lga_map.png"),
                bbox_inches='tight')

if __name__ == "__main__":
    import argparse

    logFormat = "%(asctime)s: %(funcName)s: %(message)s"
    logging.basicConfig(level='INFO', 
                        format=logFormat,
                        filename='landfallIntensity.log', 
                        filemode='w',
                        datefmt="%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(getattr(logging, 'INFO'))
    formatter = logging.Formatter(logFormat,
                                  datefmt='%H:%M:%S', )
    console.setFormatter(formatter)
    LOGGER.addHandler(console)
    LOGGER.info(f"Started {sys.argv[0]} (pid {os.getpid()})")

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path")
    args = parser.parse_args()

    basepath = args.path

    gatefile = "/g/data/w85/QFES_SWHA/hazard/input/QLD_LGA_2019.shp"
    gates = loadGateFile(gatefile)
    
    datapath = pjoin(basepath, "tracks")
    plotpath = pjoin(basepath, "plots")
    processpath = pjoin(basepath, "process")
    gatedata = loadLandfallRates(datapath, gates)
    plotIntensityDistribution(gatedata, plotpath)
    plotLandfallIntensity(gatedata, plotpath)
    plotLandfallMap(gatedata, plotpath)
    gatedata_nogeom = pd.DataFrame(gatedata.drop(columns='geometry'))
    gatedata_nogeom.to_csv(os.path.join(processpath, "simulated_landfall_rates.csv"), index=False)
    gatedata.to_file(os.path.join(processpath, "simulated_landfall_rates.shp"))
