'''
Formats track files from http://www.ral.ucar.edu/guidance/realtime/current/
which are 'BDECK' format.

inputPath = folder where the .dat file/s are
outputPath = where the output file is to be written to
config_file = location of tcrm .ini file with 'BDECK' format included
'''

import os
from os.path import join as pjoin, splitext

import numpy as np
from datetime import datetime

from Utilities.tracks2shp import tracks2point, tracks2line
from Utilities.loadData import loadTrackFile

inputPath = "B:/CHARS/B_Wind/data/derived/tc/events/bsh132016/NCAR_RAL_track"
outputPath = "B:/CHARS/B_Wind/data/derived/tc/events/bsh132016/NCAR_RAL_track"
def fmtlon(lonstr):
    if lonstr.endswith('W'):
        lon = 360. + float(lonstr.strip('NSEW')) / -10. 
    else:
        lon = float(lonstr.strip('NSEW')) / 10.
    return lon

def fmtlat(latstr):
    if latstr.endswith('S'):
        lat = float(latstr.strip('NSEW')) / -10.
    else:
        lat = float(latstr.strip('NSEW')) / 10.
    return lat

# Specify the b-deck format
# See http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/shindex.php
# for a detailed description of the format of b-deck files.

bdeck = {
    "delimiter": ",",
    "names" : ("basin", "num", "date", "lat", "lon",
               "vmax", "pcentre", "poci", "rmax", "name"),
    "dtype" : ("|U2", "i", "object", "f8", "f8", "f8",
               "f8", "f8", "f8", "|U10"),
    "usecols" : (0, 1, 2, 6, 7, 8, 9, 17, 19, 27),
    "converters" : {
                'basin': lambda s: s.strip(),
                'num': lambda s: s.strip(),
                'date': lambda s: datetime.strptime(s.strip(), "%Y%m%d%H"),
                'lat': lambda s: fmtlat(s),
                'lon': lambda s: fmtlon(s),
                'vmax': lambda s: float(s.strip()),
                'pcentre': lambda s: float(s.strip()),
                'poci': lambda s: float(s.strip()),
                'rmax': lambda s: float(s.strip()) * 1.852
            },
    "autostrip" : True
    }

source="BDECK"
config_file="B:/CHARS/B_Wind/data/derived/tc/events/bsh132016/TCDebbie.ini"
for f in os.listdir(inputPath):
    inputFile = pjoin(inputPath, f)
    data = np.genfromtxt(inputFile, **bdeck)
    print(inputFile)
    header = 'basin,num,date,lat,lon,vmax,pcentre,poci,rmax,name'
    fmt = '%s,%i,%s,%8.2f,%8.2f,%6.1f,%7.1f,%7.1f,%6.2f,%s'
    outputFile = pjoin(outputPath, f)
    np.savetxt(outputFile, data, fmt=fmt, delimiter=",", header=header)

    fname, ext = splitext(outputFile)
    pt_output_file = fname + '_pt.shp'
    line_output_file = fname + '_line.shp'
    dissolve_output_file = fname + '_dissolve.shp'
    tracks = loadTrackFile(config_file, outputFile, source,
                           calculateWindSpeed=False)

    tracks2point(tracks, pt_output_file)
    tracks2line(tracks, line_output_file)
    tracks2line(tracks, dissolve_output_file, dissolve=True)
