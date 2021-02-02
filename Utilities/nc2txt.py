"""
    Tropical Cyclone Risk Model (TCRM) - Version 2.0
    Copyright (C) 2015 Commonwealth of Australia (Geoscience Australia)

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

 Title: nc2txt
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 02/07/11 11:10:AM
 Description:

 Version :$Rev: 642 $

 $Id: nc2txt.py 642 2012-02-21 07:54:04Z nsummons $
"""
import os
import sys
import logging
import getopt

import numpy as np
import Utilities.grid as grid
import Utilities.nctools as nctools
__version__ = "%Id$"

def usage():
    """
    Short description of the program and how to run it
    """
    print("nc2txt.py:")
    print(__version__)
    print("Usage:")
    print("nc2txt.py -f <filename> -v <variable> [-r <record dimension>]")
    print("Options:")
    print("-f, --filename <input filename> - convert this file to an ascii grd format")
    print("-v, --variable <variable name> - convert this variable")
    print("-m, --missingvalue <missing value> - missing value as set in the input ncfile")
    print("-r, --record <record dimension name> - if the variable has more")
    print("      than 2 dimensions, this specifies the name of the third")
    print("      dimension")

def process(argv):
    recdim = None
    mv = -9999.
    try:
        opts, args = getopt.getopt(argv, "f:hm:r:v:",
                                   ["filename=",
                                    "help",
                                    "missingvalue=",
                                    "record=",
                                    "variable="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-f", "--filename"):
            filename = arg
        elif opt in ("-m", "--missingvalue"):
            mv = arg
        elif opt in ("-v", "--variable"):
            variable = arg
        elif opt in ("-r", "--record"):
            recdim = arg

    ncobj = nctools.ncLoadFile(filename)

    lat = nctools.ncGetDims(ncobj, 'lat')
    lon = nctools.ncGetDims(ncobj, 'lon')
    delta = lon[1] - lon[0]

    # Fix incorrectly reported corner of lower left pixel
    lon = lon - delta/2.
    lat = lat - delta/2.
    if recdim:
        recval = nctools.ncGetDims(ncobj, recdim)

    data = nctools.ncGetData(ncobj, variable)
    ncobj.close()
    if recdim:
        for i, v in enumerate(recval):
            outputfile = "%s.%s.%s"%(os.path.splitext(filename)[0],
                                     repr(recval[i]), 'txt')
            print("Saving data to %s"%outputfile)
            grid.grdSave(outputfile, np.flipud(data[i]), lon, lat,
                         delta, delimiter=' ', nodata=mv, fmt='%6.2f')
    else:
        outputfile = "%s.%s"%(os.path.splitext(filename)[0], 'txt')
        print("Saving data to %s"%outputfile)
        grid.grdSave(outputfile, np.flipud(data), lon, lat,
                     delta, delimiter=' ', nodata=mv, fmt='%6.2f')

if __name__ == '__main__':
    process(sys.argv[1:])
