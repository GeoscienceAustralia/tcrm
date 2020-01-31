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


 Title: plotGates.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2007-10-19
 Description:
 Reference:
 SeeAlso:
 Constraints:

 Version: $Rev: 556 $
 ModifiedBy:
 ModifiedDate:
 Modification:

 $Id: PlotGates.py 556 2007-10-23 23:51:09Z carthur $
"""

import Utilities.my_tool as myutils
import math
import pylab
from matplotlib.toolkits.basemap import Basemap

data = myutils.load("N:\\cyclone\\sandpits\\carthur\\landfall\\gates_50km.v0.2.dat",
                    delimiter=',')
lon = data[:,1]
lat = data[:,2]

pylab.clf()

b = 5.
# Plot the region surrounding the coastline rounded to the nearest 5 degree lat/lon.
lonmin = b*math.floor(min(lon)/b)
lonmax = b*math.ceil(max(lon)/b)
latmin = b*math.floor(min(lat)/b)
latmax = b*math.ceil(max(lat)/b)
meridian = (list(range(lonmin, lonmax+b, b)))
parallel = (list(range(latmin, latmax+b, b)))
prj = 'cyl'
res = 'i'
map = Basemap(projection=prj,
              llcrnrlon=lonmin,
              llcrnrlat=latmin,
              urcrnrlon=lonmax,
              urcrnrlat=latmax,
              resolution=res)

map.plot(lon, lat, 'r-')
map.plot(lon, lat, 'k+')

fs = 10

pylab.text(115.9,-31.9, "Perth", fontsize=fs, ha='left', va='top')
pylab.text(114.1,-21.9, "Exmouth", fontsize=fs, ha='left', va='top')
pylab.text(118.5,-20.4, "Port Hedland", fontsize=fs, ha='left', va='top')
pylab.text(130.8,-12.6, "Darwin", fontsize=fs, ha='left', va='top')
pylab.text(145.75,-16.9, "Cairns", fontsize=fs, ha='right', va='top')
pylab.text(149.2,-21.2, "Mackay", fontsize=fs, ha='right', va='top')
pylab.text(153.0,-27.5, "Brisbane", fontsize=fs, ha='right', va='top')
pylab.text(153.1,-30.25, "Coffs Harbour", fontsize=fs, ha='right', va='top')
pylab.text(151.1,-33.9, "Sydney", fontsize=fs, ha='right', va='top')
pylab.text(150.5,-35.3, "Ulladulla", fontsize=fs, ha='right', va='top')
map.fillcontinents(color='0.9')
map.drawcoastlines()
map.drawparallels(parallel, label=[1,0,0,1], fontsize=8)
map.drawmeridians(meridian, label=[1,0,0,1], fontsize=8)
pylab.show()

pylab.savefig("N:\\cyclone\\sandpits\\carthur\\landfall\\gates.v0.2.png")
