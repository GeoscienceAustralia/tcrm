#!/usr/bin/env python
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

Title: generate_localities_database.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2012-02-16
Description: Generates country/division/locality database from text file provided on the
"world-gazetteer.com" website.  The database is used by the configuration editor GUI
to determine place names and location coordinates.

Data licence details from website:

      "all data (c) 2005 by Stefan Helders
      This project is to be regarded as a free data provider.  Some requests let
      me precise the copyright information. If you use its data or images the only
      thing I ask you is to promote this site. If you would like to republish the
      data presented here, please do not change the data and use a copyright note
      as described as follows: (c) by Stefan Helders www.world-gazetteer.com
      all data (c) 2005 by Stefan Helders"
"""

import sqlite3
import codecs
import os.path
import numpy

# Settings
dBaseName = 'localities.dat'
inputDataFile = 'dataen.txt'
noDivisionStr = 'No division (select for localities with no division)'

# If no locality database found, then create one
if os.path.isfile(dBaseName):
    print("Database already exists")
else:
    print("Generating database.....")
    conn = sqlite3.connect(dBaseName)

    c = conn.cursor()

    # Create table
    c.execute('''create table localities (placeID text, placename text, placetype text, population real, lat real, lon real, parentcountry text, parentdivision text)''')

    f = codecs.open(inputDataFile, encoding='utf-8')
    for line in f:
        z = line.split('\t')
        if len(z) == 12:
            lat_str = z[6]
            lon_str = z[7]

            if (z[4] == 'locality') and (len(lat_str)==0 or len(lon_str)==0):
                # Skip if locality has no lat/lon coordinates
                pass
            elif (z[4] == 'locality') and (lat_str == '0') and (lon_str == '9999'):
                # Skip if locality has missing value codes for lat/lon coordinates
                pass
            else:
                # Convert lat/lon format
                if len(lat_str) > 0:
                    lat_str = "%.2f" % (float(lat_str)/100)
                if len(lon_str) > 0:
                    lon_str = "%.2f" % (numpy.mod(float(lon_str)/100, 360))
                if (z[4] == 'locality') and (z[9] == ''):
                    # For countries like Niue with no divisions
                    z[9] = noDivisionStr
                c.execute('insert into localities values (?,?,?,?,?,?,?,?)', (z[0], z[1], z[4], z[5], lat_str, lon_str, z[8], z[9]))
    conn.commit()

    # Remove divisions that contain no localities
    c.execute('select placename from localities where placetype=?', ('country',))
    countries = [z[0] for z in c.fetchall()]
    countries.sort()

    # Delete divisions that contain no localities
    count = 0
    for country in countries:
        count = count + 1
        if numpy.mod(count, 5) == 0:
            print('Progress:  %(#).1f' % {"#": (float(count) / len(countries)) * 100} + '%')
        c.execute('select placename from localities where parentcountry=? and placetype<>?', (country, 'locality'))
        divisions = [z[0] for z in c.fetchall()]

        for division in divisions:
            c.execute('select placename from localities where parentcountry=? and parentdivision=? and placetype=?', (country, division, 'locality'))
            locations = [z[0] for z in c.fetchall()]
            if len(locations)==0:
                c.execute('delete from localities where placename=? and parentcountry=? and placetype<>?', (division, country, 'locality'))
    conn.commit()

    # Handle cases when localities have no division (e.g. Niue)
    for country in countries:
        c.execute('select placename from localities where parentcountry=? and parentdivision=?', (country, noDivisionStr))
        noDivisionMatch = [z[0] for z in c.fetchall()]
        if len(noDivisionMatch) > 0:
            c.execute('insert into localities values (?,?,?,?,?,?,?,?)', ('', noDivisionStr, 'division', '', '', '', country, ''))
    conn.commit()
