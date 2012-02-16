#!/usr/bin/env python
"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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


Title: generate_cities_database.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2012-02-16
Description: Generates countries/cities database from text file provided on the
"world-gazetteer.com" website.  The database is used by the TCRM graphical user 
interface to determine place names and location coordinates.

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
import unicodedata
import os.path

# Settings
dbasename = 'cities.dat'
inputdatafile = 'dataen.txt'

# If no locality database found, then create one
if os.path.isfile(dbasename):
   print "Database already exists"
else:
   print "Generating database....."
   conn = sqlite3.connect(dbasename)
   
   c = conn.cursor()

   # Create table
   c.execute('''create table localities (placename text, placetype text, lat real, lon real, parentcountry text, parentdivision)''')

   f = codecs.open('dataen.txt', encoding='utf-8')

   for line in f:
      linedata = unicodedata.normalize('NFKD', line).encode('ascii', 'ignore')
      z = linedata.split('\t')
      if len(z) == 12:
         c.execute('insert into localities values (?,?,?,?,?,?)', (z[1], z[4], z[6], z[7], z[8], z[9]))
   conn.commit()
   print "Done!"