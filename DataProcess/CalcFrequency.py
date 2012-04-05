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

Title: CalcTrackDomain.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2011-08-3
Description: Calculates annual genesis frequency in track generator domain
"""

import os, sys, pdb, logging
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile
import numpy

class CalcFrequency:

    def __init__(self, configFile, autoCalc_gridLimit):

        tg_domain = cnfGetIniValue(configFile, 'TrackGenerator', 'gridLimit', '')
        if tg_domain == '':
            self.tg_domain = autoCalc_gridLimit
        else:
            self.tg_domain = eval(tg_domain)
        self.outputPath = cnfGetIniValue(configFile, 'Output', 'Path')

    def calc(self):
        origin_year = numpy.array(flLoadFile(os.path.join(self.outputPath, 'process', 'origin_year'), '%', ','), 
                                  dtype='int')
        origin_lon_lat = flLoadFile(os.path.join(self.outputPath, 'process', 'origin_lon_lat'), '%', ',')
        origin_lon = origin_lon_lat[:,0]
        origin_lat = origin_lon_lat[:,1]
        min_year = origin_year.min()
        # Skip last year from average since may contain only partial year record
        max_year = origin_year.max() - 1

        freq_count = numpy.zeros(3000)
        
        for yr in range(min_year, max_year + 1):
            freq_count[yr] = sum((origin_year==yr) & \
                                 (origin_lon > self.tg_domain['xMin']) & \
                                 (origin_lon < self.tg_domain['xMax']) & \
                                 (origin_lat > self.tg_domain['yMin']) & \
                                 (origin_lat < self.tg_domain['yMax']))        
        freq = numpy.mean(freq_count[min_year:max_year + 1])
        # Round to 2 decimal places
        freq = numpy.round(freq*100)/100
        return freq
