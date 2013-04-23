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

import os
import logging
from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile
import numpy

logger = logging.getLogger(__name__)

class CalcFrequency:
    """
    Calculate the annual mean frequency of TC events
    based on input dataset. The frequency is calculated
    for the given domain
    """

    def __init__(self, config_file, auto_calc_grid_limit):

        tg_domain = cnfGetIniValue(config_file, 'TrackGenerator',
                                   'gridLimit', '')
        if tg_domain == '':
            self.tg_domain = auto_calc_grid_limit
        else:
            self.tg_domain = eval(tg_domain)
            
        self.output_path = cnfGetIniValue(config_file, 'Output', 'Path')

    def calc(self):
        """
        Calculate the frequency of TC events in a pre-defined
        domain, based on the input dataset and the full range
        of years contained in the 'origin_year' file

        The 'origin_year' file is created in DataProcess.processData()
        and restricts the range to a user-selected range of years.
        """
        logger.info("Calculating annual frequency of TC events")
        origin_year = numpy.array(flLoadFile(os.path.join(self.output_path,
                                                          'process', 'origin_year'),
                                                          '%', ','), dtype='int')
        
        origin_lon_lat = flLoadFile(os.path.join(self.output_path,
                                                 'process', 'origin_lon_lat'),
                                                 '%', ',')
        origin_lon = origin_lon_lat[:, 0]
        origin_lat = origin_lon_lat[:, 1]
        min_year = origin_year.min()
        # Skip last year from average since may contain only partial year record
        max_year = origin_year.max() - 1

        freq_count = numpy.zeros(3000)
        
        for year in range(min_year, max_year + 1):
            freq_count[year] = sum((origin_year == year) & \
                                 (origin_lon > self.tg_domain['xMin']) & \
                                 (origin_lon < self.tg_domain['xMax']) & \
                                 (origin_lat > self.tg_domain['yMin']) & \
                                 (origin_lat < self.tg_domain['yMax']))        
        freq = numpy.mean(freq_count[min_year:max_year + 1])
        # Round to 2 decimal places
        freq = numpy.round(freq*100)/100
        return freq
