"""
:mod:`CalcFrequency` -- Calculate annual genesis frequency
==========================================================

Calculate annual genesis frequency of TCs within the track
generation domain.

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

"""

import os
import logging
import numpy

from Utilities.files import flLoadFile
from Utilities.config import ConfigParser

logger = logging.getLogger(__name__)

class CalcFrequency:
    """
    Calculate the annual mean frequency of TC events
    based on input dataset. The frequency is calculated
    for the given domain

    :type  tg_domain: :class:`dict`
    :param tg_domain: the domain where the tracks will be generated.
                      The :class:`dict` should contain the keys :attr:`xMin`,
                      :attr:`xMax`, :attr:`yMin` and :attr:`yMax`. The *x*
                      variable bounds the longitude and the *y* variable bounds
                      the latitude.
    """

    def __init__(self, configFile, auto_calc_grid_limit):
        """
        :type  configFile: string
        :param configFile: Configuration file name

        :type  auto_calc_grid_limit: :class:`dict`
        :param auto_calc_grid_limit: the domain where the frequency will be calculated.
                                     The :class:`dict` should contain the keys 
                                     :attr:`xMin`, :attr:`xMax`, :attr:`yMin` 
                                     and :attr:`yMax`. The *x*  variable bounds the 
                                     longitude and the *y* variable bounds
                                     the latitude.
        """

        config = ConfigParser()
        config.read(configFile)

        if config.has_option('TrackGenerator', 'gridLimit'):
            self.tg_domain = config.geteval('TrackGenerator', 'gridLimit')
        else:
            self.tg_domain = auto_calc_grid_limit
           
        self.outputPath = config.get('Output', 'Path')

    def calc(self):
        """
        Calculate the frequency of TC events in a pre-defined domain, based 
        on the input dataset and the full range of years contained in the 
        :attr:`origin_year` file.

        The :attr:`origin_year` file is created in 
        :class:`DataProcess.processData` and restricts the range to a 
        user-selected range of years.
        """

        logger.info("Calculating annual frequency of TC events")
        origin_year = numpy.array(flLoadFile(os.path.join(self.outputPath,
                                                          'process', 'origin_year'),
                                                          '%', ','), dtype='int')
        
        origin_lon_lat = flLoadFile(os.path.join(self.outputPath,
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
        freq = numpy.round(freq*100)/100

        return freq
