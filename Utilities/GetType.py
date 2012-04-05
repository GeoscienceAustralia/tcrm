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

Title: getType.py
Author: Nicholas Summons, nicholas.summons@ga.gov.au
CreationDate: 2010-04-30
Description:  The class returns default data type associations.  This is used by
              columns.py when a data type is not defined in the configuration file.
"""
import os, sys, pdb, logging


class GetType:


    def __init__(self):

        self.logger = logging.getLogger()

        self.type_dict = {'age':'float',
                         'bearing':'float',
                         'date':'string',
                         'day':'float',
                         'direction':'float',
                         'hour':'float',
                         'index':'int',
                         'intensity':'string',
                         'lat':'float',
                         'lon':'float',
                         'maxspeed':'float',
                         'minute':'float',
                         'month':'float',
                         'mpi':'float',
                         'name':'string',
                         'num':'int',
                         'penv':'float',
                         'pressure':'float',
                         'rmax':'float',
                         'season':'float',
                         'second':'float',
                         'shear':'float',
                         'skip':'string',
                         'speed':'float',
                         'surfcode':'string',
                         'tcserialno':'string',
                         'unit':'string',
                         'vmax':'float',
                         'year':'float'}


    def getType(self, key):

        try:
            data_type = self.type_dict[key]
        except KeyError:
            self.logger.warn("No data type association found for '" + str(key) + "'.  Defaulting to data type: 'float'.")
            data_type = 'float'

        return data_type


    def getKeys(self):

        key_list = self.type_dict.keys()

        return key_list
