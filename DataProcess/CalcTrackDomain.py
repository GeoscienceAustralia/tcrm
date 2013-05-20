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
CreationDate: 2010-03-23
Description: Determine track generator domain, ensuring it encompasses
             all tracks entering the windfield generator domain.
"""

import os, sys, pdb, logging

import math

from Utilities.config import cnfGetIniValue
from Utilities.files import flLoadFile


class CalcTrackDomain:

    def __init__(self, configFile):

        self.outputPath = cnfGetIniValue(configFile, 'Output', 'Path')
        self.wf_domain = eval(cnfGetIniValue(configFile, 'Region', 'gridLimit'))

    def _calc(self, index, lons, lats):

        tg_domain = self.wf_domain.copy()
        track_limits = {'xMin':9999,'xMax':-9999,'yMin':9999,'yMax':-9999}
        for [idx, lon, lat] in zip(index, lons, lats):

            if idx == 1:
                # Reset cyclone lon/lon limits
                track_limits = {'xMin':9999, 'xMax':-9999,
                                'yMin':9999, 'yMax':-9999}

            track_limits['xMin'] = min(track_limits['xMin'], lon)
            track_limits['xMax'] = max(track_limits['xMax'], lon)
            track_limits['yMin'] = min(track_limits['yMin'], lat)
            track_limits['yMax'] = max(track_limits['yMax'], lat)

            if (self.wf_domain['xMin'] <= lon <= self.wf_domain['xMax']) & \
               (self.wf_domain['yMin'] <= lat <= self.wf_domain['yMax']):

                tg_domain['xMin'] = min(tg_domain['xMin'],
                                        track_limits['xMin'])
                tg_domain['xMax'] = max(tg_domain['xMax'],
                                        track_limits['xMax'])
                tg_domain['yMin'] = min(tg_domain['yMin'],
                                        track_limits['yMin'])
                tg_domain['yMax'] = max(tg_domain['yMax'],
                                        track_limits['yMax'])

            # Extend domain to closest integer lat/lon value
            tg_domain['xMin'] = math.floor(tg_domain['xMin'])
            tg_domain['xMax'] = math.ceil(tg_domain['xMax'])
            tg_domain['yMin'] = math.floor(tg_domain['yMin'])
            tg_domain['yMax'] = math.ceil(tg_domain['yMax'])

        return tg_domain

    def calcDomainFromTracks(self, index, lons, lats):
        tg_domain = self._calc(index, lons, lats)
        return tg_domain

    def calcDomainFromFile(self):

        #tg_domain = self.wf_domain.copy()
        #track_limits = {'xMin':9999,'xMax':-9999,'yMin':9999,'yMax':-9999}

        # Load tracks from file
        cyclone_tracks = flLoadFile(os.path.join(self.outputPath,
                                                 'process', 'cyclone_tracks'),
                                                 '%', ',')
        tg_domain = self._calc(cyclone_tracks[0],
                               cyclone_tracks[1],
                               cyclone_tracks[2])

        return tg_domain
        """
        for idx, lon, lat in cyclone_tracks:

            if idx == 1:
                # Reset cyclone lon/lon limits
                track_limits = {'xMin':9999, 'xMax':-9999,
                                'yMin':9999, 'yMax':-9999}

            track_limits['xMin'] = min(track_limits['xMin'], lon)
            track_limits['xMax'] = max(track_limits['xMax'], lon)
            track_limits['yMin'] = min(track_limits['yMin'], lat)
            track_limits['yMax'] = max(track_limits['yMax'], lat)

            if (self.wf_domain['xMin'] <= track_lon <= self.wf_domain['xMax']) & \
               (self.wf_domain['yMin'] <= track_lat <= self.wf_domain['yMax']):

                tg_domain['xMin'] = min(tg_domain['xMin'],
                                        track_limits['xMin'])
                tg_domain['xMax'] = max(tg_domain['xMax'],
                                        track_limits['xMax'])
                tg_domain['yMin'] = min(tg_domain['yMin'],
                                        track_limits['yMin'])
                tg_domain['yMax'] = max(tg_domain['yMax'],
                                        track_limits['yMax'])

            # Extend domain to closest integer lat/lon value
            tg_domain['xMin'] = math.floor(tg_domain['xMin'])
            tg_domain['xMax'] = math.ceil(tg_domain['xMax'])
            tg_domain['yMin'] = math.floor(tg_domain['yMin'])
            tg_domain['yMax'] = math.ceil(tg_domain['yMax'])

        return tg_domain
        """