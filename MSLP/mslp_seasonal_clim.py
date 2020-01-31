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


Title: mslp_seasonal_clim.py

Author: Nicholas Summons, nicholas.summons@ga.gov.au

Last modified: 4 June 2010

Description: Utility for creating Mean Sea Level Pressure (MSLP) seasonal climatology maps.
             Uses NCEP-DOE Reanalysis 2 data averaged over date range: 1980-2007.
             This script can either be run stand alone to create a NetCDF output file or
             the class MSLPGrid can be invoked to return the MSLP seasonal average grid.

Acknowledgements:
             NCEP-DOE Reanalysis 2 data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA,
             from their Web site at http://www.esrl.noaa.gov/psd/

Input data: mslp_seasonal_clim.nc (contains monthly means averaged over 28 year period)
"""

import os
import numpy as np
import Utilities.nctools as nctools
from Utilities import pathLocator


class MSLPGrid:

    def __init__(self, selected_months, filename=''):
        if not os.path.isfile(filename):
            tcrm_dir = pathLocator.getRootDirectory()
            filename = os.path.join(tcrm_dir, 'MSLP', 'mslp_monthly_clim.nc')
            if not os.path.isfile(filename):
                error_msg = "MSLP data file not found"
                raise IOError(error_msg)
        selected_months = set(selected_months)
        ncobj = nctools.ncLoadFile(filename)
        mslp_all = nctools.ncGetData(ncobj, 'mslp')
        self.lon = nctools.ncGetDims(ncobj, 'lon')
        self.lat = nctools.ncGetDims(ncobj, 'lat')
        dim0,dim1,dim2 = np.shape(mslp_all)

        # Average over selected months
        mslp_sum = np.zeros([dim1, dim2], dtype='float32')
        for month in selected_months:
            mslp_sum = mslp_sum + mslp_all[month-1,:,:]
        self.mslp_av = np.flipud(mslp_sum / len(selected_months))

    def sampleGrid(self, lon, lat):
        """sampleGrid(self, lon, lat):
        Grab nearest value to given location.
        No interpolation performed!
        """
        indi = self.lon.searchsorted(lon)-1
        indj = self.lat.searchsorted(lat)-1

        return self.mslp_av[indj, indi]

    def returnGrid(self):
        return self.lon, self.lat, self.mslp_av


#def main(configFile):
#    selected_months_str = str(cnfGetIniValue(configFile, 'DataProcess', 'selected_months', arange(13)))
#    selected_months = set(selected_months_str.strip('[]{}() ').replace(',', ' ').split(' '))
#    selected_months.discard('')
#    if selected_months.issubset([str(k) for k in range(1,13)]):
#        selected_months = [int(k) for k in selected_months]
#        months_str = ', '.join([calendar.month_abbr[i] for i in sort(list(selected_months))])
#        print "Creating Mean Sea Level Pressure (MSLP) seasonal climatology:"
#        print "Months specified for seasonal average: " + months_str
#        print "Using NCEP Reanalysis-2 data from 1980-2007"
#
#        msp = MSLPGrid(selected_months)
#        lon, lat, mslp_av = msp.returnGrid()
#
#        #Create output file
#        output_filename = "mslp_clim_" + ''.join([calendar.month_abbr[i][0] for i in sort(list(selected_months))]) + '.nc'
#        data_title = 'MSLP (NCEP Reanalysis-2) seasonal climatology.  Averaging period: ' \
#                     +  months_str + ' ' + '1980-2007.'
#        dimensions = {0:{'name':'lat','values':lat,'dtype':'f','atts':{'long_name':'Latitude',
#                                                                       'units':'degrees_north'} },
#                      1:{'name':'lon','values':lon,'dtype':'f','atts':{'long_name':'Longitude',
#                                                                       'units':'degrees_east'} } }
#
#        variables = {0:{'name':'mslp','dims':('lat','lon'),
#                        'values':array(mslp_av),'dtype':'f',
#                        'atts':{'long_name':'Mean sea level pressure',
#                                'units':'hPa'} } }
#        nctools.ncSaveGrid( output_filename, dimensions, variables,
#                            nodata=-9999,datatitle=data_title )
#
#        print "Created output file: " + output_filename
#
#
#if __name__ == "__main__":
#    try:
#        configFile = sys.argv[1]
#    except IndexError:
#        # Try loading config file with same name as python script
#        configFile = __file__.rstrip('.py') + '.ini'
#        # If no filename is specified and default filename doesn't exist => raise error
#        if not os.path.exists(configFile):
#            error_msg = "No configuration file specified"
#            raise IOError, error_msg
#    # If config file doesn't exist => raise error
#    if not os.path.exists(configFile):
#        error_msg = "Configuration file '" + configFile +"' not found"
#        raise IOError, error_msg
#
#    main(configFile)
