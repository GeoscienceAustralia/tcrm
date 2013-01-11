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

 Title: columns.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 2009-01-28 10:30 AM
 Description: Helper functions for reading CSV formatted data.

 Version :$Rev: 642 $

 $Id: columns.py 642 2012-02-21 07:54:04Z nsummons $
"""
import os, sys, pdb, logging

from config import cnfGetIniValue
from GetType import GetType
import numpy
import csv
gT = GetType()

__version__ = '$Id: columns.py 642 2012-02-21 07:54:04Z nsummons $'
logger = logging.getLogger()

def colReadCSV(configFile, dataFile, source, nullValue=sys.maxint, configSettings=None):
    """colReadCSV(dataFile, source, nullValue=sys.maxint):
    Loads a csv file containing 'column' data into a dictionary of
    (numpy) arrays with keys labelled by 'fields'. Values of the
    dictionary are arrays containing the data in the columns in the
    file.  This is necessary as the csv.reader() function
    (and hence DictReader) returns all values as strings. No automatic
    data type conversion is performed.
    Input: dataFile - path to the csv file to read
           source - name of a data source format describing the data.
                    The source should have a corresponding section in
                    gConfigFile.
    Output: data - a dictionary of arrays, with keys given by the
                   column names specified in the configuration file.
    Example: data = colReadCSV(dataFile, source)
    """

    if configFile is not None:
        delimiter = cnfGetIniValue(configFile,source, 'FieldDelimiter', ',')
        cols = cnfGetIniValue(configFile, source, 'Columns').split(delimiter)
        fields = cnfGetIniValue(configFile, source, 'Fields', '')
        if fields == '':
            fields = []
        else:
            fields = cnfGetIniValue(configFile, source, 'Fields', '').split(delimiter)
        headingLine = cnfGetIniValue(configFile, source, 'HeadingLine', False)
        if headingLine:
            numHeadingLines = 1
        else:
            numHeadingLines = 0
        # Overwrite deprecated 'HeadingLine' setting if 'NumberOfHeadingLines' is specified
        numHeadingLines = int(cnfGetIniValue(configFile, source, 'NumberOfHeadingLines', numHeadingLines))
        logger.debug("Opening %s, using format %s"%(dataFile, source))
    else:
        # If no configFile provided => load settings from configSettings dictionary
        # [Note: this option is required by the config editor GUI]
        delimiter = configSettings['delimiter']
        cols = configSettings['Columns'].split(delimiter)
        numHeadingLines = int(configSettings['NumberOfHeadingLines'])
        fields = []


    # Add default field names from GetType.py
    fields = fields + gT.getKeys()

    # Check the data file is readable:
    try:
        fh = open(dataFile, 'r')
    except IOError:
        logger.critical("Cannot open %s"%dataFile)
        raise IOError, ("Cannot open %s"%dataFile)
    else:
        numRecords = len(fh.readlines())
        fh.close()

    # A dictionary of column names to field names:
    mappings = {}
    for colname in cols:
        if colname not in fields:
            mappings[colname] = cnfGetIniValue(configFile, 'Mappings',
                                                      colname)
        else:
            mappings[colname] = colname

    # Create a temporary dictionary of *lists*:
    temp = {}
    for colname in cols:
        temp[colname] = []

    data = csv.DictReader(file(dataFile), cols, restval=nullValue,
                          delimiter=str(delimiter))

    dtypeDict = {}

    for key in cols:
        if configFile is not None:
            dtype = cnfGetIniValue(configFile, 'Types', key, 'None')
        else:
            dtype = 'None'
        if dtype != 'None':
            dtypeDict[key] = dtype
        else:
            # If data type not specified, get default value from getType.py
            dtypeDict[key] = gT.getType(key)

    ii = 0
    for record in data:
        ii += 1
        if ii <= numHeadingLines:
            if ii == 1:
                logger.debug("Skipping text header lines")
            continue
        for key in record.keys():
            if key in temp.keys():
                if dtypeDict[key] == 'int':
                    try:
                        temp[key].append(int(record[key]))
                    except ValueError:
                        temp[key].append(nullValue)
                elif dtypeDict[key] == 'float':
                    try:
                        temp[key].append(float(record[key]))
                    except ValueError:
                        temp[key].append(nullValue)
                else:
                    temp[key].append(record[key])

    output = {}
    for key in temp.keys():
        output[mappings[key]] = numpy.array(temp[key])

    return output
