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

 Title: config.py - configuration settings.
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 26/06/2008 2:09PM
 Description: Provides functions for manipulating configuration files
 e.g. reading settings from a configuration file, reading a list of
 values from a configuration file.
 Constraints:
 SeeAlso: files.py
 Version :$Rev: 512 $

 $Id: config.py 512 2011-10-31 07:20:38Z nsummons $
"""

import os, sys, pdb, logging
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    execfile(filename)
import time
import ConfigParser
from files import flConfigFile
import inspect

__version__ = '$Id: config.py 512 2011-10-31 07:20:38Z nsummons $'
config_dict = {}
logger = logging.getLogger()

def cnfCacheIniFile(configFile=None):
    """
    A wrapper to the ConfigParser module which caches a dictionary of
    dictionaries. Each dictionary in the parent dictionary is named
    after sections in the configuration file. Keys in the
    sub-directories are named by the options in each section.
    Input: configuration file name (optional, defaults to output from
           flConfigFile())
    Output: configuration dictionary
    Example: config_dict = cnfCacheIniFile(configFile)
    """
    if configFile:
        try:
            FH = open(configFile)
        except IOError:
            logger.info("Cannot open %s" % configFile)
            return config_dict
    elif len(sys.argv) > 1:
        try:
            FH = open(sys.argv[1])
        except IOError:
            logger.info("No configuration file given at command line")
    else:
        try:
            FH = open(flConfigFile(level=len(inspect.stack())))
            logger.info("Opening default config file %s" % \
                         flConfigFile(level=len(inspect.stack())))
        except IOError:
            logger.info("Cannot open default config file %s" % \
                         flConfigFile(level=len(inspect.stack())))
            return config_dict

    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(FH)
    for sec in cp.sections():
        name = sec
        if name not in config_dict:
            config_dict[name] = {}
        for opt in cp.options(sec):
            try:
                config_dict[name][opt] = cp.getint(sec, opt)
            except ValueError:
                try:
                    config_dict[name][opt] = cp.getfloat(sec, opt)
                except ValueError:
                    try:
                        config_dict[name][opt] = cp.getboolean(sec, opt)
                    except ValueError:
                        config_dict[name][opt] = cp.get(sec, opt)
    FH.close()
    return config_dict

def cnfGetCachedIniValue(section, option):
    """
    Get the cached value of an option
    Assumes the global config_dict dictionary exists & is populated
    Input: section, option
    Output: value
    """
    if section in config_dict:
        if option in config_dict[section]:
            return config_dict[section][option]
        else:
            return None
    else:
        return None

def cnfSetCachedIniValue(section, option, value):
    """
    Set the cache value of an option
    Assumes the global config_dict dictionary exists
    Input: section, option, value
    Output: None
    """
    if section in config_dict:
        if option in config_dict[section]:
            config_dict[section][option] = value
    else:
        config_dict[section] = {option:value}
    return

def cnfGetIniValue(configFile, section, option, default=None):
    """
    Get the value of an attribute from a specified INI file
    Input: ini file name, section, option, default
    Output: value
    Example: value = cnfGetIniValue(fileName, section, option, default)

    Notes:
    The search is case-sensitive.
    If you must get the value directly from the file, use
    cnfGetIniFileValue.  If the file is already cached then the cached
    copy is used. Otherwise the ini filename is first cached, and then
    the attribute's value extracted. If that fails, return the default
    (if defined) or an undefined value.
    """
    value = cnfGetCachedIniValue(section, option)
    if value is not None:
        return value
    else:
        config_dict = cnfCacheIniFile(configFile)
    value = cnfGetCachedIniValue(section, option)
    if value is not None:
        return value
    elif default is not None:
        return default
    else:
        logger.info("No value set for section: %s, option: %s (and no default value was given)" % \
                     (section,option))
        return None


def cnfGetIniList(configFile, section, first=1, last=None):
    """
    Get a list of values for integer options in a configuration file.
    Input: ini file name, section, first value (optional, default 1)
    Output: list, or number of values in scalar context
    Example: out = cnfGetIniList(fileName, section, first=1, last=None)
    First value defaults to 1. Last is the last value which will be
    tried.  Allows gaps in the sequence, but there should not be
    duplicate values (we can't define which would be retrieved!).
    """
    if configFile is None:
        FH = open(flConfigFile())
    else:
        FH = open(configFile)
    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(FH)
    values = []
    options = cp.options(section)
    options = [int(o) for o in options]
    options.sort()
    for opt in options:
        if last:
            if (int(opt) >= first) and (int(opt) <= last):
                try:
                    values.append(cp.getint(section, str(opt)))
                except ValueError:
                    try:
                        values.append(cp.getfloat(section, str(opt)))
                    except ValueError:
                        try:
                            values.append(cp.getboolean(section, str(opt)))
                        except ValueError:
                            values.append(cp.get(section, str(opt)))
        else:
            if (int(opt) >= first):
                try:
                    values.append(cp.getint(section, str(opt)))
                except ValueError:
                    try:
                        values.append(cp.getfloat(section, str(opt)))
                    except ValueError:
                        try:
                            values.append(cp.getboolean(section, str(opt)))
                        except ValueError:
                            values.append(cp.get(section, str(opt)))
    FH.close()
    return values

def cnfGetIniFileValue(configFile, section, option, default=None):
    """
    Get a value directly from the configuration file, rather than using
    a cached value.  This shouldn't be used in anger, as it may produce
    adverse effects if you make changes to the config file while a
    program is running. Optionally takes a default value to return if
    the option is not in the configuration file.
    Input: configuration file name, section name, option name, default
    value.
    Output: value directly from the configuration file.
    Example: value = cnfGetIniFileValue(configFile, section, option,
                                        default)
    """
    if configFile is None:
        FH = open(flConfigFile())
    else:
        FH = open(configFile)
    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(FH)
    sections = cp.sections()
    if section in sections:
        options = cp.options(section)
        if option in options:
            try:
                value = cp.getint(section, option)
            except ValueError:
                try:
                    value = cp.getfloat(section, option)
                except ValueError:
                    try:
                        value = cp.getboolean(section, option)
                    except ValueError:
                        value = cp.get(section, option)
        else:
            value = default
    else:
        value = default
    FH.close()
    return value
