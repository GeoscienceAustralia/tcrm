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

 Title: config.py - configuration settings.
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: 26/06/2008 2:09PM
 Description: Provides functions for manipulating configuration files
 e.g. reading settings from a configuration file, reading a list of
 values from a configuration file.
 Constraints:
 SeeAlso: files.py
 Version :$Rev: 642 $

"""

import sys
import logging
import ConfigParser
from Utilities.files import flConfigFile
import inspect
import re

CONFIG_DICT = {}
LOG = logging.getLogger()


class ConfigError(Exception):

    """
    Exception class to handle configuration errors
    """

    def __init__(self, value):
        Exception.__init__()
        self.value = value

    def __str__(self):
        return repr(self.value)


def cnfCacheIniFile(configFile=None):
    """
    A wrapper to the ConfigParser module which caches a dictionary of
    dictionaries. Each dictionary in the parent dictionary is named
    after sections in the configuration file. Keys in the
    sub-directories are named by the options in each section.

    Input: configuration file name (optional, defaults to output from
           flConfigFile())
    Output: configuration dictionary
    Example: CONFIG_DICT = cnfCacheIniFile(configFile)
    """
    if configFile:
        try:
            fh = open(configFile)
        except IOError:
            LOG.info("Cannot open %s", configFile)
            return CONFIG_DICT
    elif len(sys.argv) > 1:
        try:
            fh = open(sys.argv[1])
        except IOError:
            LOG.info("No configuration file given at command line")
    else:
        try:
            fh = open(flConfigFile(level=len(inspect.stack())))
            LOG.info("Opening default config file %s",
                     flConfigFile(level=len(inspect.stack())))
        except IOError:
            LOG.info("Cannot open default config file %s",
                     flConfigFile(level=len(inspect.stack())))
            return CONFIG_DICT

    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(fh)
    for sec in cp.sections():
        name = sec
        if name not in CONFIG_DICT:
            CONFIG_DICT[name] = {}
        for opt in cp.options(sec):
            try:
                CONFIG_DICT[name][opt] = cp.getint(sec, opt)
            except ValueError:
                try:
                    CONFIG_DICT[name][opt] = cp.getfloat(sec, opt)
                except ValueError:
                    try:
                        CONFIG_DICT[name][opt] = cp.getboolean(sec, opt)
                    except ValueError:
                        CONFIG_DICT[name][opt] = cp.get(sec, opt)
    fh.close()
    return CONFIG_DICT


def cnfGetCachedIniValue(section, option):
    """
    Get the cached value of an option
    Assumes the global CONFIG_DICT dictionary exists & is populated

    Input: section, option
    Output: value

    Example: value = cnfGetCachedIniValue( section, option )
    """
    try:
        if section in CONFIG_DICT:
            if option in CONFIG_DICT[section]:
                return CONFIG_DICT[section][option]
            else:
                return None
        else:
            return None
    except ConfigParser.NoSectionError:
        LOG.exception("No %s section in the configuration file", section)
        raise


def cnfSetCachedIniValue(section, option, value):
    """
    Set the cached value of an option in the CONFIG_DICT object
    Assumes the global CONFIG_DICT dictionary exists.

    Input: section, option, value
    Output: None
    Example: cnfSetCachedIniValue( section, option, value )
    """

    if section in CONFIG_DICT:
        if option in CONFIG_DICT[section]:
            CONFIG_DICT[section][option] = value
    else:
        CONFIG_DICT[section] = {option: value}
    return


def __cnfGetIniValue(configFile, section, option, default=None):
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
        CONFIG_DICT = cnfCacheIniFile(configFile)
    value = cnfGetCachedIniValue(section, option)
    if value is not None:
        return value
    elif default is not None:
        return default
    else:
        raise ConfigError("No value set for section: %s, option: %s (and no \
        default value was given)", section, repr(option))
        # return None


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
        fh = open(flConfigFile())
    else:
        fh = open(configFile)
    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(fh)
    values = []
    try:
        options = cp.options(section)
    except ConfigParser.NoSectionError:
        LOG.exception("No section named %s in configuration file %s",
                      section, configFile)
        fh.close()
        raise

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
            if int(opt) >= first:
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
    fh.close()
    return values


def cnfGetIniFileValue(configFile, section, option, default=None):
    """
    Get a value directly from the configuration file, rather than
    using a cached value.  This shouldn't be used in anger, as it may
    produce adverse effects if you make changes to the config file
    while a program is running. Optionally takes a default value to
    return if the option is not in the configuration file.

    Input: configuration file name, section name, option name, default
    value.
    Output: value directly from the configuration file.
    Example: value = cnfGetIniFileValue(configFile, section, option,
                                        default)
    """
    if configFile is None:
        fh = open(flConfigFile())
    else:
        fh = open(configFile)
    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(fh)
    try:
        sections = cp.sections()
    except ConfigParser.NoSectionError:
        LOG.exception("No section named %s in configuration file %s",
                      section, configFile)
        fh.close()
        raise

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
    fh.close()
    return value


def cnfGetUnorderedList(configFile, section):
    """
    Get a list of unordered values in the given section.

    Input: ini file name, section
    Output: list, or number of values in scalar context
    Example: out = cnfGetUnorderdedList( filename, section )

    This function requires a different approach to reading
    the configuration file, as ConfigParser in Python versions
    prior to 2.7 cannot handle options with no value. The workaround
    reads the configuration file and constructs a dict similar
    to CONFIG_DICT used elsewhere in this module.
    Lists of unordered values are included in the dict with the
    key and subkey equal to the section name.
    i.e. cfgDict[ section ][ section ] = [list]


    """
    sect = None
    sect_list = []
    cfgDict = {}
    if configFile is None:
        configFile = flConfigFile()
    try:
        fh = open(configFile)
    except IOError:
        LOG.warn("Cannot open %s", configFile)
    else:
        for line in fh:
            line = line.lstrip().rstrip('\n')
            cm = re.match('^;', line)
            if cm:
                # Ignore lines that begin with a comment character
                continue
            sm = re.match(r'^\[(\w*)\]', line)
            am = re.match(r'^([^=]+)=(.+)', line)
            if sm:
                new_sect = sm.group(1)
                if sect:
                    key = sect
                    subkey = key
                    cfgDict[key][subkey] = sect_list
                sect = new_sect
                sect_list = []

            elif am:
                # Attribute/value pair
                att = am.group(1)
                val = am.group(2)
                att = att.rstrip()
                val = val.rstrip().lstrip()
                if cfgDict.has_key(sect):
                    cfgDict[sect][att] = val
                else:
                    cfgDict[sect] = {}
                    cfgDict[sect][att] = val
            elif len(line):
                sect_list.append(line)

        fh.close()

    return cfgDict[section][section]


def _cnfCacheIniFile(filename):
    """
    This version of cnfCacheIniFile doesn't rely on configParser.
    The disadvantage of this is that all values are treated as strings
    and so type checking/conversion will need to be performed in the
    calling program.

    Input: configuration file name (optional, defaults to output from
           flConfigFile())
    Output: configuration dictionary
    Example: CONFIG_DICT = _cnfCacheIniFile(configFile)
    """

    cfgDict = {}
    try:
        fh = open(filename)
    except IOError:
        LOG.warn("Cannot open %s", filename)
    else:
        for line in fh:
            line = line.rstrip('\n')
            line = line.lstrip()
            cm = re.match('^[;#]', line)
            if cm:
                # Ignore comment lines
                continue
            sm = re.match(r'^\[(\w*)\]', line)
            am = re.match(r'^([^=]+)=(.+)', line)
            if sm:
                new_sect = sm.group(1)
                if sect:  # pylint: disable=E0601
                    key = sect
                    subkey = key
                    cfgDict[key][subkey] = sect_list  # pylint: disable=E0601
                sect = new_sect
                sect_list = []

            elif am:
                # Attribute/value pair
                att = am.group(1)
                val = am.group(2)
                att = att.rstrip()
                val = val.rstrip().lstrip()
                if cfgDict.has_key(sect):
                    cfgDict[sect][att] = val
                else:
                    cfgDict[sect] = {}
                    cfgDict[sect][att] = val
            elif len(line):
                sect_list.append(line)

        fh.close()

    return cfgDict


def cnfRefreshCachedIniFile(configFile):
    """
    Reload a configuration file into the configuration
    dictionary. This will replace all values with the values from
    the file, so use with caution as the input configuration file
    may have been changed.

    Input: ini file name
    Output: None. The CONFIG_DICT object is updated
    Example: cnfRefreshCachedIniFile( configFile )
    """
    CONFIG_DICT = cnfCacheIniFile(configFile)
