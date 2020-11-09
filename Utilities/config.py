"""
:mod:`config` -- reading configuration files
============================================

.. module:: config
    :synopsis: Provides functions for manipulating configuration files
               e.g. reading setting from a configuration file.

.. moduleauthor:: Dale Roberts <dale.roberts@ga.gov.au>

"""

import io
from configparser import RawConfigParser
import os.path

#from ast import literal_eval as eval

def parseBool(txt):
    """
    Parser for boolean options

    :param str txt: String from config file to parse.

    :returns: ``True`` if the string is 'True', ``False`` otherwise.
    :rtype: boolean

    """

    return txt == 'True'

def parseList(txt):
    """
    Parse a comma-separated line into a list.

    :param str txt: String from config file to parse.

    :return: List, based on the input string.
    :rtype: list

    """

    return txt.split(',')


def formatList(lst):
    """
    Convert a list into a comma-joined string.

    :param list lst: Input list to join.

    :return: A string comprised of the list elements joined by commas.
    :rtype: str

    """

    return ','.join([str(l) for l in lst])


FORMATERS = {
    'Input_mslpgrid': formatList,
    'TCRM_columns': formatList
}

PARSERS = {
    'Actions_dataprocess': parseBool,
    'Actions_executehazard': parseBool,
    'Actions_executestat': parseBool,
    'Actions_executetrackgenerator': parseBool,
    'Actions_executewindfield': parseBool,
    'Actions_plotdata': parseBool,
    'Actions_plothazard': parseBool,
    'Actions_createdatabase': parseBool,
    'Actions_downloaddata': parseBool,
    'Actions_executeevaluate': parseBool,
    'DataProcess_inputfile': str,
    'DataProcess_source': str,
    'DataProcess_startseason': int,
    'DataProcess_filterseasons': parseBool,
    'Hazard_calculateci': parseBool,
    'Hazard_minimumrecords': int,
    'Hazard_plotspeedunits': str,
    'Hazard_years': parseList,
    'Hazard_samplesize': int,
    'Hazard_percentilerange': int,
    'Hazard_extremevaluedistribution': str,
    'Hazard_SmoothPlots': parseBool,
    'Input_landmask': str,
    'Input_locationfile': str,
    'Input_mslpgrid': parseList,
    'Input_mslpfile': str,
    'Input_mslpvariablename': str,
    'Logging_logfile': str,
    'Logging_loglevel': str,
    'Logging_progressbar': parseBool,
    'Logging_verbose': parseBool,
    'Logging_datestamp':parseBool,
    'Output_path': str,
    'Process_datfile': str,
    'Process_excludepastprocessed': parseBool,
    'RMW_getrmwdistfrominputdata': parseBool,
    'RMW_mean': float,
    'RMW_sigma': float,
    'Region_gridlimit': eval,
    'Region_gridspace': eval,
    'Region_gridinc': eval,
    'Region_localityid': int,
    'Region_localityname': str,
    'StatInterface_gridinc': eval,
    'StatInterface_gridspace': eval,
    'StatInterface_kde2dtype': str,
    'StatInterface_kdestep': float,
    'StatInterface_kdetype': str,
    'StatInterface_minsamplescell': int,
    'TCRM_columns': parseList,
    'TCRM_fielddelimiter': str,
    'TCRM_numberofheadinglines': int,
    'TCRM_pressureunits': str,
    'TCRM_speedunits': str,
    'Timeseries_Extract': parseBool,
    'Timeseries_LocationFile': str,
    'Timeseries_StationID': str,
    'Timeseries_Windfields': parseBool,
    'TrackGenerator_numsimulations': int,
    'TrackGenerator_seasonseed': int,
    'TrackGenerator_trackseed': int,
    'TrackGenerator_numtimesteps': int,
    'TrackGenerator_timestep': float,
    'WindfieldInterface_beta': float,
    'WindfieldInterface_beta1': float,
    'WindfieldInterface_beta2': float,
    'WindfieldInterface_margin': float,
    'WindfieldInterface_profiletype': str,
    'WindfieldInterface_resolution': float,
    'WindfieldInterface_domain': str,
    'WindfieldInterface_source': str,
    'WindfieldInterface_thetamax': float,
    'WindfieldInterface_trackfile': str,
    'WindfieldInterface_trackpath': str,
    'WindfieldInterface_windfieldtype': str,
    'WindfieldInterface_plotoutput': parseBool}

DEFAULTS = """
[Actions]
DataProcess=True
ExecuteStat=True
ExecuteTrackGenerator=True
ExecuteWindfield=True
ExecuteHazard=True
ExecuteEvaluate=True
PlotData=True
PlotHazard=True
CreateDatabase=True
DownloadData=True

[Region]
gridSpace={'x':1.0,'y':1.0}
gridInc={'x':1.0,'y':0.5}

[DataProcess]
StartSeason=1981
FilterSeasons=True
InputFile=Allstorms.ibtracs_wmo.v03r10.csv
Source=IBTRACS

[StatInterface]
kdeType=Gaussian
kde2DType=Gaussian
kdeStep=0.2
minSamplesCell=100

[TrackGenerator]
NumSimulations=500
NumTimeSteps=360
TimeStep=1.0
Format=csv
SeasonSeed=1
TrackSeed=1

[WindfieldInterface]
profileType=holland
windFieldType=kepert
beta=1.3
beta1=1.3
beta2=1.3
thetaMax=70.0
Margin=2
Resolution=0.05
PlotOutput=False
Domain=bounded

[Hazard]
Years=2,5,10,20,25,50,100,200,250,500,1000
MinimumRecords=50
CalculateCI=True
PercentileRange=90
SampleSize=50
PlotSpeedUnits=mps
ExtremeValueDistribution=GPD
SmoothPlots=True

[RMW]
GetRMWDistFromInputData=False

[Input]
LocationFile=input/stationlist.shp
LandMask=input/landmask.nc
MSLPFile=MSLP/slp.day.ltm.nc
MSLPVariableName=slp
Datasets=IBTRACS,LTMSLP

[Output]
Path=output
Format=txt

[Timeseries]
Extract=False
StationID=WMO
LocationFile=./input/stationlist.shp
Windfields=False

[Logging]
ProgressBar=False
LogFile=main.log
LogLevel=INFO
Verbose=False
Datestamp=False

[Source]
FieldDelimiter=,
NumberOfHeadingLines=0
SpeedUnits=mps
PressureUnits=hPa
LengthUnits=km

[IBTRACS]
URL=ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r10/wmo/csv/Allstorms.ibtracs_wmo.v03r10.csv.gz
path=input
filename=Allstorms.ibtracs_wmo.v03r10.csv
Columns=tcserialno,season,num,skip,skip,skip,date,skip,lat,lon,skip,pressure
FieldDelimiter=,
NumberOfHeadingLines=3
PressureUnits=hPa
LengthUnits=km
DateFormat=%Y-%m-%d %H:%M:%S
SpeedUnits=kph

[LTMSLP]
URL=ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/slp.day.1981-2010.ltm.nc
path=MSLP
filename=slp.day.ltm.nc

"""

class _ConfigParser(RawConfigParser):

    """
    A configuration file parser that extends
    :class:`ConfigParser.RawConfigParser` with a few helper functions
    and default options.
    """
    ignoreSubsequent = True
    def __init__(self, defaults=DEFAULTS):
        RawConfigParser.__init__(self)
        self.readfp(io.StringIO(defaults))
        self.readonce = False
        
    def geteval(self, section, option):
        """
        :param str section: Section name to evaluate.
        :param str option: Option name to evaluate.

        :return: an evaluated setting.

        """
        return self._get(section, eval, option)

    def read(self, filename):
        """
        Read a configuration file, and set the :attr:`readonce` attribute
        to ``True``.

        :param str filename: Path to the configuration file to read.

        """

        if filename is None:
            return
        if self.readonce:
            return
        if not os.path.exists(filename):
            raise ValueError("config file does not exist: {}".format(filename))
        RawConfigParser.read(self, filename)
        self.readonce = True

    def items(self, section):
        """
        Return the parsed option, value pairs for a section of the configuration.

        :param str section: Section name.

        :returns: (Option, value) tuple pairs for the given section.
        :rtype: list

        """

        raw = RawConfigParser.items(self, section)
        parsed = {}
        for name, value in raw:
            try:
                parse = PARSERS['%s_%s' % (section, name)]
                parsed[name] = parse(value)
            except KeyError:
                parsed[name] = value
        return list(parsed.items())

    def set(self, section, option, value=None):
        """
        Set the value of a specific section and option in the configuration.

        :param str section: Section to be updated.
        :param str option: Option to be updated.
        :param value: Value to set.

        """

        try:
            formatter = FORMATERS['%s_%s' % (section, option)]
            newvalue = formatter(value)
        except KeyError:
            newvalue = value
        RawConfigParser.set(self, section, option, newvalue)

singleton = _ConfigParser(defaults=DEFAULTS)
def ConfigParser():
    return singleton
def reset():
    """Re-instantiate ConfigParser (only for use in tests)"""
    global singleton
    singleton = _ConfigParser(defaults=DEFAULTS)

def cnfGetIniValue(configFile, section, option, default=None):
    """
    Helper function to interface with code that uses the
    old config parser.

    :param str configFile: path to the configuration file to read.
    :param str section: Section name to read.
    :param str option: Option name to read.
    :param default: Optional default value to use if the section/option
                    pair is not present in the file.

    :returns: Value recorded in the section/option of the config file, if
              present; the default value otherwise (if defined).

    :raises NoSectionError, NoOptionError: if the section/option is not
              defined, and no default value given.
    """
    config = ConfigParser()
    config.read(configFile)

    if not config.has_section(section):
        return default
    if not config.has_option(section, option):
        return default

    if default is None:
        try:
            res = config.geteval(section, option)
        except (NameError, SyntaxError):
            res = config.get(section, option)
        return res
    if isinstance(default, str):
        return config.get(section, option)
    if isinstance(default, bool):
        return config.getboolean(section, option)
    if isinstance(default, int):
        return config.getint(section, option)
    if isinstance(default, float):
        return config.getfloat(section, option)
    if isinstance(default, dict):
        return config.geteval(section, option)
