import io
from ConfigParser import RawConfigParser


def parseBool(txt):
    return txt == 'True'


def parseList(txt):
    return txt.split(',')


def formatList(lst):
    return ','.join(map(str, lst))


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
    'Input_landmask': str,
    'Input_mslpgrid': parseList,
    'Logging_logfile': str,
    'Logging_loglevel': str,
    'Logging_progressbar': parseBool,
    'Logging_verbose': parseBool,
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
    'TrackGenerator_numsimulations': int,
    'TrackGenerator_seasonseed': int,
    'TrackGenerator_trackseed': int,
    'TrackGenerator_yearspersimulation': int,
    'TrackGenerator_numtimesteps': int,
    'TrackGenerator_timestep': float,
    'WindfieldInterface_beta': float,
    'WindfieldInterface_beta1': float,
    'WindfieldInterface_beta2': float,
    'WindfieldInterface_margin': float,
    'WindfieldInterface_profiletype': str,
    'WindfieldInterface_resolution': float,
    'WindfieldInterface_source': str,
    'WindfieldInterface_thetamax': float,
    'WindfieldInterface_trackfile': str,
    'WindfieldInterface_trackpath': str,
    'WindfieldInterface_windfieldtype': str}

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
DownloadData=True

[Region]
gridSpace={'x':1.0,'y':1.0}
gridInc={'x':1.0,'y':0.5}

[DataProcess]
StartSeason=1981
FilterSeasons=True
InputFile=Allstorms.ibtracs_wmo.v03r05.csv
Source=IBTRACS

[StatInterface]
kdeType=Gaussian
kde2DType=Gaussian
kdeStep=0.2
minSamplesCell=100

[TrackGenerator]
NumSimulations=500
YearsPerSimulation=1
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

[Hazard]
Years=2,5,10,20,25,50,100,200,250,500,1000
MinimumRecords=50
CalculateCI=True
PercentileRange=90
SampleSize=100
PlotSpeedUnits=mps

[RMW]
GetRMWDistFromInputData=False

[Input]
LandMask=input/landmask.nc
MSLPFile=MSLP/slp.day.ltm.nc
Datasets=IBTRACS,LTMSLP

[Output]
Path=output
Format=txt

[Logging]
ProgressBar=False
LogFile=main.log
LogLevel=INFO
Verbose=False

[Source]
FieldDelimiter=,
NumberOfHeadingLines=0
SpeedUnits=mps
PressureUnits=hPa
LengthUnits=km

[IBTRACS]
URL=ftp://eclipse.ncdc.noaa.gov/pub/ibtracs/v03r05/wmo/csv/Allstorms.ibtracs_wmo.v03r05.csv.gz
path=input
filename=Allstorms.ibtracs_wmo.v03r05.csv
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


def singleton(cls):
    instances = {}

    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance


@singleton
class ConfigParser(RawConfigParser):

    """
    A configuration file parser that extends
    :class:`ConfigParser.RawConfigParser` with a few helper functions
    and default options.
    """

    def __init__(self):
        RawConfigParser.__init__(self)
        self.readfp(io.BytesIO(DEFAULTS))
        self.read_once = False

    def geteval(self, section, option):
        """
        :return: an evaluated setting.
        """
        return self._get(section, eval, option)

    def read(self, filename):
        if filename is None:
            return
        if self.read_once:
            return
        RawConfigParser.read(self, filename)
        self.read_once = True

    def items(self, section):
        raw = RawConfigParser.items(self, section)
        parsed = {}
        for name, value in raw:
            try:
                parse = PARSERS['%s_%s' % (section, name)]
                parsed[name] = parse(value)
            except KeyError:
                parsed[name] = value
        return parsed.items()

    def set(self, section, option, value):
        try:
            formatter = FORMATERS['%s_%s' % (section, option)]
            newvalue = formatter(value)
        except KeyError:
            newvalue = value
        RawConfigParser.set(self, section, option, newvalue)


def cnfGetIniValue(configFile, section, option, default=None):
    """
    Helper function to interface with code that uses the
    old config parser.
    """
    config = ConfigParser()
    config.read(configFile)
    if not config.has_option(section, option):
        return default
    if default is None:
        return config.get(section, option)
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
