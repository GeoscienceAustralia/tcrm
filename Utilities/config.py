import os, io
from ConfigParser import RawConfigParser

def parseGrid(txt):
    return eval('[' + txt + ']')

def parseColumns(txt):
    return txt.split(',')

PARSERS = {
    'Actions_dataprocess': bool,
    'Actions_executehazard': bool,
    'Actions_executestat': bool,
    'Actions_executetrackgenerator': bool,
    'Actions_executewindfield': bool,
    'Actions_plotdata': bool,
    'Actions_plothazard': bool,
    'DataProcess_inputfile': str,
    'DataProcess_source': str,
    'DataProcess_startseason': int,
    'HazardInterface_calculateci': bool,
    'HazardInterface_inputpath': str,
    'HazardInterface_minimumrecords': int,
    'HazardInterface_numsim': int,
    'HazardInterface_plotspeedunits': str,
    'HazardInterface_resolution': float,
    'HazardInterface_years': tuple,
    'HazardInterface_yearspersimulation': int,
    'Input_landmask': str,
    'Input_mslpgrid': parseGrid,
    'Logging_logfile': str,
    'Logging_loglevel': str,
    'Logging_progressbar': bool,
    'Logging_verbose': bool,
    'Output_path': str,
    'Process_datfile': str,
    'Process_excludepastprocessed': bool,
    'RMW_getrmwdistfrominputdata': bool,
    'RMW_mean': float,
    'RMW_sigma': float,
    'Region_gridlimit': eval,
    'Region_localityid': int,
    'Region_localityname': str,
    'StatInterface_gridinc': eval,
    'StatInterface_gridspace': eval,
    'StatInterface_kde2dtype': str,
    'StatInterface_kdestep': float,
    'StatInterface_kdetype': str,
    'StatInterface_minsamplescell': int,
    'TCRM_columns': parseColumns,
    'TCRM_fielddelimiter': str,
    'TCRM_numberofheadinglines': int,
    'TCRM_pressureunits': str,
    'TCRM_speedunits': str,
    'TrackGenerator_gridinc': eval,
    'TrackGenerator_gridspace': eval,
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
    'WindfieldInterface_numberoffiles': int,
    'WindfieldInterface_profiletype': str,
    'WindfieldInterface_resolution': float,
    'WindfieldInterface_source': str,
    'WindfieldInterface_thetamax': float,
    'WindfieldInterface_trackfile': str,
    'WindfieldInterface_trackpath': str,
    'WindfieldInterface_windfieldtype': str}

DEFAULTS = """
[Actions]
DataProcess=False
ExecuteStat=False
ExecuteTrackGenerator=False
ExecuteWindfield=False
ExecuteHazard=False
PlotData=False
PlotHazard=False

[DataProcess]
StartSeason=1981

[TrackGenerator]
NumSimulations=50
YearsPerSimulation=10
NumTimeSteps=360
TimeStep=1.0
Format=csv
gridSpace={'x':1.0,'y':1.0}
gridInc={'x':1.0,'y':0.5}

[WindfieldInterface]
profileType=holland
windFieldType=kepert
beta=1.3
beta1=1.3
beta2=1.3
thetaMax=70.0
Margin=2
Resolution=0.01
PlotOutput=False

[RMW]
GetRMWDistFromInputData=False

[Input]
LandMask=%(cwd)s/input/landmask.nc
MSLPFile=%(cwd)s/MSLP/mslp_daily_ltm.nc

[Output]
Path=%(cwd)s/output
Format=txt

[Logging]
ProgressBar=False
LogFile=main.log
LogLevel=INFO
Verbose=False

[StatInterface]
kdeType=Gaussian

[Source]
FieldDelimiter=,
NumberOfHeadingLines=0
SpeedUnits=mps
PressureUnits=hPa
LengthUnits=km
""" % {'cwd': os.getcwd()}

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
        if filename is None: return
        if self.read_once: return
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
