import os, io
from ConfigParser import RawConfigParser

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
