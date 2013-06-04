import os, io
from ConfigParser import RawConfigParser

DEFAULTS = """
[TrackGenerator]
NumSimulations=50
YearsPerSimulation=10
NumTimeSteps=360
TimeStep=1.0
Format=csv

[Output]
Path=%(cwd)s/output
""" % {'cwd': os.getcwd()}

class ConfigParser(RawConfigParser):
    """
    A configuration file parser that extends
    :class:`ConfigParser.RawConfigParser` with a few helper functions
    and default options.
    """

    def __init__(self, *args, **kwargs):
        RawConfigParser.__init__(self, *args, **kwargs)
        self.readfp(io.BytesIO(DEFAULTS))

    def geteval(self, section, option):
        """
        :return: an evaluated setting.
        """
        return self._get(section, eval, option)

