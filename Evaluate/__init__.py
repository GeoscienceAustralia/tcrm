from pressureDistribution import PressureDistribution
from trackDensity import TrackDensity


def run(configFile):
    PD = PressureDistribution(configFile)
    TD = TrackDensity(configFile)
    
    PD.run()
    TD.run()