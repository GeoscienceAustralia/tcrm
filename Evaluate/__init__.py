from pressureDistribution import PressureDistribution
from trackDensity import TrackDensity
from longitudeCrossing import LongitudeCrossing


def run(configFile):
    PD = PressureDistribution(configFile)
    TD = TrackDensity(configFile)
    LC = LongitudeCrossing(configFile)
    
    PD.run()
    TD.run()
    LC.run()
