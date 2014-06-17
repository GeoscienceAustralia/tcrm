from pressureDistribution import PressureDistribution
from trackDensity import TrackDensity
from longitudeCrossing import LongitudeCrossing
from landfallRates import LandfallRates


def run(configFile):
    PD = PressureDistribution(configFile)
    TD = TrackDensity(configFile)
    LC = LongitudeCrossing(configFile)
    LF = LandfallRates(configFile)
    
    PD.run()
    TD.run()
    LC.run()
    LF.run()
