"""
:mod:`Evaluate` -- run evaluation of track model
================================================
.. module:: Evaluate
    :synopsis: Provide some evaluation methods for testing
               the track generation model.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

from pressureDistribution import PressureDistribution
from trackDensity import TrackDensity
from longitudeCrossing import LongitudeCrossing
from landfallRates import LandfallRates
from genesisDensity import GenesisDensity


def run(configFile):
    """
    Run the evaluation methods for the pressure distributions, track
    density, landfall rates and longitude crossing rates.

    :param str configFile: path to the configuration file.

    """

    PD = PressureDistribution(configFile)
    TD = TrackDensity(configFile)
    LC = LongitudeCrossing(configFile)
    LF = LandfallRates(configFile)
    GD = GenesisDensity(configFile)

    PD.run()
    TD.run()
    LC.run()
    LF.run()
    GD.run()
