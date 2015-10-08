"""
:mod:`Evaluate` -- run evaluation of track model
================================================
.. module:: Evaluate
    :synopsis: Provide some evaluation methods for testing
               the track generation model.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

from Evaluate.pressureDistribution import PressureDistribution
from Evaluate.trackDensity import TrackDensity
from Evaluate.longitudeCrossing import LongitudeCrossing
from Evaluate.landfallRates import LandfallRates
from Evaluate.genesisDensity import GenesisDensity


def run(configFile):
    """
    Run the evaluation methods for the pressure distributions, track
    density, landfall rates and longitude crossing rates.

    :param str configFile: path to the configuration file.

    """

    pd = PressureDistribution(configFile)
    td = TrackDensity(configFile)
    lc = LongitudeCrossing(configFile)
    lf = LandfallRates(configFile)
    gd = GenesisDensity(configFile)

    pd.run()
    td.run()
    lc.run()
    lf.run()
    gd.run()
