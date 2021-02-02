"""
:mod:`CalcFrequency` -- Calculate annual genesis frequency
==========================================================

.. module:: CalcFrequency
   :synopsis: Calculate annual genesis frequency of TCs within the track
              generation domain.

.. moduleauthor:: Nicholas Summons, <nicholas.summons@ga.gov.au>

"""

from os.path import join as pjoin
import logging
import numpy as np

from Utilities.files import flLoadFile
from Utilities.config import ConfigParser

logger = logging.getLogger(__name__)

class CalcFrequency(object):
    """
    Calculate the annual mean frequency of TC events
    based on input dataset. The frequency is calculated
    for the given domain.

    :param dict tg_domain: the domain where the tracks will be generated.
                      The :class:`dict` should contain the keys :attr:`xMin`,
                      :attr:`xMax`, :attr:`yMin` and :attr:`yMax`. The *x*
                      variable bounds the longitude and the *y* variable bounds
                      the latitude.
    """

    def __init__(self, configFile, auto_calc_grid_limit):
        """
        :type  configFile: string
        :param configFile: Configuration file name

        :type  auto_calc_grid_limit: :class:`dict`
        :param auto_calc_grid_limit: the domain where the frequency will be calculated.
                                     The :class:`dict` should contain the keys
                                     :attr:`xMin`, :attr:`xMax`, :attr:`yMin`
                                     and :attr:`yMax`. The *x*  variable bounds the
                                     longitude and the *y* variable bounds
                                     the latitude.
        """

        config = ConfigParser()
        config.read(configFile)

        if config.has_option('TrackGenerator', 'gridLimit'):
            self.tg_domain = config.geteval('TrackGenerator', 'gridLimit')
        else:
            self.tg_domain = auto_calc_grid_limit

        self.outputPath = config.get('Output', 'Path')

    def calc(self):
        """
        Calculate the frequency of TC events in a pre-defined domain, based
        on the input dataset and the full range of years contained in the
        :attr:`origin_year` file.

        The :attr:`origin_year` file is created in
        :class:`DataProcess.processData` and restricts the range to a
        user-selected range of years.
        """

        logger.info("Calculating annual frequency of TC events")
        origin_year = np.array(flLoadFile(pjoin(self.outputPath,
                                                'process', 'origin_year'),
                                          '%', ','), dtype='int')

        origin_lon_lat = flLoadFile(pjoin(self.outputPath,
                                          'process', 'origin_lon_lat'),
                                    '%', ',')
        origin_lon = origin_lon_lat[:, 0]
        origin_lat = origin_lon_lat[:, 1]
        min_year = origin_year.min()
        # Skip last year from average since may contain only partial year record
        max_year = origin_year.max() - 1

        freq_count = np.zeros(3000)

        for year in range(min_year, max_year + 1):
            freq_count[year] = sum((origin_year == year) & \
                                 (origin_lon > self.tg_domain['xMin']) & \
                                 (origin_lon < self.tg_domain['xMax']) & \
                                 (origin_lat > self.tg_domain['yMin']) & \
                                 (origin_lat < self.tg_domain['yMax']))

        freq = np.mean(freq_count[min_year:max_year + 1])
        freq = np.round(freq*100)/100

        fname = pjoin(self.outputPath, 'process', 'region_frequency')
        data = np.array([np.arange(min_year, max_year + 1),
                         freq_count[min_year:max_year + 1]])
        header = "Year,count"
        np.savetxt(fname, data.T, fmt="%d", delimiter=",", header=header)

        return freq
