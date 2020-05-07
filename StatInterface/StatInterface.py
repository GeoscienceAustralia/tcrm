"""
:mod:`StatInterface` -- statistical analysis of input datasets
==============================================================

.. module:: StatInterface
    :synopsis: Generate cumulative distribution functions and probability
               density functions of the various parameters,
               largely using kernel density estimation methods.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>


"""

import sys
import logging as log
from . import KDEOrigin
from . import KDEParameters

from os.path import join as pjoin
from Utilities.config import cnfGetIniValue, ConfigParser
from .GenerateDistributions import GenerateDistributions
from .generateStats import GenerateStats


class StatInterface(object):
    """
    Main interface to the statistical analysis module of TCRM. This
    module generates cumulative distribution functions and probability
    density functions of the various parameters, largely using kernel
    density estimation methods.

    :param str configFile: Path to configuration file.
    :param autoCalc_gridLimit: function to calculate the extent of a domain.
    :param progressBar: a :meth:`SimpleProgressBar` object to print
                        progress to STDOUT.


    """

    def __init__(self, configFile, autoCalc_gridLimit=None,
                 progressbar=None):
        """
        Initialize the data and variables required for the interface
        """
        self.configFile = configFile
        config = ConfigParser()
        config.read(configFile)
        self.progressbar = progressbar

        log.info("Initialising StatInterface")

        self.kdeType = config.get('StatInterface', 'kdeType')
        minSamplesCell = config.getint('StatInterface', 'minSamplesCell')
        self.kdeStep = config.getfloat('StatInterface', 'kdeStep')
        self.outputPath = config.get('Output', 'Path')
        self.processPath = pjoin(self.outputPath, 'process')

        missingValue = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'MissingValue', sys.maxsize)

        gridLimitStr = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'gridLimit', '')

        if gridLimitStr is not '':
            try:
                self.gridLimit = eval(gridLimitStr)
            except SyntaxError:
                log.exception('Error! gridLimit is not a dictionary')
        else:
            self.gridLimit = autoCalc_gridLimit
            log.info('No gridLimit specified - using automatic' +
                     ' selection: ' + str(self.gridLimit))

        try:
            gridSpace = config.geteval('Region', 'gridSpace')
            gridInc = config.geteval('Region', 'gridInc')
        except SyntaxError:
            log.exception('Error! gridSpace or gridInc not dictionaries')
            raise

        self.generateDist = GenerateDistributions(self.configFile,
                                                  self.gridLimit,
                                                  gridSpace, gridInc,
                                                  self.kdeType,
                                                  minSamplesCell,
                                                  missingValue)
        self.gridSpace = gridSpace
        self.gridInc = gridInc

    def kdeOrigin(self):
        """
        Generate 2D PDFs relating to the origin of cyclones.

        """
        log.info('Generating 2D PDF of TC origins')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'origin_lon_lat'))

        kde = KDEOrigin.KDEOrigin(self.configFile,
                                  self.gridLimit, 0.1,
                                  progressbar=self.progressbar)
        kde.generateKDE(save=True, plot=True)
        kde.generateCdf()

    def kdeGenesisDate(self):
        """
        Generate CDFs relating to the genesis day-of-year of cyclones
        for each grid cell in teh model domain.

        """
        log.info('Generating CDFs for TC genesis day')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'jdays'))
        pList = pjoin(self.processPath, 'jdays')
        lonLat = pjoin(self.processPath, 'origin_lon_lat')
        kde = KDEParameters.KDEParameters(self.kdeType)
        #kde.generateGenesisDateCDF(jdays, lonLat, bw=14,
        #                           genesisKDE=pjoin(self.processPath,
        #                                            'cdfGenesisDays'))
        self.generateDist.allDistributions(lonLat, pList, 'init_day',
                                           kdeStep=0.25, periodic=365)

    def cdfCellBearing(self):
        """
        Generate CDFs relating to the bearing of cyclones for each
        grid cell in the model domain.

        """
        log.info('Generating CDFs for TC bearing')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'init_lon_lat'))
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'init_bearing'))
        log.debug('Outputting data into %s', pjoin(
            self.processPath, 'all_cell_cdf_init_bearing'))

        lonLat = pjoin(self.processPath, 'init_lon_lat')
        pList = pjoin(self.processPath, 'init_bearing')
        self.generateDist.allDistributions(
            lonLat, pList, 'init_bearing', 1, True)

    def cdfCellSpeed(self):
        """
        Generate CDFs relating to the speed of motion of cyclones for each
        grid cell in the model domain.

        """
        log.info('Generating CDFs for TC speed')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'init_lon_lat'))
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'init_speed'))
        log.debug('Outputting data into %s',
                  pjoin(self.processPath,
                        'all_cell_cdf_init_speed'))

        lonLat = pjoin(self.processPath, 'init_lon_lat')
        pList = pjoin(self.processPath, 'init_speed')
        self.generateDist.allDistributions(lonLat, pList, 'init_speed',
                                           self.kdeStep)

    def cdfCellPressure(self):
        """
        Generate CDFs relating to the pressures of cyclones
        in each grid cell in the model domain.

        """
        log.info('Generating CDFs for TC pressure')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'origin_lon_lat'))
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'init_pressure'))
        log.debug('Outputting data into %s',
                  pjoin(self.processPath,
                        'all_cell_cdf_init_pressure'))

        lonLat = pjoin(self.processPath, 'origin_lon_lat')
        pList = pjoin(self.processPath, 'init_pressure')
        self.generateDist.allDistributions(lonLat, pList, 'init_pressure',
                                           self.kdeStep)

    def cdfCellSize(self):
        """
        Generate CDFs relating to the size (radius of maximum wind)
        of cyclones in each grid cell in the model domain.

        """
        log.info('Generating CDFs for TC size')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'origin_lon_lat'))
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'init_rmax'))
        log.debug('Outputting data into %s',
                  pjoin(self.processPath,
                        'all_cell_cdf_init_rmax'))

        lonLat = pjoin(self.processPath, 'origin_lon_lat')
        pList = pjoin(self.processPath, 'init_rmax')
        self.generateDist.allDistributions(lonLat, pList, 'init_rmax',
                                           self.kdeStep)

    def calcCellStatistics(self, minSample=100):
        """
        Calculate the cell statistics for speed, bearing, pressure, and
        pressure rate of change for all the grid cells in the domain.

        The statistics calculated are mean, variance, and
        autocorrelation.

        The cell statistics are calculated on a grid defined by
        :attr:`gridLimit`, :attr:`gridSpace` and :attr:`gridInc` using
        an instance of
        :class:`StatInterface.generateStats.GenerateStats`.

        An optional :attr:`minSample` (default=100) can be given which
        sets the minimum number of observations in a given cell to
        calculate the statistics.

        """

        path = self.processPath

        def calculate(filename, angular=False):
            """
            Helper function to calculate the statistics.
            """
            return GenerateStats(
                pjoin(path, filename),
                pjoin(path, 'all_lon_lat'),
                self.gridLimit,
                self.gridSpace,
                self.gridInc,
                minSample=minSample,
                angular=angular)

        log.debug('Calculating cell statistics for speed')
        vStats = calculate('all_speed')
        vStats.plotStatistics(pjoin(self.outputPath, 'plots',
                                    'stats', 'speed_stats'))
        log.debug('Saving cell statistics for speed to netcdf file')
        vStats.save(pjoin(path, 'speed_stats.nc'), 'speed')

        dvStats = calculate('speed_rate')
        dvStats.plotStatistics(pjoin(self.outputPath, 'plots',
                                     'stats', 'speed_rate_stats'))
        log.debug('Saving cell statistics for speed rate to netcdf file')
        dvStats.save(pjoin(path, 'speed_rate_stats.nc'), 'speed_rate')

        log.debug('Calculating cell statistics for pressure')
        pStats = calculate('all_pressure')
        pStats.plotStatistics(pjoin(self.outputPath, 'plots',
                                    'stats', 'pressure_stats'))
        log.debug('Saving cell statistics for pressure to netcdf file')
        pStats.save(pjoin(path, 'pressure_stats.nc'), 'pressure')

        log.debug('Calculating cell statistics for pressure rate' +
                  ' of change')
        dpStats = calculate('pressure_rate')
        dpStats.plotStatistics(pjoin(self.outputPath, 'plots',
                                     'stats', 'pressure_rate_stats'))

        log.debug('Saving cell statistics for pressure rate to netcdf file')
        dpStats.save(pjoin(path, 'pressure_rate_stats.nc'), 'pressure_rate')

        log.debug('Calculating cell statistics for bearing')
        bStats = calculate('all_bearing', angular=True)
        bStats.plotStatistics(pjoin(self.outputPath, 'plots',
                                    'stats', 'bearing_stats'))
        log.debug('Saving cell statistics for bearing to netcdf file')
        bStats.save(pjoin(path, 'bearing_stats.nc'), 'bearing')

        log.debug('Calculating cell statistics for bearing rate of change')
        dbStats = calculate('bearing_rate', angular=True)
        dbStats.plotStatistics(pjoin(self.outputPath, 'plots',
                                     'stats', 'bearing_rate_stats'))
        log.debug('Saving cell statistics for bearing to netcdf file')
        dbStats.save(pjoin(path, 'bearing_rate_stats.nc'), 'bearing_rate')

