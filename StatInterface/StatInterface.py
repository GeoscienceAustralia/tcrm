"""
Model Calibration
"""
import sys
import logging as log
import KDEOrigin
import KDEParameters

from os.path import join as pjoin
from Utilities.config import cnfGetIniValue
from GenerateDistributions import GenerateDistributions
from generateStats import GenerateStats


class StatInterface(object):

    """
    Parameters
    ----------
    None

    Members
    -------
    lon_lat : string (file name including path)
        latitude and longitude data for cyclone origins
    init_bearing : string (file name including path)
        initial bearings of cyclones
    bearing_no_init : string (file name including path)
        bearings of cyclones with no initial bearing
    all_bearing : string (file name including path)
        all bearings of cyclones
    init_speed : string (file name including path)
        initial speeds of cyclones
    speed_no_init : string (file name including path)
        speeds of cyclones with no initial speed
    all_speed : string (file name including path)
        all speeds of cyclones
    init_pressure : string (file name including path)
        initial pressures of cyclones
    pressure_no_init : string (file name including path)
        pressures of cyclones with no initial pressure
    all_pressure : string (file name including path)
        all pressures of cyclones
    kdeType : string
        kernel density estimation type
    ns : int
        number of samples

    Methods
    -------
    kdeBearing()
        generate KDEs relating to the bearings of cyclones
    kdeSpeed()
        generate KDEs relating to the speeds of cyclones
    kdePressure()
        generate KDEs relating to the pressures of cyclones
    kdeCell(cellNo)
        generate cyclone parameter KDEs relating to a particular cell

    """

    def __init__(self, configFile, autoCalc_gridLimit=None,
                 progressbar=None):
        """
        Initialize the data and variables required for the interface
        """
        self.configFile = configFile
        self.progressbar = progressbar

        log.info("Initialising StatInterface")

        self.kdeType = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'kdeType', 'Biweight')
        self.kde2DType = cnfGetIniValue(self.configFile, 'StatInterface',
                                        'kde2DType', 'Gaussian')
        self.ns = cnfGetIniValue(self.configFile, 'StatInterface',
                                 'Samples', 50000)
        minSamplesCell = cnfGetIniValue(self.configFile, 'StatInterface',
                                        'minSamplesCell', 100)
        self.kdeStep = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'kdeStep', 0.2)
        self.outputPath = cnfGetIniValue(self.configFile, 'Output', 'Path')
        self.processPath = pjoin(self.outputPath, 'process')
        missingValue = cnfGetIniValue(self.configFile, 'StatInterface',
                                      'MissingValue', sys.maxint)

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
            gridSpace = eval(cnfGetIniValue(self.configFile, 'StatInterface',
                                            'gridSpace'))
            gridInc = eval(cnfGetIniValue(self.configFile, 'StatInterface',
                                          'gridInc'))
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
        Generate 2D KDEs relating to the origin of cyclones
        """
        log.info('Generating 2D PDF of TC origins')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'origin_lon_lat'))

        kde = KDEOrigin.KDEOrigin(self.configFile, self.kde2DType,
                                  self.gridLimit, 0.1,
                                  progressbar=self.progressbar)
        kde.generateKDE(None, save=True)
        kde.generateCdf()

    def kdeGenesisDate(self):
        """
        Generate CDFs relating to the genesis day-of-year of cyclones
        """
        log.info('Generating CDFs for TC genesis day')
        log.debug('Reading data from %s',
                  pjoin(self.processPath, 'jdays'))
        pList = pjoin(self.processPath, 'jdays')
        lonLat = pjoin(self.processPath, 'init_lon_lat')
        kde = KDEParameters.KDEParameters(self.kdeType)
        #kde.generateGenesisDateCDF(jdays, lonLat, bw=14,
        #                           genesisKDE=pjoin(self.processPath,
        #                                            'cdfGenesisDays'))
        self.generateDist.allDistributions(lonLat, pList, 'init_day', 
                                           kdeStep=0.25, periodic=365)

    def cdfCellBearing(self):
        """
        generate CDFs relating to the bearing of cyclones
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
        generate CDFs relating to the pressures of cyclones
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
        generate CDFs relating to the pressures of cyclones
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
        generate CDFs relating to the pressures of cyclones
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

        An optional :attr:`minSample` can be given which sets the
        minimum number of observations in a given cell to calculate the
        statistics.  """

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

        log.debug('Saving cell statistics for speed to netcdf file')
        vStats.save(pjoin(path, 'speed_stats.nc'), 'speed')

        log.debug('Calculating cell statistics for pressure')
        pStats = calculate('all_pressure')

        log.debug('Saving cell statistics for pressure to netcdf file')
        pStats.save(pjoin(path, 'pressure_stats.nc'), 'pressure')

        log.debug('Calculating cell statistics for bearing')
        bStats = calculate('all_bearing', angular=True)

        log.debug('Saving cell statistics for bearing to netcdf file')
        bStats.save(pjoin(path, 'bearing_stats.nc'), 'bearing')

        log.debug('Calculating cell statistics for pressure rate' +
                  ' of change')
        dpStats = calculate('pressure_rate')

        log.debug('Saving cell statistics for pressure rate to netcdf file')
        dpStats.save(pjoin(path, 'pressure_rate_stats.nc'), 'pressure rate')
