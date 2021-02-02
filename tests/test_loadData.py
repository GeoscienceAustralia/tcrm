import unittest
import numpy as np
from numpy.testing import assert_almost_equal

import os
import sys
from datetime import datetime
import pickle
from . import NumpyTestCase
try:
    from . import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

import Utilities.loadData as loadData

class TestInitialPositions(unittest.TestCase):
    """
    Test performance of getInitialPositions()
    """

    def setUp(self):

        self.inputData = pickle.load(open(os.path.join(unittest_dir,
                                                        'test_data',
                                                        'loadDataInput.pkl'), 'rb'))
        #self.indexData = dict(index=self.inputData['index'])
        self.serialData = dict(tcserialno=self.inputData['tcserialno'])
        self.seasonData = dict(season=self.inputData['season'],
                               num=self.inputData['num'])
        self.missingFields = dict(lon=self.inputData['lon'],
                                  lat=self.inputData['lat'])

        self.numData = pickle.load(open(os.path.join(unittest_dir,
                                                      'test_data',
                                                      'loadDataNumber.pkl'), 'rb'))

        self.testIndex = pickle.load(open(os.path.join(unittest_dir,
                                                        'test_data',
                                                        'loadDataIndex.pkl'), 'rb'))
        self.numIndex = pickle.load(open(os.path.join(unittest_dir,
                                                       'test_data',
                                                       'loadNumIndex.pkl'), 'rb'))

    def test_getInitPos_fromSerialNo(self):
        """Test to ensure the function returns correct values based on serial number"""
        idx = loadData.getInitialPositions(self.serialData)
        assert_almost_equal(idx, self.testIndex)

    def test_getInitPos_fromSeason(self):
        """Test to ensure the function returns correct values based on season"""
        idx = loadData.getInitialPositions(self.seasonData)
        assert_almost_equal(idx, self.testIndex)

    def test_getInitPos_fromTCNum(self):
        """Test to ensure the function returns correct values based on TC number"""
        idx = loadData.getInitialPositions(self.numData)
        assert_almost_equal(idx, self.numIndex)

    def test_getInitPos_failure(self):
        """Ensure getInitialPositions fails if insufficient data provided"""
        self.assertRaises(KeyError, loadData.getInitialPositions,
                          self.missingFields)



class TestDateParsing(unittest.TestCase):
    """
    Test performance of ParseDates()
    """

    def setUp(self):
        """ """
        input_file = open(os.path.join(unittest_dir, 'test_data',
                                       'parseDates.pkl'), 'rb')
        self.dateformat = '%Y-%m-%d %H:%M:%S'
        self.inputData = pickle.load(input_file)
        self.indicator = pickle.load(input_file)
        self.year = pickle.load(input_file)
        self.month = pickle.load(input_file)
        self.day = pickle.load(input_file)
        self.hour = pickle.load(input_file)
        self.minute = pickle.load(input_file)
        # For testing 'HHMM' formatted times:
        self.hourmin = pickle.load(input_file)

        input_file.close()
        self.input_dates = dict(date=self.inputData['date'])


    def test_dateInput(self):
        """Test parseDates returns correct values when passed date info"""
        year, month, day, hour, minute, dt = loadData.parseDates(self.input_dates,
                                                             self.indicator)
        assert_almost_equal(year, self.year)
        assert_almost_equal(month, self.month)
        assert_almost_equal(day, self.day)
        assert_almost_equal(hour, self.hour)
        assert_almost_equal(minute, self.minute)

    def test_parseDatesYMDHMInput(self):
        """Test parseDates with year, month, day, hour, minute input"""
        inputdata = dict(year=self.year,
                         month=self.month,
                         day=self.day,
                         hour=self.hour,
                         minute=self.minute)
        year, month, day, hour, minute, dt = loadData.parseDates(inputdata,
                                                             self.indicator)

        assert_almost_equal(year, self.year)
        assert_almost_equal(month, self.month)
        assert_almost_equal(day, self.day)
        assert_almost_equal(hour, self.hour)
        assert_almost_equal(minute, self.minute)

    def test_parseDatesYMDHInput(self):
        """Test parseDates with year, month, day, hourminute (HHMM) input"""
        inputdata = dict(year=self.year,
                         month=self.month,
                         day=self.day,
                         hour=self.hourmin)
        year, month, day, hour, minute, dt = loadData.parseDates(inputdata,
                                                             self.indicator)

        assert_almost_equal(year, self.year)
        assert_almost_equal(month, self.month)
        assert_almost_equal(day, self.day)
        assert_almost_equal(hour, self.hour)
        assert_almost_equal(minute, self.minute)

    def test_ParseDatesNoMinsInput(self):
        """Test parseDates with year, month, day, hour (no minutes) input"""
        inputdata = dict(year=self.year,
                         month=self.month,
                         day=self.day,
                         hour=self.hour)
        year, month, day, hour, minute, dt = loadData.parseDates(inputdata,
                                                             self.indicator)

        assert_almost_equal(year, self.year)
        assert_almost_equal(month, self.month)
        assert_almost_equal(day, self.day)
        assert_almost_equal(hour, self.hour)
        assert_almost_equal(minute, np.zeros((self.hour.size), 'i'))

class TestDateConversion(unittest.TestCase):
    def setUp(self):

        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'date2ymhd.pkl'), 'rb')
        self.goodInputDates = pickle.load(inputFile)
        self.badInputDates = pickle.load(inputFile)
        self.dateformat = '%Y-%m-%d %H:%M:%S'
        self.outputYear = pickle.load(inputFile)
        self.outputMonth = pickle.load(inputFile)
        self.outputDay = pickle.load(inputFile)
        self.outputHour = pickle.load(inputFile)
        self.outputMinute = pickle.load(inputFile)
        inputFile.close()


    def test_date2ymdh(self):
        """Test date2ymdh function"""
        year, month, day, hour, minute, dt = loadData.date2ymdh(self.goodInputDates)
        assert_almost_equal(year, self.outputYear)
        assert_almost_equal(month, self.outputMonth)
        assert_almost_equal(day, self.outputDay)
        assert_almost_equal(hour, self.outputHour)
        assert_almost_equal(minute, self.outputMinute)

    def test_date2ymdhBadFormat(self):
        """Test date2ymdh raises ValueError for poorly formatted year data"""

        datefmt = '%H:%M %m/%d/%y'
        now = datetime.now().strftime(datefmt)
        self.assertRaises(ValueError, loadData.date2ymdh, now, datefmt)

    def test_date2ymdhFormats(self):
        """Test date2ymdh with different input date formats"""

        formats = ['%Y-%m-%d %H:%M:%S',
                   '%Y%m%dT%H%M',
                   '%H:%M %d/%m/%Y',
                   '%H:%M %m/%d/%Y',
                   '%I:%M %p %d/%m/%Y']
        for fmt in formats:
            dates = []
            for d in self.goodInputDates:
                dtobj = datetime.strptime(d, self.dateformat)
                datestr = dtobj.strftime(fmt)
                dates.append(datestr)
            year, month, day, hour, minute, dt = loadData.date2ymdh(dates, fmt)
            assert_almost_equal(self.outputYear, year)
            assert_almost_equal(self.outputMonth, month)
            assert_almost_equal(self.outputDay, day)
            assert_almost_equal(self.outputHour, hour)
            assert_almost_equal(self.outputMinute, minute)

    def test_badData(self):
        """Test date2ymdh raises ValueError for dodgy input date"""
        self.assertRaises(ValueError, loadData.date2ymdh, self.badInputDates)


class TestAgeParsing(unittest.TestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                                   'parseAge.pkl'), 'rb')

        self.inputData = pickle.load(inputFile)
        self.indicator = pickle.load(inputFile)

        self.outputYear = pickle.load(inputFile)
        self.outputMonth = pickle.load(inputFile)
        self.outputDay = pickle.load(inputFile)
        self.outputHour = pickle.load(inputFile)
        self.outputMinute = pickle.load(inputFile)
        inputFile.close()

#    def test_parseAge(self):
#        """Test parseAge function"""
#        year, month, day, hour, minute = loadData.parseAge(self.inputData,
#                                                           self.indicator)
#        assert_almost_equal(self.outputYear, year)
#        assert_almost_equal(self.outputMonth, month)
#        assert_almost_equal(self.outputDay, day)
#        assert_almost_equal(self.outputHour, hour)
#        assert_almost_equal(self.outputMinute, minute)


class TestTimeDeltas(unittest.TestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'getTimeDelta.pkl'), 'rb')
        self.inputYear = pickle.load(inputFile)
        self.inputMonth = pickle.load(inputFile)
        self.inputDay = pickle.load(inputFile)
        self.inputHour = pickle.load(inputFile)
        self.inputMinute = pickle.load(inputFile)
        self.outputDT = pickle.load(inputFile)
        inputFile.close()

    def test_getTimeDelta(self):
        """Test getTimeDelta function"""
        dt = loadData.getTimeDelta(self.inputYear,
                                   self.inputMonth,
                                   self.inputDay,
                                   self.inputHour,
                                   self.inputMinute)

        assert_almost_equal(dt, self.outputDT)

    def test_getTimeDeltaBadInput(self):
        """Test getTimeDelta raises ValueError on bad input"""
        inputMonth = self.inputMonth
        inputMonth[345] = 13

        badMonthArgs = [self.inputYear, inputMonth, self.inputDay,
                        self.inputHour, self.inputMinute]

        inputYear = self.inputYear
        inputYear[126] = -1
        badYearArgs = [inputYear, self.inputMonth, self.inputDay,
                        self.inputHour, self.inputMinute]

        self.assertRaises(ValueError, loadData.getTimeDelta,
                          *badMonthArgs)
        self.assertRaises(ValueError, loadData.getTimeDelta,
                          *badYearArgs)

class TestTime(unittest.TestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'getTime.pkl'), 'rb')
        self.inputYear = pickle.load(inputFile)
        self.inputMonth = pickle.load(inputFile)
        self.inputDay = pickle.load(inputFile)
        self.inputHour = pickle.load(inputFile)
        self.inputMinute = pickle.load(inputFile)
        self.outputTime = pickle.load(inputFile)
        inputFile.close()

    def test_getTime(self):
        """Test getTime function"""
        time = loadData.getTime(self.inputYear,
                                self.inputMonth,
                                self.inputDay,
                                self.inputHour,
                                self.inputMinute)

        assert_almost_equal(time, self.outputTime)

    def test_getTimeBadInput(self):
        """Test getTime raises ValueError on bad input"""
        inputMonth = self.inputMonth
        inputMonth[345] = 13

        badMonthArgs = [self.inputYear, inputMonth, self.inputDay,
                        self.inputHour, self.inputMinute]

        inputYear = self.inputYear
        inputYear[126] = -1
        badYearArgs = [inputYear, self.inputMonth, self.inputDay,
                        self.inputHour, self.inputMinute]

        self.assertRaises(ValueError, loadData.getTime,
                          *badMonthArgs)
        self.assertRaises(ValueError, loadData.getTime,
                          *badYearArgs)

class TestJulianDays(unittest.TestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'julianDays.pkl'), 'rb')
        self.inputYear = pickle.load(inputFile)
        self.inputMonth = pickle.load(inputFile)
        self.inputDay = pickle.load(inputFile)
        self.inputHour = pickle.load(inputFile)
        self.inputMinute = pickle.load(inputFile)
        self.outputJdays = pickle.load(inputFile)
        inputFile.close()

    def test_julianDays(self):
        """Test julianDays function"""
        jday = loadData.julianDays(self.inputYear,
                                   self.inputMonth,
                                   self.inputDay,
                                   self.inputHour,
                                   self.inputMinute)

        assert_almost_equal(jday, self.outputJdays)

    def test_julianDaysBadInput(self):
        """Test julianDays raises ValueError on bad input"""
        inputMonth = self.inputMonth
        inputMonth[345] = 13

        badMonthArgs = [self.inputYear, inputMonth, self.inputDay,
                        self.inputHour, self.inputMinute]

        inputYear = self.inputYear
        inputYear[126] = -1
        badYearArgs = [inputYear, self.inputMonth, self.inputDay,
                        self.inputHour, self.inputMinute]

        self.assertRaises(ValueError, loadData.julianDays,
                          *badMonthArgs)
        self.assertRaises(ValueError, loadData.julianDays,
                          *badYearArgs)

class TestGetPoci(unittest.TestCase):

    def setUp(self):
        np.random.seed(10)
        self.penv = np.arange(1000, 1011, 1)
        self.pcentre = np.arange(900, 1001, 10, dtype=float)
        self.lat = np.arange(-24, -2, 2)
        self.jdays = np.arange(1, 365, 36)

        self.pociOutput = np.array([1016.11315678, 1014.75370624,
                                    1014.00665338, 1013.6080981,
                                    1013.23285869, 1012.61540945,
                                    1011.64863179, 1010.42365243,
                                    1009.1959527, 1008.29035302,
                                    1007.98020964])
        self.pociOutputCoeffs = np.array([])

    def test_getPociDefaults(self):
        """Test getPoci returns correct value based on defaults"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        Poci = loadData.getPoci(1000, 900, -24, 1, eps)
        assert_almost_equal(Poci, 1016.1131567762006)

    def test_getPociWrongLengths(self):
        """getPoci raises exception when inputs are different lengths"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        self.assertRaises(Exception, loadData.getPoci, self.penv[:-1],
                          self.pcentre, self.lat, self.jdays, eps)

    def test_getPociArrayInput(self):
        """Test getPoci with array input"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        Poci = loadData.getPoci(self.penv, self.pcentre, self.lat,
                                self.jdays, eps)
        assert_almost_equal(Poci, self.pociOutput)

    def test_getPociArrayMissingValues(self):
        """getPoci filters values where input data is missing"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        pcentre = self.pcentre
        pcentre[-1] = sys.maxsize
        Poci = loadData.getPoci(self.penv, pcentre, self.lat,
                                self.jdays, eps)
        PociOutput = np.array([1016.11315678, 1014.75370624,
                               1014.00665338, 1013.6080981,
                               1013.23285869, 1012.61540945,
                               1011.64863179, 1010.42365243,
                               1009.1959527, 1008.29035302,
                               np.nan])
        assert_almost_equal(Poci, PociOutput)

    def test_getPociArrayInvalidInput(self):
        """getPoci filters values where penv < pcentre"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        penv = self.penv
        penv[-1] = 999
        Poci = loadData.getPoci(self.penv, self.pcentre, self.lat,
                                self.jdays, eps)
        PociOutput = np.array([1016.11315678, 1014.75370624,
                               1014.00665338, 1013.6080981,
                               1013.23285869, 1012.61540945,
                               1011.64863179, 1010.42365243,
                               1009.1959527, 1008.29035302,
                               np.nan])
        assert_almost_equal(Poci, PociOutput)
    
    def test_getPociWithCoeffs(self):
        """getPoci with user-defined set of coefficients"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        coeffs = [2288, -0.65, -1.33, 7.0e-04, 5e-03, -1.5]
        Poci = loadData.getPoci(1000, 900, -24, 1, eps, coeffs)
        self.assertAlmostEqual(Poci, 1012.8047170899531)

    def test_getPociIncompleteCoeffs(self):
        """Test getPoci falls back to default coeffs if input incomplete"""
        eps = np.random.normal(0, scale=2.5717404300409674)
        coeffs = [-0.6496398,-1.33467, 
                  7.085303e-04, 4.87049101e-03,
                  -1.43573905]
        Poci = loadData.getPoci(self.penv, self.pcentre, self.lat,
                                self.jdays, eps, coeffs)
        assert_almost_equal(Poci, self.pociOutput)
    
        
#class TestFilterPressure(unittest.TestCase):
#
#    def setUp(self):
#        inputFile = open(os.path.join(unittest_dir, 'test_data',
#                                      'filterPressure.pkl'), 'rb')
#        self.inputdata = cPickle.load(inputFile)
#        self.outputdata = cPickle.load(inputFile)
#        inputFile.close()
#
#    def test_filterPressure(self):
#        """Test filterPressure function"""
#        result = loadData.filterPressure(self.inputdata)
#        assert_almost_equal(result, self.outputdata, decimal=5)

if __name__ == "__main__":
    unittest.main(verbosity=2)
