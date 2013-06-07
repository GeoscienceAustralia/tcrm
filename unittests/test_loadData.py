#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 Title: test_loadData.py - unittests for loadData.py
 Author: Craig Arthur, craig.arthur@ga.gov.au
 CreationDate: Mon May 20 13:24:00 2013
 Description: Unit test suite for loadData module.

 Version: $Rev$
 Id: $Id$

"""

import os
import sys
import numpy
from datetime import datetime
from unittest import TestSuite, TestLoader, TextTestRunner
import cPickle
import NumpyTestCase
try:
    import pathLocate
except:
    from unittests import pathLocate

# Add parent folder to python path
unittest_dir = pathLocate.getUnitTestDirectory()
sys.path.append(pathLocate.getRootDirectory())

import Utilities.loadData as loadData

class test_getInitialPositions(NumpyTestCase.NumpyTestCase):
    """
    Test performance of getInitialPositions()
    """

    def setUp(self):

        self.inputData = cPickle.load(open(os.path.join(unittest_dir,
                                                        'test_data',
                                                        'loadDataInput.pck')))
        #self.indexData = dict(index=self.inputData['index'])
        self.serialData = dict(tcserialno=self.inputData['tcserialno'])
        self.seasonData = dict(season=self.inputData['season'],
                               num=self.inputData['num'])
        self.missingFields = dict(lon=self.inputData['lon'],
                                  lat=self.inputData['lat'])

        self.numData = cPickle.load(open(os.path.join(unittest_dir,
                                                      'test_data',
                                                      'loadDataNumber.pck')))

        self.testIndex = cPickle.load(open(os.path.join(unittest_dir,
                                                        'test_data',
                                                        'loadDataIndex.pck')))
        self.numIndex = cPickle.load(open(os.path.join(unittest_dir,
                                                       'test_data',
                                                       'loadNumIndex.pck')))

    #def test_getInitPos_fromIndex(self):
    #    """Test to ensure the function returns correct values based on index"""
    #    idx = loadData.getInitialPositions(self.indexData)
    #    self.numpyAssertEqual(idx, self.testIndex)


    def test_getInitPos_fromSerialNo(self):
        """Test to ensure the function returns correct values based on serial number"""
        idx = loadData.getInitialPositions(self.serialData)
        self.numpyAssertEqual(idx, self.testIndex)

    def test_getInitPos_fromSeason(self):
        """Test to ensure the function returns correct values based on season"""
        idx = loadData.getInitialPositions(self.seasonData)
        self.numpyAssertEqual(idx, self.testIndex)

    def test_getInitPos_fromTCNum(self):
        """Test to ensure the function returns correct values based on TC number"""
        idx = loadData.getInitialPositions(self.numData)
        self.numpyAssertEqual(idx, self.numIndex)

    def test_getInitPos_failure(self):
        """Ensure getInitialPositions fails if insufficient data provided"""

        self.assertRaises(KeyError, loadData.getInitialPositions,
                                      self.missingFields)



class test_parseDates(NumpyTestCase.NumpyTestCase):
    """
    Test performance of ParseDates()
    """

    def setUp(self):
        """ """
        input_file = open(os.path.join(unittest_dir, 'test_data',
                                       'parseDates.pck'))
        self.dateformat = '%Y-%m-%d %H:%M:%S'
        self.inputData = cPickle.load(input_file)
        self.indicator = cPickle.load(input_file)
        self.year = cPickle.load(input_file)
        self.month = cPickle.load(input_file)
        self.day = cPickle.load(input_file)
        self.hour = cPickle.load(input_file)
        self.minute = cPickle.load(input_file)
        # For testing 'HHMM' formatted times:
        self.hourmin = cPickle.load(input_file)

        input_file.close()
        self.input_dates = dict(date=self.inputData['date'])


    def test_dateInput(self):
        """Test parseDates returns correct values when passed date info"""
        year, month, day, hour, minute = loadData.parseDates(self.input_dates,
                                                             self.indicator)
        self.numpyAssertEqual(year, self.year)
        self.numpyAssertEqual(month, self.month)
        self.numpyAssertEqual(day, self.day)
        self.numpyAssertEqual(hour, self.hour)
        self.numpyAssertEqual(minute, self.minute)

    def test_parseDatesYMDHMInput(self):
        """Test parseDates with year, month, day, hour, minute input"""
        inputdata = dict(year=self.year,
                         month=self.month,
                         day=self.day,
                         hour=self.hour,
                         minute=self.minute)
        year, month, day, hour, minute = loadData.parseDates(inputdata,
                                                             self.indicator)

        self.numpyAssertEqual(year, self.year)
        self.numpyAssertEqual(month, self.month)
        self.numpyAssertEqual(day, self.day)
        self.numpyAssertEqual(hour, self.hour)
        self.numpyAssertEqual(minute, self.minute)

    def test_parseDatesYMDHInput(self):
        """Test parseDates with year, month, day, hourminute (HHMM) input"""
        inputdata = dict(year=self.year,
                         month=self.month,
                         day=self.day,
                         hour=self.hourmin)
        year, month, day, hour, minute = loadData.parseDates(inputdata,
                                                             self.indicator)

        self.numpyAssertEqual(year, self.year)
        self.numpyAssertEqual(month, self.month)
        self.numpyAssertEqual(day, self.day)
        self.numpyAssertEqual(hour, self.hour)
        self.numpyAssertEqual(minute, self.minute)

    def test_ParseDatesNoMinsInput(self):
        """Test parseDates with year, month, day, hour (no minutes) input"""
        inputdata = dict(year=self.year,
                         month=self.month,
                         day=self.day,
                         hour=self.hour)
        year, month, day, hour, minute = loadData.parseDates(inputdata,
                                                             self.indicator)

        self.numpyAssertEqual(year, self.year)
        self.numpyAssertEqual(month, self.month)
        self.numpyAssertEqual(day, self.day)
        self.numpyAssertEqual(hour, self.hour)
        self.numpyAssertEqual(minute, numpy.zeros((self.hour.size), 'i'))

class test_date2ymdh(NumpyTestCase.NumpyTestCase):
    def setUp(self):

        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'date2ymhd.pck'))
        self.goodInputDates = cPickle.load(inputFile)
        self.badInputDates = cPickle.load(inputFile)
        self.dateformat = '%Y-%m-%d %H:%M:%S'
        self.outputYear = cPickle.load(inputFile)
        self.outputMonth = cPickle.load(inputFile)
        self.outputDay = cPickle.load(inputFile)
        self.outputHour = cPickle.load(inputFile)
        self.outputMinute = cPickle.load(inputFile)
        inputFile.close()


    def test_date2ymdh(self):
        """Test date2ymdh function"""
        year, month, day, hour, minute = loadData.date2ymdh(self.goodInputDates)
        self.numpyAssertEqual(year, self.outputYear)
        self.numpyAssertEqual(month, self.outputMonth)
        self.numpyAssertEqual(day, self.outputDay)
        self.numpyAssertEqual(hour, self.outputHour)
        self.numpyAssertEqual(minute, self.outputMinute)

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
            year, month, day, hour, minute = loadData.date2ymdh(dates, fmt)
            self.numpyAssertEqual(self.outputYear, year)
            self.numpyAssertEqual(self.outputMonth, month)
            self.numpyAssertEqual(self.outputDay, day)
            self.numpyAssertEqual(self.outputHour, hour)
            self.numpyAssertEqual(self.outputMinute, minute)

    def test_badData(self):
        """Test date2ymdh raises ValueError for dodgy input date"""
        self.assertRaises(ValueError, loadData.date2ymdh, self.badInputDates)


class test_parseAge(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                                   'parseAge.pck'))

        self.inputData = cPickle.load(inputFile)
        self.indicator = cPickle.load(inputFile)

        self.outputYear = cPickle.load(inputFile)
        self.outputMonth = cPickle.load(inputFile)
        self.outputDay = cPickle.load(inputFile)
        self.outputHour = cPickle.load(inputFile)
        self.outputMinute = cPickle.load(inputFile)
        inputFile.close()

    def test_parseAge(self):
        """Test parseAge function"""
        year, month, day, hour, minute = loadData.parseAge(self.inputData,
                                                           self.indicator)
        self.numpyAssertEqual(self.outputYear, year)
        self.numpyAssertEqual(self.outputMonth, month)
        self.numpyAssertEqual(self.outputDay, day)
        self.numpyAssertEqual(self.outputHour, hour)
        self.numpyAssertEqual(self.outputMinute, minute)


class test_getTimeDelta(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'getTimeDelta.pck'))
        self.inputYear = cPickle.load(inputFile)
        self.inputMonth = cPickle.load(inputFile)
        self.inputDay = cPickle.load(inputFile)
        self.inputHour = cPickle.load(inputFile)
        self.inputMinute = cPickle.load(inputFile)
        self.outputDT = cPickle.load(inputFile)
        inputFile.close()

    def test_getTimeDelta(self):
        """Test getTimeDelta function"""
        dt = loadData.getTimeDelta(self.inputYear,
                                   self.inputMonth,
                                   self.inputDay,
                                   self.inputHour,
                                   self.inputMinute)

        self.numpyAssertEqual(dt, self.outputDT)

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

class test_getTime(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'getTime.pck'))
        self.inputYear = cPickle.load(inputFile)
        self.inputMonth = cPickle.load(inputFile)
        self.inputDay = cPickle.load(inputFile)
        self.inputHour = cPickle.load(inputFile)
        self.inputMinute = cPickle.load(inputFile)
        self.outputTime = cPickle.load(inputFile)
        inputFile.close()

    def test_getTime(self):
        """Test getTime function"""
        time = loadData.getTime(self.inputYear,
                                self.inputMonth,
                                self.inputDay,
                                self.inputHour,
                                self.inputMinute)

        self.numpyAssertAlmostEqual(time, self.outputTime)

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

class test_julianDays(NumpyTestCase.NumpyTestCase):

    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'julianDays.pck'))
        self.inputYear = cPickle.load(inputFile)
        self.inputMonth = cPickle.load(inputFile)
        self.inputDay = cPickle.load(inputFile)
        self.inputHour = cPickle.load(inputFile)
        self.inputMinute = cPickle.load(inputFile)
        self.outputJdays = cPickle.load(inputFile)
        inputFile.close()

    def test_julianDays(self):
        """Test julianDays function"""
        jday = loadData.julianDays(self.inputYear,
                                   self.inputMonth,
                                   self.inputDay,
                                   self.inputHour,
                                   self.inputMinute)

        self.numpyAssertAlmostEqual(jday, self.outputJdays)

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

class test_loadTrackFile(NumpyTestCase.NumpyTestCase):
    def setUp(self):
        self.config_file = os.path.join(unittest_dir, 'test_data',
                                      'test.ini')
        self.track_file = os.path.join(unittest_dir, 'test_data',
                                      'test_trackset.csv')
        self.source = 'TESTSOURCE'

        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'loadTrackFile.pck'))
        self.trackData = cPickle.load(inputFile)

class test_filterPressure(NumpyTestCase.NumpyTestCase):
    def setUp(self):
        inputFile = open(os.path.join(unittest_dir, 'test_data',
                                      'filterPressure.pck'))
        self.inputdata = cPickle.load(inputFile)
        self.outputdata = cPickle.load(inputFile)
        inputFile.close()

    def test_filterPressure(self):
        """Test filterPressure function"""
        result = loadData.filterPressure(self.inputdata)
        self.numpyAssertEqual(result, self.outputdata)

if __name__ == "__main__":
    #testSuite = unittest.makeSuite(test_getInitialPositions,'test')
    initPosSuite = TestLoader().loadTestsFromTestCase(test_getInitialPositions)
    parseAgeSuite = TestLoader().loadTestsFromTestCase(test_parseAge)
    parseDateSuite = TestLoader().loadTestsFromTestCase(test_parseDates)
    date2ymdhSuite = TestLoader().loadTestsFromTestCase(test_date2ymdh)
    timeDeltaSuite = TestLoader().loadTestsFromTestCase(test_getTimeDelta)
    timeSuite = TestLoader().loadTestsFromTestCase(test_getTime)
    jdaySuite = TestLoader().loadTestsFromTestCase(test_julianDays)
    FilterPrsSuite = TestLoader().loadTestsFromTestCase(test_filterPressure)

    testSuite = TestSuite([initPosSuite,
                           parseAgeSuite,
                           parseDateSuite,
                           date2ymdhSuite,
                           timeDeltaSuite,
                           timeSuite,
                           jdaySuite,
                           FilterPrsSuite])

    TextTestRunner(verbosity=1).run(testSuite)
