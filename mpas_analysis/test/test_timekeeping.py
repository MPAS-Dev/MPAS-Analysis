"""
Unit test infrastructure for the Date class

Xylar Asay-Davis
11/02/2016
"""

import pytest
import datetime
from dateutil.relativedelta import relativedelta
from mpas_analysis.test import TestCase
from mpas_analysis.shared.timekeeping.utility import stringToDatetime, \
    stringToRelativedelta, clampToNumpyDatetime64


class TestTimekeeping(TestCase):
    def test_timekeeping(self):

        # test each possible format:
        # YYYY-MM-DD_hh:mm:ss
        # YYYY-MM-DD_hh.mm.ss
        # YYYY-MM-DD_SSSSS
        # DDD_hh:mm:ss
        # DDD_hh.mm.ss
        # DDD_SSSSS
        # hh.mm.ss
        # hh:mm:ss
        # YYYY-MM-DD
        # SSSSS

        # test datetime.datetime
        # YYYY-MM-DD_hh:mm:ss
        date1 = stringToDatetime('0001-01-01_00:00:00')
        date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                  second=0)
        self.assertEqual(date1, date2)

        # test relativedelta
        # YYYY-MM-DD_hh:mm:ss
        # gregorian_noleap
        date1 = stringToRelativedelta('0001-00-00_00:00:00',
                                      calendar='gregorian_noleap')
        date2 = relativedelta(years=1, months=0, days=0, hours=0,
                              minutes=0, seconds=0, leapdays=False)
        self.assertEqual(date1, date2)

        # gregorian
        date1 = stringToRelativedelta('0001-00-00_00:00:00',
                                      calendar='gregorian')
        date2 = relativedelta(years=1, months=0, days=0, hours=0,
                              minutes=0, seconds=0, leapdays=True)
        self.assertEqual(date1, date2)

        # YYYY-MM-DD_hh.mm.ss
        date1 = stringToDatetime('0001-01-01_00.00.00')
        date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                  second=0)
        self.assertEqual(date1, date2)

        # YYYY-MM-DD_SSSSS
        date1 = stringToDatetime('0001-01-01_00002')
        date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                  second=2)
        self.assertEqual(date1, date2)

        # DDD_hh:mm:ss
        date1 = stringToRelativedelta('0001_00:00:01',
                                      calendar='gregorian_noleap')
        date2 = relativedelta(years=0, months=0, days=1, hours=0,
                              minutes=0, seconds=1, leapdays=False)
        self.assertEqual(date1, date2)

        # DDD_hh.mm.ss
        date1 = stringToRelativedelta('0002_01.00.01',
                                      calendar='gregorian_noleap')
        date2 = relativedelta(years=0, months=0, days=2, hours=1,
                              minutes=0, seconds=1, leapdays=False)
        self.assertEqual(date1, date2)

        # DDD_SSSSS
        date1 = stringToRelativedelta('0002_00003',
                                      calendar='gregorian_noleap')
        date2 = relativedelta(years=0, months=0, days=2, hours=0,
                              minutes=0, seconds=3, leapdays=False)
        self.assertEqual(date1, date2)

        # hh:mm:ss
        date1 = stringToDatetime('00:00:01')
        date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                  second=1)
        self.assertEqual(date1, date2)

        # hh.mm.ss
        date1 = stringToRelativedelta('00.00.01',
                                      calendar='gregorian')
        date2 = relativedelta(years=0, months=0, days=0, hours=0,
                              minutes=0, seconds=1, leapdays=True)
        self.assertEqual(date1, date2)

        # YYYY-MM-DD
        date1 = stringToDatetime('0001-01-01')
        date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                  second=0)
        self.assertEqual(date1, date2)

        # SSSSS
        date1 = stringToRelativedelta('00005',
                                      calendar='gregorian')
        date2 = relativedelta(years=0, months=0, days=0, hours=0,
                              minutes=0, seconds=5, leapdays=True)
        self.assertEqual(date1, date2)

        date1 = stringToDatetime('1996-01-15')
        date2 = stringToRelativedelta('0005-00-00',
                                      calendar='gregorian')
        diff = date1-date2
        self.assertEqual(diff, stringToDatetime('1991-01-15'))

        date1 = stringToDatetime('1996-01-15')
        date2 = stringToRelativedelta('0000-02-00',
                                      calendar='gregorian')
        diff = date1-date2
        self.assertEqual(diff, stringToDatetime('1995-11-15'))

        date1 = stringToDatetime('1996-01-15')
        date2 = stringToRelativedelta('0000-00-20',
                                      calendar='gregorian')
        diff = date1-date2
        self.assertEqual(diff, stringToDatetime('1995-12-26'))

        # since pandas and xarray use the numpy type 'datetime[ns]`, which
        # has a limited range of dates, the date 0001-01-01 gets increased to
        # the minimum allowed year boundary, 1678-01-01 to avoid invalid
        # dates.
        date1 = clampToNumpyDatetime64(stringToDatetime('0001-01-01'),
                                       yearOffset=0)
        date2 = datetime.datetime(year=1678, month=1, day=1)
        self.assertEqual(date1, date2)

        date1 = clampToNumpyDatetime64(stringToDatetime('0001-01-01'),
                                       yearOffset=1849)
        date2 = datetime.datetime(year=1850, month=1, day=1)
        self.assertEqual(date1, date2)

        # since pandas and xarray use the numpy type 'datetime[ns]`, which
        # has a limited range of dates, the date 9999-01-01 gets decreased to
        # the maximum allowed year boundary, 2262-01-01 to avoid invalid
        # dates.
        date = clampToNumpyDatetime64(stringToDatetime('9999-01-01'),
                                      yearOffset=0)
        datetime2 = datetime.datetime(year=2262, month=1, day=1)
        self.assertEqual(date1, date2)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
