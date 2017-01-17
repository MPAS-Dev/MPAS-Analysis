"""
Unit test infrastructure for the Date class

Xylar Asay-Davis
11/02/2016
"""

import pytest
import datetime
from mpas_analysis.test import TestCase
from mpas_analysis.shared.timekeeping.Date import Date

class TestDate(TestCase):
    def test_date(self):

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

        # test with isInterval == False
        # YYYY-MM-DD_hh:mm:ss
        date1 = Date(dateString='0001-01-01_00:00:00', isInterval=False)
        date2 = Date(years=1, months=0, days=0, hours=0, minutes=0, seconds=0,
                     isInterval=False)
        self.assertEqual(date1, date2)

        # test with isInterval == True
        # YYYY-MM-DD_hh:mm:ss
        date1 = Date(dateString='0001-00-00_00:00:00', isInterval=True)
        date2 = Date(years=1, months=0, days=0, hours=0, minutes=0, seconds=0,
                     isInterval=True)
        self.assertEqual(date1, date2)

        # YYYY-MM-DD_hh.mm.ss
        date1 = Date(dateString='0001-01-02_00.01.00', isInterval=False)
        date2 = Date(years=1, months=0, days=1, hours=0, minutes=1, seconds=0,
                     isInterval=False)
        self.assertEqual(date1, date2)

        # YYYY-MM-DD_SSSSS
        date1 = Date(dateString='0001-01-01_00002', isInterval=False)
        date2 = Date(years=1, months=0, days=0, hours=0, minutes=0, seconds=2,
                     isInterval=False)
        self.assertEqual(date1, date2)

        # DDD_hh:mm:ss
        date1 = Date(dateString='0001_00:00:01', isInterval=True)
        date2 = Date(years=0, months=0, days=1, hours=0, minutes=0, seconds=1,
                     isInterval=True)
        self.assertEqual(date1, date2)

        # DDD_hh.mm.ss
        date1 = Date(dateString='0002_01.00.01', isInterval=True)
        date2 = Date(years=0, months=0, days=2, hours=1, minutes=0, seconds=1,
                     isInterval=True)
        self.assertEqual(date1, date2)

        # DDD_SSSSS
        date1 = Date(dateString='0002_00003', isInterval=True)
        date2 = Date(years=0, months=0, days=2, hours=0, minutes=0, seconds=3,
                     isInterval=True)
        self.assertEqual(date1, date2)

        # hh:mm:ss
        date1 = Date(dateString='00:00:01', isInterval=False)
        date2 = Date(years=0, months=0, days=0, hours=0, minutes=0, seconds=1,
                     isInterval=False)
        self.assertEqual(date1, date2)

        # hh.mm.ss
        date1 = Date(dateString='00.00.01', isInterval=True)
        date2 = Date(years=0, months=0, days=0, hours=0, minutes=0, seconds=1,
                     isInterval=True)
        self.assertEqual(date1, date2)

        # YYYY-MM-DD
        date1 = Date(dateString='0001-01-01', isInterval=False)
        date2 = Date(years=1, months=0, days=0, hours=0, minutes=0, seconds=0,
                     isInterval=False)
        self.assertEqual(date1, date2)

        # SSSSS
        date1 = Date(dateString='00005', isInterval=True)
        date2 = Date(years=0, months=0, days=0, hours=0, minutes=0, seconds=5,
                     isInterval=True)
        self.assertEqual(date1, date2)


        # test operators
        date1 = Date(dateString='1992-02-01', isInterval=False)
        date2 = Date(dateString='1991-03-01', isInterval=False)
        diff = date1-date2
        self.assertEqual(diff, Date(dateString='0000-11-00', isInterval=True))
        self.assertEqual(date1 < date2, False)
        self.assertEqual(date2 < date1, True)
        self.assertEqual(date1 < date1, False)

        date1 = Date(dateString='1996-01-15', isInterval=False)
        date2 = Date(dateString='0005-00-00', isInterval=True)
        diff = date1-date2
        self.assertEqual(diff, Date(dateString='1991-01-15', isInterval=False))

        date1 = Date(dateString='1996-01-15', isInterval=False)
        date2 = Date(dateString='0000-02-00', isInterval=True)
        diff = date1-date2
        self.assertEqual(diff, Date(dateString='1995-11-15', isInterval=False))

        date1 = Date(dateString='1996-01-15', isInterval=False)
        date2 = Date(dateString='0000-00-20', isInterval=True)
        diff = date1-date2
        self.assertEqual(diff, Date(dateString='1995-12-26', isInterval=False))

        date = Date(dateString='1996-01-15', isInterval=False)
        datetime1 = date.to_datetime(yearOffset=0)
        datetime2 = datetime.datetime(year=1996, month=1, day=15)
        self.assertEqual(datetime1, datetime2)

        date = Date(dateString='0000-00-20', isInterval=True)
        timedelta1 = date.to_timedelta()
        timedelta2 = datetime.timedelta(days=20)
        self.assertEqual(timedelta1, timedelta2)

        # since pandas and xarray use the numpy type 'datetime[ns]`, which
        # has a limited range of dates, the date 0001-01-01 gets increased to
        # the minimum allowed year boundary, 1678-01-01 to avoid invalid
        # dates.
        date = Date(dateString='0001-01-01', isInterval=False)
        datetime1 = date.to_datetime(yearOffset=0)
        datetime2 = datetime.datetime(year=1678, month=1, day=1)
        self.assertEqual(datetime1, datetime2)

        date = Date(dateString='0001-01-01', isInterval=False)
        datetime1 = date.to_datetime(yearOffset=1849)
        datetime2 = datetime.datetime(year=1850, month=1, day=1)
        self.assertEqual(datetime1, datetime2)

        # since pandas and xarray use the numpy type 'datetime[ns]`, which
        # has a limited range of dates, the date 9999-01-01 gets decreased to
        # the maximum allowed year boundary, 2262-01-01 to avoid invalid
        # dates.
        date = Date(dateString='9999-01-01', isInterval=False)
        datetime1 = date.to_datetime(yearOffset=0)
        datetime2 = datetime.datetime(year=2262, month=1, day=1)
        self.assertEqual(datetime1, datetime2)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
