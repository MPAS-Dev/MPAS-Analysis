"""
Unit test infrastructure for the Date class

Author
------
Xylar Asay-Davis

Last Modified
-------------
02/17/2017
"""

import pytest
import datetime
from mpas_analysis.shared.timekeeping.MpasRelativeDelta \
    import MpasRelativeDelta
from mpas_analysis.test import TestCase
from mpas_analysis.shared.timekeeping.utility import string_to_datetime, \
    string_to_relative_delta, string_to_days_since_date, days_to_datetime, \
    datetime_to_days, date_to_days


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

        for calendar in ['gregorian', 'gregorian_noleap']:
            # test datetime.datetime
            # YYYY-MM-DD_hh:mm:ss
            date1 = string_to_datetime('0001-01-01_00:00:00')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=0)
            self.assertEqual(date1, date2)

            delta1 = string_to_relative_delta('0001-00-00_00:00:00',
                                              calendar=calendar)
            delta2 = MpasRelativeDelta(years=1, months=0, days=0, hours=0,
                                       minutes=0, seconds=0, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # YYYY-MM-DD_hh.mm.ss
            date1 = string_to_datetime('0001-01-01_00.00.00')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=0)
            self.assertEqual(date1, date2)

            # YYYY-MM-DD_SSSSS
            date1 = string_to_datetime('0001-01-01_00002')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=2)
            self.assertEqual(date1, date2)

            # DDD_hh:mm:ss
            delta1 = string_to_relative_delta('0001_00:00:01',
                                              calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=1, hours=0,
                                       minutes=0, seconds=1, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # DDD_hh.mm.ss
            delta1 = string_to_relative_delta('0002_01.00.01',
                                              calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=2, hours=1,
                                       minutes=0, seconds=1, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # DDD_SSSSS
            delta1 = string_to_relative_delta('0002_00003',
                                              calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=2, hours=0,
                                       minutes=0, seconds=3, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # hh:mm:ss
            date1 = string_to_datetime('00:00:01')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=1)
            self.assertEqual(date1, date2)

            # hh.mm.ss
            delta1 = string_to_relative_delta('00.00.01',
                                              calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=0, hours=0,
                                       minutes=0, seconds=1, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # YYYY-MM-DD
            date1 = string_to_datetime('0001-01-01')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=0)
            self.assertEqual(date1, date2)

            # SSSSS
            delta1 = string_to_relative_delta('00005',
                                              calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=0, hours=0,
                                       minutes=0, seconds=5, calendar=calendar)
            self.assertEqual(delta1, delta2)

            date1 = string_to_datetime('1996-01-15')
            delta = string_to_relative_delta('0005-00-00',
                                             calendar=calendar)
            date2 = date1-delta
            self.assertEqual(date2, string_to_datetime('1991-01-15'))

            date1 = string_to_datetime('1996-01-15')
            delta = string_to_relative_delta('0000-02-00',
                                             calendar=calendar)
            date2 = date1-delta
            self.assertEqual(date2, string_to_datetime('1995-11-15'))

            date1 = string_to_datetime('1996-01-15')
            delta = string_to_relative_delta('0000-00-20',
                                             calendar=calendar)
            date2 = date1-delta
            self.assertEqual(date2, string_to_datetime('1995-12-26'))

    def test_MpasRelativeDeltaOps(self):
        # test if the calendars behave as they should close to leap day
        # also, test addition and subtraction of the form
        # datetime.datetime +/- MpasRelativeDelta above
        # both calendars with adding one day
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-03-01']):
            self.assertEqual(string_to_datetime('2016-02-28') +
                             string_to_relative_delta('0000-00-01',
                                                      calendar=calendar),
                             string_to_datetime(expected))

        # both calendars with subtracting one day
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-02-28']):
            self.assertEqual(string_to_datetime('2016-03-01') -
                             string_to_relative_delta('0000-00-01',
                                                      calendar=calendar),
                             string_to_datetime(expected))

        # both calendars with adding one month
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-02-28']):
            self.assertEqual(string_to_datetime('2016-01-31') +
                             string_to_relative_delta('0000-01-00',
                                                      calendar=calendar),
                             string_to_datetime(expected))

        # both calendars with subtracting one month
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-02-28']):
            self.assertEqual(string_to_datetime('2016-03-31') -
                             string_to_relative_delta('0000-01-00',
                                                      calendar=calendar),
                             string_to_datetime(expected))

        for calendar in ['gregorian', 'gregorian_noleap']:

            delta1 = string_to_relative_delta('0000-01-00',  calendar=calendar)
            delta2 = string_to_relative_delta('0000-00-01',  calendar=calendar)
            deltaSum = string_to_relative_delta('0000-01-01',
                                                calendar=calendar)
            # test MpasRelativeDelta + MpasRelativeDelta
            self.assertEqual(delta1 + delta2, deltaSum)
            # test MpasRelativeDelta - MpasRelativeDelta
            self.assertEqual(deltaSum - delta2, delta1)

            # test MpasRelativeDelta(date1, date2)
            date1 = string_to_datetime('0002-02-02')
            date2 = string_to_datetime('0001-01-01')
            delta = string_to_relative_delta('0001-01-01',  calendar=calendar)
            self.assertEqual(MpasRelativeDelta(dt1=date1, dt2=date2,
                                               calendar=calendar),
                             delta)

            # test MpasRelativeDelta + datetime.datetime (an odd order but
            # it's allowed...)
            date1 = string_to_datetime('0001-01-01')
            delta = string_to_relative_delta('0001-01-01',  calendar=calendar)
            date2 = string_to_datetime('0002-02-02')
            self.assertEqual(delta + date1, date2)

            # test multiplication/division by scalars
            delta1 = string_to_relative_delta('0001-01-01',  calendar=calendar)
            delta2 = string_to_relative_delta('0002-02-02',  calendar=calendar)
            self.assertEqual(2*delta1, delta2)
            self.assertEqual(delta2/2, delta1)

        # make sure there's an error when we try to add MpasRelativeDeltas
        # with different calendars
        with self.assertRaisesRegexp(ValueError,
                                     'MpasRelativeDelta objects can only be '
                                     'added if their calendars match.'):
            delta1 = string_to_relative_delta('0000-01-00',
                                              calendar='gregorian')
            delta2 = string_to_relative_delta('0000-00-01',
                                              calendar='gregorian_noleap')
            deltaSum = delta1 + delta2

    def test_string_to_days_since_date(self):
        referenceDate = '0001-01-01'
        for calendar in ['gregorian', 'gregorian_noleap']:
            for dateString, expected_days in [('0001-01-01', 0.),
                                              ('0001-01-02', 1.),
                                              ('0001-02-01', 31.),
                                              ('0002-01-01', 365.)]:
                days = string_to_days_since_date(dateString=dateString,
                                                 calendar=calendar,
                                                 referenceDate=referenceDate)
                self.assertEqual(days, expected_days)

        referenceDate = '2016-01-01'
        for calendar, expected_days in [('gregorian', 366.),
                                        ('gregorian_noleap', 365.)]:
            days = string_to_days_since_date(dateString='2017-01-01',
                                             calendar=calendar,
                                             referenceDate=referenceDate)
            self.assertEqual(days, expected_days)

    def test_days_to_datetime(self):
        referenceDate = '0001-01-01'
        for calendar in ['gregorian', 'gregorian_noleap']:
            for dateString, days in [('0001-01-01', 0.),
                                     ('0001-01-02', 1.),
                                     ('0001-02-01', 31.),
                                     ('0002-01-01', 365.)]:
                datetime = days_to_datetime(days=days,
                                            calendar=calendar,
                                            referenceDate=referenceDate)
                self.assertEqual(datetime, string_to_datetime(dateString))

        referenceDate = '2016-01-01'
        for calendar, days in [('gregorian', 366.),
                               ('gregorian_noleap', 365.)]:
            datetime = days_to_datetime(days=days,
                                        calendar=calendar,
                                        referenceDate=referenceDate)
            self.assertEqual(datetime, string_to_datetime('2017-01-01'))

    def test_datetime_to_days(self):
        referenceDate = '0001-01-01'
        for calendar in ['gregorian', 'gregorian_noleap']:
            for dateString, expected_days in [('0001-01-01', 0.),
                                              ('0001-01-02', 1.),
                                              ('0001-02-01', 31.),
                                              ('0002-01-01', 365.)]:
                days = datetime_to_days(dates=string_to_datetime(dateString),
                                        calendar=calendar,
                                        referenceDate=referenceDate)
                self.assertEqual(days, expected_days)

        referenceDate = '2016-01-01'
        for calendar, expected_days in [('gregorian', 366.),
                                        ('gregorian_noleap', 365.)]:
            days = datetime_to_days(dates=string_to_datetime('2017-01-01'),
                                    calendar=calendar,
                                    referenceDate=referenceDate)
            self.assertEqual(days, expected_days)

    def test_date_to_days(self):
        referenceDate = '0001-01-01'
        for calendar in ['gregorian', 'gregorian_noleap']:
            days = date_to_days(year=1, month=1, day=1, calendar=calendar,
                                referenceDate=referenceDate)
            self.assertEqual(days, 0.)
            days = date_to_days(year=1, month=1, day=2, calendar=calendar,
                                referenceDate=referenceDate)
            self.assertEqual(days, 1.)
            days = date_to_days(year=1, month=2, day=1, calendar=calendar,
                                referenceDate=referenceDate)
            self.assertEqual(days, 31.)
            days = date_to_days(year=2, month=1, day=1, calendar=calendar,
                                referenceDate=referenceDate)
            self.assertEqual(days, 365.)

        referenceDate = '2016-01-01'
        for calendar, expected_days in [('gregorian', 366.),
                                        ('gregorian_noleap', 365.)]:
            days = date_to_days(year=2017, month=1, day=1,
                                calendar=calendar,
                                referenceDate=referenceDate)
            self.assertEqual(days, expected_days)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
