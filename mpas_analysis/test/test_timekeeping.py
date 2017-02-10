"""
Unit test infrastructure for the Date class

Author
------
Xylar Asay-Davis

Last Modified
-------------
02/09/2017
"""

import pytest
import datetime
from mpas_analysis.shared.timekeeping.MpasRelativeDelta \
    import MpasRelativeDelta
from mpas_analysis.test import TestCase
from mpas_analysis.shared.timekeeping.utility import stringToDatetime, \
    stringToRelativeDelta, clampToNumpyDatetime64


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
            date1 = stringToDatetime('0001-01-01_00:00:00')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=0)
            self.assertEqual(date1, date2)

            delta1 = stringToRelativeDelta('0001-00-00_00:00:00',
                                           calendar=calendar)
            delta2 = MpasRelativeDelta(years=1, months=0, days=0, hours=0,
                                       minutes=0, seconds=0, calendar=calendar)
            self.assertEqual(delta1, delta2)

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
            delta1 = stringToRelativeDelta('0001_00:00:01',
                                           calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=1, hours=0,
                                       minutes=0, seconds=1, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # DDD_hh.mm.ss
            delta1 = stringToRelativeDelta('0002_01.00.01',
                                           calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=2, hours=1,
                                       minutes=0, seconds=1, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # DDD_SSSSS
            delta1 = stringToRelativeDelta('0002_00003',
                                           calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=2, hours=0,
                                       minutes=0, seconds=3, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # hh:mm:ss
            date1 = stringToDatetime('00:00:01')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=1)
            self.assertEqual(date1, date2)

            # hh.mm.ss
            delta1 = stringToRelativeDelta('00.00.01',
                                           calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=0, hours=0,
                                       minutes=0, seconds=1, calendar=calendar)
            self.assertEqual(delta1, delta2)

            # YYYY-MM-DD
            date1 = stringToDatetime('0001-01-01')
            date2 = datetime.datetime(year=1, month=1, day=1, hour=0, minute=0,
                                      second=0)
            self.assertEqual(date1, date2)

            # SSSSS
            delta1 = stringToRelativeDelta('00005',
                                           calendar=calendar)
            delta2 = MpasRelativeDelta(years=0, months=0, days=0, hours=0,
                                       minutes=0, seconds=5, calendar=calendar)
            self.assertEqual(delta1, delta2)

            date1 = stringToDatetime('1996-01-15')
            delta = stringToRelativeDelta('0005-00-00',
                                          calendar=calendar)
            date2 = date1-delta
            self.assertEqual(date2, stringToDatetime('1991-01-15'))

            date1 = stringToDatetime('1996-01-15')
            delta = stringToRelativeDelta('0000-02-00',
                                          calendar=calendar)
            date2 = date1-delta
            self.assertEqual(date2, stringToDatetime('1995-11-15'))

            date1 = stringToDatetime('1996-01-15')
            delta = stringToRelativeDelta('0000-00-20',
                                          calendar=calendar)
            date2 = date1-delta
            self.assertEqual(date2, stringToDatetime('1995-12-26'))

        # since pandas and xarray use the numpy type 'datetime[ns]`, which
        # has a limited range of dates, the date 0001-01-01 gets increased
        # to the minimum allowed year boundary, 1678-01-01 to avoid invalid
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
        # has a limited range of dates, the date 9999-01-01 gets decreased
        # to the maximum allowed year boundary, 2262-01-01 to avoid invalid
        # dates.
        date1 = clampToNumpyDatetime64(stringToDatetime('9999-01-01'),
                                       yearOffset=0)
        date2 = datetime.datetime(year=2262, month=1, day=1)
        self.assertEqual(date1, date2)

    def test_MpasRelativeDeltaOps(self):
        # test if the calendars behave as they should close to leap day
        # also, test addition and subtraction of the form
        # datetime.datetime +/- MpasRelativeDelta above
        # both calendars with adding one day
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-03-01']):
            self.assertEqual(stringToDatetime('2016-02-28') +
                             stringToRelativeDelta('0000-00-01',
                                                   calendar=calendar),
                             stringToDatetime(expected))

        # both calendars with subtracting one day
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-02-28']):
            self.assertEqual(stringToDatetime('2016-03-01') -
                             stringToRelativeDelta('0000-00-01',
                                                   calendar=calendar),
                             stringToDatetime(expected))

        # both calendars with adding one month
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-02-28']):
            self.assertEqual(stringToDatetime('2016-01-31') +
                             stringToRelativeDelta('0000-01-00',
                                                   calendar=calendar),
                             stringToDatetime(expected))

        # both calendars with subtracting one month
        for calendar, expected in zip(['gregorian', 'gregorian_noleap'],
                                      ['2016-02-29', '2016-02-28']):
            self.assertEqual(stringToDatetime('2016-03-31') -
                             stringToRelativeDelta('0000-01-00',
                                                   calendar=calendar),
                             stringToDatetime(expected))

        for calendar in ['gregorian', 'gregorian_noleap']:

            delta1 = stringToRelativeDelta('0000-01-00',  calendar=calendar)
            delta2 = stringToRelativeDelta('0000-00-01',  calendar=calendar)
            deltaSum = stringToRelativeDelta('0000-01-01',  calendar=calendar)
            # test MpasRelativeDelta + MpasRelativeDelta
            self.assertEqual(delta1 + delta2, deltaSum)
            # test MpasRelativeDelta - MpasRelativeDelta
            self.assertEqual(deltaSum - delta2, delta1)

            # test MpasRelativeDelta(date1, date2)
            date1 = stringToDatetime('0002-02-02')
            date2 = stringToDatetime('0001-01-01')
            delta = stringToRelativeDelta('0001-01-01',  calendar=calendar)
            self.assertEqual(MpasRelativeDelta(dt1=date1, dt2=date2,
                                               calendar=calendar),
                             delta)

            # test MpasRelativeDelta + datetime.datetime (an odd order but
            # it's allowed...)
            date1 = stringToDatetime('0001-01-01')
            delta = stringToRelativeDelta('0001-01-01',  calendar=calendar)
            date2 = stringToDatetime('0002-02-02')
            self.assertEqual(delta + date1, date2)

            # test multiplication/division by scalars
            delta1 = stringToRelativeDelta('0001-01-01',  calendar=calendar)
            delta2 = stringToRelativeDelta('0002-02-02',  calendar=calendar)
            self.assertEqual(2*delta1, delta2)
            self.assertEqual(delta2/2, delta1)


        # make sure there's an error when we try to add MpasRelativeDeltas
        # with different calendars
        with self.assertRaisesRegexp(ValueError,
                                     'MpasRelativeDelta objects can only be '
                                     'added if their calendars match.'):
            delta1 = stringToRelativeDelta('0000-01-00',
                                           calendar='gregorian')
            delta2 = stringToRelativeDelta('0000-00-01',
                                           calendar='gregorian_noleap')
            deltaSum = delta1 + delta2


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
