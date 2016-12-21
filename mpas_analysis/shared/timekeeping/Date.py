"""
    Module for the Date class used to parse and compare dates and times

    Xylar Asay-Davis
    Last modified: 11/02/2016
"""

import functools
import numpy
import datetime

@functools.total_ordering
class Date(object):
    """
    Class for representing dates on a 365-day calendar.
    Date objects can be created either from a formatted string or
    from a number of seconds (mostly intended for internal use).
    Date objects can be added to or subtracted from one another and
    can be compared with one another.
    """

    # constructor
    def __init__(self, dateString=None, isInterval=False, totalSeconds=None,
                 years=None, months=None, days=None,
                 hours=None, minutes=None, seconds=None):
        """
        creates a new Date object.  If the dateString is supplied, it should
        have one of the following formats:
        YYYY-MM-DD_hh:mm:ss
        YYYY-MM-DD_hh.mm.ss
        YYYY-MM-DD_SSSSS
        DDD_hh:mm:ss
        DDD_hh.mm.ss
        DDD_SSSSS
        hh.mm.ss
        hh:mm:ss
        YYYY-MM-DD
        SSSSS

        isInterval indicates whether the date is an interval (difference
        between dates) or a normal (non-interval) date.  Intervals mean that
        the month and day start with 0, while strings representing non-interval
        dates have day and months starting with 1.

        If a dateString is not supplied, totalSeconds can be used to supply
        the date as a number of seconds (as a 64-bit integer).

        If neither dateString nor totalSeconds is given, all of years, months,
        days, hours, minutes and seconds are required to represent the date.
        These argument are intended mostly for internal use.
        """

        self.isInterval = isInterval
        if dateString is not None:
            self._parseDate(dateString)
        elif totalSeconds is not None:
            self._secondsToDate(totalSeconds)
        else:
            if years is None:
                raise ValueError('years must be set')
            self.years = numpy.int64(years)
            if months is None:
                raise ValueError('months must be set')
            self.months = numpy.int64(months)
            if days is None:
                raise ValueError('days must be set')
            self.days = numpy.int64(days)
            if hours is None:
                raise ValueError('hours must be set')
            self.hours = numpy.int64(hours)
            if minutes is None:
                raise ValueError('minutes must be set')
            self.minutes = numpy.int64(minutes)
            if seconds is None:
                raise ValueError('seconds must be set')
            self.seconds = numpy.int64(seconds)
            self._setTotalSeconds()

    def to_datetime(self, yearOffset=0):
        """
        Converts the date object to a datetime object.
        The yearOffset is added to this date's year, and
        the resulting date is clamped to the range supported by
        numpy's datetime64[ns], used internally by xarray an
        pandas

        Last modified: 11/28/2016
        Author: Xylar Asay-Davis
        """
        if self.isInterval:
            raise ValueError("self.isInterval == True. Use to_timedelta "
                             "instead of to_datetime")

        year = numpy.maximum(datetime.MINYEAR,
                             numpy.minimum(datetime.MAXYEAR,
                                           self.years+yearOffset))
        outDate =  datetime.datetime(year=year, month=self.months+1,
                                     day=self.days+1, hour=self.hours,
                                     minute=self.minutes, second=self.seconds)

        minDate = datetime.datetime(year=1678, month=1, day=1,
                                    hour=0, minute=0, second=0)
        maxDate = datetime.datetime(year=2262, month=1, day=1,
                                    hour=0, minute=0, second=0)
        outDate = max(minDate, min(maxDate, outDate))
        return outDate

    def to_timedelta(self):
        """
        Converts the date object to a timedelta object

        Last modified: 11/28/2016
        Author: Xylar Asay-Davis
        """
        if not self.isInterval:
            raise ValueError("self.isInterval == False. Use to_datetime "
                             "instead of to_timedelta")

        days = 365*self.years + self._monthsToDays(self.months) + self.days
        return datetime.timedelta(days=self.days, hours=self.hours,
                                  minutes=self.minutes, seconds=self.seconds)

    def __lt__(self, other):
        if self.isInterval != other.isInterval:
            raise ValueError('Comparing interval with non-interval Date '
                             'object')
        return self.totalSeconds < other.totalSeconds

    def __eq__(self, other):
        if self.isInterval != other.isInterval:
            raise ValueError('Comparing interval with non-interval Date '
                             'object')
        return self.totalSeconds == other.totalSeconds

    def __add__(self, other):
        if self.isInterval:
            raise ValueError('Attempting to add to an interval Date object')
        if not other.isInterval:
            raise ValueError('Attempting to add a non-interval Date object')

        seconds = self.seconds + other.seconds
        minutes = self.minutes + other.minutes + seconds/60
        seconds %= 60
        hours = self.hours + other.hours + minutes/60
        minutes %= 60
        months = self.months + other.months
        years = self.years + other.years + months/12
        months %= 12
        days = (self._monthsToDays(months) + self.days + other.days + hours/24)
        years += days/365
        days %= 365
        (months, days) = self._daysToMonthsAndDays(days)
        return Date(isInterval=False, years=years, months=months, days=days,
                    hours=hours, minutes=minutes, seconds=seconds)

    def __sub__(self, other):
        if self.isInterval:
            raise ValueError('Attempting to subtract from an interval Date '
                             'object')

        isInterval = not other.isInterval
        seconds = self.seconds - other.seconds
        minutes = self.minutes - other.minutes + seconds/60
        seconds %= 60
        hours = self.hours - other.hours + minutes/60
        minutes %= 60
        months = self.months - other.months
        years = self.years - other.years + months/12
        months %= 12
        days = (self._monthsToDays(months) + self.days - other.days + hours/24)
        years += days/365
        days %= 365
        (months, days) = self._daysToMonthsAndDays(days)

        return Date(isInterval=isInterval, years=years, months=months,
                    days=days, hours=hours, minutes=minutes, seconds=seconds)

    def __str__(self):
        if self.isInterval:
            offset = 0
        else:
            offset = 1
        return '{:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}'.format(
            self.years, self.months+offset, self.days+offset,
            self.hours, self.minutes, self.seconds)

    def _diffSeconds(self, other):
        return

    def _setTotalSeconds(self):
        days = self.years*365 + self._monthsToDays(self.months) + self.days
        self.totalSeconds = (((days*24 + self.hours)*60 + self.minutes)*60 +
                             self.seconds)

    def _monthsToDays(self, months):
        daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        days = numpy.int64(0)
        for month in range(months):
            days += daysInMonth[month]
        return days

    def _daysToMonthsAndDays(self, days):
        daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        assert(days < 365)
        months = numpy.int64(0)
        while days > daysInMonth[months]:
            days -= daysInMonth[months]
            months += 1
        days = numpy.int64(days)
        return (months, days)

    def _secondsToDate(self, seconds):
        self.totalSeconds = seconds
        self.years = numpy.int64(seconds / 31536000)
        seconds %= 31536000
        days = numpy.int64(seconds / 86400)
        (self.months, self.days) = self._daysToMonthsAndDays(days)
        seconds %= 86400
        self.hours = numpy.int64(seconds / 3600)
        seconds %= 3600
        self.minutes = numpy.int64(seconds / 60)
        seconds %= 60
        self.seconds = seconds

    def _parseDate(self, dateString):
        """
        parses a dateString in one of the following formats into
        a Date object:
        YYYY-MM-DD_hh:mm:ss
        YYYY-MM-DD_hh.mm.ss
        YYYY-MM-DD_SSSSS
        DDD_hh:mm:ss
        DDD_hh.mm.ss
        DDD_SSSSS
        hh.mm.ss
        hh:mm:ss
        YYYY-MM-DD
        YYYY-MM
        SSSSS
        """
        if self.isInterval:
            offset = numpy.int64(0)
        else:
            offset = numpy.int64(1)

        if '_' in dateString:
            ymd, hms = dateString.split('_')
        else:
            if '-' in dateString:
                ymd = dateString
                # error can result if dateString = '1990-01'
                # assume this means '1990-01-01'
                if len(ymd.split('-')) == 2:
                    ymd += '-01'
                hms = '00:00:00'
            else:
                if self.isInterval:
                    ymd = '0000-00-00'
                else:
                    ymd = '0000-01-01'
                hms = dateString

        if '.' in hms:
            hms = hms.replace('.', ':')

        if '-' in ymd:
            (self.years, self.months, self.days) \
                = [numpy.int64(sub) for sub in ymd.split('-')]
            self.months -= offset
            self.days -= offset
        else:
            self.days = numpy.int64(ymd) - offset
            self.years = numpy.int64(0)
            self.months = numpy.int64(0)

        if ':' in hms:
            (self.hours, self.minutes, self.seconds) \
                = [numpy.int64(sub) for sub in hms.split(':')]
        else:
            self.seconds = numpy.int64(hms)
            self.minutes = numpy.int64(0)
            self.hours = numpy.int64(0)
        self._setTotalSeconds()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
