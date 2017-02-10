import datetime
from dateutil.relativedelta import relativedelta
from calendar import monthrange, isleap


class MpasRelativeDelta(relativedelta):
    """
    MpasRelativeDelta is a subclass of dateutil.relativedelta for relative time
    intervals with different MPAS calendars.

    Only relative intervals (years, months, etc.) are supported and not the
    absolute date specifications (year, month, etc.).  Addition/subtraction
    of datetime.datetime objects or other MpasRelativeDelta (but currently not
    datetime.date, datetime.timedelta or other related objects) is supported.

    Author
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/09/2017
    """

    def __init__(self, dt1=None, dt2=None, years=0, months=0, days=0,
                 hours=0, minutes=0, seconds=0, calendar='gregorian'):
        if calendar not in ['gregorian', 'gregorian_noleap']:
            raise ValueError('Unsupported MPAs calendar {}'.format(calendar))
        self.calendar = calendar
        super(MpasRelativeDelta, self).__init__(dt1=dt1, dt2=dt2, years=years,
                                                months=months, days=days,
                                                hours=hours, minutes=minutes,
                                                seconds=seconds)

    def __add__(self, other):
        if not isinstance(other, (datetime.datetime, MpasRelativeDelta)):
            return NotImplemented

        if isinstance(other, MpasRelativeDelta):
            if self.calendar != other.calendar:
                raise ValueError('MpasRelativeDelta objects can only be added '
                                 'if their calendars match.')
            years = self.years + other.years
            months = self.months + other.months
            if months > 12:
                years += 1
                months -= 12
            elif months < 1:
                years -= 1
                months += 12

            return self.__class__(years=years,
                                  months=months,
                                  days=self.days + other.days,
                                  hours=self.hours + other.hours,
                                  minutes=self.minutes + other.minutes,
                                  seconds=self.seconds + other.seconds,
                                  calendar=self.calendar)

        year = other.year+self.years

        month = other.month
        if self.months != 0:
            assert 1 <= abs(self.months) <= 12
            month += self.months
            if month > 12:
                year += 1
                month -= 12
            elif month < 1:
                year -= 1
                month += 12

        if self.calendar == 'gregorian':
            daysInMonth = monthrange(year, month)[1]
        elif self.calendar == 'gregorian_noleap':
            # use year 0001, which is not a leapyear
            daysInMonth = monthrange(1, month)[1]

        day = min(daysInMonth, other.day)
        repl = {"year": year, "month": month, "day": day}

        days = self.days
        if self.calendar == 'gregorian_noleap' and isleap(year):
            if month == 2 and day+days >= 29:
                # skip forward over the leap day
                days += 1
            elif month == 3 and day+days <= 0:
                # skip backward over the leap day
                days -= 1

        return (other.replace(**repl) +
                datetime.timedelta(days=days,
                                   hours=self.hours,
                                   minutes=self.minutes,
                                   seconds=self.seconds))

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __sub__(self, other):
        if not isinstance(other, MpasRelativeDelta):
            return NotImplemented
        return self.__add__(other.__neg__())

    def __neg__(self):
        return self.__class__(years=-self.years,
                              months=-self.months,
                              days=-self.days,
                              hours=-self.hours,
                              minutes=-self.minutes,
                              seconds=-self.seconds,
                              calendar=self.calendar)

    def __mul__(self, other):
        try:
            f = float(other)
        except TypeError:
            return NotImplemented

        return self.__class__(years=int(self.years * f),
                              months=int(self.months * f),
                              days=int(self.days * f),
                              hours=int(self.hours * f),
                              minutes=int(self.minutes * f),
                              seconds=int(self.seconds * f),
                              calendar=self.calendar)

    __rmul__ = __mul__

    def __div__(self, other):
        try:
            reciprocal = 1 / float(other)
        except TypeError:
            return NotImplemented

        return self.__mul__(reciprocal)

    __truediv__ = __div__

    def __repr__(self):
        outList = []
        for attr in ["years", "months", "days", "leapdays",
                     "hours", "minutes", "seconds", "microseconds"]:
            value = getattr(self, attr)
            if value:
                outList.append("{attr}={value:+g}".format(attr=attr,
                                                          value=value))
        outList.append("calendar='{}'".format(self.calendar))
        return "{classname}({attrs})".format(classname=self.__class__.__name__,
                                             attrs=", ".join(outList))

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
