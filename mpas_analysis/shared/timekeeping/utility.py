"""
Time keeping utility functions

Author
------
Xylar Asay-Davis

Last Modified
-------------
02/06/2017
"""

import datetime
from dateutil.relativedelta import relativedelta


def stringToDatetime(dateString):  # {{{
    """
    Given a date string and a calendar, returns a `datetime.datetime`

    Parameters
    ----------
    dateString : string
        A date and time in one of the following formats:
        - YYYY-MM-DD hh:mm:ss
        - YYYY-MM-DD hh.mm.ss
        - YYYY-MM-DD SSSSS
        - DDD hh:mm:ss
        - DDD hh.mm.ss
        - DDD SSSSS
        - hh.mm.ss
        - hh:mm:ss
        - YYYY-MM-DD
        - YYYY-MM
        - SSSSS

        Note: either underscores or spaces can be used to separate the date
        from the time portion of the string.

    Returns
    -------
    datetime : A `datetime.datetime` object

    Raises
    ------
    ValueError
        If an invalid `dateString` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/04/2017
    """

    (year, month, day, hour, minute, second) = \
        _parseDateString(dateString, isInterval=False)

    return datetime.datetime(year=year, month=month, day=day, hour=hour,
                             minute=minute, second=second)  # }}}


def stringToRelativedelta(dateString, calendar='gregorian'):  # {{{
    """
    Given a date string and a calendar, returns an instance of
    `dateutil.relativedelta.relativedelta`

    Parameters
    ----------
    dateString : string
        A date and time in one of the following formats:
        - YYYY-MM-DD hh:mm:ss
        - YYYY-MM-DD hh.mm.ss
        - YYYY-MM-DD SSSSS
        - DDD hh:mm:ss
        - DDD hh.mm.ss
        - DDD SSSSS
        - hh.mm.ss
        - hh:mm:ss
        - YYYY-MM-DD
        - YYYY-MM
        - SSSSS

        Note: either underscores or spaces can be used to separate the date
        from the time portion of the string.

    calendar: {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    Returns
    -------
    relativedelta : A `dateutil.relativedelta.relativedelta` object

    Raises
    ------
    ValueError
        If an invalid `dateString` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/04/2017
    """

    (years, months, days, hours, minutes, seconds) = \
        _parseDateString(dateString, isInterval=True)

    if calendar == 'gregorian':
        leapdays = True
    elif calendar == 'gregorian_noleap':
        leapdays = False
    else:
        raise ValueError('Unsupported calendar {}'.format(calendar))

    return relativedelta(years=years, months=months, days=days, hours=hours,
                         minutes=minutes, seconds=seconds, leapdays=leapdays)
    # }}}


def clampToNumpyDatetime64(date, yearOffset):
    """
    Temporary function for adding an offset year and clamping a datetime to
    range supported by `numpy.datetime64`.
    """

    year = date.year + yearOffset
    if year < 1678:
        return datetime.datetime(year=1678, month=1, day=1, hour=0,
                                 minute=0, second=0)

    if year >= 2262:
        return datetime.datetime(year=2262, month=1, day=1, hour=0,
                                 minute=0, second=0)

    return datetime.datetime(year, date.month, date.day,
                             date.hour, date.minute, date.second)


def _parseDateString(dateString, isInterval=False):  # {{{
    """
    Given a string containing a date, returns a tuple defining a date of the
    form (year, month, day, hour, minute, second) appropriate for constructing
    a datetime or timedelta

    Parameters
    ----------
    dateString : string
        A date and time in one of the followingformats:
        - YYYY-MM-DD hh:mm:ss
        - YYYY-MM-DD hh.mm.ss
        - YYYY-MM-DD SSSSS
        - DDD hh:mm:ss
        - DDD hh.mm.ss
        - DDD SSSSS
        - hh.mm.ss
        - hh:mm:ss
        - YYYY-MM-DD
        - YYYY-MM
        - SSSSS

        Note: either underscores or spaces can be used to separate the date
        from the time portion of the string.

    isInterval : bool, optional
        If ``isInterval=True``, the result is appropriate for constructing
        a `datetime.timedelta` object rather than a `datetime`.

    Returns
    -------
    date : A tuple of (year, month, day, hour, minute, second)

    Raises
    ------
    ValueError
        If an invalid `dateString` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/04/2017
    """
    if isInterval:
        offset = 0
    else:
        offset = 1

    # change underscores to spaces so both can be supported
    dateString = dateString.replace('_', ' ')
    if ' ' in dateString:
        ymd, hms = dateString.split(' ')
    else:
        if '-' in dateString:
            ymd = dateString
            # error can result if dateString = '1990-01'
            # assume this means '1990-01-01'
            if len(ymd.split('-')) == 2:
                ymd += '-01'
            hms = '00:00:00'
        else:
            if isInterval:
                ymd = '0000-00-00'
            else:
                ymd = '0001-01-01'
            hms = dateString

    if '.' in hms:
        hms = hms.replace('.', ':')

    if '-' in ymd:
        (year, month, day) \
            = [int(sub) for sub in ymd.split('-')]
    else:
        day = int(ymd)
        year = 0
        month = offset

    if ':' in hms:
        (hour, minute, second) \
            = [int(sub) for sub in hms.split(':')]
    else:
        second = int(hms)
        minute = 0
        hour = 0
    return (year, month, day, hour, minute, second)  # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
