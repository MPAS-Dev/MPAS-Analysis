"""
Time keeping utility functions

Author
------
Xylar Asay-Davis

Last Modified
-------------
02/11/2017
"""

import datetime
import netCDF4
import numpy

from .MpasRelativeDelta import MpasRelativeDelta


def get_simulation_start_time(streams):
    """
    Given a StreamsFile object, returns the simulation start time parsed from
    a restart file.

    Parameters
    ----------
    steams : StreamsFile object
        For parsing an MPAS streams file

    Returns
    -------
    simulation_start_time : string
        The start date of the simulation parsed from a restart file identified
        by the contents of `streams`.

    Raises
    ------
    IOError
        If no restart file can be found.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/11/2017
    """

    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS restart file found: need at least one '
                      'restart file for analysis to work correctly')

    ncFile = netCDF4.Dataset(restartFile, mode='r')
    simulationStartTime = ncFile.variables['simulationStartTime'][:]
    # convert from character array to str
    simulationStartTime = ''.join(simulationStartTime).strip()
    # replace underscores so it works as a CF-compliant reference date
    simulationStartTime = simulationStartTime.replace('_', ' ')
    ncFile.close()

    return simulationStartTime


def string_to_datetime(dateString):  # {{{
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
        _parse_date_string(dateString, isInterval=False)

    return datetime.datetime(year=year, month=month, day=day, hour=hour,
                             minute=minute, second=second)  # }}}


def string_to_relative_delta(dateString, calendar='gregorian'):  # {{{
    """
    Given a date string and a calendar, returns an instance of
    `MpasRelativeDelta`

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
    relativedelta : An `MpasRelativeDelta` object

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
        _parse_date_string(dateString, isInterval=True)

    return MpasRelativeDelta(years=years, months=months, days=days,
                             hours=hours, minutes=minutes, seconds=seconds,
                             calendar=calendar)
    # }}}


def string_to_days_since_date(dateString, calendar='gregorian',
                              referenceDate='0001-01-01'):
    """
    Given a date string or an array-like of date strings, a reference date
    string, and a calendar, returns the number of days (as a float or
    numpy.array of floats) since the reference date

    Parameters
    ----------
    dateStrings : str or array-like of str
        A date and time (or array of date/times) in one of the following
        formats:
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

    referenceDate : str, optional
        A reference date of the form:
            - 0001-01-01
            - 0001-01-01 00:00:00

    Returns
    -------
    days : float or numpy.array of floats
        The number of days since `referenceDate` for each date in dateString

    Raises
    ------
    ValueError
        If an invalid `dateString` or `calendar` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/04/2017
    """

    isSingleString = isinstance(dateString, str)

    if isSingleString:
        dateString = [dateString]

    dates = [string_to_datetime(string) for string in dateString]
    days = datetime_to_days(dates, calendar=calendar,
                            referenceDate=referenceDate)

    if isSingleString:
        days = days[0]
    else:
        days = numpy.array(days)
    return days


def days_to_datetime(days, calendar='gregorian', referenceDate='0001-01-01'):
    """
    Covert days to `datetime.datetime` objects given a reference date and an
    MPAS calendar (either 'gregorian' or 'gregorian_noleap').

    Parameters
    ----------
    days : float or array-like of floats
        The number of days since the reference date.

    calendar : {'gregorian', 'gregorian_noleap'}, optinal
        A calendar to be used to convert days to a `datetime.datetime` object.

    referenceDate : str, optional
        A reference date of the form:
            - 0001-01-01
            - 0001-01-01 00:00:00

    Returns
    -------
    datetime : An instance of `datetime.datetime` (or array-like of datetimes)
        The days since `referenceDate` on the given `calendar`.

    Raises
    ------
    ValueError
        If an invalid `days`, `referenceDate` or `calendar` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/04/2017
    """

    datetimes = netCDF4.num2date(days,
                                 'days since {}'.format(referenceDate),
                                 calendar=_mpas_to_netcdf_calendar(calendar))

    # convert to datetime.datetime
    if isinstance(datetimes, numpy.ndarray):
        newDateTimes = []
        for date in datetimes.flat:
            newDateTimes.append(_round_datetime(date))
        if len(newDateTimes) > 0:
            datetimes = numpy.reshape(numpy.array(newDateTimes),
                                      datetimes.shape)

    else:
        datetimes = _round_datetime(datetimes)

    return datetimes


def datetime_to_days(dates, calendar='gregorian', referenceDate='0001-01-01'):
    """
    Given date(s), a calendar and a reference date, returns the days since
    the reference date, either as a single float or an array of floats.

    Parameters
    ----------
    datetime : instance or array-like of datetime.datetime
        The date(s) to be converted to days since `referenceDate` on the
        given `calendar`.

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        A calendar to be used to convert days to a `datetime.datetime` object.

    referenceDate : str, optional
        A reference date of the form:
            - 0001-01-01
            - 0001-01-01 00:00:00

    Returns
    -------
    days : float or array of floats
        The days since `referenceDate` on the given `calendar`.

    Raises
    ------
    ValueError
        If an invalid `datetimes`, `referenceDate` or `calendar` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/11/2017
    """

    isSingleDate = False
    if isinstance(dates, datetime.datetime):
        dates = [dates]
        isSingleDate = True

    days = netCDF4.date2num(dates, 'days since {}'.format(referenceDate),
                            calendar=_mpas_to_netcdf_calendar(calendar))

    if isSingleDate:
        days = days[0]

    return days


def date_to_days(year=1, month=1, day=1, hour=0, minute=0, second=0,
                 calendar='gregorian', referenceDate='0001-01-01'):
    """
    Given a date in the form of year, month, day, etc.; a calendar; and a
    reference date, returns the days since the reference date.

    Parameters
    ----------
    year, month, day, hour, minute, second : int, optional
        The date to be converted to days since `referenceDate` on the
        given `calendar`.

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        A calendar to be used to convert days to a `datetime.datetime` object.

    referenceDate : str, optional
        A reference date of the form:
            - 0001-01-01
            - 0001-01-01 00:00:00

    Returns
    -------
    days : float
        The days since `referenceDate` on the given `calendar`.

    Raises
    ------
    ValueError
        If an invalid `referenceDate` or `calendar` is supplied.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/11/2017
    """

    calendar = _mpas_to_netcdf_calendar(calendar)

    date = datetime.datetime(year, month, day, hour, minute, second)

    return netCDF4.date2num(date, 'days since {}'.format(referenceDate),
                            calendar=calendar)


def _parse_date_string(dateString, isInterval=False):  # {{{
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


def _mpas_to_netcdf_calendar(calendar):
    """
    Convert from MPAS calendar to NetCDF4 calendar names.
    """

    if calendar == 'gregorian_noleap':
        calendar = 'noleap'
    elif calendar != 'gregorian':
        raise ValueError('Unsupported calendar {}'.format(calendar))
    return calendar


def _round_datetime(date):
    """Round a datetime object to nearest second
    date : datetime.datetime or similar objet object.
    """
    (year, month, day, hour, minute, second, microsecond) = \
        (date.year, date.month, date.day, date.hour, date.minute, date.second,
         date.microsecond)

    date = datetime.datetime(year=year, month=month, day=day,
                             hour=hour, minute=minute,
                             second=second)

    add_seconds = int(1e-6*microsecond+0.5)

    return date + datetime.timedelta(0, add_seconds)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
