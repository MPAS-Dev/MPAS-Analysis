# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
"""
Time keeping utility functions
"""
# Authors
# -------
# Xylar Asay-Davis

import datetime
import netCDF4
import xarray
import numpy

from mpas_analysis.shared.timekeeping.MpasRelativeDelta import \
    MpasRelativeDelta
from mpas_analysis.shared.io.utility import decode_strings


def get_simulation_start_time(streams):
    """
    Given a ``StreamsFile`` object, returns the simulation start time parsed
    from a restart file.

    Parameters
    ----------
    steams : ``StreamsFile`` object
        For parsing an MPAS streams file

    Returns
    -------
    simulation_start_time : str
        The start date of the simulation parsed from a restart file identified
        by the contents of ``streams``.

    Raises
    ------
    IOError
        If no restart file can be found.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS restart file found: need at least one '
                      'restart file for analysis to work correctly')

    ds = xarray.open_dataset(restartFile)
    da = ds.simulationStartTime
    if da.dtype.type is numpy.bytes_:
        simulationStartTime = bytes.decode(da.values.tobytes())
    else:
        simulationStartTime = da.values.tobytes()
    # replace underscores so it works as a CF-compliant reference date
    simulationStartTime = simulationStartTime.rstrip('\x00').replace('_', ' ')

    return simulationStartTime


def string_to_datetime(dateString):
    """
    Given a date string and a calendar, returns a ``datetime.datetime``

    Parameters
    ----------
    dateString : string
        A date and time in one of the following formats::

            YYYY-MM-DD hh:mm:ss
            YYYY-MM-DD hh.mm.ss
            YYYY-MM-DD SSSSS
            DDD hh:mm:ss
            DDD hh.mm.ss
            DDD SSSSS
            hh.mm.ss
            hh:mm:ss
            YYYY-MM-DD
            YYYY-MM
            SSSSS

        Note: either underscores or spaces can be used to separate the date
        from the time portion of the string.

    Returns
    -------
    datetime : A ``datetime.datetime`` object

    Raises
    ------
    ValueError
        If an invalid ``dateString`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    (year, month, day, hour, minute, second) = \
        _parse_date_string(dateString, isInterval=False)

    if year == 0:
        year = 1
    return datetime.datetime(year=year, month=month, day=day, hour=hour,
                             minute=minute, second=second)


def string_to_relative_delta(dateString, calendar='gregorian'):
    """
    Given a date string and a calendar, returns an instance of
    ``MpasRelativeDelta``

    Parameters
    ----------
    dateString : str
        A date and time in one of the following formats::

            YYYY-MM-DD hh:mm:ss
            YYYY-MM-DD hh.mm.ss
            YYYY-MM-DD SSSSS
            DDD hh:mm:ss
            DDD hh.mm.ss
            DDD SSSSS
            hh.mm.ss
            hh:mm:ss
            YYYY-MM-DD
            YYYY-MM
            SSSSS

        Note: either underscores or spaces can be used to separate the date
        from the time portion of the string.

    calendar: {'gregorian', 'noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    Returns
    -------
    relativedelta : An ``MpasRelativeDelta`` object

    Raises
    ------
    ValueError
        If an invalid ``dateString`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    (years, months, days, hours, minutes, seconds) = \
        _parse_date_string(dateString, isInterval=True)

    return MpasRelativeDelta(years=years, months=months, days=days,
                             hours=hours, minutes=minutes, seconds=seconds,
                             calendar=calendar)


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
        formats::

            YYYY-MM-DD hh:mm:ss
            YYYY-MM-DD hh.mm.ss
            YYYY-MM-DD SSSSS
            DDD hh:mm:ss
            DDD hh.mm.ss
            DDD SSSSS
            hh.mm.ss
            hh:mm:ss
            YYYY-MM-DD
            YYYY-MM
            SSSSS

        Note: either underscores or spaces can be used to separate the date
        from the time portion of the string.

    calendar: {'gregorian', 'noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    referenceDate : str, optional
        A reference date of the form::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    days : float or numpy.array of floats
        The number of days since ``referenceDate`` for each date in
        ``dateString``

    Raises
    ------
    ValueError
        If an invalid ``dateString`` or ``calendar`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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
    Covert days to ``datetime.datetime`` objects given a reference date and an
    MPAS calendar (either 'gregorian' or 'noleap').

    Parameters
    ----------
    days : float or array-like of floats
        The number of days since the reference date.

    calendar : {'gregorian', 'noleap'}, optional
        A calendar to be used to convert days to a ``datetime.datetime``
        object.

    referenceDate : str, optional
        A reference date of the form::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    datetime : `datetime.datetime` (or array-like of datetimes)
        The days since ``referenceDate`` on the given ``calendar``.

    Raises
    ------
    ValueError
        If an invalid ``days``, ``referenceDate`` or ``calendar`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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
        The date(s) to be converted to days since ``referenceDate`` on the
        given ``calendar``.

    calendar : {'gregorian', 'noleap'}, optional
        A calendar to be used to convert days to a ``datetime.datetime`` object.

    referenceDate : str, optional
        A reference date of the form::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    days : float or array of floats
        The days since ``referenceDate`` on the given ``calendar``.

    Raises
    ------
    ValueError
        If an invalid ``datetimes``, ``referenceDate`` or ``calendar`` is
        supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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
    Convert a date to days since the reference date.

    Parameters
    ----------
    year, month, day, hour, minute, second : int, optional
        The date to be converted to days since ``referenceDate`` on the
        given ``calendar``.

    calendar : {'gregorian', 'noleap'}, optional
        A calendar to be used to convert days to a ``datetime.datetime``
        object.

    referenceDate : str, optional
        A reference date of the form::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    days : float
        The days since ``referenceDate`` on the given ``calendar``.

    Raises
    ------
    ValueError
        If an invalid ``referenceDate`` or ``calendar`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    calendar = _mpas_to_netcdf_calendar(calendar)

    date = datetime.datetime(year, month, day, hour, minute, second)

    return netCDF4.date2num(date, 'days since {}'.format(referenceDate),
                            calendar=calendar)


def _parse_date_string(dateString, isInterval=False):
    """
    Given a string containing a date, returns a tuple defining a date of the
    form (year, month, day, hour, minute, second) appropriate for constructing
    a datetime or timedelta

    Parameters
    ----------
    dateString : string
        A date and time in one of the followingformats::

            YYYY-MM-DD hh:mm:ss
            YYYY-MM-DD hh.mm.ss
            YYYY-MM-DD SSSSS
            DDD hh:mm:ss
            DDD hh.mm.ss
            DDD SSSSS
            hh.mm.ss
            hh:mm:ss
            YYYY-MM-DD
            YYYY-MM
            SSSSS

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
        If an invalid ``dateString`` is supplied.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    if isInterval:
        offset = 0
    else:
        offset = 1

    # change underscores to spaces so both can be supported
    dateString = dateString.rstrip('\x00').replace('_', ' ').strip()
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
    return (year, month, day, hour, minute, second)


def _mpas_to_netcdf_calendar(calendar):
    """
    Convert from MPAS calendar to NetCDF4 calendar names.
    """

    if calendar == 'gregorian_noleap':
        calendar = 'noleap'
    if calendar not in ['gregorian', 'noleap']:
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

    add_seconds = int(1e-6 * microsecond + 0.5)

    return date + datetime.timedelta(0, add_seconds)
