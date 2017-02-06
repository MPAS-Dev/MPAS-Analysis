<h1> Reorganize Timekeeping <br>
Xylar Asay-Davis <br>
date: 2017/02/06 <br>
</h1>
<h2> Summary </h2>
Currently, the `Date` class is used to parse a date object from a date string
(e.g. '0001-01-01_00:00:00') taken from MPAS namelists, streams files or time
variables (e.g. `xtime`).  However, this class assumes a 365-day calendar and
cannot easily be adapted to the Gregorian calendar also supported by MPAS
components (`config_calendar_type = 'gregorian'`).  Furthermore, existing 
routines exist to handle most of the capabilites
of the `Date` class.  The proposed reorganization would eliminate the `Date` class
in favor of a numer of helper functions that can be used to convert between various
date formats: date strings, days since a reference date, `datetime.datetime` objects
and `relativedelta` objects (see below).  The success of this reorganization will be
demonstrated when the existing analysis can be performed successfully with the new
utility functions with both MPAS calendars, the `'gregorian_noleap'` (365-day) calendar
used by most existing ACME and MPAS runs and the `'gregorian'` calendar also supported
in MPAS components.


<h1> Requirements </h1>

<h2> Requirement: Date string parsing supports both MPAS calendars <br>
Date last modified: 2017/02/06 <br>
Contributors: Xylar Asay-Davis
</h2>

There must be a way to parse dates from MPAS that is aware of the appropriate calendar 
stored in the `config_calendar_type` namelist option, either `'gregorian'` or 
`'gregorian_noleap'`.

<h2> Requirement: Capability of incrementing dates by a number of years and/or months <br>
Date last modified: 2017/02/06 <br>
Contributors: Xylar Asay-Davis
</h2>

The analysis requires a way of incrementing a given date by an interval specified in
not only days, hours, minutes and seconds but also months and years.  The standard
`datetime.timedelta` does not support increments by years and months because they are
not fixed periods of time.  The existing `Date` class in MPAS-Analysis supports
increments in months and years, but only for the `'gregorian_noleap'` (365-day) calendar.
A method must exist to increment dates on either calendar by a given number of years
and/or months (in addition to days, hours, etc.).


<h1> Design and Implementation </h1>

<h2> Implementation: Date string parsing supports both MPAS calendars <br>
Date last modified: 2017/02/06 <br>
Contributors: Xylar Asay-Davis
</h2>

The implementation is on the branch:
https://github.com/xylar/MPAS-Analysis/tree/timekeeping_reorg
and in PR #102

The function for converting a date string to a `datetime.datetime` is documented as follows:
```python
def stringToDatetime(dateString):
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
```

As long as `relativedelta` objects rather than `datetime.timedelta` objects are used to increment
`datetime.datetime` objects, `datetime.datetime` can be used to represent dates on either the Gregorian
or the 365-day calendar.

<h2> Implementation: Capability of incrementing dates by a number of years and/or months <br>
Date last modified: 2017/02/09 <br>
Contributors: Xylar Asay-Davis
</h2>

The implementation is on the branch:
https://github.com/xylar/MPAS-Analysis/tree/timekeeping_reorg
and in PR #102

The proposed implementation adds a new class MpasRelativeDelta derived from 
`dateutil.relativedelta.relativedelta` to compute the expected 
increments in years and months (as well as days, hours, minutes and seconds, as needed).
The class is documented as follows
```python
class MpasRelativeDelta(relativedelta):
    """
    MpasRelativeDelta is a subclass of dateutil.relativedelta for relative time
    intervals with different MPAS calendars.

    Only relative intervals (years, months, etc.) are supported and not the
    absolute date specifications (year, month, etc.).  Addition/subtraction
    of datetime.datetime objects (but not other MpasRelativeDelta,
    datetime.timedelta or other related objects) is supported.

    Author
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/09/2017
```

The function for converting a date string to a `MpasRelativeDelta` is documented as follows:
```python
from dateutil.relativedelta import relativedelta
...
def stringToRelativedelta(dateString, calendar='gregorian'):
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
```

<h1> Testing </h1>

<h2> Testing and Validation: Date string parsing supports both MPAS calendars <br>
Date last modified: 2017/02/08 <br>
Contributors: Xylar Asay-Davis
</h2>
Analysis will be run on Edison with all available configurations found in `configs/edison`.  As there
are currently no plans to run with the `gregorian` calendar option, we do not have test runs that use this
calendar.  If this situation changes in the future, we'll test at that time.

Regression tests previously for `Date` has been modified to test the new utility functions. New tests
have been added to test that dates with both `gregorian` and `gregorian_noleap` calendars behave as
expected, particularly around the leap day.

<h1> Testing </h1>
<h2> Testing and Validation: Capability of incrementing dates by a number of years and/or months <br>
Date last modified: 2017/02/06 <br>
Contributors: Xylar Asay-Davis
</h2>

Same as above.
