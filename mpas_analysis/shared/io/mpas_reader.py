"""
Utility functions for reading a single MPAS file into xarray and for removing
all but a given list of variables from a data set.

Authors
-------
Xylar Asay-Davis
"""

import xarray

from ..timekeeping.utility import string_to_days_since_date, days_to_datetime


def open_mpas_dataset(fileName, calendar,
                      timeVariableNames=['xtime_startMonthly',
                                         'xtime_endMonthly'],
                      variableList=None, startDate=None, endDate=None):  # {{{
    """
    Opens and returns an xarray data set given file name(s) and the MPAS
    calendar name.

    Parameters
    ----------
    fileName : str
        File path to read

    calendar : {``'gregorian'``, ``'gregorian_noleap'``}, optional
        The name of one of the calendars supported by MPAS cores

    timeVariableName : string, optional
        The name of the time variable (typically ``'xtime'``
        or ``['xtime_startMonthly', 'xtime_endMonthly']``)

    variableList : list of strings, optional
        If present, a list of variables to be included in the data set

    startDate, endDate : string or datetime.datetime, optional
        If present, the first and last dates to be used in the data set.  The
        time variable is sliced to only include dates within this range.

    Returns
    -------
    ds : ``xarray.Dataset``

    Raises
    ------
    TypeError
        If the time variable has an unsupported type (not a date string,
        a floating-pont number of days since the start of the simulation
        or a ``numpy.datatime64`` object).

    ValueError
        If the time variable is not found in the data set

    Authors
    -------
    Xylar Asay-Davis
    """

    ds = xarray.open_dataset(fileName, decode_cf=True, decode_times=False,
                             lock=False)

    ds = _parse_dataset_time(ds, timeVariableNames, calendar)

    if startDate is not None and endDate is not None:
        if isinstance(startDate, str):
            startDate = string_to_days_since_date(dateString=startDate,
                                                  calendar=calendar)
        if isinstance(endDate, str):
            endDate = string_to_days_since_date(dateString=endDate,
                                                calendar=calendar)

        # select only the data in the specified range of dates
        ds = ds.sel(Time=slice(startDate, endDate))

    if ds.dims['Time'] == 0:
        raise ValueError('The data set contains no Time entries between '
                         'dates {} and {}.'.format(
                             days_to_datetime(startDate, calendar=calendar),
                             days_to_datetime(endDate, calendar=calendar)))
    if variableList is not None:
        ds = subset_variables(ds, variableList)

    return ds  # }}}


def subset_variables(ds, variableList):  # {{{
    """
    Given a data set and a list of variable names, returns a new data set that
    contains only variables with those names.

    Parameters
    ----------
    ds : ``xarray.DataSet`` object
        The data set from which a subset of variables is to be extracted.

    variableList : string or list of strings
        The names of the variables to be extracted.

    Returns
    -------
    ds : ``xarray.DataSet`` object
        A copy of the original data set with only the variables in
        variableList.

    Raises
    ------
    ValueError
        If the resulting data set is empty.

    Authors
    -------
    Phillip J. Wolfram, Xylar Asay-Davis
    """

    allvars = ds.data_vars.keys()

    # get set of variables to drop (all ds variables not in vlist)
    dropvars = set(allvars) - set(variableList)

    # drop spurious variables
    ds = ds.drop(dropvars)

    # must also drop all coordinates that are not associated with the variables
    coords = set()
    for avar in ds.data_vars.keys():
        coords |= set(ds[avar].coords.keys())
    dropcoords = set(ds.coords.keys()) - coords

    # drop spurious coordinates
    ds = ds.drop(dropcoords)

    if len(ds.data_vars.keys()) == 0:
        raise ValueError(
                'Empty dataset is returned.\n'
                'Variables {}\n'
                'are not found within the dataset '
                'variables: {}.'.format(variableList, allvars))

    return ds  # }}}


def _parse_dataset_time(ds, inTimeVariableName, calendar,
                        outTimeVariableName='Time',
                        referenceDate='0001-01-01'):  # {{{
    """
    A helper function for computing a time coordinate from an MPAS time
    variable.  Given a data set and a time variable name (or tuple of 2
    time names), returns a new data set with time coordinate
    `outTimeVariableName` filled with days since `referenceDate`

    Parameters
    ----------
    ds : ``xarray.DataSet``
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate.

    inTimeVariableName : str or tuple or list of str
        The name of the time variable in the MPAS data set that will be
        used to build the 'Time' coordinate.  The array(s) named by
        inTimeVariableName should contain date strings.  Typically,
        inTimeVariableName is ``'xtime'``. If a list of two variable
        names is provided, times from the two are averaged together to
        determine the value of the time coordinate.  In such cases,
        inTimeVariableName is typically
        ``['xtime_startMonthly', 'xtime_endMonthly']``.

    calendar : {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

    outTimeVariableName : str
        The name of the coordinate to assign times to, typically 'Time'.

    referenceDate : str, optional
        The reference date for the time variable, typically '0001-01-01',
        taking one of the following forms::

            0001-01-01
            0001-01-01 00:00:00

    Returns
    -------
    dsOut : ``xarray.DataSet``
        A copy of the input data set with the `outTimeVariableName`
        coordinate containing the time coordinate parsed from
        `inTimeVariableName`.

    Raises
    ------
    TypeError
        If the time variable has an unsupported type (not a date string
        or a floating-pont number of days since the start of the simulatio).

    Authors
    -------
    Xylar Asay-Davis
    """

    if isinstance(inTimeVariableName, (tuple, list)):
        # we want to average the two
        assert(len(inTimeVariableName) == 2)

        dsStart = _parse_dataset_time(
            ds=ds,
            inTimeVariableName=inTimeVariableName[0],
            calendar=calendar,
            outTimeVariableName=outTimeVariableName,
            referenceDate=referenceDate)
        dsEnd = _parse_dataset_time(
            ds=ds,
            inTimeVariableName=inTimeVariableName[1],
            calendar=calendar,
            outTimeVariableName=outTimeVariableName,
            referenceDate=referenceDate)
        starts = dsStart[outTimeVariableName].values
        ends = dsEnd[outTimeVariableName].values

        # replace the time in starts with the mean of starts and ends
        dsOut = dsStart.copy()

        dsOut.coords['startTime'] = (outTimeVariableName, starts)
        dsOut.coords['endTime'] = (outTimeVariableName, ends)

        dsOut.coords[outTimeVariableName] = (outTimeVariableName,
                                             [starts[i] +
                                                 (ends[i] - starts[i])/2
                                                 for i in range(len(starts))])

    else:
        # there is just one time variable (either because we're recursively
        # calling the function or because we're not averaging).

        timeVar = ds[inTimeVariableName]

        if timeVar.dtype != '|S64':
            raise TypeError("timeVar of unsupported type {}.  String variable "
                            "expected.".format(timeVar.dtype))

        # this is an array of date strings like 'xtime'
        # convert to string
        timeStrings = [''.join(xtime).strip() for xtime in timeVar.values]
        days = string_to_days_since_date(dateString=timeStrings,
                                         referenceDate=referenceDate,
                                         calendar=calendar)

        dsOut = ds.copy()
        dsOut.coords[outTimeVariableName] = (outTimeVariableName, days)

    return dsOut  # }}}


# vim: ai ts=4 sts=4 et sw=4 ft=python
