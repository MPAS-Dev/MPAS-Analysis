import numpy as np
import xarray
from functools import partial

from ..timekeeping.utility import string_to_days_since_date, \
    string_to_datetime, days_to_datetime, datetime_to_days

"""
Utility functions for importing MPAS files into xarray.

open_multifile_dataset : open an xarray data set from MPAS data files
subset_variables : Keep only a subset of variables in a dataset
preprocess : preprocess a single file of an xarray dataset
remove_repeated_time_index : remove redundant indices in the 'Time' coordinate

Authors
-------
Phillip J. Wolfram, Xylar Asay-Davis

Last modified
-------------
02/22/2017
"""


def open_multifile_dataset(fileNames, calendar,
                           simulationStartTime=None,
                           timeVariableName='xtime',
                           variableList=None, selValues=None,
                           iselValues=None):  # {{{
    """
    Opens and returns an xarray data set given file name(s) and the MPAS
    calendar name.

    Parameters
    ----------
    fileNames : list of strings
        A lsit of file paths to read

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    simulationStartTime : string, optional
        The start date of the simulation, used to convert from time variables
        expressed as days since the start of the simulation to days since the
        reference date. `simulationStartTime` takes one of the following
        forms:
            0001-01-01

            0001-01-01 00:00:00

        simulationStartTime is only required if the MPAS time variable
        (identified by timeVariableName) is a number of days since the
        start of the simulation.

    timeVariableName : string, optional
        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    variableList : list of strings, optional
        If present, a list of variables to be included in the data set

    selectCorrdValues : dict, optional
        A dictionary of coordinate names (keys) and values or arrays of
        values used to slice the variales in the data set.  See
        xarray.dataset.sel() for details on how this dictonary is used.
        An example:
            selectCorrdValues = {'cellLon': 180.0}

    iselValues : dict, optional
        A dictionary of coordinate names (keys) and indices, slices or
        arrays of indices used to slice the variales in the data set.  See
        xarray.dataset.isel() for details on how this dictonary is used.
        An example:
            iselValues = {'nVertLevels': slice(0, 3),
                           'nCells': cellIDs}

    Returns
    -------
    ds : An xarray data set.

    Raises
    ------
    TypeError
        If the time variable has an unsupported type (not a date string or
        a floating-pont number of days since the start of the simulation).

    ValueError
        If the time variable is not found in the data set or if the time
        variable is a number of days since the start of the simulation but
        simulationStartTime is None.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/17/2017
    """

    preprocess_partial = partial(preprocess,
                                 calendar=calendar,
                                 simulationStartTime=simulationStartTime,
                                 timeVariableName=timeVariableName,
                                 variableList=variableList,
                                 selValues=selValues,
                                 iselValues=iselValues)

    ds = xarray.open_mfdataset(fileNames,
                               preprocess=preprocess_partial,
                               decode_times=False, concat_dim='Time')

    ds = remove_repeated_time_index(ds)

    return ds  # }}}


def subset_variables(ds, variableList):  # {{{
    """
    Given a data set and a list of variable names, returns a new data set that
    contains only variables with those names.

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set from which a subset of variables is to be extracted.

    variableList : string or list of strings
        The names of the variables to be extracted.

    Returns
    -------
    ds : xarray.DataSet object
        A copy of the original data set with only the variables in
        variableList.

    Raises
    ------
    ValueError
        If the resulting data set is empty.

    Authors
    -------
    Phillip J. Wolfram, Xylar Asay-Davis

    Last modified
    -------------
    02/16/2017
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


def preprocess(ds, calendar, simulationStartTime, timeVariableName,
               variableList, selValues, iselValues, maxChunkSize=1000):  # {{{
    """
    Builds correct time specification for MPAS, allowing a date offset
    because the time must be between 1678 and 2262 based on the xarray
    library.  Also, if slicing information (`selValues` and/or
    `iselValues`) was provided in `openMultifileDataSet`, this
    function performs the appropriate slicing on the data set.

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate.

    calendar : {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

    simulationStartTime : string, optinal
        The start date of the simulation, used to convert from time
        variables expressed as days since the start of the simulation to
        days since the reference date. `simulationStartTime` takes one
        of the following forms:
            0001-01-01

            0001-01-01 00:00:00

        simulationStartTime is only required if the MPAS time variable
        (identified by timeVariableName) is a number of days since the
        start of the simulation.

    timeVariableName : string
        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    variableList : list of strings
        If present, a list of variables to be included in the data set

    selectCorrdValues : dict
        A dictionary of coordinate names (keys) and values or arrays of
        values used to slice the variales in the data set.  See
        xarray.DataSet.sel() for details on how this dictonary is used.
        An example:
            selectCorrdValues = {'cellLon': 180.0}

    iselValues : dict
        A dictionary of coordinate names (keys) and indices, slices or
        arrays of indices used to slice the variales in the data set.  See
        xarray.DataSet.isel() for details on how this dictonary is used.
        An example:
            iselValues = {'nVertLevels': slice(0, 3),
                           'nCells': cellIDs}

    maxChunkSize : int
       Specifies the maximum chunk size to limit chunks used by dask to
       prevent out of memory errors for large datasets.

    Returns
    -------
    ds : xarray.DataSet object
        A copy of the data set with the time coordinate set and which
        has been sliced.

    Authors
    -------
    Phillip J. Wolfram, Milena Veneziani, Luke van Roekel
    and Xylar Asay-Davis

    Last modified
    -------------
    03/30/2017
    """

    ds = _parse_dataset_time(ds=ds,
                             inTimeVariableName=timeVariableName,
                             calendar=calendar,
                             simulationStartTime=simulationStartTime,
                             outTimeVariableName='Time',
                             referenceDate='0001-01-01')

    if variableList is not None:
        ds = subset_variables(ds,
                              _ensure_list(variableList))

    _assert_valid_selections(ds, selValues,
                             iselValues)

    if selValues is not None:
        ds = ds.sel(**selValues)

    if iselValues is not None:
        ds = ds.isel(**iselValues)

    chunks = {}
    for name in ds.chunks.keys():
        chunklim = np.asarray(ds.chunks[name]).max()
        chunks[name] = np.minimum(maxChunkSize, chunklim)

    ds = ds.chunk(chunks)

    return ds  # }}}


def remove_repeated_time_index(ds):  # {{{
    """
    Remove repeated times from xarray dataset.

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set potentially containing repeated time indices.

    Returns
    -------
    ds : xarray.DataSet object
        A copy of the original data set with any repeated time indices removed.

    Authors
    -------
    Phillip J. Wolfram, Xylar Asay-Davis

    Last modified
    -------------
    02/11/2017
    """
    # get repeated indices
    times = ds.Time.values
    indices = range(len(times))
    uniqueTimes = set()
    remove = []
    for timeIndex, time in enumerate(times):
        if time not in uniqueTimes:
            uniqueTimes.add(time)
        else:
            remove.append(timeIndex)

    # remove repeaded indices, working backwards from the last
    remove.reverse()
    for timeIndex in remove:
        indices.pop(timeIndex)

    # remove repeated indices
    ds = ds.isel(Time=indices)

    return ds  # }}}


def _assert_valid_selections(ds, selvals, iselvals):  # {{{
    """
    Ensure that dataset selections are compatable.

    It is possible selVals and iselVals may conflict, e.g., selVals
    restricts the dataset to a point where iselvals is unable to be
    satisfied, hence a check is needed to make sure that keys in selvals
    and iselvals are unique.  Additionally, keys for selvals and iselvals
    are tested to make sure they are dataset dimensions that can be used
    for selection.

    Authors
    -------
    Phillip J. Wolfram, Xylar Asay-Davis

    Last modified
    -------------
    02/10/2017
    """
    def test_vals_in_ds(vals, dims):
        if vals is not None:
            for val in vals.keys():
                assert val in dims, \
                    '{} is not a dimension in the dataset ' \
                    'that can be used for selection.'.format(val)

    if (selvals is not None) and (iselvals is not None):
        duplicatedkeys = len(np.intersect1d(selvals.keys(),
                                            iselvals.keys()))
        assert len(duplicatedkeys) == 0, \
            'Duplicated selection of variables {} was found!  ' \
            'Selection is ambiguous.'.format(duplicatedkeys)

    test_vals_in_ds(selvals, ds.dims)
    test_vals_in_ds(iselvals, ds.dims)

    return  # }}}


def _ensure_list(alist):  # {{{
    """
    Ensure that variables used as a list are actually lists.

    Authors
    -------
    Phillip J. Wolfram, Xylar Asay-Davis

    Last modified
    -------------
    02/10/2017
    """

    if isinstance(alist, str):
        # print 'Warning, converting %s to a list'%(alist)
        alist = [alist]

    return alist  # }}}


def _parse_dataset_time(ds, inTimeVariableName, calendar,
                        simulationStartTime, outTimeVariableName,
                        referenceDate):  # {{{
    """
    A helper function for computing a time coordinate from an MPAS time
    variable.  Given a data set and a time variable name (or tuple of 2
    time names), returns a new data set with time coordinate
    `outTimeVariableName` filled with days since `referenceDate`

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate.

    inTimeVariableName : string or tuple or list of strings
        The name of the time variable in the MPAS data set that will be
        used to build the 'Time' coordinate.  The array(s) named by
        inTimeVariableName should contain date strings or the number of
        days since the start of the simulation. Typically,
        inTimeVariableName is one of {'daysSinceStartOfSim','xtime'}.
        If a list of two variable
        names is provided, times from the two are averaged together to
        determine the value of the time coordinate.  In such cases,
        inTimeVariableName is typically {['xtime_start', 'xtime_end']}.

    calendar : {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores


    simulationStartTime : string
        The start date of the simulation, used to convert from time variables
        expressed as days since the start of the simulation to days since the
        reference date. `simulationStartTime` takes one of the following
        forms:
            0001-01-01

            0001-01-01 00:00:00

        simulationStartTime is only required if the MPAS time variable
        (identified by timeVariableName) is a number of days since the
        start of the simulation.

    outTimeVariableName : string
        The name of the coordinate to assign times to, typically 'Time'.

    referenceDate : string
        The reference date for the time variable, typically '0001-01-01',
        taking one of the following forms:
            0001-01-01

            0001-01-01 00:00:00

    Returns
    -------
    dataset : xarray.dataset object
        A copy of the input data set with the `outTimeVariableName`
        coordinate containing the time coordinate parsed from
        `inTimeVariableName`.

    Raises
    ------
    TypeError
        If the time variable has an unsupported type (not a date string
        or a floating-pont number of days since the start of the simulatio).
    ValueError
        If  the time variable is a number of days since the start of the
        simulation but simulationStartTime is None.

    Authors
    -------
    Xylar Asay-Davis

    Last modified
    -------------
    02/16/2017
    """

    if isinstance(inTimeVariableName, (tuple, list)):
        # we want to average the two
        assert(len(inTimeVariableName) == 2)

        dsStart = _parse_dataset_time(
            ds=ds,
            inTimeVariableName=inTimeVariableName[0],
            calendar=calendar,
            simulationStartTime=simulationStartTime,
            outTimeVariableName=outTimeVariableName,
            referenceDate=referenceDate)
        dsEnd = _parse_dataset_time(
            ds=ds,
            inTimeVariableName=inTimeVariableName[1],
            calendar=calendar,
            simulationStartTime=simulationStartTime,
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

        # The contents of the time variable is expected to be either a string
        # (|S64) or a float (meaning days since start of the simulation).

        timeVar = ds[inTimeVariableName]

        if timeVar.dtype == '|S64':
            # this is an array of date strings like 'xtime'
            # convert to string
            timeStrings = [''.join(xtime).strip() for xtime in timeVar.values]
            days = string_to_days_since_date(dateString=timeStrings,
                                             referenceDate=referenceDate,
                                             calendar=calendar)

        elif timeVar.dtype == 'float64':
            # this array contains floating-point days like
            # 'daysSinceStartOfSim'

            if simulationStartTime is None:
                raise ValueError('MPAS time variable {} appears to be a number of days since start \n'
                                 'of sim but simulationStartTime was not  supplied.'.format(inTimeVariableName))

            if (string_to_datetime(referenceDate) ==
                    string_to_datetime(simulationStartTime)):
                days = timeVar.values
            else:
                # a conversion may be required
                dates = days_to_datetime(days=timeVar.values,
                                         referenceDate=simulationStartTime,
                                         calendar=calendar)
                days = datetime_to_days(dates=dates,
                                        referenceDate=referenceDate,
                                        calendar=calendar)

        elif timeVar.dtype == 'timedelta64[ns]':
            raise TypeError('timeVar of unsupported type {}.  This is likely because xarray.open_dataset \n'
                            'was called with decode_times=True, which can mangle MPAS times.'.format(timeVar.dtype))
        else:
            raise TypeError("timeVar of unsupported type {}".format(
                timeVar.dtype))

        dsOut = ds.copy()
        dsOut.coords[outTimeVariableName] = (outTimeVariableName, days)

    return dsOut  # }}}


# vim: ai ts=4 sts=4 et sw=4 ft=python
