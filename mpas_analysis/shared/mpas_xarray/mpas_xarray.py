import datetime
import numpy as np
import pandas as pd
import xarray
from functools import partial

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
02/16/2017
"""


def open_multifile_dataset(fileNames, timeVariableName='Time',
                           variableList=None, selValues=None,
                           iselValues=None, yearOffset=0):  # {{{
    """
    Opens and returns an xarray data set given file name(s) and the MPAS
    calendar name.

    Parameters
    ----------
    fileNames : list of strings
        A lsit of file paths to read

    timeVariableName : string, optional
        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    variableList : list of strings, optional
        If present, a list of variables to be included in the data set

    selectCorrdValues : dict, optional
        A dictionary of coordinate names (keys) and values or arrays of
        values used to slice the variales in the data set.  See
        xarray.DataSet.sel() for details on how this dictonary is used.
        An example:
            selectCorrdValues = {'cellLon': 180.0}

    iselValues : dict, optional
        A dictionary of coordinate names (keys) and indices, slices or
        arrays of indices used to slice the variales in the data set.  See
        xarray.DataSet.isel() for details on how this dictonary is used.
        An example:
            iselValues = {'nVertLevels': slice(0, 3),
                           'nCells': cellIDs}

    yearOffset : float, optional
        An offset used to convert an MPAS date to a date in the range supported
        by xarray (numpy.datetime64).  Resulting dates must be between 1678
        and 2622.

    Returns
    -------
    ds : An xarray data set.

    Raises
    ------
    TypeError
        If the time variable has an unsupported type (not a date string,
        a floating-pont number of days since the start of the simulation
        or a numpy.datatime64 object).

    ValueError
        If the time variable is not found in the data set.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/16/2017
    """

    preprocess_partial = partial(preprocess,
                                 timeVariableName=timeVariableName,
                                 variableList=variableList,
                                 selValues=selValues,
                                 iselValues=iselValues,
                                 yearOffset=yearOffset)

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


def preprocess(ds, timeVariableName, variableList, selValues,
               iselValues, yearOffset):  # {{{
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

    yearOffset : float
        An offset used to convert an MPAS date to a date in the range supported
        by xarray (numpy.datetime64).  Resulting dates must be between 1678
        and 2622.

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
    02/16/2017
    """

    datetimes = _get_datetimes(ds, timeVariableName,
                               yearOffset)

    _assert_valid_datetimes(datetimes, yearOffset)

    # append the corret time information
    ds.coords['Time'] = datetimes
    # record the yroffset
    ds.attrs.__setitem__('time_yearoffset', str(yearOffset))

    if variableList is not None:
        ds = subset_variables(ds,
                              _ensure_list(variableList))

    _assert_valid_selections(ds, selValues,
                             iselValues)

    if selValues is not None:
        ds = ds.sel(**selValues)

    if iselValues is not None:
        ds = ds.isel(**iselValues)

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
    02/10/2017
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


def _assert_valid_datetimes(datetimes, yearOffset):  # {{{
    """
    Ensure that datatimes are compatable with xarray

    Authors
    -------
    Phillip J. Wolfram, Xylar Asay-Davis

    Last modified
    -------------
    02/16/2017
    """
    assert datetimes[0].year > 1678, \
        'ERROR: yearOffset={}'.format(yearOffset) + \
        ' must be large enough to ensure datetimes larger than year 1678'
    assert datetimes[-1].year < 2262, \
        'ERROR: yearOffset={}'.format(yearOffset) + \
        ' must be small enough to ensure datetimes smaller than year 2262'

    return  # }}}


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


def _get_datetimes(ds, timeVariableName, yearOffset):  # {{{
    """
    A helper function for computing a time coordinate from an MPAS time
    variable.  Given a data set and a time variable name (or tuple of 2
    time names), returns a list of objects representing the time coordinate

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate.

    timeVariableName : string or tuple or list of strings, optional
        The name of the time variable in the MPAS data set that will be
        used to build the 'Time' coordinate.  The array(s) named by
        timeVariableName should contain date strings or the number of
        days since the start of the simulation. Typically,
        timeVariableName is one of {'daysSinceStartOfSim','xtime'}.
        If a list of two variable
        names is provided, times from the two are averaged together to
        determine the value of the time coordinate.  In such cases,
        timeVariableName is typically {['xtime_start', 'xtime_end']}.

    yearOffset : int, optional
        An offset to be added to all years in the resulting array of dates.
        This offset is typically required to convert MPAS dates
        starting with year 0001 into dates in the range supported by xarray
        (1678 <= year <= 2262).

    Returns
    -------
    datetimes : list of datetime.datetime objects
        An array that represents the xarray time coordinate for the
        data set, based on the given MPAs time variable name.

    Raises
    ------
    TypeError
        If the time variable has an unsupported type (not a date string,
        a floating-pont number of days since the start of the simulation
        or a numpy.datatime64 object).

    Authors
    -------
    Xylar Asay-Davis

    Last modified
    -------------
    02/16/2017
    """

    if isinstance(timeVariableName, (tuple, list)):
        # we want to average the two
        assert(len(timeVariableName) == 2)
        starts = _get_datetimes(ds, timeVariableName[0], yearOffset)
        ends = _get_datetimes(ds, timeVariableName[1], yearOffset)
        datetimes = [starts[i] + (ends[i] - starts[i])/2
                     for i in range(len(starts))]
        return datetimes

    timeVariable = ds[timeVariableName]

    if timeVariable.dtype == '|S64':
        # this is a variable like date strings like 'xtime'
        time = [''.join(atime).strip() for atime in timeVariable.values]
        datetimes = [datetime.datetime(yearOffset + int(x[:4]),
                                       int(x[5:7]), int(x[8:10]),
                                       int(x[11:13]), int(x[14:16]),
                                       int(x[17:19]))
                     for x in time]
    elif timeVariable.dtype == 'float64':
        # this array contains floating-point days like
        # 'daysSinceStartOfSim'
        start = datetime.datetime(year=yearOffset+1, month=1, day=1)
        datetimes = [start + datetime.timedelta(x)
                     for x in timeVariable.values]
    elif timeVariable.dtype == 'timedelta64[ns]':
        # this array contains a variable like 'daysSinceStartOfSim' as a
        # timedelta64
        start = datetime.datetime(year=yearOffset+1, month=1, day=1)
        datetimes = [start + x for x in
                     pd.to_timedelta(timeVariable.values, unit='ns')]
    else:
        raise TypeError("timeVariable of unsupported type {}".format(
            timeVariable.dtype))

    return datetimes  # }}}

# vim: ai ts=4 sts=4 et sw=4 ft=python
