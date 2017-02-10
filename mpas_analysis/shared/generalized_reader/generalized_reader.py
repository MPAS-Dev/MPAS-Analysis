"""
Utility functions for importing MPAS files into xarray. These functions extend
the capabilities of mpas_xarray to include mapping variable names from MPAS
names to MPAS-Analysis generalized names and support for slicing to given
start and end dates.

open_multifile_dataset : opens a data set, maps variable names, preprocess
    the data set removes repeated time indices, and slices the time coordinate
    to lie between desired start and end dates.

Authors
-------
Xylar Asay-Davis

Last modified
-------------
02/16/2017
"""

import xarray
from functools import partial

from ..mpas_xarray import mpas_xarray
from ..timekeeping.utility import stringToDatetime, clampToNumpyDatetime64


def open_multifile_dataset(fileNames, calendar, timeVariableName='Time',
                           variableList=None, selValues=None,
                           iselValues=None, variableMap=None,
                           startDate=None, endDate=None,
                           yearOffset=0):  # {{{
    """
    Opens and returns an xarray data set given file name(s) and the MPAS
    calendar name.

    Parameters
    ----------
    fileNames : list of strings
        A lsit of file paths to read

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    timeVariableName : string, optional
        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    variableList : list of strings, optional
        If present, a list of variables to be included in the data set

    selValues : dict, optional
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

    variableMap : dict, optional
        A dictionary with keys that are variable names used by
        MPAS-Analysis and values that are lists of possible names for the same
        variable in the MPAS dycore that produced the data set (which may
        differ between versions).

    startDate, endDate : string or datetime.datetime, optional
        If present, the first and last dates to be used in the data set.  The
        time variable is sliced to only include dates within this range.

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

    preprocess_partial = partial(_preprocess,
                                 calendar=calendar,
                                 timeVariableName=timeVariableName,
                                 variableList=variableList,
                                 selValues=selValues,
                                 iselValues=iselValues,
                                 variableMap=variableMap,
                                 startDate=startDate,
                                 endDate=endDate,
                                 yearOffset=yearOffset)

    ds = xarray.open_mfdataset(fileNames,
                               preprocess=preprocess_partial,
                               decode_times=False, concat_dim='Time')

    ds = mpas_xarray.remove_repeated_time_index(ds)

    if startDate is not None and endDate is not None:
        if isinstance(startDate, str):
            startDate = clampToNumpyDatetime64(stringToDatetime(
                startDate), yearOffset)
        if isinstance(endDate, str):
            endDate = clampToNumpyDatetime64(stringToDatetime(
                endDate), yearOffset)

    # select only the data in the specified range of dates
    ds = ds.sel(Time=slice(startDate, endDate))

    return ds  # }}}


def _preprocess(ds, calendar, timeVariableName, variableList, selValues,
                iselValues, variableMap, startDate, endDate,
                yearOffset):  # {{{
    """
    Performs variable remapping, then calls mpas_xarray.preprocess, to
    perform the remainder of preprocessing.

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate and with variable names to be
        substituted.

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    timeVariableName : string
        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    variableList : list of strings
        If present, a list of variables to be included in the data set

    selValues : dict
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

    variableMap : dict
        A dictionary with keys that are variable names used by
        MPAS-Analysis and values that are lists of possible names for the same
        variable in the MPAS dycore that produced the data set (which may
        differ between versions).

    startDate, endDate : string or datetime.datetime
        If present, the first and last dates to be used in the data set.  The
        time variable is sliced to only include dates within this range.

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
    Xylar Asay-Davis

    Last modified
    -------------
    02/16/2017
    """

    submap = variableMap

    # time_variable_names is a special case so we take it out of the map
    # and handle it manually (adding a new variable rather than renaming
    # an existing one)
    if variableMap is not None and timeVariableName in variableMap:
        # make a copy of variableMap and remove timeVariableName
        submap = variableMap.copy()
        submap.pop(timeVariableName, None)
        # temporarily change the time variable name
        timeVariableName = \
            _map_variable_name(timeVariableName,
                               ds,
                               variableMap)

    if submap is not None:
        ds = _rename_variables(ds, submap)

    # now that the variables are mapped, do the normal preprocessing in
    # mpas_xarray
    ds = mpas_xarray.preprocess(ds,
                                timeVariableName=timeVariableName,
                                variableList=variableList,
                                selValues=selValues,
                                iselValues=iselValues,
                                yearOffset=yearOffset)

    return ds  # }}}


def _map_variable_name(variableName, ds, variableMap):  # {{{
    """
    Given a `variableName` in a `variableMap` and an xarray `ds`,
    return the name of the the first variable in `variableMap[variableName]`
    that is found in ds.

    variableMap is a dictionary with keys that are variable names used by
    MPAS-Analysis and values that are lists of possible names for the same
    variable in the MPAS dycore that produced the data set (which may differ
    between versions).

    Parameters
    ----------
    variableName : string
        Name of a variable in `varriableMap`

    ds : `xarray.DataSet` object
        A data set in which the mapped variable name should be found

    variableMap : dict
        A dictionary with keys that are variable names used by
        MPAS-Analysis and values that are lists of possible names for the same
        variable in the MPAS dycore that produced the data set (which may
        differ between versions).

    Returns
    -------
    mappedVariableName : The corresponding variable name to `variableName`
        found in `ds`.

    Raises
    ------
    ValueError
        If none of the possible variable names in `variableMap[variableName]`
        can be found in `ds`.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/08/2017
    """
    possibleVariables = variableMap[variableName]
    for variable in possibleVariables:
        if isinstance(variable, (list, tuple)):
            allFound = True
            for subvariable in variable:
                if subvariable not in ds.data_vars.keys():
                    allFound = False
                    break
            if allFound:
                return variable

        elif variable in ds.data_vars.keys():
            return variable

    raise ValueError('Variable {} could not be mapped. None of the '
                     'possible mapping variables {}\n match any of the '
                     'variables in {}.'.format(
                         variableName, possibleVariables,
                         ds.data_vars.keys()))
    # }}}


def _rename_variables(ds, variableMap):  # {{{
    """
    Given an `xarray.DataSet` object `ds` and a dictionary mapping
    variable names `variableMap`, returns a new data set in which variables
    from `ds` with names equal to values in `variableMap` are renamed
    to the corresponding key in `variableMap`.

    Parameters
    ----------
    ds : `xarray.DataSet` object
        A data set in which the mapped variable names should be renamed

    variableMap : dict
        A dictionary with keys that are variable names used by
        MPAS-Analysis and values that are lists of possible names for the same
        variable in the MPAS dycore that produced the data set (which may
        differ between versions).

    Returns
    -------
    outDataSEt : A new `xarray.DataSet` object with the variable renamed.

    Author
    ------
    Xylar Asay-Davis

    Last modified
    -------------
    02/08/2017
    """

    renameDict = {}
    for datasetVariable in ds.data_vars:
        for mapVariable in variableMap:
            renameList = variableMap[mapVariable]
            if datasetVariable in renameList:
                renameDict[datasetVariable] = mapVariable
                break

    return ds.rename(renameDict)  # }}}

# vim: ai ts=4 sts=4 et sw=4 ft=python
