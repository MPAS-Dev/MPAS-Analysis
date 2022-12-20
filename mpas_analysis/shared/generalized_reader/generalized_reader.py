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
Utility functions for importing MPAS files into xarray. These functions extend
the capabilities of mpas_xarray to include mapping variable names from MPAS
names to MPAS-Analysis generalized names and support for slicing to given
start and end dates.

open_multifile_dataset : opens a data set, maps variable names, preprocess
    the data set removes repeated time indices, and slices the time coordinate
    to lie between desired start and end dates.
"""
# Authors
# -------
# Xylar Asay-Davis

import xarray
from functools import partial
import resource

from mpas_analysis.shared.mpas_xarray import mpas_xarray
from mpas_analysis.shared.timekeeping.utility import \
    string_to_days_since_date, days_to_datetime


def open_multifile_dataset(fileNames, calendar, config,
                           simulationStartTime=None,
                           timeVariableName='Time',
                           variableList=None, selValues=None,
                           iselValues=None, variableMap=None,
                           startDate=None, endDate=None,
                           chunking=None):
    """
    Opens and returns an xarray data set given file name(s) and the MPAS
    calendar name.

    Parameters
    ----------
    fileNames : list of strings
        A lsit of file paths to read

    calendar : {``'gregorian'``, ``'noleap'``}, optional
        The name of one of the calendars supported by MPAS cores

    config : mpas_tools.config.MpasConfigParser
        Contains configuration options

    simulationStartTime : string, optional
        The start date of the simulation, used to convert from time variables
        expressed as days since the start of the simulation to days since the
        reference date. ``simulationStartTime`` takes one of the following
        forms::

            0001-01-01
            0001-01-01 00:00:00

        ``simulationStartTime`` is only required if the MPAS time variable
        (identified by ``timeVariableName``) is a number of days since the
        start of the simulation.

    timeVariableName : string, optional
        The name of the time variable (typically ``'Time'`` if using a
        ``variableMap`` or ``'xtime'`` if not using a ``variableMap``)

    variableList : list of strings, optional
        If present, a list of variables to be included in the data set

    selValues : dict, optional
        A dictionary of coordinate names (keys) and values or arrays of
        values used to slice the variales in the data set.  See
        ``xarray.DataSet.sel()`` for details on how this dictonary is used.
        An example::

            selectCorrdValues = {'cellLon': 180.0}

    iselValues : dict, optional
        A dictionary of coordinate names (keys) and indices, slices or
        arrays of indices used to slice the variales in the data set.  See
        ``xarray.DataSet.isel()`` for details on how this dictonary is used.
        An example::

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

    chunking : None, int, True, dict, optional
        If integer is present, applies maximum chunk size from config file
        value ``maxChunkSize``, otherwise if None do not perform chunking.  If
        True, use automated chunking using default config value
        ``maxChunkSize``. If chunking is a dict use dictionary values for
        chunking.

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
        If the time variable is not found in the data set or if the time
        variable is a number of days since the start of the simulation but
        simulationStartTime is None.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Phillip J. Wolfram

    preprocess_partial = partial(_preprocess,
                                 calendar=calendar,
                                 simulationStartTime=simulationStartTime,
                                 timeVariableName=timeVariableName,
                                 variableList=variableList,
                                 selValues=selValues,
                                 iselValues=iselValues,
                                 variableMap=variableMap,
                                 startDate=startDate,
                                 endDate=endDate)

    ds = xarray.open_mfdataset(fileNames,
                               preprocess=preprocess_partial,
                               combine='nested',
                               concat_dim='Time',
                               decode_times=False)

    ds = mpas_xarray.remove_repeated_time_index(ds)

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
    # process chunking
    if chunking is True:
        # limit chunk size to prevent memory error
        chunking = config.getint('input', 'maxChunkSize')

    ds = mpas_xarray.process_chunking(ds, chunking)

    return ds


def _preprocess(ds, calendar, simulationStartTime, timeVariableName,
                variableList, selValues, iselValues, variableMap,
                startDate, endDate):
    """
    Performs variable remapping, then calls mpas_xarray.preprocess, to
    perform the remainder of preprocessing.

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate and with variable names to be
        substituted.

    calendar : {'gregorian', 'noleap'}
        The name of one of the calendars supported by MPAS cores

        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    simulationStartTime : string
        The start date of the simulation, used to convert from time variables
        expressed as days since the start of the simulation to days since the
        reference date. `simulationStartTime` takes one of the following
        forms::

            0001-01-01
            0001-01-01 00:00:00

        simulationStartTime is only required if the MPAS time variable
        (identified by time_variable_name) is a number of days since the
        start of the simulation.

    timeVariableName : string
        The name of the time variable (typically 'Time' if using a variable_map
        or 'xtime' if not using a variable_map)

    variableList : list of strings
        If present, a list of variables to be included in the data set

    selValues : dict
        A dictionary of coordinate names (keys) and values or arrays of
        values used to slice the variales in the data set.  See
        xarray.DataSet.sel() for details on how this dictonary is used.
        An example::

            selectCorrdValues = {'cellLon': 180.0}

    iselValues : dict
        A dictionary of coordinate names (keys) and indices, slices or
        arrays of indices used to slice the variales in the data set.  See
        xarray.DataSet.isel() for details on how this dictonary is used.
        An example::

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

    Returns
    -------
    ds : xarray.DataSet object
        A copy of the data set with the time coordinate set and which
        has been sliced.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Phillip J. Wolfram

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
                                calendar=calendar,
                                simulationStartTime=simulationStartTime,
                                timeVariableName=timeVariableName,
                                variableList=variableList,
                                selValues=selValues,
                                iselValues=iselValues)

    return ds


def _map_variable_name(variableName, ds, variableMap):
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
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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


def _rename_variables(ds, variableMap):
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
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    renameDict = {}
    for datasetVariable in ds.data_vars:
        for mapVariable in variableMap:
            renameList = variableMap[mapVariable]
            if datasetVariable in renameList:
                renameDict[datasetVariable] = mapVariable
                break

    return ds.rename(renameDict)


# vim: ai ts=4 sts=4 et sw=4 ft=python
