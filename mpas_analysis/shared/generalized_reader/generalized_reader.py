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
02/23/2017
"""

import xarray
from functools import partial
import resource

from ..mpas_xarray import mpas_xarray
from ..timekeeping.utility import string_to_days_since_date


def open_multifile_dataset(fileNames, calendar, config,
                           simulationStartTime=None,
                           timeVariableName='Time',
                           variableList=None, selValues=None,
                           iselValues=None, variableMap=None,
                           startDate=None, endDate=None):  # {{{
    """
    Opens and returns an xarray data set given file name(s) and the MPAS
    calendar name.

    Parameters
    ----------
    fileNames : list of strings
        A lsit of file paths to read

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores

    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    simulationStartTime : string, optional
        The start date of the simulation, used to convert from time variables
        expressed as days since the start of the simulation to days since the
        reference date. `simulationStartTime` takes one of the following
        forms:
            0001-01-01

            0001-01-01 00:00:00

        simulationStartTime is only required if the MPAS time variable
        (identified by time_variable_name) is a number of days since the
        start of the simulation.

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
        If the time variable is not found in the data set or if the time
        variable is a number of days since the start of the simulation but
        simulationStartTime is None.

    Author
    ------
    Xylar Asay-Davis, Phillip J. Wolfram

    Last modified
    -------------
    03/29/2017
    """

    # limit chunk size to prevent memory error
    maxChunkSize = config.getint('input', 'maxChunkSize')

    preprocess_partial = partial(_preprocess,
                                 calendar=calendar,
                                 simulationStartTime=simulationStartTime,
                                 timeVariableName=timeVariableName,
                                 variableList=variableList,
                                 selValues=selValues,
                                 iselValues=iselValues,
                                 variableMap=variableMap,
                                 startDate=startDate,
                                 endDate=endDate,
                                 maxChunkSize=maxChunkSize)

    kwargs = {'decode_times': False,
              'concat_dim': 'Time'}

    autocloseFileLimitFraction = config.getfloat('input',
                                                 'autocloseFileLimitFraction')

    # get the number of files that can be open at the same time.  We want the
    # "soft" limit because we'll get a crash if we exceed it.
    softLimit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]

    # use autoclose if we will use more than autocloseFileLimitFraction (50%
    # by default) of the soft limit of open files
    autoclose = len(fileNames) > softLimit*autocloseFileLimitFraction

    try:
        ds = xarray.open_mfdataset(fileNames,
                                   preprocess=preprocess_partial,
                                   autoclose=autoclose, **kwargs)
    except TypeError as e:
        if 'autoclose' in str(e):
            if autoclose:
                # This indicates that xarray version doesn't support autoclose
                print 'Warning: open_multifile_dataset is trying to use autoclose=True but\n' \
                      'it appears your xarray version doesn\'t support this argument. Will\n' \
                      'try again without autoclose argument.'

            ds = xarray.open_mfdataset(fileNames,
                                       preprocess=preprocess_partial,
                                       **kwargs)
        else:
            raise e

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

    # private record of autoclose use
    ds.attrs['_autoclose'] = autoclose

    return ds  # }}}


def _preprocess(ds, calendar, simulationStartTime, timeVariableName,
                variableList, selValues, iselValues, variableMap,
                startDate, endDate, maxChunkSize):  # {{{
    """
    Performs variable remapping, then calls mpas_xarray.preprocess, to
    perform the remainder of preprocessing.

    Parameters
    ----------
    ds : xarray.DataSet object
        The data set containing an MPAS time variable to be used to build
        an xarray time coordinate and with variable names to be
        substituted.

    calendar : {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

        The name of the time variable (typically 'Time' if using a variableMap
        or 'xtime' if not using a variableMap)

    simulationStartTime : string
        The start date of the simulation, used to convert from time variables
        expressed as days since the start of the simulation to days since the
        reference date. `simulationStartTime` takes one of the following
        forms:
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
    Xylar Asay-Davis, Phillip J. Wolfram

    Last modified
    -------------
    03/30/2017
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
                                calendar=calendar,
                                simulationStartTime=simulationStartTime,
                                timeVariableName=timeVariableName,
                                variableList=variableList,
                                selValues=selValues,
                                iselValues=iselValues,
                                maxChunkSize=maxChunkSize)

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
