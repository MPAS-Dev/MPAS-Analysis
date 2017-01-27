#!/usr/bin/env python
"""
mpas_xarray.py
==============================================================
Wrapper to handle importing MPAS files into xarray.

 Module:
 1. converts MPAS time in various formats to xarray time.  The MPAS time
    variable is provided via
    `preprocess_mpas(..., timestr='xtime', ...)`.
    `timestr` can either be a single variable name or a pair of variable
    names.  In the latter case, each time variable is converted to an
    xarray time and the mean of the two times is used as the final xarray
    time.  Each variable name in `timestr` can refer either to a float
    array containing the the number of days since the start of the
    simulation (e.g. `daysSinceStartOfSim`) or a string variable with the
    date and time (e.g. `xtime`) in the usual MPAS format:
    YYYY-MM-DD_hh:mm:ss
 2. provides capability to remove redundant time entries from reading of
    multiple netCDF datasets via `remove_repeated_time_index`.
 3. provides capability to build a variable map between MPAS dycore variable
    names and those used in mpas_analysis.  This aids in supporting multiple
    versions of MPAS dycores.  The function `map_variable(...)` can be used
    to find the associated MPAS dycore variable name in a dataset given a
    variable name as used in mpas_analysis.  The function
    `rename_variables(...)` can be used to rename all variables in a variable
    map from their MPAS dycore names to the corresponding mpas_analysis names.

 Example Usage:

>>> from mpas_xarray import preprocess_mpas, remove_repeated_time_index
>>>
>>> ds = xarray.open_mfdataset('globalStats*nc', preprocess=preprocess_mpas)
>>> ds = remove_repeated_time_index(ds)

Phillip J. Wolfram, Xylar Asay-Davis
Last modified: 12/07/2016
"""

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


def subset_variables(ds, vlist):  # {{{
    """
    Reduces an xarray dataset ds to only contain the variables in vlist.

    Phillip J. Wolfram
    01/10/2017
    """

    allvars = ds.data_vars.keys()

    # get set of variables to drop (all ds variables not in vlist)
    dropvars = set(allvars) - set(vlist)

    # drop spurious variables
    ds = ds.drop(dropvars)

    # must also drop all coordinates that are not associated with the variables
    coords = set()
    for avar in ds.data_vars.keys():
        coords |= set(ds[avar].coords.keys())
    dropcoords = set(ds.coords.keys()) - coords

    # drop spurious coordinates
    ds = ds.drop(dropcoords)

    assert len(ds.variables.keys()) > 0, \
        'MPAS_XARRAY ERROR: Empty dataset is returned.\n' \
        'Variables {}\nare not found within the dataset ' \
        'variables: {}.'.format(vlist, allvars)

    return ds  # }}}


def assert_valid_datetimes(datetimes, yearoffset):  # {{{
    """
    Ensure that datatimes are compatable with xarray

    Phillip J. Wolfram
    04/20/2016
    """
    assert datetimes[0].year > 1678, \
        'ERROR: yearoffset={}'.format(yearoffset) + \
        ' must be large enough to ensure datetimes larger than year 1678'
    assert datetimes[-1].year < 2262, \
        'ERROR: yearoffset={}'.format(yearoffset) + \
        ' must be small enough to ensure datetimes smaller than year 2262'

    return  # }}}


def assert_valid_selections(ds, selvals, iselvals):  # {{{
    """
    Ensure that dataset selections are compatable.

    It is possible selVals and iselVals may conflict, e.g., selVals restricts
    the dataset to a point where iselvals is unable to be satisfied, hence a
    check is needed to make sure that keys in selvals and iselvals are unique.
    Additionally, keys for selvals and iselvals are tested to make sure they
    are dataset dimensions that can be used for selection.

    Phillip J. Wolfram
    Last modified: 12/07/2016
    """

    if (selvals is not None) and (iselvals is not None):
        duplicatedkeys = len(np.intersect1d(selvals.keys(), iselvals.keys()))
        assert len(duplicatedkeys) == 0, \
            'Duplicated selection of variables {} was found!  ' \
            'Selection is ambiguous.'.format(duplicatedkeys)

    def test_vals_in_ds(vals, dims):
        if vals is not None:
            for val in vals.keys():
                assert val in dims, \
                    '{} is not a dimension in the dataset ' \
                    'that can be used for selection.'.format(val)

    test_vals_in_ds(selvals, ds.dims)
    test_vals_in_ds(iselvals, ds.dims)

    return  # }}}


def ensure_list(alist):  # {{{
    """
    Ensure that variables used as a list are actually lists.

    Phillip J. Wolfram
    09/08/2016
    """

    if isinstance(alist, str):
        # print 'Warning, converting %s to a list'%(alist)
        alist = [alist]

    return alist  # }}}


def get_datetimes(ds, timestr, yearoffset):  # {{{
    """
    Computes a list of datetimes from the time variable in the dataset ds with
    variable name (or a list of 2 names) given by timestr, typically one of
    'daysSinceStartOfSim', 'xtime', or ['xtime_start', 'xtime_end'].

    The variable(s) pointed to by timestr should contain time information as a
    date string, a floating-point number of days or a number of days
    represented as a pandas timedelta (in ns).  The result is a list of
    datetimes corresponding to the input dates offset as appropriate by the
    yearoffset.

    Xylar Asay-Davis
    Last modified: 12/05/2016
    """

    if isinstance(timestr, (tuple, list)):
        # we want to average the two
        assert(len(timestr) == 2)
        starts = get_datetimes(ds, timestr[0], yearoffset)
        ends = get_datetimes(ds, timestr[1], yearoffset)
        datetimes = [starts[i] + (ends[i] - starts[i])/2
                     for i in range(len(starts))]
        return datetimes

    time_var = ds[timestr]

    if time_var.dtype == '|S64':
        # this is a variable like date strings like 'xtime'
        time = [''.join(atime).strip() for atime in time_var.values]
        datetimes = [datetime.datetime(yearoffset + int(x[:4]), int(x[5:7]),
                                       int(x[8:10]), int(x[11:13]),
                                       int(x[14:16]), int(x[17:19]))
                     for x in time]
    elif time_var.dtype == 'float64':
        # this array contains floating-point days like 'daysSinceStartOfSim'
        start = datetime.datetime(year=yearoffset+1, month=1, day=1)
        datetimes = [start + datetime.timedelta(x)
                     for x in time_var.values]
    elif time_var.dtype == 'timedelta64[ns]':
        # this array contains a variable like 'daysSinceStartOfSim' as a
        # timedelta64
        start = datetime.datetime(year=yearoffset+1, month=1, day=1)
        datetimes = [start + x for x in
                     pd.to_timedelta(time_var.values, unit='ns')]
    else:
        raise TypeError("time_var of unsupported type {}".format(
            time_var.dtype))

    return datetimes  # }}}


def map_variable(variable_name, ds, varmap):  # {{{
    """
    Find the variable (or list of variables) in dataset ds that map to the
    mpas_analysis variable given by variable_name.

    varmap is a dictionary with keys that are variable names used by
    MPAS-Analysis and values that are lists of possible names for the same
    variable in the MPAS dycore that produced the data set (which may differ
    between versions).

    Xylar Asay-Davis
    12/04/2016
    """
    possible_variables = varmap[variable_name]
    for var in possible_variables:
        if isinstance(var, (list, tuple)):
            allFound = True
            for subvar in var:
                if subvar not in ds.data_vars.keys():
                    allFound = False
                    break
            if allFound:
                return var

        elif var in ds.data_vars.keys():
            return var

    raise ValueError('Variable {} could not be mapped. None of the '
                     'possible mapping variables {}\n match any of the '
                     'variables in {}.'.format(
                         variable_name, possible_variables,
                         ds.data_vars.keys()))
    # }}}


def rename_variables(ds, varmap, timestr):  # {{{
    """
    Rename all variables in ds based on which are found in varmap.

    varmap is a dictionary with keys that are variable names used by
    MPAS-Analysis and values that are lists of possible names for the same
    variable in the MPAS dycore that produced the data set (which may differ
    between versions).

    timestr is points to the time variable(s), which are treated as a special
    case since they may need to be averaged.

    Returns a new timestr after mapping in timestr is in varmap, otherwise
    returns timestr unchanged.

    Xylar Asay-Davis
    12/08/2016
    """

    submap = varmap
    if timestr in varmap:
        # make a copy of varmap and remove timestr
        submap = varmap.copy()
        submap.pop(timestr, None)

    rename_dict = {}
    for ds_var in ds.data_vars:
        for map_var in submap:
            rename_list = varmap[map_var]
            if ds_var in rename_list:
                rename_dict[ds_var] = map_var
                break

    ds.rename(rename_dict, inplace=True)

    if timestr in varmap:
        timestr = map_variable(timestr, ds, varmap)

    return timestr  # }}}


def preprocess_mpas(ds, onlyvars=None, selvals=None, iselvals=None,
                    timestr='xtime', yearoffset=1849,
                    varmap=None):  # {{{
    """
    Builds correct time specification for MPAS, allowing a date offset because
    the time must be between 1678 and 2262 based on the xarray library.

    The time specification is relevant for so-called time-slice model
    experiments, in which CO2 and greenhouse gas conditions are kept
    constant over the entire model simulation. Typical time-slice experiments
    are run with 1850 (pre-industrial) conditions and 2000 (present-day)
    conditions. Hence, a default date offset is chosen to be yearoffset=1849,
    (year 0001 of an 1850 run will correspond with Jan 1st, 1850).

    The data set is assumed to have an array of date strings with variable
    name (or list of 2 names) given by timestr, typically one of
    'daysSinceStartOfSim', 'xtime', or ['xtime_start', 'xtime_end'].

    The onlyvars option reduces the dataset to only include variables in the
    onlyvars list. If onlyvars=None, include all dataset variables.

    iselvals and selvals provide index and value-based slicing operations for
    individual datasets prior to their merge via xarray.
    iselvals is a dictionary, e.g. iselvals = {'nVertLevels': slice(0, 3),
                                               'nCells': cellIDs}
    selvals is a dictionary, e.g. selvals = {'cellLon': 180.0}

    varmap is an optional dictionary that can be used to rename
    variables in the data set to standard names expected by mpas_analysis.
    If timestr is present in varmap, the values of varmap[timestr]
    will be used to determine the associated time variable in ds. However, the
    variable(s) associated with timestr in ds will not be renamed.  This is
    because there may be more than one variable in ds that maps to timestr
    (e.g. xtime_start and xtime_end), so that a one-to-one mapping is not
    possible for this variable.

    Phillip J. Wolfram, Milena Veneziani, Luke van Roekel and Xylar Asay-Davis
    Last modified: 12/05/2016
    """

    if varmap is not None:
        timestr = rename_variables(ds, varmap, timestr)

    datetimes = get_datetimes(ds, timestr, yearoffset)

    assert_valid_datetimes(datetimes, yearoffset)

    # append the corret time information
    ds.coords['Time'] = datetimes
    # record the yroffset
    ds.attrs.__setitem__('time_yearoffset', str(yearoffset))

    if onlyvars is not None:
        ds = subset_variables(ds, ensure_list(onlyvars))

    assert_valid_selections(ds, selvals, iselvals)

    if selvals is not None:
        ds = ds.sel(**selvals)

    if iselvals is not None:
        ds = ds.isel(**iselvals)

    return ds  # }}}


def remove_repeated_time_index(ds):  # {{{
    """
    Remove repeated times from xarray dataset.

    Phillip J. Wolfram
    12/01/2015
    """
    # get repeated indices
    time = ds.Time.values
    index = range(len(time))
    uniquetime = set()
    remove = []
    for tid, atime in enumerate(time):
        if atime not in uniquetime:
            uniquetime.add(atime)
        else:
            remove.append(tid)

    remove.reverse()
    for tid in remove:
        index.pop(tid)

    # remove repeated indices
    ds = ds.isel(Time=index)

    return ds  # }}}


def test_load_mpas_xarray_datasets(path):  # {{{
    ds = xr.open_mfdataset(path, preprocess=lambda x:
                           preprocess_mpas(x, yearoffset=1850))
    ds = remove_repeated_time_index(ds)

    # make a simple plot from the data
    ds.Time.plot()
    plt.show()

    return  # }}}


def test_load_mpas_xarray_timeSeriesStats_datasets(path):  # {{{
    timestr = 'timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1'
    ds = xr.open_mfdataset(path, preprocess=lambda x:
                           preprocess_mpas(x,
                                           timeSeriesStats=True,
                                           timestr=timestr))
    ds = remove_repeated_time_index(ds)
    ds2 = xr.open_mfdataset(path, preprocess=lambda x:
                            preprocess_mpas(x, yearoffset=1850))
    ds2 = remove_repeated_time_index(ds2)

    # make a simple plot from the data
    def plot_data(ds):
        var = ds["timeSeriesStatsMonthly_avg_iceAreaCell_1"]
        return var.where(var > 0).mean('nCells').plot()

    plot_data(ds)
    plot_data(ds2)
    plt.title("Curve centered around right times (b) \n " +
              "Curve shifted towards end of avg period (g)")
    plt.show()

    return  # }}}


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="files to be opened with xarray, could be of form "
                      "'output*.nc'",
                      metavar="FILE")
    parser.add_option("--istimeavg", dest="istimeavg",
                      help="option to use the preprocess for "
                      "timeSeriesStatsAM fields")

    options, args = parser.parse_args()
    if not options.inputfilename:
        parser.error("Input filename or expression ('-f') is a required"
                     "input, e.g. -f 'output*.npz'")

    if not options.istimeavg:
        test_load_mpas_xarray_datasets(options.inputfilename)
    else:
        test_load_mpas_xarray_timeSeriesStats_datasets(options.inputfilename)

# vim: ai ts=4 sts=4 et sw=4 ft=python
