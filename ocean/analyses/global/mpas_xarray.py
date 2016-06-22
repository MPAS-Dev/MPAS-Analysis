#!/usr/bin/env python
"""
mpas_xarray.py
==============================================================
Wrapper to handle importing MPAS files into xarray.

 Module does
 1. Converts MPAS "xtime" to xarray time.  Time dimension is assigned via
    `preprocess_mpas`.
 2. Converts MPAS "timeSinceStartOfSim" to xarray time for MPAS fields coming from the
    timeSeriesStatsAM.  Time dimension is assigned via `preprocess_mpas_timeSeriesStats`.
 3. Provides capability to remove redundant time entries from reading of multiple netCDF
    datasets via `remove_repeated_time_index`.

 Example Usage:

>>> from mpas_xarray import preprocess_mpas, remove_repeated_time_index
>>>
>>> ds = xarray.open_mfdataset('globalStats*nc', preprocess=preprocess_mpas)
>>> ds = remove_repeated_time_index(ds)

Phillip J. Wolfram
12/01/2015
"""

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

def subset_variables(ds, vlist): #{{{
    """
    Reduces an xarray dataset ds to only contain the variables in vlist.

    Phillip J. Wolfram
    05/05/2016
    """

    # get set of variables to drop (all ds variables not in vlist)
    dropvars = set(ds.data_vars.keys()) - set(vlist)

    # drop spurious variables
    ds = ds.drop(dropvars)

    # must also drop all coordinates that are not associated with the variables
    coords = set()
    for avar in ds.data_vars.keys():
        coords |= set(ds[avar].coords.keys())
    dropcoords = set(ds.coords.keys()) - coords

    # drop spurious coordinates
    ds = ds.drop(dropcoords)

    return ds #}}}

def assert_valid_datetimes(datetimes, yearoffset): #{{{
    """
    Ensure that datatimes are compatable with xarray

    Phillip J. Wolfram
    04/20/2016
    """
    assert datetimes[0].year > 1678, 'ERROR: yearoffset=%s'%(yearoffset) + \
            ' must be large enough to ensure datetimes larger than year 1678'
    assert datetimes[-1].year < 2262, 'ERROR: yearoffset=%s'%(yearoffset) + \
            ' must be large enough to ensure datetimes larger than year 2262'

    return #}}}

def general_processing(ds, datetimes, yearoffset, onlyvars): #{{{
    """
    Simple commonly used preprocessing commands general to multiple
    preprocessing routines.

    Phillip J. Wolfram
    05/05/2016
    """

    assert_valid_datetimes(datetimes, yearoffset)

    # append the corret time information
    ds.coords['Time'] = pd.to_datetime(datetimes)
    # record the yroffset
    ds.attrs.__setitem__('time_yearoffset', str(yearoffset))

    if onlyvars is not None:
        ds = subset_variables(ds, onlyvars)

    return ds #}}}

def preprocess_mpas(ds, yearoffset=1850, onlyvars=None): #{{{
    """
    Builds corret time specification for MPAS, allowing a year offset because the
    time must be between 1678 and 2262 based on the xarray library.

    The time specification is relevant for so-called time-slice model
    experiments, in which CO2 and greenhouse gas conditions are kept
    constant over the entire model simulation. Typical time-slice experiments
    are run with 1850 (pre-industrial) conditions and 2000 (present-day)
    conditions. Hence, a default date offset is chosen to be yearoffset=1849,
    monthoffset=12, dayoffset=31 (day 1 of an 1850 run will be seen as
    Jan 1st, 1850).

    The onlyvars option reduces the dataset to only include variables in the onlyvars list.
    If onlyvars=None, include all dataset variables.

    Phillip J. Wolfram
    05/05/2016
    """

    # compute shifted datetimes
    time = np.array([''.join(atime).strip() for atime in ds.xtime.values])
    datetimes = [datetime.datetime(yearoffset + int(x[:4]), int(x[5:7]), \
            int(x[8:10]), int(x[11:13]), int(x[14:16]), int(x[17:19])) for x in time]

    ds = general_processing(ds, datetimes, yearoffset, onlyvars)

    return ds #}}}

def preprocess_mpas_timeSeriesStats(ds,
        timestr='timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1',
        yearoffset=1849, monthoffset=12, dayoffset=31, onlyvars=None): #{{{
    """
    Builds corret time specification for MPAS timeSeriesStats analysis member fields,
    allowing a date offset because the time must be between 1678 and 2262
    based on the xarray library.

    This time specification is relevant for so-called time-slice model
    experiments, in which CO2 and greenhouse gas conditions are kept
    constant over the entire model simulation. Typical time-slice experiments
    are run with 1850 (pre-industrial) conditions and 2000 (present-day)
    conditions. Hence, a default date offset is chosen to be yearoffset=1849,
    monthoffset=12, dayoffset=31 (day 1 of an 1850 run will be seen as
    Jan 1st, 1850).

    The onlyvars option reduces the dataset to only include variables in the onlyvars list.
    If onlyvars=None, include all dataset variables.

    Milena Veneziani and Phillip J. Wolfram
    05/05/2016
    """

    # compute shifted datetimes
    daysSinceStart = ds[timestr]
    if (timestr == 'timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1'): 
    	datetimes = [datetime.datetime(yearoffset, monthoffset, dayoffset) + datetime.timedelta(x) for x in daysSinceStart.values]
    else:
    	daysSinceStart = pd.to_timedelta(daysSinceStart.values,unit='ns')
    	datetimes = [datetime.datetime(yearoffset, monthoffset, dayoffset) + x for x in daysSinceStart]

    ds = general_processing(ds, datetimes, yearoffset, onlyvars)

    return ds #}}}

def remove_repeated_time_index(ds): #{{{
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

    return ds #}}}

def test_load_mpas_xarray_datasets(path): #{{{
    ds = xr.open_mfdataset(path, preprocess=preprocess_mpas)
    ds = remove_repeated_time_index(ds)

    # make a simple plot from the data
    ds.Time.plot()
    plt.show()

    return #}}}

def test_load_mpas_xarray_timeSeriesStats_datasets(path): #{{{
    ds = xr.open_mfdataset(path, preprocess=preprocess_mpas_timeSeriesStats)
    ds = remove_repeated_time_index(ds)
    ds2 = xr.open_mfdataset(path, preprocess=preprocess_mpas)
    ds2 = remove_repeated_time_index(ds2)

    # make a simple plot from the data
    def plot_data(ds):
        var = ds["timeSeriesStatsMonthly_avg_iceAreaCell_1"]
        return var.where(var > 0).mean('nCells').plot()

    plot_data(ds)
    plot_data(ds2)
    plt.title("Curve centered around right times (b) \n "+\
              "Curve shifted towards end of avg period (g)")
    plt.show()

    return #}}}


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="inputfilename",
                      help="files to be opened with xarray, could be of form 'output*.nc'", \
                      metavar="FILE")
    parser.add_option("--istimeavg", dest="istimeavg",
                      help="option to use the preprocess for timeSeriesStatsAM fields")

    options, args = parser.parse_args()
    if not options.inputfilename:
        parser.error("Input filename or expression ('-f') is a required input..."+\
                " e.g., -f 'output*.npz'")

    if not options.istimeavg:
        test_load_mpas_xarray_datasets(options.inputfilename)
    else:
        test_load_mpas_xarray_timeSeriesStats_datasets(options.inputfilename)

# vim: ai ts=4 sts=4 et sw=4 ft=python
