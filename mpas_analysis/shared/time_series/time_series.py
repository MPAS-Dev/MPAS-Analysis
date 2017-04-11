"""
Utility functions related to time-series data sets

Authors
-------
Xylar Asay-Davis

Last Modified
-------------
04/08/2017
"""

import xarray as xr
import numpy
import os
import warnings

from ..timekeeping.utility import days_to_datetime


def cache_time_series(timesInDataSet, timeSeriesCalcFunction, cacheFileName,
                      calendar, yearsPerCacheUpdate=1,
                      printProgress=False):  # {{{
    '''
    Create or update a NetCDF file ``cacheFileName`` containing the given time
    series, calculated with ``timeSeriesCalcFunction`` over the given times,
    start and end year, and time frequency with which results are cached.

    Note: only works with climatologies where the mask (locations of ``NaN``
    values) doesn't vary with time.

    Parameters
    ----------
    timesInDataSet : array-like
        Times at which the time series is to be calculated, typically taken
        from ``ds.Times.values`` for a data set from which the time series
        will be extracted or computed.

    timeSeriesCalcFunction : function
        A function with arguments ``timeIndices``, indicating the entries in
        ``timesInDataSet`` to be computed, and ``firstCall``, indicating
        whether this is the first call to the funciton (useful for printing
        progress information).

    cacheFileName :  str
        The absolute path to the cache file where the times series will be
        stored

    calendar : ``{'gregorian', 'gregorian_noleap'}``
        The name of one of the calendars supported by MPAS cores, used to
        determine ``year`` and ``month`` from ``Time`` coordinate

    yearsPerCacheUpdate : int, optional
        The frequency with which the cache file is updated as the computation
        progresses.  If the computation is expensive, it may be useful to
        output the file frequently.  If not, there will be needless overhead
        in caching the file too frequently.

    printProgress: bool, optional
        Whether progress messages should be printed as the climatology is
        computed

    Returns
    -------
    climatology : object of same type as ``ds``
        A data set without the ``'Time'`` coordinate containing the mean
        of ds over all months in monthValues, weighted by the number of days
        in each month.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/08/2017

    '''

    timesProcessed = numpy.zeros(len(timesInDataSet), bool)
    # figure out which files to load and which years go in each file
    continueOutput = os.path.exists(cacheFileName)
    cacheDataSetExists = False
    if continueOutput:
        if printProgress:
            print '   Read in previously computed time series'
        # read in what we have so far

        try:
            dsCache = xr.open_dataset(cacheFileName, decode_times=False)
            cacheDataSetExists = True
        except IOError:
        # assuming the cache file is corrupt, so deleting it.
            message = 'Deleting cache file {}, which appears to have ' \
                      'been corrupted.'.format(cacheFileName)
            warnings.warn(message)
            os.remove(cacheFileName)

        if cacheDataSetExists:
            # force loading and then close so we can overwrite the file later
            dsCache.load()
            dsCache.close()
            for time in dsCache.Time.values:
                timesProcessed[timesInDataSet == time] = True

    datetimes = days_to_datetime(timesInDataSet, calendar=calendar)
    yearsInDataSet = numpy.array([date.year for date in datetimes])

    startYear = yearsInDataSet[0]
    endYear = yearsInDataSet[-1]

    firstProcessed = True
    for firstYear in range(startYear, endYear+1, yearsPerCacheUpdate):
        years = range(firstYear, numpy.minimum(endYear+1,
                                               firstYear+yearsPerCacheUpdate))

        mask = numpy.zeros(len(yearsInDataSet), bool)
        for year in years:
            mask = numpy.logical_or(mask, yearsInDataSet == year)
        mask = numpy.logical_and(mask, numpy.logical_not(timesProcessed))

        timeIndices = numpy.nonzero(mask)[0]

        if len(timeIndices) == 0:
            # no unprocessed time entries in this data range
            continue

        if printProgress:
            if firstProcessed:
                print '   Process and save time series'
            if yearsPerCacheUpdate == 1:
                print '     {:04d}'.format(years[0])
            else:
                print '     {:04d}-{:04d}'.format(years[0], years[-1])

        ds = timeSeriesCalcFunction(timeIndices, firstProcessed)
        firstProcessed = False

        if cacheDataSetExists:
            dsCache = xr.concat([dsCache, ds], dim='Time')
        else:
            dsCache = ds
            cacheDataSetExists = True

        dsCache.to_netcdf(cacheFileName)

    return dsCache.sel(Time=slice(timesInDataSet[0],timesInDataSet[-1]))  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
