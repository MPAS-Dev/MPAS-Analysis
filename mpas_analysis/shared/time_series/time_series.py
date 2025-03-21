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
Utility functions related to time-series data sets
"""
# Authors
# -------
# Xylar Asay-Davis

import xarray as xr
import numpy
import os
import shutil
import glob
import subprocess

from mpas_analysis.shared.timekeeping.utility import days_to_datetime


def combine_time_series_with_ncrcat(inFileNames, outFileName,
                                    variableList=None, logger=None):
    """
    Uses ncrcat to extact time series from a series of files

    inFileNames : str or list of str
        A file name with wildcard(s) or a list of input files from which to
        extract the time series.

    outFileName : str
        The output NetCDF file where the time series should be written.

    variableList : list of str, optional
        A list of varibles to include.  All variables are included by default

    logger : `logging.Logger``, optional
        A logger to which ncclimo output should be redirected

    Raises
    ------
    OSError
        If ``ncrcat`` is not in the system path.

    Author
    ------
    Xylar Asay-Davis
    """

    if shutil.which('ncrcat') is None:
        raise OSError('ncrcat not found. Make sure the latest nco '
                      'package is installed: \n'
                      'conda install nco\n'
                      'Note: this presumes use of the conda-forge '
                      'channel.')

    if os.path.exists(outFileName):
        return

    if isinstance(inFileNames, str):
        inFileNames = sorted(glob.glob(inFileNames))

    args = ['ncrcat', '-4', '--record_append', '--no_tmp_fl']

    if variableList is not None:
        args.extend(['-v', ','.join(variableList)])

    printCommand = '{} {} ... {} {}'.format(' '.join(args), inFileNames[0],
                                            inFileNames[-1],
                                            outFileName)
    args.extend(inFileNames)
    args.append(outFileName)

    if logger is None:
        print('running: {}'.format(printCommand))
    else:
        logger.info('running: {}'.format(printCommand))
        for handler in logger.handlers:
            handler.flush()
    process = subprocess.Popen(args, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if stdout:
        stdout = stdout.decode('utf-8')
        for line in stdout.split('\n'):
            if logger is None:
                print(line)
            else:
                logger.info(line)
    if stderr:
        stderr = stderr.decode('utf-8')
        for line in stderr.split('\n'):
            if logger is None:
                print(line)
            else:
                logger.error(line)

    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode,
                                            ' '.join(args))


def cache_time_series(timesInDataSet, timeSeriesCalcFunction, cacheFileName,
                      calendar, yearsPerCacheUpdate=1,
                      logger=None):
    """
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

    calendar : {'gregorian', 'noleap'}
        The name of one of the calendars supported by MPAS cores, used to
        determine ``year`` and ``month`` from ``Time`` coordinate

    yearsPerCacheUpdate : int, optional
        The frequency with which the cache file is updated as the computation
        progresses.  If the computation is expensive, it may be useful to
        output the file frequently.  If not, there will be needless overhead
        in caching the file too frequently.

    logger : ``logging.Logger``, optional
        A logger to which to write output as the time series is computed

    Returns
    -------
    climatology : object of same type as ``ds``
        A data set without the ``'Time'`` coordinate containing the mean
        of ds over all months in monthValues, weighted by the number of days
        in each month.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    timesProcessed = numpy.zeros(len(timesInDataSet), bool)
    # figure out which files to load and which years go in each file
    continueOutput = os.path.exists(cacheFileName)
    cacheDataSetExists = False
    if continueOutput:
        if logger is not None:
            logger.info('   Read in previously computed time series')
        # read in what we have so far

        try:
            dsCache = xr.open_dataset(cacheFileName, decode_times=False)
            cacheDataSetExists = True
        except IOError:
            # assuming the cache file is corrupt, so deleting it.
            message = 'Deleting cache file {}, which appears to have ' \
                      'been corrupted.'.format(cacheFileName)
            if logger is None:
                print('Warning: {}'.format(message))
            else:
                logger.warning(message)
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
    for firstYear in range(startYear, endYear + 1, yearsPerCacheUpdate):
        years = range(firstYear, numpy.minimum(endYear + 1,
                                               firstYear + yearsPerCacheUpdate))

        mask = numpy.zeros(len(yearsInDataSet), bool)
        for year in years:
            mask = numpy.logical_or(mask, yearsInDataSet == year)
        mask = numpy.logical_and(mask, numpy.logical_not(timesProcessed))

        timeIndices = numpy.nonzero(mask)[0]

        if len(timeIndices) == 0:
            # no unprocessed time entries in this data range
            continue

        if logger is not None:
            if firstProcessed:
                logger.info('   Process and save time series')
            if yearsPerCacheUpdate == 1:
                logger.info('     {:04d}'.format(years[0]))
            else:
                logger.info('     {:04d}-{:04d}'.format(years[0], years[-1]))

        ds = timeSeriesCalcFunction(timeIndices, firstProcessed)
        firstProcessed = False

        if cacheDataSetExists:
            dsCache = xr.concat([dsCache, ds], dim='Time')
            # now sort the Time dimension:
            dsCache = dsCache.loc[{'Time': sorted(dsCache.Time.values)}]
        else:
            dsCache = ds
            cacheDataSetExists = True

        dsCache.to_netcdf(cacheFileName)

    return dsCache.sel(Time=slice(timesInDataSet[0], timesInDataSet[-1]))
