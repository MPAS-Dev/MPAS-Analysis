# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
Functions for creating climatologies from monthly time series data
"""
# Authors
# -------
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import os
import numpy
from tempfile import TemporaryDirectory

from pyremap import Remapper, LatLonGridDescriptor, ProjectionGridDescriptor

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.timekeeping.utility import days_to_datetime

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, fingerprint_generator
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.climatology.comparison_descriptors import \
    get_comparison_descriptor, known_comparison_grids


def get_remapper(config, sourceDescriptor, comparisonDescriptor,
                 mappingFilePrefix, method, logger=None):  # {{{
    """
    Given config options and descriptions of the source and comparison grids,
    returns a ``pyremap.Remapper`` object that can be used to remap from source
    files or data sets to corresponding data sets on the comparison grid.

    If necessary, creates the mapping file containing weights and indices
    needed to perform remapping.

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        Contains configuration options

    sourceDescriptor : ``MeshDescriptor`` subclass object
        A description of the source mesh or grid

    comparisonDescriptor : ``MeshDescriptor`` subclass object
        A description of the comparison grid

    mappingFilePrefix : str
        A prefix to be prepended to the mapping file name

    method : {'bilinear', 'neareststod', 'conserve'}
        The method of interpolation used.

    logger : ``logging.Logger``, optional
        A logger to which ncclimo output should be redirected

    Returns
    -------
    remapper : ``pyremap.Remapper`` object
        A remapper that can be used to remap files or data sets from the source
        grid or mesh to the comparison grid.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    mappingFileName = None

    if not _matches_comparison(sourceDescriptor, comparisonDescriptor):
        # we need to remap because the grids don't match

        mappingBaseName = '{}_{}_to_{}_{}.nc'.format(
            mappingFilePrefix,
            sourceDescriptor.meshName,
            comparisonDescriptor.meshName,
            method)

        tryCustom = config.get('diagnostics', 'customDirectory') != 'none'
        if tryCustom:
            # first see if mapping files are in the custom directory
            mappingSubdirectory = build_config_full_path(
                config, 'diagnostics', 'mappingSubdirectory',
                baseDirectoryOption='customDirectory')

            mappingFileName = '{}/{}'.format(mappingSubdirectory,
                                         mappingBaseName)
        if not tryCustom or not os.path.exists(mappingFileName):
            # second see if mapping files are in the base directory

            mappingSubdirectory = build_config_full_path(
                config, 'diagnostics', 'mappingSubdirectory',
                baseDirectoryOption='base_path')

            mappingFileName = '{}/{}'.format(mappingSubdirectory,
                                             mappingBaseName)

        if not os.path.exists(mappingFileName):
            # we don't have a mapping file yet, so get ready to create one
            # in the output subfolder if needed
            mappingSubdirectory = \
                build_config_full_path(config, 'output',
                                       'mappingSubdirectory')
            make_directories(mappingSubdirectory)
            mappingFileName = '{}/{}'.format(mappingSubdirectory,
                                             mappingBaseName)

    remapper = Remapper(sourceDescriptor, comparisonDescriptor,
                        mappingFileName)

    mpiTasks = config.getWithDefault('execute', 'mapMpiTasks', 1)
    esmf_parallel_exec = config.get('execute', 'mapParallelExec')
    if esmf_parallel_exec == 'None':
        esmf_parallel_exec = None

    mappingSubdirectory = \
        build_config_full_path(config, 'output',
                               'mappingSubdirectory')
    make_directories(mappingSubdirectory)
    with TemporaryDirectory(dir=mappingSubdirectory) as tempdir:
        remapper.build_mapping_file(method=method, logger=logger,
                                    mpiTasks=mpiTasks, tempdir=tempdir,
                                    esmf_parallel_exec=esmf_parallel_exec)

    return remapper  # }}}


def compute_monthly_climatology(ds, calendar=None, maskVaries=True):  # {{{
    """
    Compute monthly climatologies from a data set.  The mean is weighted but
    the number of days in each month of the data set, ignoring values masked
    out with NaNs.  If the month coordinate is not present, a data array
    ``month`` will be added based on ``Time`` and the provided calendar.

    Parameters
    ----------
    ds : xarray.Dataset or xarray.DataArray
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01 or ``month`` coordinate

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate, so must be supplied if
        ``ds`` does not already have a ``month`` coordinate or data array

    maskVaries : bool, optional
        If the mask (where variables in ``ds`` are ``NaN``) varies with time.
        If not, the weighted average does not need make extra effort to account
        for the mask.  Most MPAS fields will have masks that don't vary in
        time, whereas observations may sometimes be present only at some
        times and not at others, requiring ``maskVaries = True``.

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

    def compute_one_month_climatology(ds):
        monthValues = list(ds.month.values)
        return compute_climatology(ds, monthValues, calendar, maskVaries)

    ds = add_years_months_days_in_month(ds, calendar)

    monthlyClimatology = \
        ds.groupby('month').map(compute_one_month_climatology)

    return monthlyClimatology  # }}}


def compute_climatology(ds, monthValues, calendar=None,
                        maskVaries=True):  # {{{
    """
    Compute a monthly, seasonal or annual climatology data set from a data
    set.  The mean is weighted but the number of days in each month of
    the data set, ignoring values masked out with NaNs.  If the month
    coordinate is not present, a data array ``month`` will be added based
    on ``Time`` and the provided calendar.

    Parameters
    ----------
    ds : xarray.Dataset or xarray.DataArray
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01 or ``month`` coordinate

    monthValues : int or array-like of ints
        A single month or an array of months to be averaged together

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate, so must be supplied if
        ``ds`` does not already have a ``month`` coordinate or data array

    maskVaries : bool, optional
        If the mask (where variables in ``ds`` are ``NaN``) varies with time.
        If not, the weighted average does not need make extra effort to account
        for the mask.  Most MPAS fields will have masks that don't vary in
        time, whereas observations may sometimes be present only at some
        times and not at others, requiring ``maskVaries = True``.

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

    ds = add_years_months_days_in_month(ds, calendar)

    mask = xr.zeros_like(ds.month, bool)

    for month in monthValues:
        mask = numpy.logical_or(mask, ds.month == month)

    climatologyMonths = ds.where(mask, drop=True)

    climatology = _compute_masked_mean(climatologyMonths, maskVaries)

    return climatology  # }}}


def add_years_months_days_in_month(ds, calendar=None):  # {{{
    '''
    Add ``year``, ``month`` and ``daysInMonth`` as data arrays in ``ds``.
    The number of days in each month of ``ds`` is computed either using the
    ``startTime`` and ``endTime`` if available or assuming ``gregorian_noleap``
    calendar and ignoring leap years.  ``year`` and ``month`` are computed
    accounting correctly for the the calendar.

    Parameters
    ----------
    ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01

    calendar : {'gregorian', 'gregorian_noleap'}, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``year`` and ``month`` from ``Time`` coordinate

    Returns
    -------
    ds : object of same type as ``ds``
        The data set with ``year``, ``month`` and ``daysInMonth`` data arrays
        added (if not already present)
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    if ('year' in ds.coords and 'month' in ds.coords and
            'daysInMonth' in ds.coords):
        return ds

    ds = ds.copy()

    if 'year' not in ds.coords or 'month' not in ds.coords:
        if calendar is None:
            raise ValueError('calendar must be provided if month and year '
                             'coordinate is not in ds.')
        datetimes = days_to_datetime(ds.Time, calendar=calendar)

    if 'year' not in ds.coords:
        ds.coords['year'] = ('Time', [date.year for date in datetimes])

    if 'month' not in ds.coords:
        ds.coords['month'] = ('Time', [date.month for date in datetimes])

    if 'daysInMonth' not in ds.coords:
        if 'startTime' in ds.coords and 'endTime' in ds.coords:
            ds.coords['daysInMonth'] = ds.endTime - ds.startTime
        else:
            if calendar == 'gregorian':
                print('Warning: The MPAS run used the Gregorian calendar '
                      'but does not appear to have\n'
                      'supplied start and end times.  Climatologies '
                      'will be computed with\n'
                      'month durations ignoring leap years.')

            daysInMonth = numpy.array(
                [constants.daysInMonth[int(month) - 1] for
                 month in ds.month.values], float)
            ds.coords['daysInMonth'] = ('Time', daysInMonth)

    return ds  # }}}


def remap_and_write_climatology(config, climatologyDataSet,
                                climatologyFileName, remappedFileName,
                                remapper, logger=None):  # {{{
    """
    Given a field in a climatology data set, use the ``remapper`` to remap
    horizontal dimensions of all fields, write the results to an output file,
    and return the remapped data set.

    Note that ``climatologyFileName`` and ``remappedFileName`` will be
    overwritten if they exist, so if this behavior is not desired, the calling
    code should skip this call if the files exist and simply load the contents
    of ``remappedFileName``.

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        Contains configuration options

    climatologyDataSet : ``xarray.DataSet`` or ``xarray.DataArray`` object
        A data set containing a climatology

    fieldName : str
        A field within the climatology to be remapped

    climatologyFileName : str
        The name of the output file to which the data set should be written
        before remapping (if using ncremap).

    remappedFileName : str
        The name of the output file to which the remapped data set should
        be written.

    remapper : ``pyremap.Remapper`` object
        A remapper that can be used to remap files or data sets to a
        comparison grid.

    logger : ``logging.Logger``, optional
        A logger to which ncclimo output should be redirected

    Returns
    -------
    remappedClimatology : ``xarray.DataSet`` or ``xarray.DataArray`` object
        A data set containing the remapped climatology
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    useNcremap = config.getboolean('climatology', 'useNcremap')

    if remapper.mappingFileName is None:
        # no remapping is needed
        remappedClimatology = climatologyDataSet
    else:
        renormalizationThreshold = config.getfloat(
            'climatology', 'renormalizationThreshold')
        parallel_exec = config.get(
            'execute', 'ncremapParallelExec')
        if parallel_exec == 'None':
            parallel_exec = None

        if useNcremap:
            if not os.path.exists(climatologyFileName):
                write_netcdf(climatologyDataSet, climatologyFileName)
            remapper.remap_file(inFileName=climatologyFileName,
                                outFileName=remappedFileName,
                                overwrite=True,
                                renormalize=renormalizationThreshold,
                                logger=logger,
                                parallel_exec=parallel_exec)
            remappedClimatology = xr.open_dataset(remappedFileName)
        else:

            remappedClimatology = remapper.remap(climatologyDataSet,
                                                 renormalizationThreshold)
            write_netcdf(remappedClimatology, remappedFileName)
    return remappedClimatology  # }}}


def get_unmasked_mpas_climatology_directory(config, op='avg'):  # {{{
    """
    Get the directory for an unmasked MPAS climatology produced by ncclimo,
    making the directory if it doesn't already exist

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        configuration options

    op : {'avg', 'min', 'max'}
         operator for monthly stats
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    climatologyOpDirectory = get_climatology_op_directory(config, op)

    mpasMeshName = config.get('input', 'mpasMeshName')

    directory = '{}/unmasked_{}'.format(climatologyOpDirectory,
                                        mpasMeshName)

    make_directories(directory)
    return directory  # }}}


def get_unmasked_mpas_climatology_file_name(config, season, componentName,
                                            op='avg'):
    # {{{
    """
    Get the file name for an unmasked MPAS climatology produced by ncclimo

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        configuration options

    season : str
        One of the seasons in ``constants.monthDictionary``

    componentName : {'ocean', 'seaIce'}
        The MPAS component for which the climatology is being computed

    op : {'avg', 'min', 'max'}
         operator for monthly stats
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')

    if componentName == 'ocean':
        ncclimoModel = 'mpaso'
    elif componentName == 'seaIce':
        ncclimoModel = 'mpascice'
    else:
        raise ValueError('component {} is not supported by ncclimo.\n'
                         'Check with Charlie Zender and Xylar Asay-Davis\n'
                         'about getting it added'.format(componentName))

    directory = get_unmasked_mpas_climatology_directory(config, op)

    make_directories(directory)
    monthValues = sorted(constants.monthDictionary[season])
    startMonth = monthValues[0]
    endMonth = monthValues[-1]

    suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
        startYear, startMonth, endYear, endMonth)

    if season in constants.abrevMonthNames:
        season = '{:02d}'.format(monthValues[0])
    fileName = '{}/{}_{}_{}.nc'.format(directory, ncclimoModel,
                                       season, suffix)
    return fileName  # }}}


def get_masked_mpas_climatology_file_name(config, season, componentName,
                                          climatologyName, op='avg'):  # {{{
    """
    Get the file name for a masked MPAS climatology

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        Configuration options

    season : str
        One of the seasons in ``constants.monthDictionary``

    componentName : {'ocean', 'seaIce'}
        The MPAS component for which the climatology is being computed

    climatologyName : str
        The name of the climatology (typically the name of a field to mask
        and later remap)

    op : {'avg', 'min', 'max'}
         operator for monthly stats
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')
    mpasMeshName = config.get('input', 'mpasMeshName')

    if componentName == 'ocean':
        ncclimoModel = 'mpaso'
    elif componentName == 'seaIce':
        ncclimoModel = 'mpascice'
    else:
        raise ValueError('component {} is not supported by ncclimo.\n'
                         'Check with Charlie Zender and Xylar Asay-Davis\n'
                         'about getting it added'.format(componentName))

    climatologyOpDirectory = get_climatology_op_directory(config, op)

    stageDirectory = '{}/masked'.format(climatologyOpDirectory)

    directory = '{}/{}_{}'.format(
        stageDirectory, climatologyName,
        mpasMeshName)

    make_directories(directory)

    monthValues = sorted(constants.monthDictionary[season])
    startMonth = monthValues[0]
    endMonth = monthValues[-1]

    suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
        startYear, startMonth, endYear, endMonth)

    if season in constants.abrevMonthNames:
        season = '{:02d}'.format(monthValues[0])
    fileName = '{}/{}_{}_{}.nc'.format(
        directory, ncclimoModel, season, suffix)

    return fileName  # }}}


def get_remapped_mpas_climatology_file_name(config, season, componentName,
                                            climatologyName,
                                            comparisonGridName,
                                            op='avg'):  # {{{
    """
    Get the file name for a masked MPAS climatology

    Parameters
    ----------
    config : mpas_tools.config.MpasConfigParser
        Configuration options

    season : str
        One of the seasons in ``constants.monthDictionary``

    componentName : {'ocean', 'seaIce'}
        The MPAS component for which the climatology is being computed

    climatologyName : str
        The name of the climatology (typically the name of a field to mask
        and later remap)

    comparisonGridName : str
        The name of the comparison grid to use for remapping.  If it is one
        of the known comparison grid names, the full grid name is looked up via
        :py:func:`mpas_analysis.shared.climatology.get_comparison_descriptor()`

    op : {'avg', 'min', 'max'}
         operator for monthly stats
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')
    mpasMeshName = config.get('input', 'mpasMeshName')

    if componentName == 'ocean':
        ncclimoModel = 'mpaso'
    elif componentName == 'seaIce':
        ncclimoModel = 'mpascice'
    else:
        raise ValueError('component {} is not supported by ncclimo.\n'
                         'Check with Charlie Zender and Xylar Asay-Davis\n'
                         'about getting it added'.format(componentName))

    climatologyOpDirectory = get_climatology_op_directory(config, op)

    if comparisonGridName in known_comparison_grids:
        comparisonDescriptor = get_comparison_descriptor(config,
                                                         comparisonGridName)
        comparisonFullMeshName = comparisonDescriptor.meshName
    else:
        comparisonFullMeshName = comparisonGridName.replace(' ', '_')

    stageDirectory = '{}/remapped'.format(climatologyOpDirectory)

    directory = '{}/{}_{}_to_{}'.format(stageDirectory, climatologyName,
                                        mpasMeshName, comparisonFullMeshName)

    make_directories(directory)

    monthValues = sorted(constants.monthDictionary[season])
    startMonth = monthValues[0]
    endMonth = monthValues[-1]

    suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
        startYear, startMonth, endYear, endMonth)

    if season in constants.abrevMonthNames:
        season = '{:02d}'.format(monthValues[0])
    fileName = '{}/{}_{}_{}.nc'.format(
        directory, ncclimoModel, season, suffix)

    return fileName  # }}}


def get_climatology_op_directory(config, op='avg'):
    '''
    Get the output directory for MPAS climatologies from output with the given
    monthly operator: avg, min or max
    '''
    climatologyBaseDirectory = build_config_full_path(
        config, 'output', 'mpasClimatologySubdirectory')

    return '{}/{}'.format(climatologyBaseDirectory, op)


def _compute_masked_mean(ds, maskVaries):  # {{{
    '''
    Compute the time average of data set, masked out where the variables in ds
    are NaN and, if ``maskVaries == True``, weighting by the number of days
    used to compute each monthly mean time in ds.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def ds_to_weights(ds):
        # make an identical data set to ds but replacing all data arrays with
        # nonnull applied to that data array
        weights = ds.copy(deep=True)
        if isinstance(ds, xr.core.dataarray.DataArray):
            weights = ds.notnull()
        elif isinstance(ds, xr.core.dataset.Dataset):
            for var in ds.data_vars:
                weights[var] = ds[var].notnull()
        else:
            raise TypeError('ds must be an instance of either xarray.Dataset '
                            'or xarray.DataArray.')

        return weights

    if maskVaries:
        dsWeightedSum = (ds * ds.daysInMonth).sum(dim='Time', keep_attrs=True)

        weights = ds_to_weights(ds)

        weightSum = (weights * ds.daysInMonth).sum(dim='Time')

        timeMean = dsWeightedSum / weightSum.where(weightSum > 0.)
    else:
        days = ds.daysInMonth.sum(dim='Time')

        dsWeightedSum = (ds * ds.daysInMonth).sum(dim='Time', keep_attrs=True)

        timeMean = dsWeightedSum / days.where(days > 0.)

    return timeMean  # }}}


def _matches_comparison(obsDescriptor, comparisonDescriptor):  # {{{
    '''
    Determine if the two meshes are the same
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    if isinstance(obsDescriptor, ProjectionGridDescriptor) and \
            isinstance(comparisonDescriptor, ProjectionGridDescriptor):
        # pretty hard to determine if projections are the same, so we'll rely
        # on the grid names
        match = obsDescriptor.meshName == comparisonDescriptor.meshName and \
            len(obsDescriptor.x) == len(comparisonDescriptor.x) and \
            len(obsDescriptor.y) == len(comparisonDescriptor.y) and \
            numpy.all(numpy.isclose(obsDescriptor.x,
                                    comparisonDescriptor.x)) and \
            numpy.all(numpy.isclose(obsDescriptor.y,
                                    comparisonDescriptor.y))
    elif isinstance(obsDescriptor, LatLonGridDescriptor) and \
            isinstance(comparisonDescriptor, LatLonGridDescriptor):
        match = ((('degree' in obsDescriptor.units and
                   'degree' in comparisonDescriptor.units) or
                  ('radian' in obsDescriptor.units and
                   'radian' in comparisonDescriptor.units)) and
                 len(obsDescriptor.lat) == len(comparisonDescriptor.lat) and
                 len(obsDescriptor.lon) == len(comparisonDescriptor.lon) and
                 numpy.all(numpy.isclose(obsDescriptor.lat,
                                         comparisonDescriptor.lat)) and
                 numpy.all(numpy.isclose(obsDescriptor.lon,
                                         comparisonDescriptor.lon)))
    else:
        match = False

    return match  # }}}


def _setup_climatology_caching(ds, startYearClimo, endYearClimo,
                               yearsPerCacheFile, cachePrefix,
                               monthValues):  # {{{
    '''
    Determine which cache files already exist, which are incomplete and which
    years are present in each cache file (whether existing or to be created).
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    cacheInfo = []

    cacheIndices = -1 * numpy.ones(ds.dims['Time'], int)
    monthsInDs = ds.month.values
    yearsInDs = ds.year.values

    # figure out which files to load and which years go in each file
    for firstYear in range(startYearClimo, endYearClimo + 1,
                           yearsPerCacheFile):
        years = range(firstYear, firstYear + yearsPerCacheFile)

        yearString, fileSuffix = _get_year_string(years[0], years[-1])
        outputFileClimo = '{}_{}.nc'.format(cachePrefix, fileSuffix)

        done = False
        if os.path.exists(outputFileClimo):
            # already cached
            dsCached = None
            try:
                dsCached = xr.open_dataset(outputFileClimo)
            except IOError:
                # assuming the cache file is corrupt, so deleting it.
                print('Warning: Deleting cache file {}, which appears to '
                      'have been corrupted.'.format(outputFileClimo))

                os.remove(outputFileClimo)

            monthsIfDone = len(monthValues) * len(years)
            if ((dsCached is not None) and
                    (dsCached.attrs['totalMonths'] == monthsIfDone)):
                # also complete, so we can move on
                done = True
            if dsCached is not None:
                dsCached.close()

        cacheIndex = len(cacheInfo)
        for year in years:
            for month in monthValues:
                mask = numpy.logical_and(yearsInDs == year,
                                         monthsInDs == month)
                cacheIndices[mask] = cacheIndex

        if numpy.count_nonzero(cacheIndices == cacheIndex) == 0:
            continue

        cacheInfo.append((outputFileClimo, done, yearString))

    ds = ds.copy()
    ds.coords['cacheIndices'] = ('Time', cacheIndices)

    return cacheInfo, cacheIndices  # }}}


def _cache_individual_climatologies(ds, cacheInfo, printProgress,
                                    yearsPerCacheFile, monthValues,
                                    calendar):  # {{{
    '''
    Cache individual climatologies for later aggregation.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    for cacheIndex, info in enumerate(cacheInfo):
        outputFileClimo, done, yearString = info
        if done:
            continue
        dsYear = ds.where(ds.cacheIndices == cacheIndex, drop=True)

        if printProgress:
            print('     {}'.format(yearString))

        totalDays = dsYear.daysInMonth.sum(dim='Time').values

        monthCount = dsYear.dims['Time']

        climatology = compute_climatology(dsYear, monthValues, calendar,
                                          maskVaries=False)

        climatology.attrs['totalDays'] = totalDays
        climatology.attrs['totalMonths'] = monthCount
        climatology.attrs['fingerprintClimo'] = fingerprint_generator()

        write_netcdf(climatology, outputFileClimo)
        climatology.close()

    # }}}


def _cache_aggregated_climatology(startYearClimo, endYearClimo, cachePrefix,
                                  printProgress, monthValues,
                                  cacheInfo):  # {{{
    '''
    Cache aggregated climatology from individual climatologies.
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    yearString, fileSuffix = _get_year_string(startYearClimo, endYearClimo)
    outputFileClimo = '{}_{}.nc'.format(cachePrefix, fileSuffix)

    done = False
    if len(cacheInfo) == 0:
        climatology = None
        done = True

    if os.path.exists(outputFileClimo):
        # already cached
        climatology = None
        try:
            climatology = xr.open_dataset(outputFileClimo)

        except IOError:
            # assuming the cache file is corrupt, so deleting it.
            print('Warning: Deleting cache file {}, which appears to have '
                  'been corrupted.'.format(outputFileClimo))
            os.remove(outputFileClimo)

        if len(cacheInfo) == 1 and outputFileClimo == cacheInfo[0][0]:
            # theres only one cache file and it already has the same name
            # as the aggregated file so no need to aggregate
            done = True

        elif climatology is not None:
            monthsIfDone = (
                endYearClimo - startYearClimo + 1) * len(monthValues)
            if climatology.attrs['totalMonths'] == monthsIfDone:
                # also complete, so we can move on
                done = True
            else:
                climatology.close()

    if not done:
        if printProgress:
            print('   Computing aggregated climatology '
                  '{}...'.format(yearString))

        first = True
        for cacheIndex, info in enumerate(cacheInfo):
            inFileClimo = info[0]
            ds = xr.open_dataset(inFileClimo)
            days = ds.attrs['totalDays']
            months = ds.attrs['totalMonths']
            if first:
                totalDays = days
                totalMonths = months
                climatology = ds * days
                first = False
            else:
                totalDays += days
                totalMonths += months
                climatology = climatology + ds * days

            ds.close()
        climatology = climatology / totalDays

        climatology.attrs['totalDays'] = totalDays
        climatology.attrs['totalMonths'] = totalMonths
        climatology.attrs['fingerprintClimo'] = fingerprint_generator()

        write_netcdf(climatology, outputFileClimo)

    return climatology  # }}}


def _get_year_string(startYear, endYear):
    if startYear == endYear:
        yearString = '{:04d}'.format(startYear)
        fileSuffix = 'year{}'.format(yearString)
    else:
        yearString = '{:04d}-{:04d}'.format(startYear, endYear)
        fileSuffix = 'years{}'.format(yearString)

    return yearString, fileSuffix


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
