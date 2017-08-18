"""
Functions for creating climatologies from monthly time series data

Authors
-------
Xylar Asay-Davis
"""

import xarray as xr
import os
import numpy
import warnings

from ..constants import constants

from ..timekeeping.utility import days_to_datetime

from ..io.utility import build_config_full_path, make_directories, \
    fingerprint_generator
from ..io import write_netcdf

from ..interpolation import Remapper
from ..grid import LatLonGridDescriptor, ProjectionGridDescriptor


def get_remapper(config, sourceDescriptor, comparisonDescriptor,
                 mappingFilePrefix, method):  # {{{
    """
    Given config options and descriptions of the source and comparison grids,
    returns a ``Remapper`` object that can be used to remap from source files
    or data sets to corresponding data sets on the comparison grid.

    If necessary, creates the mapping file containing weights and indices
    needed to perform remapping.

    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
        Contains configuration options

    sourceDescriptor : ``MeshDescriptor`` subclass object
        A description of the source mesh or grid

    comparisonDescriptor : ``MeshDescriptor`` subclass object
        A description of the comparison grid

    mappingFilePrefix : str
        A prefix to be prepended to the mapping file name

    method : {'bilinear', 'neareststod', 'conserve'}
        The method of interpolation used.

    Returns
    -------
    remapper : ``Remapper`` object
        A remapper that can be used to remap files or data sets from the source
        grid or mesh to the comparison grid.

    Authors
    -------
    Xylar Asay-Davis
    """

    mappingFileName = None

    if not _matches_comparison(sourceDescriptor, comparisonDescriptor):
        # we need to remap because the grids don't match

        mappingBaseName = '{}_{}_to_{}_{}.nc'.format(
                mappingFilePrefix,
                sourceDescriptor.meshName,
                comparisonDescriptor.meshName,
                method)

        if config.has_option('input', 'mappingDirectory'):
            # a mapping directory was supplied, so we'll see if there's
            # a mapping file there that we can use
            mappingSubdirectory = config.get('input', 'mappingDirectory')
            mappingFileName = '{}/{}'.format(mappingSubdirectory,
                                             mappingBaseName)
            if not os.path.exists(mappingFileName):
                # nope, looks like we still need to make a mapping file
                mappingFileName = None

        if mappingFileName is None:
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

    remapper.build_mapping_file(method=method)

    return remapper  # }}}


def get_observation_climatology_file_names(config, fieldName, monthNames,
                                           componentName, remapper):  # {{{
    """
    Given config options, the name of a field and a string identifying the
    months in a seasonal climatology, returns the full path for observation
    climatology files before and after remapping.

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    fieldName : str
        Name of the field being mapped, used as a prefix for the climatology
        file name.

    monthNames : str
        A string identifying the months in a seasonal climatology (e.g. 'JFM')

    remapper : ``Remapper`` object
        A remapper that used to remap files or data sets from the
        observation grid to a comparison grid

    Returns
    -------
    climatologyFileName : str
        The absolute path to a file where the climatology should be stored
        before remapping.

    remappedFileName : str
        The absolute path to a file where the climatology should be stored
        after remapping.

    Authors
    -------
    Xylar Asay-Davis
    """

    obsSection = '{}Observations'.format(componentName)

    climatologyDirectory = build_config_full_path(
        config=config, section='output',
        relativePathOption='climatologySubdirectory',
        relativePathSection=obsSection)

    remappedDirectory = build_config_full_path(
        config=config, section='output',
        relativePathOption='remappedClimSubdirectory',
        relativePathSection=obsSection)

    obsGridName = remapper.sourceDescriptor.meshName
    comparisonGridName = remapper.destinationDescriptor.meshName

    climatologyFileName = '{}/{}_{}_{}.nc'.format(
        climatologyDirectory, fieldName, obsGridName, monthNames)
    remappedFileName = '{}/{}_{}_to_{}_{}.nc'.format(
        remappedDirectory, fieldName, obsGridName, comparisonGridName,
        monthNames)

    make_directories(climatologyDirectory)

    if not _matches_comparison(remapper.sourceDescriptor,
                               remapper.destinationDescriptor):
        make_directories(remappedDirectory)

    return (climatologyFileName, remappedFileName)  # }}}


def compute_monthly_climatology(ds, calendar=None, maskVaries=True):  # {{{
    """
    Compute monthly climatologies from a data set.  The mean is weighted but
    the number of days in each month of the data set, ignoring values masked
    out with NaNs.  If the month coordinate is not present, a data array
    ``month`` will be added based on ``Time`` and the provided calendar.

    Parameters
    ----------
    ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01 or ``month`` coordinate

    calendar : ``{'gregorian', 'gregorian_noleap'}``, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate, so must be supplied if
        ``ds`` does not already have a ``month`` coordinate or data array

    maskVaries: bool, optional
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

    Authors
    -------
    Xylar Asay-Davis
    """

    def compute_one_month_climatology(ds):
        monthValues = list(ds.month.values)
        return compute_climatology(ds, monthValues, calendar, maskVaries)

    ds = add_years_months_days_in_month(ds, calendar)

    monthlyClimatology = \
        ds.groupby('month').apply(compute_one_month_climatology)

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
    ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01 or ``month`` coordinate

    monthValues : int or array-like of ints
        A single month or an array of months to be averaged together

    calendar : ``{'gregorian', 'gregorian_noleap'}``, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate, so must be supplied if
        ``ds`` does not already have a ``month`` coordinate or data array

    maskVaries: bool, optional
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

    Authors
    -------
    Xylar Asay-Davis
    """

    ds = add_years_months_days_in_month(ds, calendar)

    mask = xr.zeros_like(ds.month, bool)

    for month in monthValues:
        mask = xr.ufuncs.logical_or(mask, ds.month == month)

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

    calendar : ``{'gregorian', 'gregorian_noleap'}``, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``year`` and ``month`` from ``Time`` coordinate

    Returns
    -------
    ds : object of same type as ``ds``
        The data set with ``year``, ``month`` and ``daysInMonth`` data arrays
        added (if not already present)

    Authors
    -------
    Xylar Asay-Davis
    '''

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
                message = 'The MPAS run used the Gregorian calendar but ' \
                          'does not appear to have\n' \
                          'supplied start and end times.  Climatologies ' \
                          'will be computed with\n' \
                          'month durations ignoring leap years.'
                warnings.warn(message)

            daysInMonth = numpy.array([constants.daysInMonth[month-1] for
                                       month in ds.month.values], float)
            ds.coords['daysInMonth'] = ('Time', daysInMonth)

    return ds  # }}}


def remap_and_write_climatology(config, climatologyDataSet,
                                climatologyFileName, remappedFileName,
                                remapper):  # {{{
    """
    Given a field in a climatology data set, use the ``remapper`` to regrid
    horizontal dimensions of all fields, write the results to an output file,
    and return the remapped data set.

    Note that ``climatologyFileName`` and ``remappedFileName`` will be
    overwritten if they exist, so if this behavior is not desired, the calling
    code should skip this call if the files exist and simply load the contents
    of ``remappedFileName``.

    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
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

    remapper : ``Remapper`` object
        A remapper that can be used to remap files or data sets to a
        comparison grid.

    Returns
    -------
    remappedClimatology : ``xarray.DataSet`` or ``xarray.DataArray`` object
        A data set containing the remapped climatology

    Authors
    -------
    Xylar Asay-Davis
    """
    useNcremap = config.getboolean('climatology', 'useNcremap')

    if remapper.mappingFileName is None:
        # no remapping is needed
        remappedClimatology = climatologyDataSet
    else:
        if useNcremap:
            if not os.path.exists(climatologyFileName):
                write_netcdf(climatologyDataSet, climatologyFileName)
            remapper.remap_file(inFileName=climatologyFileName,
                                outFileName=remappedFileName,
                                overwrite=True)
            remappedClimatology = xr.open_dataset(remappedFileName)
        else:
            renormalizationThreshold = config.getfloat(
                'climatology', 'renormalizationThreshold')

            remappedClimatology = remapper.remap(climatologyDataSet,
                                                 renormalizationThreshold)
            write_netcdf(remappedClimatology, remappedFileName)
    return remappedClimatology  # }}}


def _compute_masked_mean(ds, maskVaries):  # {{{
    '''
    Compute the time average of data set, masked out where the variables in ds
    are NaN and, if ``maskVaries == True``, weighting by the number of days
    used to compute each monthly mean time in ds.

    Authors
    -------
    Xylar Asay-Davis
    '''
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

    Authors
    -------
    Xylar Asay-Davis
    '''

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

    Authors
    -------
    Xylar Asay-Davis
    '''

    cacheInfo = []

    cacheIndices = -1*numpy.ones(ds.dims['Time'], int)
    monthsInDs = ds.month.values
    yearsInDs = ds.year.values

    # figure out which files to load and which years go in each file
    for firstYear in range(startYearClimo, endYearClimo+1, yearsPerCacheFile):
        years = range(firstYear, firstYear+yearsPerCacheFile)

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
                message = 'Deleting cache file {}, which appears to have ' \
                          'been corrupted.'.format(outputFileClimo)
                warnings.warn(message)
                os.remove(outputFileClimo)

            monthsIfDone = len(monthValues)*len(years)
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

    Authors
    -------
    Xylar Asay-Davis
    '''

    for cacheIndex, info in enumerate(cacheInfo):
        outputFileClimo, done, yearString = info
        if done:
            continue
        dsYear = ds.where(ds.cacheIndices == cacheIndex, drop=True)

        if printProgress:
            print '     {}'.format(yearString)

        totalDays = dsYear.daysInMonth.sum(dim='Time').values

        monthCount = dsYear.dims['Time']

        climatology = compute_climatology(dsYear,  monthValues, calendar,
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

    Authors
    -------
    Xylar Asay-Davis
    '''

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
            message = 'Deleting cache file {}, which appears to have ' \
                      'been corrupted.'.format(outputFileClimo)
            warnings.warn(message)
            os.remove(outputFileClimo)

        if len(cacheInfo) == 1 and outputFileClimo == cacheInfo[0][0]:
            # theres only one cache file and it already has the same name
            # as the aggregated file so no need to aggregate
            done = True

        elif climatology is not None:
            monthsIfDone = (endYearClimo-startYearClimo+1)*len(monthValues)
            if climatology.attrs['totalMonths'] == monthsIfDone:
                # also complete, so we can move on
                done = True
            else:
                climatology.close()

    if not done:
        if printProgress:
            print '   Computing aggregated climatology ' \
                  '{}...'.format(yearString)

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
