"""
Functions for creating climatologies from monthly time series data

Authors
-------
Xylar Asay-Davis

Last Modified
-------------
04/08/2017
"""

import xarray as xr
import os
import numpy
import netCDF4
import warnings

from ..constants import constants

from ..timekeeping.utility import days_to_datetime

from ..io.utility import build_config_full_path, make_directories

from ..interpolation import Remapper
from ..grid import MpasMeshDescriptor, LatLonGridDescriptor


def get_mpas_remapper(config, meshFileName):  # {{{
    """
    Given config options and the name of an MPAS mesh file, returns a
    ``Remapper`` object that can be used to remap from MPAS files or data sets
    to corresponding data sets on the comparison grid.

    If necessary, creates the mapping file containing weights and indices
    needed to perform remapping.

    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
        Contains configuration options

    meshFileName : str
        The path of the file containing the source MPAS mesh

    Returns
    -------
    remapper : ``Remapper`` object
        A remapper that can be used to remap files or data sets from the MPAS
        mesh to the comparison grid.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/06/2017
    """

    climSection = 'climatology'

    method = config.getWithDefault(climSection, 'mpasInterpolationMethod',
                                   'bilinear')

    lat, lon, comparisonLatRes, comparisonLonRes = \
        _get_comparison_grid(config)

    mappingFileOption = 'mpasMappingFile'
    if config.has_option(climSection, mappingFileOption):
        # a mapping file was supplied, so we'll use that name
        mpasMappingFileName = config.get(climSection, mappingFileOption)
    else:
        # we need to build the path to the mapping file and an appropriate
        # file name
        mappingSubdirectory = build_config_full_path(config, 'output',
                                                     'mappingSubdirectory')

        make_directories(mappingSubdirectory)

        meshName = config.get('input', 'mpasMeshName')

        mpasMappingFileName = '{}/map_{}_to_{}x{}degree_{}.nc'.format(
            mappingSubdirectory, meshName, comparisonLatRes, comparisonLonRes,
            method)

        config.set(climSection, mappingFileOption, mpasMappingFileName)

    sourceDescriptor = MpasMeshDescriptor(meshFileName)
    destinationDescriptor = LatLonGridDescriptor()
    destinationDescriptor.create(lat, lon, units='degrees')

    remapper = Remapper(sourceDescriptor, destinationDescriptor,
                        mpasMappingFileName)

    remapper.build_mapping_file(method=method)

    return remapper  # }}}


def get_observations_remapper(config, componentName, fieldName,
                              gridFileName, latVarName='lat',
                              lonVarName='lon'):  # {{{
    """
    Given config options, the name of the component being analyzed and a
    grid file containing 1D lat and lon arrays, returns a ``Remapper`` object
    that can be used to remap from observations files or data sets to
    corresponding data sets on the comparison grid.

    If necessary, creates the mapping file containing weights and indices
    needed to perform remapping.

    If no remapping is needed, returns ``None``


    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
        Contains configuration options

    componentName : ``{'ocean', 'seaIce'}``
        Name of the component, used to look up climatology and observation
        options

    fieldName : str
        Name of the field being mapped, used to give each set of
        observation weights a unique name.

    gridFileName : str
        The path of the file containing the source lat-lon grid

    latVarName, lonVarName : str, optional
        The name of the latitude and longitude variables in the source grid
        file

    Returns
    -------
    remapper : ``Remapper`` object
        A remapper that can be used to remap files or data sets from the
        observations grid to the comparison grid.  Returns ``None`` if no
        mapping is required.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/06/2017
    """

    mappingFileOption = '{}ClimatologyMappingFile'.format(fieldName)
    obsSection = '{}Observations'.format(componentName)

    outLat, outLon, comparisonLatRes, comparisonLonRes = \
        _get_comparison_grid(config)

    method = config.getWithDefault(obsSection, 'interpolationMethod',
                                   'bilinear')

    if config.has_option(obsSection, mappingFileOption):
        obsMappingFileName = config.get(obsSection, mappingFileOption)
    else:

        (gridName, matchesComparison) = _get_grid_name(gridFileName,
                                                       latVarName,
                                                       lonVarName,
                                                       comparisonLatRes,
                                                       comparisonLonRes)

        if matchesComparison:
            # no need to remap the observations
            obsMappingFileName = None
        else:
            mappingSubdirectory = build_config_full_path(config, 'output',
                                                         'mappingSubdirectory')

            make_directories(mappingSubdirectory)

            obsMappingFileName = \
                '{}/map_obs_{}_{}_to_{}x{}degree_{}.nc'.format(
                    mappingSubdirectory, fieldName, gridName,
                    comparisonLatRes, comparisonLonRes, method)

            config.set(obsSection, mappingFileOption, obsMappingFileName)

    if obsMappingFileName is None:
        remapper = None
    else:
        sourceDescriptor = LatLonGridDescriptor()
        sourceDescriptor.read(gridFileName, latVarName=latVarName,
                              lonVarName=lonVarName)
        destinationDescriptor = LatLonGridDescriptor()
        destinationDescriptor.create(outLat, outLon, units='degrees')

        remapper = Remapper(sourceDescriptor, destinationDescriptor,
                            obsMappingFileName)
        remapper.build_mapping_file(method=method)

    return remapper  # }}}


def get_mpas_climatology_file_names(config, fieldName, monthNames):  # {{{
    """
    Given config options, the name of a field and a string identifying the
    months in a seasonal climatology, returns the full path for MPAS
    climatology files before and after regridding.

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    fieldName : str
        Name of the field being mapped, used as a prefix for the climatology
        file name.

    monthNames : str
        A string identifying the months in a seasonal climatology (e.g. 'JFM')

    Returns
    -------
    climatologyFileName : str
        The absolute path to a file where the climatology should be stored
        before regridding.

    climatologyPrefix : str
        The prfix including absolute path for climatology cache files before
        regridding.

    regriddedFileName : str
        The absolute path to a file where the climatology should be stored
        after regridding.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/03/2017
    """

    climSection = 'climatology'
    startYear = config.getint(climSection, 'startYear')
    endYear = config.getint(climSection, 'endYear')

    meshName = config.get('input', 'mpasMeshName')

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    climatologyDirectory = build_config_full_path(
        config, 'output', 'mpasClimatologySubdirectory')

    regriddedDirectory = build_config_full_path(
        config, 'output', 'mpasRegriddedClimSubdirectory')

    make_directories(regriddedDirectory)
    make_directories(climatologyDirectory)

    climatologyPrefix = '{}/{}_{}_{}'.format(climatologyDirectory, fieldName,
                                             meshName, monthNames)

    yearString, fileSuffix = _get_year_string(startYear, endYear)
    climatologyFileName = '{}_{}.nc'.format(climatologyPrefix, fileSuffix)
    regriddedFileName = \
        '{}/{}_{}_to_{}x{}degree_{}_{}.nc'.format(
            regriddedDirectory, fieldName, meshName, comparisonLatRes,
            comparisonLonRes, monthNames, fileSuffix)

    return (climatologyFileName, climatologyPrefix, regriddedFileName)  # }}}


def get_observation_climatology_file_names(config, fieldName, monthNames,
                                           componentName, gridFileName,
                                           latVarName='lat',
                                           lonVarName='lon'):  # {{{
    """
    Given config options, the name of a field and a string identifying the
    months in a seasonal climatology, returns the full path for observation
    climatology files before and after regridding.

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    fieldName : str
        Name of the field being mapped, used as a prefix for the climatology
        file name.

    monthNames : str
        A string identifying the months in a seasonal climatology (e.g. 'JFM')

    gridFileName : str
        The path of the file containing the source lat-lon grid

    latVarName, lonVarName : str, optional
        The name of the latitude and longitude variables in the source grid
        file

    Returns
    -------
    climatologyFileName : str
        The absolute path to a file where the climatology should be stored
        before regridding.

    regriddedFileName : str
        The absolute path to a file where the climatology should be stored
        after regridding.


    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/03/2017
    """

    climSection = 'climatology'
    obsSection = '{}Observations'.format(componentName)

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    climatologyDirectory = build_config_full_path(
        config=config, section='output',
        relativePathOption='climatologySubdirectory',
        relativePathSection=obsSection)

    regriddedDirectory = build_config_full_path(
        config=config, section='output',
        relativePathOption='regriddedClimSubdirectory',
        relativePathSection=obsSection)

    (gridName, matchesComparison) = _get_grid_name(gridFileName,
                                                   latVarName,
                                                   lonVarName,
                                                   comparisonLatRes,
                                                   comparisonLonRes)

    climatologyFileName = '{}/{}_{}_{}.nc'.format(
        climatologyDirectory, fieldName, gridName, monthNames)
    regriddedFileName = '{}/{}_{}_to_{}x{}degree_{}.nc'.format(
        regriddedDirectory, fieldName, gridName, comparisonLatRes,
        comparisonLonRes, monthNames)

    make_directories(climatologyDirectory)

    if not matchesComparison:
        make_directories(regriddedDirectory)

    return (climatologyFileName, regriddedFileName)  # }}}


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

    Last Modified
    -------------
    04/08/2017
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

    Last Modified
    -------------
    04/08/2017
    """

    ds = add_years_months_days_in_month(ds, calendar)

    mask = xr.zeros_like(ds.month, bool)

    for month in monthValues:
        mask = xr.ufuncs.logical_or(mask, ds.month == month)

    climatologyMonths = ds.where(mask, drop=True)

    climatology = _compute_masked_mean(climatologyMonths, maskVaries)

    return climatology  # }}}


def cache_climatologies(ds, monthValues, config, cachePrefix, calendar,
                        printProgress=False):  # {{{
    '''
    Cache NetCDF files for each year of an annual climatology, and then use
    the cached files to compute a climatology for the full range of years.
    The start and end years of the climatology are taken from ``config``, and
    are updated in ``config`` if the data set ``ds`` doesn't contain this
    full range.

    Note: only works with climatologies where the mask (locations of ``NaN``
    values) doesn't vary with time.

    Parameters
    ----------
    ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01

    monthValues : int or array-like of ints
        A single month or an array of months to be averaged together

    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    cachePrefix :  str
        The file prefix (including path) to which the year (or years) will be
        appended as cache files are stored

    calendar : ``{'gregorian', 'gregorian_noleap'}``
        The name of one of the calendars supported by MPAS cores, used to
        determine ``year`` and ``month`` from ``Time`` coordinate

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
    04/11/2017
    '''
    startYearClimo = config.getint('climatology', 'startYear')
    endYearClimo = config.getint('climatology', 'endYear')
    yearsPerCacheFile = config.getint('climatology', 'yearsPerCacheFile')

    if printProgress:
        print '   Computing and caching climatologies covering {}-year ' \
              'spans...'.format(yearsPerCacheFile)

    ds = add_years_months_days_in_month(ds, calendar)

    cacheInfo, cacheIndices = _setup_climatology_caching(ds, startYearClimo,
                                                         endYearClimo,
                                                         yearsPerCacheFile,
                                                         cachePrefix,
                                                         monthValues)

    ds = ds.copy()
    ds.coords['cacheIndices'] = ('Time', cacheIndices)

    # compute and store each cache file with interval yearsPerCacheFile
    _cache_individual_climatologies(ds, cacheInfo, printProgress,
                                    yearsPerCacheFile, monthValues,
                                    calendar)

    # compute the aggregate climatology
    climatology = _cache_aggregated_climatology(startYearClimo, endYearClimo,
                                                cachePrefix, printProgress,
                                                monthValues, cacheInfo)

    return climatology  # }}}


def update_start_end_year(ds, config, calendar):  # {{{
    """
    Given a monthly climatology, compute a seasonal climatology weighted by
    the number of days in each month (on the no-leap-year calendar).

    Parameters
    ----------
    ds : instance of xarray.Dataset
        A data set from which start and end years will be determined

    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    calendar : {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

    Returns
    -------
    changed : bool
        Whether the start and end years were changed

    startYear, endYear : int
        The start and end years of the data set

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/25/2017
    """
    requestedStartYear = config.getint('climatology', 'startYear')
    requestedEndYear = config.getint('climatology', 'endYear')

    startYear = days_to_datetime(ds.Time.min().values, calendar=calendar).year
    endYear = days_to_datetime(ds.Time.max().values,  calendar=calendar).year
    changed = False
    if startYear != requestedStartYear or endYear != requestedEndYear:
        message = "climatology start and/or end year different from " \
                  "requested\n" \
                  "requestd: {:04d}-{:04d}\n" \
                  "actual:   {:04d}-{:04d}\n".format(requestedStartYear,
                                                     requestedEndYear,
                                                     startYear,
                                                     endYear)
        warnings.warn(message)
        config.set('climatology', 'startYear', str(startYear))
        config.set('climatology', 'endYear', str(endYear))
        changed = True

    return changed, startYear, endYear  # }}}


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

    Last Modified
    -------------
    04/08/2017
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
                                climatologyFileName, regriddedFileName,
                                remapper):  # {{{
    """
    Given a field in a climatology data set, use the ``remapper`` to regrid
    horizontal dimensions of all fields, write the results to an output file,
    and return the regridded data set.

    Note that ``climatologyFileName`` and ``regriddedFileName`` will be
    overwritten if they exist, so if this behavior is not desired, the calling
    code should skip this call if the files exist and simply load the contents
    of ``regriddedFileName``.

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
        before regridding (if using ncremap).

    regriddedFileName : str
        The name of the output file to which the regridded data set should
        be written.

    remapper : ``Remapper`` object
        A remapper that can be used to remap files or data sets from the
        observations grid to the comparison grid.  Returns ``None`` if no
        mapping is required.

    Returns
    -------
    remappedClimatology : ``xarray.DataSet`` or ``xarray.DataArray`` object
        A data set containing the remapped climatology

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/06/2017
    """
    useNcremap = config.getboolean('climatology', 'useNcremap')

    if useNcremap:
        useNcremap = config.getboolean('climatology', 'useNcremap')
        if not os.path.exists(climatologyFileName):
            climatologyDataSet.to_netcdf(climatologyFileName)
        remapper.remap_file(inFileName=climatologyFileName,
                            outFileName=regriddedFileName,
                            overwrite=True)
        remappedClimatology = xr.open_dataset(regriddedFileName)
    else:
        renormalizationThreshold = config.getfloat('climatology',
                                                   'renormalizationThreshold')

        remappedClimatology = remapper.remap(climatologyDataSet,
                                             renormalizationThreshold)
    return remappedClimatology  # }}}


def _compute_masked_mean(ds, maskVaries):  # {{{
    '''
    Compute the time average of data set, masked out where the variables in ds
    are NaN and, if ``maskVaries == True``, weighting by the number of days
    used to compute each monthly mean time in ds.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/08/2017
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
        dsWeightedSum.compute()

        weights = ds_to_weights(ds)
        weights.compute()

        weightSum = (weights * ds.daysInMonth).sum(dim='Time')
        weightSum.compute()

        timeMean = dsWeightedSum / weightSum.where(weightSum > 0.)
        timeMean.compute()
    else:
        days = ds.daysInMonth.sum(dim='Time')
        days.compute()

        dsWeightedSum = (ds * ds.daysInMonth).sum(dim='Time', keep_attrs=True)
        dsWeightedSum.compute()

        timeMean = dsWeightedSum / days.where(days > 0.)
        timeMean.compute()

    return timeMean  # }}}


def _get_comparison_grid(config):  # {{{
    '''
    Returns the latitude and longitude arrays defining the comparison grid

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    Returns
    -------
    lon, lat : 1D numpy arrays
        latitude and longitude arrays defining the comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/08/2017
    '''

    climSection = 'climatology'

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    nLat = int((constants.latmax-constants.latmin)/comparisonLatRes)+1
    nLon = int((constants.lonmax-constants.lonmin)/comparisonLonRes)+1
    lat = numpy.linspace(constants.latmin, constants.latmax, nLat)
    lon = numpy.linspace(constants.lonmin, constants.lonmax, nLon)

    return (lat, lon, comparisonLatRes, comparisonLonRes)  # }}}


def _get_grid_name(gridFileName, latVarName, lonVarName, comparisonLatRes,
                   comparisonLonRes):  # {{{
    '''
    Given a grid file with given lat and lon variable names, finds the
    resolution of the grid and generates a grid name.  Given comparison
    lat and lon resolution, determines if the grid matches the comparison grid.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/08/2017
    '''
    inFile = netCDF4.Dataset(gridFileName, 'r')

    # Get info from input file
    inLat = numpy.array(inFile.variables[latVarName][:], float)
    inLon = numpy.array(inFile.variables[lonVarName][:], float)
    inDLat = inLat[1]-inLat[0]
    inDLon = inLon[1]-inLon[0]
    if 'degree' in inFile.variables[latVarName].units:
        inUnits = 'degree'
    else:
        inUnits = 'radian'
    inFile.close()
    gridName = '{}x{}{}'.format(abs(inDLat), abs(inDLon), inUnits)

    matchesComparison = ((inUnits == 'degree') and
                         (comparisonLatRes == inDLat) and
                         (comparisonLonRes == inDLon) and
                         (inLat[0]-0.5*inDLat == constants.latmin) and
                         (inLon[0]-0.5*inDLon == constants.lonmin))

    return (gridName, matchesComparison)  # }}}


def _setup_climatology_caching(ds, startYearClimo, endYearClimo,
                               yearsPerCacheFile, cachePrefix,
                               monthValues):  # {{{
    '''
    Determine which cache files already exist, which are incomplete and which
    years are present in each cache file (whether existing or to be created).

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/08/2017
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

    Last Modified
    -------------
    04/19/2017
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

        climatology.to_netcdf(outputFileClimo)
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

    Last Modified
    -------------
    04/19/2017
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

        print 'DEBUG:', outputFileClimo
        climatology.to_netcdf(outputFileClimo)

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
