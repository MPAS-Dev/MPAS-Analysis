"""
Functions for creating climatologies from monthly time series data

Authors
-------
Xylar Asay-Davis

Last Modified
-------------
03/04/2017
"""

import xarray as xr
import os
import numpy
import netCDF4
import warnings

from ..constants import constants

from ..timekeeping.utility import days_to_datetime

from ..io.utility import build_config_full_path

from ..interpolation import interpolate


def write_mpas_mapping_file(config, meshFileName):
    """
    Given config options, the name of the component being analyzed and an
    MPAS mesh file, either finds an existing MPAS-to-comparison-grid mapping
    file or creates a new mapping file.

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    meshFileName : str
        The path of the file containing the source MPAS mesh

    Returns
    -------
    mpasMappingFileName : str
        The absolute path to an existing mapping file or the location
        at which one was created.  The mapping file can be used to
        interpolate between MPAS meshes and the climatology comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/04/2017
    """

    climSection = 'climatology'

    method = config.getWithDefault(climSection, 'mpasInterpolationMethod',
                                   'bilinear')

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    (lat, lon) = _get_comparison_lat_lon(comparisonLatRes, comparisonLonRes)

    overwriteMapping = config.getWithDefault(climSection,
                                             'overwriteMapping',
                                             False)

    mappingFileOption = 'mpasMappingFile'
    if config.has_option(climSection, mappingFileOption):
        # a mapping file was supplied, so we'll use that name
        mpasMappingFileName = config.get(climSection, mappingFileOption)
    else:
        # we need to build the path to the mapping file and an appropriate
        # file name
        mappingSubdirectory = build_config_full_path(config, 'output',
                                                     'mappingSubdirectory')
        try:
            os.makedirs(mappingSubdirectory)
        except OSError:
            pass

        meshName = config.get('input', 'mpasMeshName')

        mpasMappingFileName = '{}/map_{}_to_{}x{}degree_{}.nc'.format(
            mappingSubdirectory, meshName, comparisonLatRes, comparisonLonRes,
            method)

        config.set(climSection, mappingFileOption, mpasMappingFileName)

    interpolate.build_remap_weights(sourceFileName=meshFileName,
                                    outWeightFileName=mpasMappingFileName,
                                    sourceFileType='mpas',
                                    destinationLat=lat,
                                    destinationLon=lon,
                                    desitnationUnits='degrees',
                                    method=method,
                                    overwrite=overwriteMapping)

    return mpasMappingFileName


def write_observations_mapping_file(config, componentName, fieldName,
                                    gridFileName, latVarName='lat',
                                    lonVarName='lon'):
    """
    Given config options, the name of the component being analyzed and a
    grid file containing 1D lat and lon arrays, either finds an existing
    obs-grid-to-comparison-grid mapping file or creates a new mapping file.

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    componentName : {'ocean', 'seaIce'}
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
    obsMappingFileName : str
        The absolute path to a mapping file (or the location
        at which to create one) for interpolation between MPAS meshes and
        the climatology comparison grid

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/04/2017
    """

    mappingFileOption = '{}ClimatologyMappingFile'.format(fieldName)
    climSection = 'climatology'
    obsSection = '{}Observations'.format(componentName)

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    (outLat, outLon) = _get_comparison_lat_lon(comparisonLatRes,
                                               comparisonLonRes)

    method = config.getWithDefault(obsSection, 'interpolationMethod',
                                   'bilinear')

    overwriteMapping = config.getWithDefault(climSection,
                                             'overwriteMapping',
                                             False)

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

            try:
                os.makedirs(mappingSubdirectory)
            except OSError:
                pass

            obsMappingFileName = \
                '{}/map_obs_{}_{}_to_{}x{}degree_{}.nc'.format(
                    mappingSubdirectory, fieldName, gridName,
                    comparisonLatRes, comparisonLonRes, method)

            config.set(obsSection, mappingFileOption, obsMappingFileName)

    if obsMappingFileName is not None:
        interpolate.build_remap_weights(
            sourceFileName=gridFileName,
            outWeightFileName=obsMappingFileName,
            sourceFileType='latlon',
            sourceLatVarName=latVarName,
            sourceLonVarName=lonVarName,
            destinationLat=outLat,
            destinationLon=outLon,
            desitnationUnits='degrees',
            method=method,
            overwrite=overwriteMapping)

    return obsMappingFileName


def get_mpas_climatology_file_names(config, fieldName, monthNames):
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
    try:
        os.makedirs(regriddedDirectory)
    except OSError:
        pass
    try:
        os.makedirs(climatologyDirectory)
    except OSError:
        pass

    climatologyFileName = '{}/{}_{}_{}_years{:04d}-{:04d}.nc'.format(
        climatologyDirectory, fieldName, meshName, monthNames, startYear,
        endYear)
    regriddedFileName = \
        '{}/{}_{}_to_{}x{}degree_{}_years{:04d}-{:04d}.nc'.format(
            regriddedDirectory, fieldName, meshName, comparisonLatRes,
            comparisonLonRes, monthNames, startYear, endYear)

    return (climatologyFileName, regriddedFileName)


def get_observation_climatology_file_names(config, fieldName, monthNames,
                                           componentName, gridFileName,
                                           latVarName='lat', lonVarName='lon'):
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

    try:
        os.makedirs(climatologyDirectory)
    except OSError:
        pass

    if not matchesComparison:
        try:
            os.makedirs(regriddedDirectory)
        except OSError:
            pass

    return (climatologyFileName, regriddedFileName)


def compute_monthly_climatology(ds, calendar=None):
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

    calendar: ``{'gregorian', 'gregorian_noleap'}``, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate, so must be supplied if
        ``ds`` does not already have a ``month`` coordinate or data array

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
    03/30/2017
    """

    def compute_one_month_climatology(ds):
        monthValues = list(ds.month.values)
        return compute_climatology(ds, monthValues, calendar)

    if 'month' not in ds.coords or 'daysInMonth' not in ds.coords:
        ds = add_months_and_days_in_month(ds, calendar)

    monthlyClimatology = \
        ds.groupby('month').apply(compute_one_month_climatology)

    return monthlyClimatology


def compute_climatology(ds, monthValues, calendar=None):
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

    calendar: ``{'gregorian', 'gregorian_noleap'}``, optional
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate, so must be supplied if
        ``ds`` does not already have a ``month`` coordinate or data array

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
    03/30/2017
    """

    if ('month' not in ds.coords or 'daysInMonth' not in ds.coords):
        ds = add_months_and_days_in_month(ds, calendar)

    mask = xr.zeros_like(ds.month, bool)

    for month in monthValues:
        mask = xr.ufuncs.logical_or(mask, ds.month == month)

    climatologyMonths = ds.where(mask, drop=True)

    climatology = _compute_masked_mean(climatologyMonths)

    return climatology


def update_start_end_year(ds, config, calendar):
    """
    Given a monthly climatology, compute a seasonal climatology weighted by
    the number of days in each month (on the no-leap-year calendar).

    Parameters
    ----------
    ds : instance of xarray.Dataset
        A data set from which start and end years will be determined

    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    calendar: {'gregorian', 'gregorian_noleap'}
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
        warnings.warn("climatology start and/or end year different from requested\n"
                      "requestd: {:04d}-{:04d}\n"
                      "actual:   {:04d}-{:04d}\n".format(requestedStartYear, requestedEndYear, startYear, endYear))
        config.set('climatology', 'startYear', str(startYear))
        config.set('climatology', 'endYear', str(endYear))
        changed = True

    return changed, startYear, endYear


def add_months_and_days_in_month(ds, calendar):
    '''
    Add ``months`` and ``daysInMonth`` as data arrays in ``ds``.  The number
    of days in each month of ``ds`` is computed either using the ``startTime``
    and ``endTime`` if available or assuming ``gregorian_noleap`` calendar and
    ignoring leap years.

    Parameters
    ----------
    ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set with a ``Time`` coordinate expressed as days since
        0001-01-01

    calendar: ``{'gregorian', 'gregorian_noleap'}``
        The name of one of the calendars supported by MPAS cores, used to
        determine ``month`` from ``Time`` coordinate

    Returns
    -------
    ds : object of same type as ``ds``
        The data set with ``month`` and ``daysInMonth`` data arrays added (if
        not already present)

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/29/2017
    """    '''

    ds = ds.copy()

    if 'month' not in ds.coords:
        if calendar is None:
            raise ValueError('calendar must be provided if month coordinate is not in ds')
        months = [date.month for date in days_to_datetime(ds.Time,
                                                          calendar=calendar)]

        ds.coords['month'] = ('Time', months)

    if 'daysInMonth' not in ds.coords:
        if 'startTime' in ds.coords and 'endTime' in ds.coords:
            ds.coords['daysInMonth'] = ds.endTime - ds.startTime
        else:
            if calendar == 'gregorian':
                warnings.warn('The MPAS run used the Gregorian calendar but does not appear to have\n'
                              'supplied start and end times.  Climatologies will be computed with\n'
                              'month durations ignoring leap years.')
            # TODO: support leap years if calendar is 'gregorian'
            daysInMonth = numpy.array([constants.daysInMonth[month-1] for
                                       month in ds.month.values], float)
            ds.coords['daysInMonth'] = ('Time', daysInMonth)

    return ds


def _compute_masked_mean(ds):
    '''
    Compute the time average of data set, masked out where the variables in ds
    are NaN and weighting by the number of days used to compute each monthly
    mean time in ds.
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
            raise TypeError('ds must be an instance of either xarray.Dataset or xarray.DataArray.')

        return weights

    dsWeightedSum = (ds * ds.daysInMonth).sum(dim='Time', keep_attrs=True)

    weights = ds_to_weights(ds)

    weightSum = (weights * ds.daysInMonth).sum(dim='Time')

    timeMean = dsWeightedSum / weightSum.where(weightSum > 0.)

    return timeMean


def _get_comparison_lat_lon(comparisonLatRes, comparisonLonRes):
    '''
    Returns the lat and lon arrays defining the corners of the comparison
    grid.
    '''
    nLat = int((constants.latmax-constants.latmin)/comparisonLatRes)+1
    nLon = int((constants.lonmax-constants.lonmin)/comparisonLonRes)+1
    lat = numpy.linspace(constants.latmin, constants.latmax, nLat)
    lon = numpy.linspace(constants.lonmin, constants.lonmax, nLon)

    return (lat, lon)


def _get_grid_name(gridFileName, latVarName, lonVarName, comparisonLatRes,
                   comparisonLonRes):
    '''
    Given a grid file with given lat and lon variable names, finds the
    resolution of the grid and generates a grid name.  Given comparison
    lat and lon resolution, determines if the grid matches the comparison grid.
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

    return (gridName, matchesComparison)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
