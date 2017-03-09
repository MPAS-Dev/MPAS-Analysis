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

from ..mpas_xarray import mpas_xarray
from ..constants import constants

from ..timekeeping.utility import days_to_datetime

from ..io.utility import buildConfigFullPath

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
        mappingSubdirectory = buildConfigFullPath(config, 'output',
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
            mappingSubdirectory = buildConfigFullPath(config, 'output',
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

    climatologyDirectory = buildConfigFullPath(config, 'output',
                                               'mpasClimatologySubdirectory')

    regriddedDirectory = buildConfigFullPath(config, 'output',
                                             'mpasRegriddedClimSubdirectory')
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

    climatologyDirectory = buildConfigFullPath(
        config=config, section='output',
        relativePathOption='climatologySubdirectory',
        relativePathSection=obsSection)

    regriddedDirectory = buildConfigFullPath(
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


def compute_monthly_climatology(ds, calendar):
    """
    Compute a monthly climatology data set from a data set with Time expressed
    as days since 0001-01-01 with the given calendar.

    Parameters
    ----------
    ds : instance of xarray.DataSet
        A data set with a 'Time' coordinate expressed as days since
        0001-01-01

    calendar: {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

    Returns
    -------
    monthlyClimatology : instance of xarray.DataSet
        A data set with a new 'month' coordinate,
        containing monthly climatologies of all variables in ds

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/25/2017
    """
    months = [date.month for date in days_to_datetime(ds.Time,
                                                      calendar=calendar)]

    ds.coords['month'] = ('Time', months)
    monthlyClimatology = ds.groupby('month').mean('Time')
    return monthlyClimatology


def compute_seasonal_climatology(monthlyClimatology, monthValues,
                                 variableName):
    """
    Given a monthly climatology, compute a seasonal climatology weighted by
    the number of days in each month (on the no-leap-year calendar).

    Parameters
    ----------
    monthlyClimatology : instance of xarray.DataSet
        A data set containing a monthly climatology

    monthValues : int or array-like of ints
        A single month or an array of months to be averaged together
        before interpolation.

    Returns
    -------
    seasonalClimatology : instance of xarray.DataSet
        A data set containing the seasonal climatology

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/27/2017
    """

    # select only the desired field from the obs. data set
    monthlyClimatology = mpas_xarray.subset_variables(
        monthlyClimatology, variableList=[variableName])

    # Make a DataArray with the number of days in each month
    daysInMonth = xr.DataArray(numpy.array(constants.daysInMonth, float),
                               coords=[monthlyClimatology.month],
                               name='daysInMonth')

    daysInMonth = daysInMonth.sel(month=monthValues)
    monthlyClimatology = monthlyClimatology.sel(month=monthValues)

    varArray = monthlyClimatology[variableName]
    mask = xr.DataArray(~numpy.isnan(varArray.values),
                        coords=varArray.coords,
                        name='mask')
    seasonalClimatology = (monthlyClimatology * daysInMonth).sum(
        dim='month', keep_attrs=True)

    days = (mask * daysInMonth).sum(dim='month')
    seasonalClimatology /= days.where(days > 0.)

    return seasonalClimatology


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
