'''
General comparison of 2-d model fields against data.  Currently only supports
sea ice concentration (sic) and sea ice thickness (sit)

Authors
-------
Xylar Asay-Davis, Milena Veneziani

Last Modified
-------------
04/08/2017
'''

import os
import os.path

import numpy.ma as ma
import numpy as np

import netCDF4
import warnings
import xarray as xr

from ..shared.constants import constants

from ..shared.climatology import climatology

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io.utility import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from .utility import setup_sea_ice_task


def seaice_modelvsobs(config, streamMap=None, variableMap=None):
    '''
    Performs analysis of sea-ice properties by comparing with
    previous model results and/or observations.

    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    field : {'sst', 'sss', 'mld'}
        The name of a field to be analyized

    streamMap : dict, optional
        A dictionary of MPAS-O stream names that map to their mpas_analysis
        counterparts.

    variableMap : dict, optional
        A dictionary of MPAS-O variable names that map to their mpas_analysis
        counterparts.

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    04/08/2017
    '''

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, namelistMap, \
        streamMap, variableMap, plotsDirectory, simulationStartTime, \
        restartFileName = setup_sea_ice_task(config)

    # get a list of timeSeriesStatsMonthly output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    fileNames = historyStreams.readpath(streamName, startDate=startDate,
                                        endDate=endDate, calendar=calendar)
    print '\n  Reading files:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(fileNames[0]),
              os.path.basename(fileNames[-1]))
    # Load data
    print '  Load sea-ice data...'
    ds = open_multifile_dataset(fileNames=fileNames,
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=['iceAreaCell',
                                              'iceVolumeCell'],
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    # Compute climatologies (first motnhly and then seasonally)
    print '  Compute seasonal climatologies...'

    changed, startYear, endYear = \
        climatology.update_start_end_year(ds, config, calendar)

    mpasRemapper = climatology.get_mpas_remapper(
        config=config, meshFileName=restartFileName)

    _compute_and_plot_concentration(config, ds, mpasRemapper, calendar)

    _compute_and_plot_thickness(config, ds, mpasRemapper, calendar)


def _compute_and_plot_concentration(config, ds, mpasRemapper, calendar):
    '''
    Given a config file, monthly climatology on the mpas grid, and the data
    necessary to perform horizontal interpolation to a comparison grid,
    computes seasonal climatologies and plots model results, observations
    and biases in sea-ice concentration.

    Parameters
    ----------
    config : an instance of MpasConfigParser

    ds : ``xarray.Dataset`` object
        an xarray data set from which to compute climatologies

    mpasRemapper : ``Remapper`` object
       for remapping from the MPAS mesh to the comparison grid

    calendar: ``{'gregorian', 'gregorian_noleap'}``
        The name of one of the calendars supported by MPAS cores

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    04/08/2017
    '''

    print '  Make ice concentration plots...'

    plotsDirectory = build_config_full_path(config, 'output',
                                            'plotsSubdirectory')
    mainRunName = config.get('runs', 'mainRunName')
    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')
    overwriteMpasClimatology = config.getWithDefault(
        'climatology', 'overwriteMpasClimatology', False)

    overwriteObsClimatology = config.getWithDefault(
        'seaIceObservations', 'overwriteObsClimatology', False)

    subtitle = 'Ice concentration'

    hemisphereSeasons = {'JFM': ('NH', 'Winter'),
                         'JAS': ('NH', 'Summer'),
                         'DJF': ('SH', 'Winter'),
                         'JJA': ('SH', 'Summer')}

    obsFileNames = {}
    regriddedObsFileNames = {}

    buildObsClimatologies = overwriteObsClimatology
    for months in hemisphereSeasons:
        hemisphere, season = hemisphereSeasons[months]
        climFieldName = 'iceConcentration'
        for obsName in ['NASATeam', 'Bootstrap']:
            key = (months, obsName)
            obsFileName = build_config_full_path(
                config, 'seaIceObservations',
                'concentration{}{}_{}'.format(obsName, hemisphere, months))
            obsFieldName = '{}_{}_{}'.format(climFieldName, hemisphere,
                                             obsName)

            if not os.path.isfile(obsFileName):
                raise OSError('Obs file {} not found.'.format(
                    obsFileName))

            (climatologyFileName, regriddedFileName) = \
                climatology.get_observation_climatology_file_names(
                    config=config, fieldName=obsFieldName, monthNames=months,
                    componentName='seaIce', gridFileName=obsFileName,
                    latVarName='t_lat', lonVarName='t_lon')

            obsFileNames[key] = obsFileName
            regriddedObsFileNames[key] = regriddedFileName

            if not os.path.exists(regriddedFileName):
                buildObsClimatologies = True

    for months in hemisphereSeasons:
        hemisphere, season = hemisphereSeasons[months]
        monthValues = constants.monthDictionary[months]
        field = 'iceAreaCell'
        climFieldName = 'iceConcentration'

        # interpolate the model results
        (climatologyFileName, climatologyPrefix, regriddedFileName) = \
            climatology.get_mpas_climatology_file_names(
                    config=config, fieldName=climFieldName,
                    monthNames=months)

        if overwriteMpasClimatology or not os.path.exists(regriddedFileName):
            seasonalClimatology = climatology.cache_climatologies(
                ds, monthValues, config, climatologyPrefix, calendar,
                printProgress=True)
            if seasonalClimatology is None:
                # apparently, there was no data available to create the
                # climatology
                warnings.warn('no data to create sea ice concentration '
                              'climatology for {}'.format(months))
                continue

            remappedClimatology = climatology.remap_and_write_climatology(
                config, seasonalClimatology, climatologyFileName,
                regriddedFileName, mpasRemapper)

        else:

            remappedClimatology = xr.open_dataset(regriddedFileName)

        iceConcentration = remappedClimatology[field].values
        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config,
            'regriddedSeaIceConcThick',
            suffix='ConcResult{}'.format(season))
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config,
            'regriddedSeaIceConcThick',
            suffix='ConcDifference{}'.format(season))

        referenceLongitude = config.getfloat(
            'regriddedSeaIceConcThick',
            'referenceLongitude{}'.format(hemisphere))
        minimumLatitude = config.getfloat(
            'regriddedSeaIceConcThick',
            'minimumLatitude{}'.format(hemisphere))

        # ice concentrations from NASATeam (or Bootstrap) algorithm
        for obsName in ['NASATeam', 'Bootstrap']:
            obsFieldName = 'AICE'

            key = (months, obsName)
            regriddedFileName = regriddedObsFileNames[key]

            if buildObsClimatologies:
                obsFileName = obsFileNames[key]
                obsRemapper = \
                    climatology.get_observations_remapper(
                        config=config, componentName='seaIce',
                        fieldName='seaIce', gridFileName=obsFileName,
                        latVarName='t_lat', lonVarName='t_lon')

                seasonalClimatology = xr.open_dataset(obsFileName)

                if obsRemapper is None:
                    remappedClimatology = seasonalClimatology
                else:
                    remappedClimatology = \
                        climatology.remap_and_write_climatology(
                            config, seasonalClimatology,  climatologyFileName,
                            regriddedFileName, obsRemapper)

            obsIceConcentration = remappedClimatology[obsFieldName].values

            difference = iceConcentration - obsIceConcentration

            title = '{} ({}, years {:04d}-{:04d})'.format(
                subtitle, months, startYear, endYear)
            fileout = '{}/iceconc{}{}_{}_{}_years{:04d}-{:04d}.png'.format(
                plotsDirectory, obsName, hemisphere, mainRunName,
                months, startYear, endYear)
            plot_polar_comparison(
                config,
                lonTarg,
                latTarg,
                iceConcentration,
                obsIceConcentration,
                difference,
                colormapResult,
                colorbarLevelsResult,
                colormapDifference,
                colorbarLevelsDifference,
                title=title,
                fileout=fileout,
                plotProjection=plotProjection,
                latmin=minimumLatitude,
                lon0=referenceLongitude,
                modelTitle=mainRunName,
                obsTitle='Observations (SSM/I {})'.format(obsName),
                diffTitle='Model-Observations',
                cbarlabel='fraction')


def _compute_and_plot_thickness(config, ds,  mpasRemapper, calendar):
    '''
    Given a config file, monthly climatology on the mpas grid, and the data
    necessary to perform horizontal interpolation to a comparison grid,
    computes seasonal climatologies and plots model results, observations
    and biases in sea-ice thickness.

    Parameters
    ----------
    config : an instance of MpasConfigParser

    ds : ``xarray.Dataset`` object
        an xarray data set from which to compute climatologies

    mpasRemapper : ``Remapper`` object
       for remapping from the MPAS mesh to the comparison grid

    calendar: ``{'gregorian', 'gregorian_noleap'}``
        The name of one of the calendars supported by MPAS cores

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    04/08/2017
    '''

    print '  Make ice thickness plots...'

    subtitle = 'Ice thickness'

    plotsDirectory = build_config_full_path(config, 'output',
                                            'plotsSubdirectory')
    mainRunName = config.get('runs', 'mainRunName')
    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')
    overwriteMpasClimatology = config.getWithDefault(
        'climatology', 'overwriteMpasClimatology', False)

    overwriteObsClimatology = config.getWithDefault(
        'seaIceObservations', 'overwriteObsClimatology', False)

    obsFileNames = {}
    regriddedObsFileNames = {}

    # build a list of regridded observations files
    buildObsClimatologies = overwriteObsClimatology
    for months in ['FM', 'ON']:
        climFieldName = 'iceThickness'
        for hemisphere in ['NH', 'SH']:
            key = (months, hemisphere)
            obsFileName = build_config_full_path(
                config, 'seaIceObservations',
                'thickness{}_{}'.format(hemisphere, months))
            if not os.path.isfile(obsFileName):
                raise OSError('Obs file {} not found.'.format(
                    obsFileName))

            obsFieldName = '{}_{}'.format(climFieldName, hemisphere)
            (climatologyFileName, regriddedFileName) = \
                climatology.get_observation_climatology_file_names(
                    config=config, fieldName=obsFieldName, monthNames=months,
                    componentName='seaIce', gridFileName=obsFileName,
                    latVarName='t_lat', lonVarName='t_lon')

            obsFileNames[key] = obsFileName
            regriddedObsFileNames[key] = regriddedFileName

            if not os.path.exists(regriddedFileName):
                buildObsClimatologies = True

    for months in ['FM', 'ON']:
        monthValues = constants.monthDictionary[months]
        field = 'iceVolumeCell'
        climFieldName = 'iceThickness'

        # interpolate the model results
        (climatologyFileName, climatologyPrefix, regriddedFileName) = \
            climatology.get_mpas_climatology_file_names(
                    config=config, fieldName=climFieldName,
                    monthNames=months)

        if overwriteMpasClimatology or not os.path.exists(climatologyFileName):
            seasonalClimatology = climatology.cache_climatologies(
                ds, monthValues, config, climatologyPrefix, calendar,
                printProgress=True)
            if seasonalClimatology is None:
                # apparently, there was no data available to create the
                # climatology
                warnings.warn('no data to create sea ice thickness '
                              'climatology for {}'.format(months))
                continue

            remappedClimatology = climatology.remap_and_write_climatology(
                config, seasonalClimatology, climatologyFileName,
                regriddedFileName, mpasRemapper)

        else:

            remappedClimatology = xr.open_dataset(regriddedFileName)

        iceThickness = remappedClimatology[field].values
        iceThickness = ma.masked_values(iceThickness, 0)
        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        for hemisphere in ['NH', 'SH']:
            obsFieldName = 'HI'

            (colormapResult, colorbarLevelsResult) = setup_colormap(
                config,
                'regriddedSeaIceConcThick',
                suffix='ThickResult{}'.format(hemisphere))
            (colormapDifference, colorbarLevelsDifference) = setup_colormap(
                config,
                'regriddedSeaIceConcThick',
                suffix='ThickDifference{}'.format(hemisphere))

            referenceLongitude = config.getfloat(
                'regriddedSeaIceConcThick',
                'referenceLongitude{}'.format(hemisphere))
            minimumLatitude = config.getfloat(
                'regriddedSeaIceConcThick',
                'minimumLatitude{}'.format(hemisphere))

            # now the observations
            key = (months, hemisphere)
            regriddedFileName = regriddedObsFileNames[key]

            if buildObsClimatologies:
                obsFileName = obsFileNames[key]
                obsRemapper = \
                    climatology.get_observations_remapper(
                        config=config, componentName='seaIce',
                        fieldName='seaIce', gridFileName=obsFileName,
                        latVarName='t_lat', lonVarName='t_lon')

                seasonalClimatology = xr.open_dataset(obsFileName)

                if obsRemapper is None:
                    remappedClimatology = seasonalClimatology
                else:
                    remappedClimatology = \
                        climatology.remap_and_write_climatology(
                            config, seasonalClimatology, climatologyFileName,
                            regriddedFileName, obsRemapper)

            obsIceThickness = remappedClimatology[obsFieldName].values

            # Mask thickness fields
            obsIceThickness = ma.masked_values(obsIceThickness, 0)
            if hemisphere == 'NH':
                # Obs thickness should be nan above 86 (ICESat data)
                obsIceThickness[latTarg > 86] = ma.masked
                plotProjection = 'npstere'
            else:
                plotProjection = 'spstere'

            difference = iceThickness - obsIceThickness

            title = '{} ({}, years {:04d}-{:04d})'.format(subtitle, months,
                                                          startYear, endYear)
            fileout = '{}/icethick{}_{}_{}_years{:04d}-{:04d}.png'.format(
                plotsDirectory, hemisphere, mainRunName, months, startYear,
                endYear)
            plot_polar_comparison(
                config,
                lonTarg,
                latTarg,
                iceThickness,
                obsIceThickness,
                difference,
                colormapResult,
                colorbarLevelsResult,
                colormapDifference,
                colorbarLevelsDifference,
                title=title,
                fileout=fileout,
                plotProjection=plotProjection,
                latmin=minimumLatitude,
                lon0=referenceLongitude,
                modelTitle=mainRunName,
                obsTitle='Observations (ICESat)',
                diffTitle='Model-Observations',
                cbarlabel='m')


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
