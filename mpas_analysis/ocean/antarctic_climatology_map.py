'''
General comparison of Antarctic climatology maps against data.

Authors
-------
Xylar Asay-Davis

Last Modified
-------------
04/18/2017
'''

import xarray as xr
import numpy as np
import os

from ..shared.plot.plotting import plot_polar_projection_comparison
from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.climatology \
    import get_Antarctic_stereographic_comparison_descriptor, \
    get_remapper, get_mpas_climatology_file_names, \
    get_observation_climatology_file_names, \
    compute_climatology, cache_climatologies, update_start_end_year, \
    remap_and_write_climatology

from ..shared.grid import MpasMeshDescriptor, ProjectionGridDescriptor

from ..shared.analysis_task import setup_task

from ..shared.mpas_xarray import mpas_xarray


def antarctic_climatology_map(config, field):

    '''
    Plots a comparison of various Antarctic fields from ACME/MPAS output to
    observations and reanalysis

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    field : {'melt', 'botTemp', 'botSalin'}
        The name of the field to be analyzed

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/18/2017
    '''

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, namelistMap, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    simulationStartTime = get_simulation_start_time(runStreams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    inputFiles = historyStreams.readpath(streamName, startDate=startDate,
                                         endDate=endDate, calendar=calendar)
    print '\n  Reading files:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(inputFiles[0]),
              os.path.basename(inputFiles[-1]))

    observationsDirectory = build_config_full_path(
        config, 'oceanObservations', '{}Subdirectory'.format(field))
    mainRunName = config.get('runs', 'mainRunName')

    overwriteMpasClimatology = config.getWithDefault(
        'climatology', 'overwriteMpasClimatology', False)

    overwriteObsClimatology = config.getWithDefault(
        'oceanObservations', 'overwriteObsClimatology', False)

    try:
        restartFileName = runStreams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for ocn_modelvsobs calculation')

    if field == 'melt':
        setup = _setup_melt
        setup_obs = _setup_melt_obs
    elif field == 'botTemp':
        setup = _setup_bot_temp
        setup_obs = _setup_bot_temp_obs
    elif field == 'botSalin':
        setup = _setup_bot_salin
        setup_obs = _setup_bot_salin_obs
    else:
        raise ValueError('Field {} not supported.'.format(field))

    try:
        mpasFieldName, mpasMaskName, verticalIndexName, modelScaling,  \
            fieldInTitle, sectionName, obsFileName, obsFieldName, \
            observationTitleLabel, outFileLabel, unitsLabel,  obsDescriptor, \
            comparisonDescriptor = setup(config, namelist,
                                         observationsDirectory)
    except ValueError as e:
        print e.message, "Skipping analysis."
        return

    outputTimes = config.getExpression(sectionName, 'comparisonTimes')

    obsRemapper = get_remapper(
        config=config, sourceDescriptor=obsDescriptor,
        comparisonDescriptor=comparisonDescriptor,
        mappingFileSection='oceanObservations',
        mappingFileOption='{}ClimatologyMappingFile'.format(field),
        mappingFilePrefix='map_obs_{}'.format(field),
        method=config.get('oceanObservations', 'interpolationMethod'))

    buildObsClimatologies = _check_cached_obs_climatologies(
        config, overwriteObsClimatology, field, outputTimes, obsRemapper)

    if buildObsClimatologies:
        dsObs = setup_obs(config, obsFileName, obsFieldName)

    ds = open_multifile_dataset(fileNames=inputFiles,
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=[mpasFieldName],
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    dsRestart = xr.open_dataset(restartFileName)
    varList = ['maxLevelCell']
    if verticalIndexName is not None and verticalIndexName not in varList:
        varList.append(verticalIndexName)
    if mpasMaskName is not None:
        varList.append(mpasMaskName)
    dsRestart = mpas_xarray.subset_variables(dsRestart, varList)
    if 'Time' in dsRestart.dims:
        dsRestart = dsRestart.isel(Time=0)

    oceanMask = xr.DataArray(np.ones(dsRestart.dims['nCells']),
                             dims=['nCells'])
    oceanMask = oceanMask.where(dsRestart.maxLevelCell > 0)

    if verticalIndexName is not None:
        verticalIndex = dsRestart[verticalIndexName] - 1
        # make a copy so we can add a vertical coordinate that we need only
        # temporarily
        dsCopy = ds.copy()
        dsCopy.coords['verticalIndex'] = ('nVertLevels',
                                          np.arange(ds.dims['nVertLevels']))
        # mask only the values with the right vertical index
        da = ds[mpasFieldName].where(dsCopy.verticalIndex == verticalIndex)
        # Each vertical layer has at most one non-NaN value so the "sum" over
        # the vertical is used to collapse the array in the vertical dimension
        ds[field] = da.sum(dim='nVertLevels')
        if mpasFieldName != field:
            ds = ds.drop(mpasFieldName)
    elif mpasFieldName != field:
        # we expect a variable with the field name later on
        ds.rename(mpasFieldName, field)

    changed, startYear, endYear = \
        update_start_end_year(ds, config, calendar)

    mpasDescriptor = MpasMeshDescriptor(
        restartFileName, meshName=config.get('input', 'mpasMeshName'))

    mpasRemapper = get_remapper(
        config=config, sourceDescriptor=mpasDescriptor,
        comparisonDescriptor=comparisonDescriptor,
        mappingFileSection='climatology',
        mappingFileOption='mpasToAntarcticStereoMappingFile',
        mappingFilePrefix='map', method=config.get('climatology',
                                                   'mpasInterpolationMethod'))

    oceanMask = mpasRemapper.remap(oceanMask)
    landMask = np.ma.masked_array(
        np.ones(oceanMask.values.shape),
        mask=np.logical_not(np.isnan(oceanMask.values)))

    # Interpolate and compute biases
    for months in outputTimes:
        monthValues = constants.monthDictionary[months]

        (climatologyFileName, climatologyPrefix, regriddedFileName) = \
            get_mpas_climatology_file_names(config=config,
                                            fieldName=field,
                                            monthNames=months,
                                            remapper=mpasRemapper)

        if overwriteMpasClimatology or not os.path.exists(regriddedFileName):
            seasonalClimatology = cache_climatologies(
                ds, monthValues, config, climatologyPrefix, calendar,
                printProgress=True)

            if mpasMaskName is not None:
                seasonalClimatology[field] = \
                    seasonalClimatology[field].where(
                        dsRestart[mpasMaskName] > 0.)

            remappedModelClimatology = remap_and_write_climatology(
                config, seasonalClimatology, climatologyFileName,
                regriddedFileName, mpasRemapper, useNcremap=False)

        else:

            remappedModelClimatology = xr.open_dataset(regriddedFileName)

        modelOutput = np.ma.masked_array(
            modelScaling*remappedModelClimatology[field].values,
            mask=np.isnan(remappedModelClimatology[field].values))

        # now the observations
        (climatologyFileName, regriddedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=field, monthNames=months,
                componentName='ocean', remapper=obsRemapper)

        if buildObsClimatologies and (overwriteObsClimatology or
                                      not os.path.exists(regriddedFileName)):

            seasonalClimatology = compute_climatology(
                dsObs, monthValues, maskVaries=True)

            remappedClimatology = \
                remap_and_write_climatology(
                    config, seasonalClimatology, climatologyFileName,
                    regriddedFileName, obsRemapper)

        else:

            remappedClimatology = xr.open_dataset(regriddedFileName)

        observations = np.ma.masked_array(
            remappedClimatology[obsFieldName].values,
            mask=np.isnan(remappedClimatology[obsFieldName].values))

        bias = modelOutput - observations

        outFileName = '{}/{}_{}_{}_years{:04d}-{:04d}.png'.format(
                plotsDirectory, outFileLabel, mainRunName,
                months, startYear, endYear)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                fieldInTitle, months, startYear, endYear)

        x = comparisonDescriptor.xCorner
        y = comparisonDescriptor.yCorner
        plot_polar_projection_comparison(config,
                                         x,
                                         y,
                                         landMask,
                                         modelOutput,
                                         observations,
                                         bias,
                                         sectionName,
                                         fileout=outFileName,
                                         title=title,
                                         modelTitle='{}'.format(mainRunName),
                                         obsTitle=observationTitleLabel,
                                         diffTitle='Model - Observations',
                                         cbarlabel=unitsLabel)


def _check_cached_obs_climatologies(config, overwriteObsClimatology, field,
                                    outputTimes, obsRemapper):  # {{{
    '''
    Check if we need to build climatologies for any observations (i.e. cached
    files don't already exist for all climatologies)

    Authors
    -------
    Xylar Asay-Davis
    '''

    buildObsClimatologies = overwriteObsClimatology
    for months in outputTimes:
        (climatologyFileName, regriddedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=field, monthNames=months,
                componentName='ocean', remapper=obsRemapper)
        if not os.path.exists(regriddedFileName):
            buildObsClimatologies = True
            break

    return buildObsClimatologies  # }}}


def _setup_melt(config, namelist, observationsDirectory):  # {{{

    try:
        flux_mode = namelist.get('config_land_ice_flux_mode')
    except ValueError:
        flux_mode = 'off'

    if flux_mode not in ['standalone', 'coupled']:
        raise ValueError('sub-ice-shelf melting disabled.')

    verticalIndexName = None

    mpasFieldName = 'melt'
    mpasMaskName = 'landIceMask'

    # model output is in kg/m^2/s, but we want m/yr
    modelScaling = constants.sec_per_year/constants.rho_fw

    fieldInTitle = 'Melt Rate'
    sectionName = 'regriddedAntarcticMelt'

    obsFileName = '{}/Rignot_2013_melt_rates_6000.0x6000.0km_10.0km_' \
                  'Antarctic_stereo.nc'.format(observationsDirectory)

    obsFieldName = 'meltRate'

    # Set appropriate melt rate figure labels
    observationTitleLabel = 'Observations (Rignot et al, 2013)'
    outFileLabel = 'meltRignot'
    unitsLabel = 'm a$^{-1}$'

    comparisonDescriptor = \
        get_Antarctic_stereographic_comparison_descriptor(config)

    obsDescriptor = ProjectionGridDescriptor(
        comparisonDescriptor.projection)

    # the mesh name will be read from the file
    obsDescriptor.read(fileName=obsFileName, xVarName='x', yVarName='y')

    return mpasFieldName, mpasMaskName, verticalIndexName, modelScaling,  \
        fieldInTitle, sectionName, obsFileName, obsFieldName, \
        observationTitleLabel, outFileLabel, unitsLabel,  obsDescriptor, \
        comparisonDescriptor  # }}}


def _setup_melt_obs(config, obsFileName, obsFieldName):  # {{{
    dsObs = xr.open_mfdataset(obsFileName)
    # since this is already a climatology, add a fake month and year to
    # "average" over
    dsObs.coords['month'] = ('Time', np.ones(1, int))
    dsObs.coords['year'] = ('Time', np.ones(1, int))
    return dsObs  # }}}


def _setup_bot_temp(config, namelist, observationsDirectory):  # {{{

    verticalIndexName = 'maxLevelCell'

    mpasFieldName = 'temperature'
    mpasMaskName = None

    # model output is already what we want (degrees C)
    modelScaling = 1.

    fieldInTitle = 'Ocean Bot. Temp.'
    sectionName = 'regriddedAntarcticBotTemp'

    obsFileName = '{}/SOSE_2005-2010_pot_temp_6000.0x6000.0km_10.0km_' \
                  'Antarctic_stereo.nc'.format(observationsDirectory)

    obsFieldName = 'botTheta'

    # Set appropriate melt rate figure labels
    observationTitleLabel = 'Observations (SOSE)'
    outFileLabel = 'botTemp_SOSE'
    unitsLabel = '$^\circ$C'

    comparisonDescriptor = \
        get_Antarctic_stereographic_comparison_descriptor(config)

    obsDescriptor = ProjectionGridDescriptor(
        comparisonDescriptor.projection)

    # the mesh name will be read from the file
    obsDescriptor.read(fileName=obsFileName, xVarName='x', yVarName='y')

    return mpasFieldName, mpasMaskName, verticalIndexName, modelScaling,  \
        fieldInTitle, sectionName, obsFileName, obsFieldName, \
        observationTitleLabel, outFileLabel, unitsLabel,  obsDescriptor, \
        comparisonDescriptor  # }}}


def _setup_bot_temp_obs(config, obsFileName, obsFieldName):  # {{{
    dsObs = xr.open_mfdataset(obsFileName)
    # since this is already a climatology, add a fake month and year to
    # "average" over
    dsObs.coords['month'] = ('Time', np.ones(1, int))
    dsObs.coords['year'] = ('Time', np.ones(1, int))
    return dsObs  # }}}


def _setup_bot_salin(config, namelist, observationsDirectory):  # {{{

    verticalIndexName = 'maxLevelCell'

    mpasFieldName = 'salinity'
    mpasMaskName = None

    # model output is already what we want (degrees C)
    modelScaling = 1.

    fieldInTitle = 'Ocean Bot. Salin.'
    sectionName = 'regriddedAntarcticBotSalin'

    obsFileName = '{}/SOSE_2005-2010_salinity_6000.0x6000.0km_10.0km_' \
                  'Antarctic_stereo.nc'.format(observationsDirectory)

    obsFieldName = 'botSalinity'

    # Set appropriate melt rate figure labels
    observationTitleLabel = 'Observations (SOSE)'
    outFileLabel = 'botSalin_SOSE'
    unitsLabel = 'PSU'

    comparisonDescriptor = \
        get_Antarctic_stereographic_comparison_descriptor(config)

    obsDescriptor = ProjectionGridDescriptor(
        comparisonDescriptor.projection)

    # the mesh name will be read from the file
    obsDescriptor.read(fileName=obsFileName, xVarName='x', yVarName='y')

    return mpasFieldName, mpasMaskName, verticalIndexName, modelScaling,  \
        fieldInTitle, sectionName, obsFileName, obsFieldName, \
        observationTitleLabel, outFileLabel, unitsLabel,  obsDescriptor, \
        comparisonDescriptor  # }}}


def _setup_bot_salin_obs(config, obsFileName, obsFieldName):  # {{{
    dsObs = xr.open_mfdataset(obsFileName)
    # since this is already a climatology, add a fake month and year to
    # "average" over
    dsObs.coords['month'] = ('Time', np.ones(1, int))
    dsObs.coords['year'] = ('Time', np.ones(1, int))
    return dsObs  # }}}



# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
