'''
General comparison of 2-d model fields against data.  Currently only supports
sea surface temperature (sst), sea surface salinity (sss) and mixed layer
depth (mld)

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

Last Modified
-------------
04/08/2017
'''

import xarray as xr
import datetime
import numpy as np
import os
import warnings

from ..shared.plot.plotting import plot_global_comparison, \
    setup_colormap
from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.climatology import get_lat_lon_comparison_descriptor, \
    get_remapper, get_mpas_climatology_file_names, \
    get_observation_climatology_file_names, \
    compute_climatology, cache_climatologies, update_start_end_year, \
    remap_and_write_climatology

from ..shared.grid import MpasMeshDescriptor, LatLonGridDescriptor

from ..shared.analysis_task import setup_task

from ..shared.mpas_xarray import mpas_xarray


def ocn_modelvsobs(config, field):

    '''
    Plots a comparison of ACME/MPAS output to SST or MLD observations

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    field : {'sst', 'sss', 'mld'}
        The name of a field to be analyized

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    04/08/2017
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

    sectionName = 'regridded{}'.format(field.upper())
    outputTimes = config.getExpression(sectionName, 'comparisonTimes')

    # get a list of regridded observations files and check if they exist.  If
    # they are all there, we don't have to do anything else with the
    # observations
    obsFileNames = \
        {'mld': '{}/holtetalley_mld_climatology.nc'.format(
                observationsDirectory),
         'sst': '{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc'.format(
                observationsDirectory),
         'sss': '{}/Aquarius_V3_SSS_Monthly.nc'.format(
                observationsDirectory)}

    obsFileName = obsFileNames[field]

    obsDescriptor = LatLonGridDescriptor()
    obsDescriptor.read(fileName=obsFileName, latVarName='lat',
                       lonVarName='lon')

    comparisonDescriptor = get_lat_lon_comparison_descriptor(config)

    obsRemapper = get_remapper(
        config=config, sourceDescriptor=obsDescriptor,
        comparisonDescriptor=comparisonDescriptor,
        mappingFileSection='oceanObservations',
        mappingFileOption='{}ClimatologyMappingFile'.format(field),
        mappingFilePrefix='map_obs_{}'.format(field),
        method=config.get('oceanObservations', 'interpolationMethod'))

    buildObsClimatologies = overwriteObsClimatology
    for months in outputTimes:
        (climatologyFileName, regriddedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=field, monthNames=months,
                componentName='ocean', remapper=obsRemapper)
        if not os.path.exists(regriddedFileName):
            buildObsClimatologies = True
            break

    varList = [field]

    if field == 'mld':

        iselvals = None

        obsFieldName = 'mld_dt_mean'

        if buildObsClimatologies:
            # Load MLD observational data
            dsObs = xr.open_mfdataset(obsFileName)

            # Increment month value to be consistent with the model output
            dsObs.iMONTH.values += 1
            # Rename the dimensions to be consistent with other obs. data sets
            dsObs.rename({'month': 'calmonth', 'lat': 'latCoord',
                          'lon': 'lonCoord'}, inplace=True)
            dsObs.rename({'iMONTH': 'Time', 'iLAT': 'lat', 'iLON': 'lon'},
                         inplace=True)

            # set the coordinates now that the dimensions have the same names
            dsObs.coords['lat'] = dsObs['latCoord']
            dsObs.coords['lon'] = dsObs['lonCoord']
            dsObs.coords['Time'] = dsObs['calmonth']
            dsObs.coords['month'] = ('Time', np.array(dsObs['calmonth'], int))

            # no meaningful year since this is already a climatology
            dsObs.coords['year'] = ('Time', np.ones(dsObs.dims['Time'], int))

            dsObs = mpas_xarray.subset_variables(dsObs, [obsFieldName,
                                                         'month'])

        # Set appropriate MLD figure labels
        observationTitleLabel = \
            'Observations (HolteTalley density threshold MLD)'
        outFileLabel = 'mldHolteTalleyARGO'
        unitsLabel = 'm'

    elif field == 'sst':

        iselvals = {'nVertLevels': 0}

        climStartYear = config.getint('oceanObservations',
                                      'sstClimatologyStartYear')
        climEndYear = config.getint('oceanObservations',
                                    'sstClimatologyEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        if climStartYear < 1925:
            period = 'pre-industrial'
        else:
            period = 'present-day'

        if buildObsClimatologies:
            dsObs = xr.open_mfdataset(obsFileName)
            dsObs.rename({'time': 'Time'}, inplace=True)
            dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
            dsObs.coords['month'] = dsObs['Time.month']
            dsObs.coords['year'] = dsObs['Time.year']

        obsFieldName = 'SST'

        # Set appropriate figure labels for SST
        observationTitleLabel = \
            'Observations (Hadley/OI, {} {:04d}-{:04d})'.format(period,
                                                                climStartYear,
                                                                climEndYear)
        outFileLabel = 'sstHADOI'
        unitsLabel = r'$^o$C'

    elif field == 'sss':

        iselvals = {'nVertLevels': 0}

        timeStart = datetime.datetime(2011, 8, 1)
        timeEnd = datetime.datetime(2014, 12, 31)

        if buildObsClimatologies:
            dsObs = xr.open_mfdataset(obsFileName)
            dsObs.rename({'time': 'Time'}, inplace=True)
            dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
            dsObs.coords['month'] = dsObs['Time.month']
            dsObs.coords['year'] = dsObs['Time.year']

        obsFieldName = 'SSS'

        observationTitleLabel = 'Observations (Aquarius, 2011-2014)'
        outFileLabel = 'sssAquarius'
        unitsLabel = 'PSU'

    ds = open_multifile_dataset(fileNames=inputFiles,
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=varList,
                                iselValues=iselvals,
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    changed, startYear, endYear = update_start_end_year(ds, config, calendar)

    mpasDescriptor = MpasMeshDescriptor(
        restartFileName, meshName=config.get('input', 'mpasMeshName'))

    mpasRemapper = get_remapper(
        config=config, sourceDescriptor=mpasDescriptor,
        comparisonDescriptor=comparisonDescriptor,
        mappingFileSection='climatology', mappingFileOption='mpasMappingFile',
        mappingFilePrefix='map', method=config.get('climatology',
                                                   'mpasInterpolationMethod'))

    (colormapResult, colorbarLevelsResult) = setup_colormap(
        config, sectionName, suffix='Result')
    (colormapDifference, colorbarLevelsDifference) = setup_colormap(
        config, sectionName, suffix='Difference')

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

            if seasonalClimatology is None:
                # apparently, there was no data available to create the
                # climatology
                warnings.warn('no data to create {} climatology for {}'.format(
                    field, months))
                continue

            remappedClimatology = remap_and_write_climatology(
                config, seasonalClimatology, climatologyFileName,
                regriddedFileName, mpasRemapper)

        else:

            remappedClimatology = xr.open_dataset(regriddedFileName)

        modelOutput = remappedClimatology[field].values
        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        # now the observations
        (climatologyFileName, regriddedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=field, monthNames=months,
                componentName='ocean', remapper=obsRemapper)

        if buildObsClimatologies and (overwriteObsClimatology or
                                      not os.path.exists(regriddedFileName)):

            seasonalClimatology = compute_climatology(
                dsObs, monthValues, maskVaries=True)

            if obsRemapper is None:
                # no need to remap because the observations are on the
                # comparison grid already
                remappedClimatology = seasonalClimatology
            else:
                remappedClimatology = \
                    remap_and_write_climatology(
                        config, seasonalClimatology, climatologyFileName,
                        regriddedFileName, obsRemapper)

        else:

            remappedClimatology = xr.open_dataset(regriddedFileName)
        observations = remappedClimatology[obsFieldName].values

        bias = modelOutput - observations

        outFileName = '{}/{}_{}_{}_years{:04d}-{:04d}.png'.format(
                plotsDirectory, outFileLabel, mainRunName,
                months, startYear, endYear)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                field.upper(), months, startYear, endYear)
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput,
                               observations,
                               bias,
                               colormapResult,
                               colorbarLevelsResult,
                               colormapDifference,
                               colorbarLevelsDifference,
                               fileout=outFileName,
                               title=title,
                               modelTitle='{}'.format(mainRunName),
                               obsTitle=observationTitleLabel,
                               diffTitle='Model-Observations',
                               cbarlabel=unitsLabel)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
