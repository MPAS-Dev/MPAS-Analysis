"""
General comparison of 2-d model fields against data.  Currently only supports
sea surface temperature (sst), sea surface salinity (sss) and mixed layer
depth (mld)

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

Last Modified
-------------
03/29/2017
"""

import xarray as xr
import datetime
import numpy as np
import netCDF4
import os

from ..shared.interpolation import interpolate

from ..shared.plot.plotting import plot_global_comparison, \
    setup_colormap
from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.climatology import climatology

from ..shared.analysis_task import setup_task

from ..shared.mpas_xarray import mpas_xarray


def ocn_modelvsobs(config, field):

    """
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
    03/29/2017
    """

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    simulationStartTime = get_simulation_start_time(runStreams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    inputFiles = historyStreams.readpath(streamName, startDate=startDate,
                                         endDate=endDate, calendar=calendar)
    print 'Reading files {} through {}'.format(inputFiles[0], inputFiles[-1])

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
        {'mld': "{}/holtetalley_mld_climatology.nc".format(
                observationsDirectory),
         'sst': "{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc".format(
                observationsDirectory),
         'sss': "{}/Aquarius_V3_SSS_Monthly.nc".format(
                observationsDirectory)}

    obsFileName = obsFileNames[field]

    buildObsClimatologies = overwriteObsClimatology
    for months in outputTimes:
        (climatologyFileName, regriddedFileName) = \
            climatology.get_observation_climatology_file_names(
                config=config, fieldName=field, monthNames=months,
                componentName='ocean', gridFileName=obsFileName,
                latVarName='lat', lonVarName='lon')
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

            dsObs = mpas_xarray.subset_variables(dsObs, [obsFieldName,
                                                         'month'])
            # Reorder dataset for consistence with other obs. data sets
            dsObs = dsObs.transpose('Time', 'lat', 'lon')

        # Set appropriate MLD figure labels
        observationTitleLabel = \
            "Observations (HolteTalley density threshold MLD)"
        outFileLabel = "mldHolteTalleyARGO"
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
            dsObs = dsObs.transpose('Time', 'lat', 'lon')
            dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
            dsObs.coords['month'] = dsObs['Time.month']

        obsFieldName = 'SST'

        # Set appropriate figure labels for SST
        observationTitleLabel = \
            "Observations (Hadley/OI, {} {:04d}-{:04d})".format(period,
                                                                climStartYear,
                                                                climEndYear)
        outFileLabel = "sstHADOI"
        unitsLabel = r'$^o$C'

    elif field == 'sss':

        iselvals = {'nVertLevels': 0}

        timeStart = datetime.datetime(2011, 8, 1)
        timeEnd = datetime.datetime(2014, 12, 31)

        if buildObsClimatologies:
            dsObs = xr.open_mfdataset(obsFileName)
            dsObs.rename({'time': 'Time'}, inplace=True)
            dsObs = dsObs.transpose('Time', 'lat', 'lon')
            dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
            dsObs.coords['month'] = dsObs['Time.month']

        obsFieldName = 'SSS'

        observationTitleLabel = "Observations (Aquarius, 2011-2014)"
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

    changed, startYear, endYear = \
        climatology.update_start_end_year(ds, config, calendar)

    mpasMappingFileName = climatology.write_mpas_mapping_file(
        config=config, meshFileName=restartFileName)

    if buildObsClimatologies:
        obsMappingFileName = \
            climatology.write_observations_mapping_file(
                config=config, componentName='ocean', fieldName=field,
                gridFileName=obsFileName, latVarName='lat',  lonVarName='lon')
    else:
        obsMappingFileName = None

    (colormapResult, colorbarLevelsResult) = setup_colormap(
        config, sectionName, suffix='Result')
    (colormapDifference, colorbarLevelsDifference) = setup_colormap(
        config, sectionName, suffix='Difference')

    # Interpolate and compute biases
    for months in outputTimes:
        monthValues = constants.monthDictionary[months]

        (climatologyFileName, regriddedFileName) = \
            climatology.get_mpas_climatology_file_names(config=config,
                                                        fieldName=field,
                                                        monthNames=months)

        if overwriteMpasClimatology or not os.path.exists(climatologyFileName):
            seasonalClimatology = climatology.compute_climatology(
                ds, monthValues, calendar)
            # write out the climatology so we can interpolate it with
            # interpolate.remap
            seasonalClimatology.to_netcdf(climatologyFileName)

        interpolate.remap(inFileName=climatologyFileName,
                          outFileName=regriddedFileName,
                          inWeightFileName=mpasMappingFileName,
                          sourceFileType='mpas',
                          overwrite=overwriteMpasClimatology)

        ncFile = netCDF4.Dataset(regriddedFileName, mode='r')
        modelOutput = ncFile.variables[field][:]
        lons = ncFile.variables["lon"][:]
        lats = ncFile.variables["lat"][:]
        ncFile.close()
        lonTarg, latTarg = np.meshgrid(lons, lats)

        # now the observations
        (climatologyFileName, regriddedFileName) = \
            climatology.get_observation_climatology_file_names(
                config=config, fieldName=field, monthNames=months,
                componentName='ocean', gridFileName=obsFileName,
                latVarName='lat', lonVarName='lon')

        if buildObsClimatologies:
            if (overwriteObsClimatology or
                    (not os.path.exists(climatologyFileName) and
                     not os.path.exists(regriddedFileName))):
                seasonalClimatology = climatology.compute_climatology(
                    dsObs, monthValues)
                # Either we want to overwite files or neither the climatology
                # nor its regridded counterpart exist. Write out the
                # climatology so we can interpolate it with interpolate.remap
                seasonalClimatology.to_netcdf(climatologyFileName)

            if obsMappingFileName is None:
                # no remapping is needed
                regriddedFileName = climatologyFileName
            else:
                interpolate.remap(inFileName=climatologyFileName,
                                  outFileName=regriddedFileName,
                                  inWeightFileName=obsMappingFileName,
                                  sourceFileType='latlon',
                                  overwrite=overwriteObsClimatology)

        # read in the results from the remapped files
        ncFile = netCDF4.Dataset(regriddedFileName, mode='r')
        observations = ncFile.variables[obsFieldName][:]
        ncFile.close()

        bias = modelOutput - observations

        outFileName = "{}/{}_{}_{}_years{:04d}-{:04d}.png".format(
                plotsDirectory, outFileLabel, mainRunName,
                months, startYear, endYear)
        title = "{} ({}, years {:04d}-{:04d})".format(
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
                               modelTitle="{}".format(mainRunName),
                               obsTitle=observationTitleLabel,
                               diffTitle="Model-Observations",
                               cbarlabel=unitsLabel)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
