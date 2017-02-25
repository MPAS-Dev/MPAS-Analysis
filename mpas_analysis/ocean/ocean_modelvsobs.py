#!/usr/bin/env python
"""
General comparison of 2-d model fields against data.  Currently only supports
sea surface temperature (sst), sea surface salinity (sss) and mixed layer
depth (mld)

Author: Luke Van Roekel, Milena Veneziani, Xylar Asay-Davis
Last Modified: 02/25/2017
"""

import matplotlib.pyplot as plt
import matplotlib.colors as cols

import xarray as xr
import datetime

from ..shared.plot.plotting import plot_global_comparison
from ..shared.constants import constants

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.climatology import climatology


def ocn_modelvsobs(config, field, streamMap=None, variableMap=None):

    """
    Plots a comparison of ACME/MPAS output to SST or MLD observations

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    field is the name of a field to be analyize (currently one of 'sst', 'sss'
    or 'mld')

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Authors: Luke Van Roekel, Milena Veneziani, Xylar Asay-Davis
    Last Modified: 02/02/2017
    """

    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    calendar = namelist.get('config_calendar_type')
    simulationStartTime = get_simulation_start_time(streams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    inputFiles = streams.readpath(streamName, startDate=startDate,
                                  endDate=endDate, calendar=calendar)
    print 'Reading files {} through {}'.format(inputFiles[0], inputFiles[-1])

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')
    observationsDirectory = buildConfigFullPath(config, 'oceanObservations',
                                                '{}Subdirectory'.format(field))
    mainRunName = config.get('runs', 'mainRunName')

    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for ocn_modelvsobs calculation')

    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')

    sectionName = 'regridded{}'.format(field.upper())
    outputTimes = config.getExpression(sectionName, 'comparisonTimes')

    varList = [field]

    if field == 'mld':

        iselvals = None

        # Load MLD observational data
        obsFileName = "{}/holtetalley_mld_climatology.nc".format(
                observationsDirectory)
        dsObs = xr.open_mfdataset(obsFileName)

        # Increment month value to be consistent with the model output
        dsObs.iMONTH.values += 1

        # Rename the time dimension to be consistent with the SST dataset
        dsObs.rename({'month': 'calmonth'}, inplace=True)
        dsObs.rename({'iMONTH': 'month'}, inplace=True)
        dsObs.coords['month'] = dsObs['calmonth']

        obsFieldName = 'mld_dt_mean'

        # Reorder dataset for consistence
        dsObs = dsObs.transpose('month', 'iLON', 'iLAT')

        # Set appropriate MLD figure labels
        observationTitleLabel = \
            "Observations (HolteTalley density threshold MLD)"
        outFileLabel = "mldHolteTalleyARGO"
        unitsLabel = 'm'

    elif field == 'sst':

        iselvals = {'nVertLevels': 0}

        obsFileName = \
            "{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc".format(
                    observationsDirectory)
        dsObs = xr.open_mfdataset(obsFileName)

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

        dsTimeSlice = dsObs.sel(time=slice(timeStart, timeEnd))
        monthlyClimatology = dsTimeSlice.groupby('time.month').mean('time')

        dsObs = monthlyClimatology.transpose('month', 'lon', 'lat')
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

        obsFileName = "{}/Aquarius_V3_SSS_Monthly.nc".format(
                observationsDirectory)
        dsObs = xr.open_mfdataset(obsFileName)

        timeStart = datetime.datetime(2011, 8, 1)
        timeEnd = datetime.datetime(2014, 12, 31)

        dsTimeSlice = dsObs.sel(time=slice(timeStart, timeEnd))

        # The following line converts from DASK to numpy to supress an odd
        # warning that doesn't influence the figure output
        dsTimeSlice.SSS.values

        monthlyClimatology = dsTimeSlice.groupby('time.month').mean('time')

        # Rename the observation data for code compactness
        dsObs = monthlyClimatology.transpose('month', 'lon', 'lat')
        obsFieldName = 'SSS'

        observationTitleLabel = "Observations (Aquarius, 2011-2014)"
        outFileLabel = 'sssAquarius'
        unitsLabel = 'PSU'

    ds = open_multifile_dataset(fileNames=inputFiles,
                                calendar=calendar,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=varList,
                                iselValues=iselvals,
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    monthlyClimatology = climatology.compute_monthly_climatology(ds, calendar)

    mpasInterolationData = climatology.init_model_interpolation(restartFile)
    lonTarg = mpasInterolationData[2]
    latTarg = mpasInterolationData[3]
    obsInterolationData = climatology.init_observations_interpolation(dsObs)

    resultContourValues = config.getExpression(sectionName,
                                               'resultContourValues')
    resultColormap = plt.get_cmap(config.get(sectionName, 'resultColormap'))
    resultColormapIndices = config.getExpression(sectionName,
                                                 'resultColormapIndices')
    resultColormap = cols.ListedColormap(resultColormap(resultColormapIndices),
                                         "resultColormap")
    differenceContourValues = config.getExpression(sectionName,
                                                   'differenceContourValues')
    differenceColormap = plt.get_cmap(config.get(sectionName,
                                                 'differenceColormap'))
    differenceColormapIndices = config.getExpression(
            sectionName, 'differenceColormapIndices')
    differenceColormap = cols.ListedColormap(
            differenceColormap(differenceColormapIndices),
            "differenceColormap")

    # Interpolate and compute biases
    for timeIndex, timestring in enumerate(outputTimes):
        monthsValue = constants.monthDictionary[timestring]

        modelOutput = climatology.interpolate_model_climatology(
            mpasInterolationData, monthlyClimatology, field, monthsValue)
        observations = climatology.interpolate_observation_climatology(
            obsInterolationData, dsObs, obsFieldName, monthsValue)

        bias = modelOutput - observations

        outFileName = "{}/{}_{}_{}_years{:04d}-{:04d}.png".format(
                plotsDirectory, outFileLabel, mainRunName,
                outputTimes[timeIndex], startYear, endYear)
        title = "{} ({}, years {:04d}-{:04d})".format(
                field.upper(), outputTimes[timeIndex], startYear, endYear)
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput,
                               observations,
                               bias,
                               resultColormap,
                               resultContourValues,
                               differenceColormap,
                               differenceContourValues,
                               fileout=outFileName,
                               title=title,
                               modelTitle="{}".format(mainRunName),
                               obsTitle=observationTitleLabel,
                               diffTitle="Model-Observations",
                               cbarlabel=unitsLabel)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
