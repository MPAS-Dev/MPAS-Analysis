import numpy as np
import xarray as xr
import pandas as pd
import datetime

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, \
    remove_repeated_time_index, subset_variables

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.timekeeping.utility import stringToDatetime, \
    clampToNumpyDatetime64


def seaice_timeseries(config, streamMap=None, variableMap=None):
    """
    Performs analysis of time series of sea-ice properties.

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Author: Xylar Asay-Davis, Milena Veneziani
    Last Modified: 02/02/2017
    """

    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'seaIceNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'seaIceStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    calendar = namelist.get('config_calendar_type')

    # get a list of timeSeriesStatsMonthly output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('timeSeries', 'startDate')
    endDate = config.get('timeSeries', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    inFiles = streams.readpath(streamName, startDate=startDate,
                               endDate=endDate,  calendar=calendar)
    print 'Reading files {} through {}'.format(inFiles[0], inFiles[-1])

    variableNames = ['iceAreaCell', 'iceVolumeCell']

    plotTitles = {'iceAreaCell': 'Sea-ice area',
                  'iceVolumeCell': 'Sea-ice volume',
                  'iceThickness': 'Sea-ice thickness'}

    unitsDictionary = {'iceAreaCell': '[km$^2$]',
                       'iceVolumeCell': '[10$^3$ km$^3$]',
                       'iceThickness': '[m]'}

    obsFileNames = {
        'iceAreaCell': [buildConfigFullPath(config, 'seaIceObservations',
                                            subdir)
                        for subdir in ['areaNH', 'areaSH']],
        'iceVolumeCell': [buildConfigFullPath(config, 'seaIceObservations',
                                              subdir)
                          for subdir in ['volNH', 'volSH']]}

    # Some plotting rules
    titleFontSize = config.get('timeSeriesSeaIceAreaVol', 'titleFontSize')

    mainRunName = config.get('runs', 'mainRunName')
    preprocessedReferenceRunName = config.get('runs',
                                              'preprocessedReferenceRunName')
    preprocessedReferenceDirectory = config.get('seaIcePreprocessedReference',
                                                'baseDirectory')

    compareWithObservations = config.getboolean('timeSeriesSeaIceAreaVol',
                                                'compareWithObservations')

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')

    yearOffset = config.getint('time', 'yearOffset')

    movingAveragePoints = config.getint('timeSeriesSeaIceAreaVol',
                                        'movingAveragePoints')

    # first, check for a sea-ice restart file
    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        # get an ocean restart file, since no sea-ice restart exists
        oceanStreamsFileName = config.get('input', 'oceanStreamsFileName')
        oceanStreams = StreamsFile(oceanStreamsFileName,
                                   streamsdir=inDirectory)
        try:
            restartFile = oceanStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O or MPAS-Seaice restart file found: need '
                          'at least one restart file for seaice_timeseries '
                          'calculation')

    print '  Load sea-ice data...'
    # Load mesh
    dsMesh = xr.open_dataset(restartFile)
    dsMesh = subset_variables(dsMesh, vlist=['lonCell', 'latCell', 'areaCell'])

    # Load data
    ds = xr.open_mfdataset(
        inFiles,
        preprocess=lambda x: preprocess_mpas(x, yearoffset=yearOffset,
                                             timestr='Time',
                                             onlyvars=['iceAreaCell',
                                                       'iceVolumeCell'],
                                             varmap=variableMap))
    ds = remove_repeated_time_index(ds)

    timeStart = clampToNumpyDatetime64(stringToDatetime(startDate), yearOffset)
    timeEnd = clampToNumpyDatetime64(stringToDatetime(endDate), yearOffset)
    # select only the data in the specified range of years
    ds = ds.sel(Time=slice(timeStart, timeEnd))

    # handle the case where the "mesh" file has a spurious time dimension
    if 'Time' in dsMesh.keys():
        dsMesh = dsMesh.drop('Time')
    ds = ds.merge(dsMesh)

    yearStart = (pd.to_datetime(ds.Time.min().values)).year
    yearEnd = (pd.to_datetime(ds.Time.max().values)).year
    timeStart = datetime.datetime(yearStart, 1, 1)
    timeEnd = datetime.datetime(yearEnd, 12, 31)

    if preprocessedReferenceRunName != 'None':
        inFilesPreprocessed = '{}/icevol.{}.year*.nc'.format(
            preprocessedReferenceDirectory, preprocessedReferenceRunName)
        dsPreprocessed = xr.open_mfdataset(
            inFilesPreprocessed,
            preprocess=lambda x: preprocess_mpas(x, yearoffset=yearOffset))
        preprocessedYearEnd = (pd.to_datetime(
            dsPreprocessed.Time.max().values)).year
        if yearStart <= preprocessedYearEnd:
            dsPreprocessedTimeSlice = dsPreprocessed.sel(Time=slice(timeStart,
                                                                    timeEnd))
        else:
            print '   Warning: Preprocessed time series ends before the ' \
                'timeSeries startYear and will not be plotted.'
            preprocessedReferenceRunName = 'None'

    # Make Northern and Southern Hemisphere partition:
    areaCell = ds.areaCell
    maskNH = ds.latCell > 0
    maskSH = ds.latCell < 0
    areaCellNH = areaCell.where(maskNH)
    areaCellSH = areaCell.where(maskSH)

    for variableName in variableNames:
        obsFileNameNH = obsFileNames[variableName][0]
        obsFileNameSH = obsFileNames[variableName][1]
        plotTitle = plotTitles[variableName]
        units = unitsDictionary[variableName]

        print '  Compute NH and SH time series of {}...'.format(variableName)
        if variableName == 'iceThickCell':
            variableNamefull = 'iceVolumeCell'
        else:
            variableNamefull = variableName
        var = ds[variableNamefull]

        varNH = var.where(maskNH)*areaCellNH
        varSH = var.where(maskSH)*areaCellSH

        maskIceExtent = var > 0.15
        varNHIceExtent = varNH.where(maskIceExtent)
        varSHIceExtent = varSH.where(maskIceExtent)

        if variableName == 'iceAreaCell':
            varNH = varNH.sum('nCells')
            varSH = varSH.sum('nCells')
            varNH = 1e-6*varNH  # m^2 to km^2
            varSH = 1e-6*varSH  # m^2 to km^2
            varNHIceExtent = 1e-6*varNHIceExtent.sum('nCells')
            varSHIceExtent = 1e-6*varSHIceExtent.sum('nCells')
        elif variableName == 'iceVolumeCell':
            varNH = varNH.sum('nCells')
            varSH = varSH.sum('nCells')
            varNH = 1e-3*1e-9*varNH  # m^3 to 10^3 km^3
            varSH = 1e-3*1e-9*varSH  # m^3 to 10^3 km^3
        else:
            varNH = varNH.mean('nCells')/areaCellNH.mean('nCells')
            varSH = varSH.mean('nCells')/areaCellSH.mean('nCells')

        print '  Make plots...'

        xLabel = 'Time [years]'

        if preprocessedReferenceRunName != 'None':
            figureNameNH = '{}/{}NH_{}_{}.png'.format(
                plotsDirectory, variableName, mainRunName,
                preprocessedReferenceRunName)
            figureNameSH = '{}/{}SH_{}_{}.png'.format(
                plotsDirectory, variableName, mainRunName,
                preprocessedReferenceRunName)
        else:
            figureNameNH = '{}/{}NH_{}.png'.format(plotsDirectory,
                                                   variableName,
                                                   mainRunName)
            figureNameSH = '{}/{}SH_{}.png'.format(plotsDirectory,
                                                   variableName,
                                                   mainRunName)

        titleNH = '{} (NH), {} (r)'.format(plotTitle, mainRunName)
        titleSH = '{} (SH), {} (r)'.format(plotTitle, mainRunName)

        if compareWithObservations:
            if variableName == 'iceAreaCell':
                titleNH = \
                    '{}\nSSM/I observations, annual cycle (k)'.format(titleNH)
                titleSH = \
                    '{}\nSSM/I observations, annual cycle (k)'.format(titleSH)
            elif variableName == 'iceVolumeCell':
                titleNH = '{}\nPIOMAS, annual cycle (k)'.format(titleNH)
                titleSH = '{}\n'.format(titleSH)

        if preprocessedReferenceRunName != 'None':
            titleNH = '{}\n {} (b)'.format(titleNH,
                                           preprocessedReferenceRunName)
            titleSH = '{}\n {} (b)'.format(titleSH,
                                           preprocessedReferenceRunName)

        if variableName == 'iceAreaCell':

            if compareWithObservations:
                dsObs = xr.open_mfdataset(
                    obsFileNameNH,
                    preprocess=lambda x: preprocess_mpas(
                        x, yearoffset=yearOffset))
                dsObs = remove_repeated_time_index(dsObs)
                varNHObs = dsObs.IceArea
                varNHObs = replicate_cycle(varNH, varNHObs)

                dsObs = xr.open_mfdataset(
                    obsFileNameSH,
                    preprocess=lambda x: preprocess_mpas(
                        x, yearoffset=yearOffset))
                dsObs = remove_repeated_time_index(dsObs)
                varSHObs = dsObs.IceArea
                varSHObs = replicate_cycle(varSH, varSHObs)

            if preprocessedReferenceRunName != 'None':
                inFilesPreprocessed = '{}/icearea.{}.year*.nc'.format(
                    preprocessedReferenceDirectory,
                    preprocessedReferenceRunName)
                dsPreprocessed = xr.open_mfdataset(
                    inFilesPreprocessed,
                    preprocess=lambda x: preprocess_mpas(
                        x, yearoffset=yearOffset))
                dsPreprocessedTimeSlice = dsPreprocessed.sel(
                    Time=slice(timeStart, timeEnd))
                varNHPreprocessed = dsPreprocessedTimeSlice.icearea_nh
                varSHPreprocessed = dsPreprocessedTimeSlice.icearea_sh

        elif variableName == 'iceVolumeCell':

            if compareWithObservations:
                dsObs = xr.open_mfdataset(
                    obsFileNameNH,
                    preprocess=lambda x: preprocess_mpas(
                        x, yearoffset=yearOffset))
                dsObs = remove_repeated_time_index(dsObs)
                varNHObs = dsObs.IceVol
                varNHObs = replicate_cycle(varNH, varNHObs)

                varSHObs = None

            if preprocessedReferenceRunName != 'None':
                inFilesPreprocessed = '{}/icevol.{}.year*.nc'.format(
                    preprocessedReferenceDirectory,
                    preprocessedReferenceRunName)
                dsPreprocessed = xr.open_mfdataset(
                    inFilesPreprocessed,
                    preprocess=lambda x: preprocess_mpas(
                        x, yearoffset=yearOffset))
                dsPreprocessedTimeSlice = dsPreprocessed.sel(
                    Time=slice(timeStart, timeEnd))
                varNHPreprocessed = dsPreprocessedTimeSlice.icevolume_nh
                varSHPreprocessed = dsPreprocessedTimeSlice.icevolume_sh

        if variableName in ['iceAreaCell', 'iceVolumeCell']:
            if compareWithObservations:
                if preprocessedReferenceRunName != 'None':
                    varsNH = [varNH, varNHObs, varNHPreprocessed]
                    varsSH = [varSH, varSHObs, varSHPreprocessed]
                    lineStyles = ['r-', 'k-', 'b-']
                    lineWidths = [1.2, 1.2, 1.2]
                else:
                    # just v1 model and obs
                    varsNH = [varNH, varNHObs]
                    varsSH = [varSH, varSHObs]
                    lineStyles = ['r-', 'k-']
                    lineWidths = [1.2, 1.2]
            elif preprocessedReferenceRunName != 'None':
                # just v1 and v0 models
                varsNH = [varNH, varNHPreprocessed]
                varsSH = [varSH, varSHPreprocessed]
                lineStyles = ['r-', 'b-']
                lineWidths = [1.2, 1.2]

            if (compareWithObservations or
                    preprocessedReferenceRunName != 'None'):
                # separate plots for nothern and southern hemispheres
                timeseries_analysis_plot(config, varsNH, movingAveragePoints,
                                         titleNH,
                                         xLabel, units, figureNameNH,
                                         lineStyles=lineStyles,
                                         lineWidths=lineWidths,
                                         titleFontSize=titleFontSize)
                timeseries_analysis_plot(config, varsSH, movingAveragePoints,
                                         titleSH,
                                         xLabel, units, figureNameSH,
                                         lineStyles=lineStyles,
                                         lineWidths=lineWidths,
                                         titleFontSize=titleFontSize)
            else:
                # we will combine north and south onto a single graph
                figureName = '{}/{}.{}.png'.format(plotsDirectory, mainRunName,
                                                   variableName)
                title = '{}, NH (r), SH (k)\n{}'.format(plotTitle, mainRunName)
                timeseries_analysis_plot(config, [varNH, varSH],
                                         movingAveragePoints,
                                         title, xLabel, units, figureName,
                                         lineStyles=['r-', 'k-'],
                                         lineWidths=[1.2, 1.2],
                                         titleFontSize=titleFontSize)

        elif variableName == 'iceThickCell':

            figureName = '{}/{}.{}.png'.format(plotsDirectory, mainRunName,
                                               variableName)
            title = '{} NH (r), SH (k)\n{}'.format(plotTitle, mainRunName)
            timeseries_analysis_plot(config, [varNH, varSH],
                                     movingAveragePoints, title,
                                     xLabel, units, figureName,
                                     lineStyles=['r-', 'k-'],
                                     lineWidths=[1.2, 1.2],
                                     titleFontSize=titleFontSize)

        else:
            raise ValueError(
                'variableName variable {} not supported for plotting'.format(
                    variableName))


def replicate_cycle(ds, dsToReplicate):
    dsShift = dsToReplicate.copy()
    shiftT = ((dsShift.Time.max() - dsShift.Time.min()) +
              (dsShift.Time.isel(Time=1) - dsShift.Time.isel(Time=0)))
    startIndex = int(np.floor((ds.Time.min()-dsToReplicate.Time.min())/shiftT))
    endIndex = int(np.ceil((ds.Time.max()-dsToReplicate.Time.min())/shiftT))
    dsShift['Time'] = dsShift['Time'] + startIndex*shiftT

    # replicate cycle:
    for cycleIndex in range(startIndex, endIndex):
        dsNew = dsToReplicate.copy()
        dsNew['Time'] = dsNew['Time'] + (cycleIndex+1)*shiftT
        dsShift = xr.concat([dsShift, dsNew], dim='Time')
    # constrict replicated dsSHort to same time dimension as ds_long:
    dsShift = dsShift.sel(Time=ds.Time.values, method='nearest')
    return dsShift
