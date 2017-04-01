import xarray as xr
import os

from ..shared.plot.plotting import timeseries_analysis_plot, \
    timeseries_analysis_plot_polar

from ..shared.io.utility import build_config_full_path, make_directories

from ..shared.timekeeping.utility import date_to_days, days_to_datetime, \
    datetime_to_days
from ..shared.timekeeping.MpasRelativeDelta import MpasRelativeDelta

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset
from ..shared.mpas_xarray.mpas_xarray import subset_variables

from .utility import setup_sea_ice_task

from ..shared.time_series import time_series


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
    Last Modified: 04/08/2017
    """
    def compute_area_vol_part(timeIndices, firstCall):
        dsLocal = ds.isel(Time=timeIndices)

        if hemisphere == 'NH':
            mask = dsMesh.latCell > 0
        else:
            mask = dsMesh.latCell < 0
        dsLocal = dsLocal.where(mask)

        dsAreaSum = (dsLocal*dsMesh.areaCell).sum('nCells')
        dsAreaSum = dsAreaSum.rename({'iceAreaCell': 'iceArea',
                                      'iceVolumeCell': 'iceVolume'})
        dsAreaSum['iceThickness'] = (dsAreaSum.iceVolume /
                                     dsMesh.areaCell.sum('nCells'))

        dsAreaSum['iceArea'].attrs['units'] = 'm$^2$'
        dsAreaSum['iceArea'].attrs['description'] = \
            'Total {} sea ice area'.format(hemisphere)
        dsAreaSum['iceVolume'].attrs['units'] = 'm$^3$'
        dsAreaSum['iceVolume'].attrs['description'] = \
            'Total {} sea ice volume'.format(hemisphere)
        dsAreaSum['iceThickness'].attrs['units'] = 'm'
        dsAreaSum['iceThickness'].attrs['description'] = \
            'Mean {} sea ice volume'.format(hemisphere)

        return dsAreaSum

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, namelistMap, \
        streamMap, variableMap, plotsDirectory, simulationStartTime, \
        restartFileName = setup_sea_ice_task(config)

    # get a list of timeSeriesStatsMonthly output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('timeSeries', 'startDate')
    endDate = config.get('timeSeries', 'endDate')
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    fileNames = historyStreams.readpath(streamName, startDate=startDate,
                                        endDate=endDate,  calendar=calendar)
    print '\n  Reading files:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(fileNames[0]),
              os.path.basename(fileNames[-1]))

    plotTitles = {'iceArea': 'Sea-ice area',
                  'iceVolume': 'Sea-ice volume',
                  'iceThickness': 'Sea-ice mean thickness'}

    units = {'iceArea': '[km$^2$]',
             'iceVolume': '[10$^3$ km$^3$]',
             'iceThickness': '[m]'}

    obsFileNames = {
        'iceArea': {'NH': build_config_full_path(config, 'seaIceObservations',
                                                 'areaNH'),
                    'SH': build_config_full_path(config, 'seaIceObservations',
                                                 'areaSH')},
        'iceVolume': {'NH': build_config_full_path(config,
                                                   'seaIceObservations',
                                                   'volNH'),
                      'SH': build_config_full_path(config,
                                                   'seaIceObservations',
                                                   'volSH')}}

    # Some plotting rules
    titleFontSize = config.get('timeSeriesSeaIceAreaVol', 'titleFontSize')

    mainRunName = config.get('runs', 'mainRunName')
    preprocessedReferenceRunName = config.get('runs',
                                              'preprocessedReferenceRunName')
    preprocessedReferenceDirectory = config.get('seaIcePreprocessedReference',
                                                'baseDirectory')

    compareWithObservations = config.getboolean('timeSeriesSeaIceAreaVol',
                                                'compareWithObservations')

    movingAveragePoints = config.getint('timeSeriesSeaIceAreaVol',
                                        'movingAveragePoints')

    polarPlot = config.getboolean('timeSeriesSeaIceAreaVol', 'polarPlot')

    outputDirectory = build_config_full_path(config, 'output',
                                             'timeseriesSubdirectory')

    make_directories(outputDirectory)

    print '  Load sea-ice data...'
    # Load mesh
    dsMesh = xr.open_dataset(restartFileName)
    dsMesh = subset_variables(dsMesh,
                              variableList=['lonCell', 'latCell', 'areaCell'])

    # Load data
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

    yearStart = days_to_datetime(ds.Time.min(), calendar=calendar).year
    yearEnd = days_to_datetime(ds.Time.max(), calendar=calendar).year
    timeStart = date_to_days(year=yearStart, month=1, day=1,
                             calendar=calendar)
    timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                           calendar=calendar)

    if preprocessedReferenceRunName != 'None':
        inFilesPreprocessed = '{}/icevol.{}.year*.nc'.format(
            preprocessedReferenceDirectory, preprocessedReferenceRunName)
        dsPreprocessed = open_multifile_dataset(
            fileNames=inFilesPreprocessed,
            calendar=calendar,
            config=config,
            timeVariableName='xtime')
        preprocessedYearEnd = days_to_datetime(dsPreprocessed.Time.max(),
                                               calendar=calendar).year
        if yearStart <= preprocessedYearEnd:
            dsPreprocessedTimeSlice = dsPreprocessed.sel(Time=slice(timeStart,
                                                                    timeEnd))
        else:
            print '   Warning: Preprocessed time series ends before the ' \
                'timeSeries startYear and will not be plotted.'
            preprocessedReferenceRunName = 'None'

    norm = {'iceArea': 1e-6,  # m^2 to km^2
            'iceVolume': 1e-12,  # m^3 to 10^3 km^3
            'iceThickness': 1.}

    xLabel = 'Time [years]'

    dsTimeSeries = {}
    obs = {}
    preprocessed = {}
    figureNameStd = {}
    figureNamePolar = {}
    title = {}
    plotVars = {}

    for hemisphere in ['NH', 'SH']:
        print '   Caching {} data'.format(hemisphere)
        cacheFileName = '{}/seaIceAreaVolumeTimeSeries_{}.nc'.format(
            outputDirectory, hemisphere)

        dsTimeSeries[hemisphere] = time_series.cache_time_series(
            ds.Time.values, compute_area_vol_part, cacheFileName, calendar,
            yearsPerCacheUpdate=10, printProgress=True)

        print '  Make {} plots...'.format(hemisphere)

        for variableName in ['iceArea', 'iceVolume']:
            key = (hemisphere, variableName)

            # apply the norm to each variable
            plotVars[key] = (norm[variableName] *
                             dsTimeSeries[hemisphere][variableName])

            prefix = '{}/{}{}_{}'.format(plotsDirectory,
                                         variableName,
                                         hemisphere,
                                         mainRunName)

            figureNameStd[key] = '{}.png'.format(prefix)
            figureNamePolar[key] = '{}_polar.png'.format(prefix)

            title[key] = '{} ({}), {} (r)'.format(
                plotTitles[variableName], hemisphere, mainRunName)

        if compareWithObservations:
            key = (hemisphere, 'iceArea')
            title[key] = '{}\nSSM/I observations, annual cycle (k)'.format(
                title[key])
            if hemisphere == 'NH':
                key = (hemisphere, 'iceVolume')
                title[key] = '{}\nPIOMAS, annual cycle (k)'.format(title[key])

        if preprocessedReferenceRunName != 'None':
            for variableName in ['iceArea', 'iceVolume']:
                key = (hemisphere, variableName)
                title[key] = '{}\n {} (b)'.format(
                    title[key], preprocessedReferenceRunName)

        if compareWithObservations:
            dsObs = open_multifile_dataset(
                fileNames=obsFileNames['iceArea'][hemisphere],
                calendar=calendar,
                config=config,
                timeVariableName='xtime')
            key = (hemisphere, 'iceArea')
            obs[key] = replicate_cycle(plotVars[key], dsObs.IceArea, calendar)

            key = (hemisphere, 'iceVolume')
            if hemisphere == 'NH':
                dsObs = open_multifile_dataset(
                    fileNames=obsFileNames['iceVolume'][hemisphere],
                    calendar=calendar,
                    config=config,
                    timeVariableName='xtime')
                obs[key] = replicate_cycle(plotVars[key], dsObs.IceVol,
                                           calendar)
            else:
                obs[key] = None

        if preprocessedReferenceRunName != 'None':
            inFilesPreprocessed = '{}/icearea.{}.year*.nc'.format(
                preprocessedReferenceDirectory,
                preprocessedReferenceRunName)
            dsPreprocessed = open_multifile_dataset(
                fileNames=inFilesPreprocessed,
                calendar=calendar,
                config=config,
                timeVariableName='xtime')
            dsPreprocessedTimeSlice = dsPreprocessed.sel(
                Time=slice(timeStart, timeEnd))
            key = (hemisphere, 'iceArea')
            preprocessed[key] = dsPreprocessedTimeSlice[
                'icearea_{}'.format(hemisphere.lower())]

            inFilesPreprocessed = '{}/icevol.{}.year*.nc'.format(
                preprocessedReferenceDirectory,
                preprocessedReferenceRunName)
            dsPreprocessed = open_multifile_dataset(
                fileNames=inFilesPreprocessed,
                calendar=calendar,
                config=config,
                timeVariableName='xtime')
            dsPreprocessedTimeSlice = dsPreprocessed.sel(
                Time=slice(timeStart, timeEnd))
            key = (hemisphere, 'iceVolume')
            preprocessed[key] = dsPreprocessedTimeSlice[
                'icevolume_{}'.format(hemisphere.lower())]

        for variableName in ['iceArea', 'iceVolume']:
            key = (hemisphere, variableName)
            if compareWithObservations:
                if preprocessedReferenceRunName != 'None':
                    plotVars[key] = [plotVars[key], obs[key],
                                     preprocessed[key]]
                    lineStyles = ['r-', 'k-', 'b-']
                    lineWidths = [1.2, 1.2, 1.2]
                else:
                    # just v1 model and obs
                    plotVars[key] = [plotVars[key], obs[key]]
                    lineStyles = ['r-', 'k-']
                    lineWidths = [1.2, 1.2]
            elif preprocessedReferenceRunName != 'None':
                # just v1 and v0 models
                plotVars[key] = [plotVars[key], preprocessed[key]]
                lineStyles = ['r-', 'b-']
                lineWidths = [1.2, 1.2]

            if (compareWithObservations or
                    preprocessedReferenceRunName != 'None'):
                # separate plots for nothern and southern hemispheres
                timeseries_analysis_plot(config, plotVars[key],
                                         movingAveragePoints,
                                         title[key], xLabel,
                                         units[variableName],
                                         figureNameStd[key],
                                         lineStyles=lineStyles,
                                         lineWidths=lineWidths,
                                         titleFontSize=titleFontSize,
                                         calendar=calendar)
                if (polarPlot):
                    timeseries_analysis_plot_polar(config, plotVars[key],
                                                   movingAveragePoints,
                                                   title[key],
                                                   figureNamePolar[key],
                                                   lineStyles=lineStyles,
                                                   lineWidths=lineWidths,
                                                   titleFontSize=titleFontSize,
                                                   calendar=calendar)
    if (not compareWithObservations and
            preprocessedReferenceRunName == 'None'):
        for variableName in ['iceArea', 'iceVolume']:
            # we will combine north and south onto a single graph
            figureNameStd = '{}/{}.{}.png'.format(plotsDirectory,
                                                  mainRunName,
                                                  variableName)
            figureNamePolar = '{}/{}.{}_polar.png'.format(plotsDirectory,
                                                          mainRunName,
                                                          variableName)
            title = '{}, NH (r), SH (k)\n{}'.format(plotTitles[variableName],
                                                    mainRunName)
            varList = [plotVars[('NH', variableName)],
                       plotVars[('SH', variableName)]]
            timeseries_analysis_plot(config, varList,
                                     movingAveragePoints,
                                     title, xLabel, units[variableName],
                                     figureNameStd,
                                     lineStyles=['r-', 'k-'],
                                     lineWidths=[1.2, 1.2],
                                     titleFontSize=titleFontSize,
                                     calendar=calendar)
            if (polarPlot):
                timeseries_analysis_plot_polar(config, varList,
                                               movingAveragePoints,
                                               title, figureNamePolar,
                                               lineStyles=['r-', 'k-'],
                                               lineWidths=[1.2, 1.2],
                                               titleFontSize=titleFontSize,
                                               calendar=calendar)


def replicate_cycle(ds, dsToReplicate, calendar):
    """
    Replicates a periodic time series `dsToReplicate` to cover the timeframe
    of the dataset `ds`.

    Parameters
    ----------
    ds : dataset used to find the start and end time of the replicated cycle

    dsToReplicate : dataset to replicate.  The period of the cycle is the
        length of dsToReplicate plus the time between the first two time
        values (typically one year total).

    calendar : {'gregorian', 'gregorian_noleap'}
        The name of one of the calendars supported by MPAS cores

    Returns:
    --------
    dsShift : a cyclicly repeated version of `dsToReplicte` covering the range
        of time of `ds`.

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    02/22/2017
    """
    dsStartTime = days_to_datetime(ds.Time.min(), calendar=calendar)
    dsEndTime = days_to_datetime(ds.Time.max(), calendar=calendar)
    repStartTime = days_to_datetime(dsToReplicate.Time.min(),
                                    calendar=calendar)
    repEndTime = days_to_datetime(dsToReplicate.Time.max(),
                                  calendar=calendar)

    repSecondTime = days_to_datetime(dsToReplicate.Time.isel(Time=1),
                                     calendar=calendar)

    period = (MpasRelativeDelta(repEndTime, repStartTime) +
              MpasRelativeDelta(repSecondTime, repStartTime))

    startIndex = 0
    while(dsStartTime > repStartTime + (startIndex+1)*period):
        startIndex += 1

    endIndex = 0
    while(dsEndTime > repEndTime + endIndex*period):
        endIndex += 1

    dsShift = dsToReplicate.copy()

    times = days_to_datetime(dsShift.Time, calendar=calendar)
    dsShift.coords['Time'] = ('Time',
                              datetime_to_days(times + startIndex*period,
                                               calendar=calendar))
    # replicate cycle:
    for cycleIndex in range(startIndex, endIndex):
        dsNew = dsToReplicate.copy()
        dsNew.coords['Time'] = ('Time',
                                datetime_to_days(times + (cycleIndex+1)*period,
                                                 calendar=calendar))
        dsShift = xr.concat([dsShift, dsNew], dim='Time')

    # clip dsShift to the range of ds
    dsStartTime = dsShift.Time.sel(Time=ds.Time.min(), method='nearest').values
    dsEndTime = dsShift.Time.sel(Time=ds.Time.max(), method='nearest').values
    dsShift = dsShift.sel(Time=slice(dsStartTime, dsEndTime))

    return dsShift
