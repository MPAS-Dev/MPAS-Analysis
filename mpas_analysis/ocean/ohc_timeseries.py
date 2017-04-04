import numpy as np
import netCDF4

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime, string_to_datetime

from ..shared.analysis_task import setup_task


def ohc_timeseries(config, streamMap=None, variableMap=None):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.
    config is an instance of an MpasAnalysisConfigParser containing
    configuration options.

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Author: Xylar Asay-Davis, Milena Veneziani
    Last Modified: 03/23/2017
    """

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    simulationStartTime = get_simulation_start_time(runStreams)

    # read parameters from config file
    mainRunName = config.get('runs', 'mainRunName')
    preprocessedReferenceRunName = config.get('runs',
                                              'preprocessedReferenceRunName')
    preprocessedInputDirectory = config.get('oceanPreprocessedReference',
                                            'baseDirectory')

    compareWithObservations = config.getboolean('timeSeriesOHC',
                                                'compareWithObservations')

    movingAveragePoints = config.getint('timeSeriesOHC', 'movingAveragePoints')

    regions = config.getExpression('regions', 'regions')
    plotTitles = config.getExpression('regions', 'plotTitles')
    regionIndicesToPlot = config.getExpression('timeSeriesOHC',
                                               'regionIndicesToPlot')

    # Note: input file, not a mesh file because we need dycore specific fields
    # such as refBottomDepth and namelist fields such as config_density0, as
    # well as simulationStartTime, that are not guaranteed to be in the mesh
    # file.
    try:
        restartFile = runStreams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for OHC calculation')

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('timeSeries', 'startDate')
    endDate = config.get('timeSeries', 'endDate')
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    fileNames = historyStreams.readpath(streamName, startDate=startDate,
                                        endDate=endDate, calendar=calendar)
    print 'Reading files {} through {}'.format(fileNames[0], fileNames[-1])

    # Define/read in general variables
    print '  Read in depth and compute specific depth indexes...'
    ncFile = netCDF4.Dataset(restartFile, mode='r')
    # reference depth [m]
    depth = ncFile.variables['refBottomDepth'][:]
    ncFile.close()
    # specific heat [J/(kg*degC)]
    cp = namelist.getfloat('config_specific_heat_sea_water')
    # [kg/m3]
    rho = namelist.getfloat('config_density0')
    factor = 1e-22*rho*cp

    k700m = np.where(depth > 700.)[0][0] - 1
    k2000m = np.where(depth > 2000.)[0][0] - 1

    kbtm = len(depth)-1

    # Load data
    print '  Load ocean data...'
    variableList = ['avgLayerTemperature',
                    'sumLayerMaskValue',
                    'avgLayerArea',
                    'avgLayerThickness']
    ds = open_multifile_dataset(fileNames=fileNames,
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=variableList,
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    timeStart = string_to_datetime(startDate)
    timeEnd = string_to_datetime(endDate)

    # Select year-1 data and average it (for later computing anomalies)
    timeStartFirstYear = string_to_datetime(simulationStartTime)
    if timeStartFirstYear < timeStart:
        startDateFirstYear = simulationStartTime
        firstYear = int(startDateFirstYear[0:4])
        endDateFirstYear = '{:04d}-12-31_23:59:59'.format(firstYear)
        filesFirstYear = historyStreams.readpath(streamName,
                                                 startDate=startDateFirstYear,
                                                 endDate=endDateFirstYear,
                                                 calendar=calendar)
        dsFirstYear = open_multifile_dataset(
            fileNames=filesFirstYear,
            calendar=calendar,
            config=config,
            simulationStartTime=simulationStartTime,
            timeVariableName='Time',
            variableList=variableList,
            variableMap=variableMap,
            startDate=startDateFirstYear,
            endDate=endDateFirstYear)
    else:
        dsFirstYear = ds
        firstYear = timeStart.year

    timeStartFirstYear = date_to_days(year=firstYear, month=1, day=1,
                                      calendar=calendar)
    timeEndFirstYear = date_to_days(year=firstYear, month=12, day=31,
                                    hour=23, minute=59, second=59,
                                    calendar=calendar)

    dsFirstYear = dsFirstYear.sel(Time=slice(timeStartFirstYear,
                                             timeEndFirstYear))

    meanFirstYear = dsFirstYear.mean('Time')

    print '  Compute temperature anomalies...'
    avgLayerTemperature = ds.avgLayerTemperature
    avgLayerTemperatureFirstYear = meanFirstYear.avgLayerTemperature

    avgLayTemperatureAnomaly = (avgLayerTemperature -
                                avgLayerTemperatureFirstYear)

    yearStart = days_to_datetime(ds.Time.min(), calendar=calendar).year
    yearEnd = days_to_datetime(ds.Time.max(), calendar=calendar).year
    timeStart = date_to_days(year=yearStart, month=1, day=1,
                             calendar=calendar)
    timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                           calendar=calendar)

    if preprocessedReferenceRunName != 'None':
        print '  Load in OHC from preprocessed reference run...'
        inFilesPreprocessed = '{}/OHC.{}.year*.nc'.format(
            preprocessedInputDirectory, preprocessedReferenceRunName)
        dsPreprocessed = open_multifile_dataset(
            fileNames=inFilesPreprocessed,
            calendar=calendar,
            config=config,
            simulationStartTime=simulationStartTime,
            timeVariableName='xtime')
        yearEndPreprocessed = days_to_datetime(dsPreprocessed.Time.max(),
                                               calendar=calendar).year
        if yearStart <= yearEndPreprocessed:
            dsPreprocessedTimeSlice = dsPreprocessed.sel(Time=slice(timeStart,
                                                                    timeEnd))
        else:
            print '   Warning: Preprocessed time series ends before the ' \
                'timeSeries startYear and will not be plotted.'
            preprocessedReferenceRunName = 'None'

    sumLayerMaskValue = ds.sumLayerMaskValue
    avgLayerArea = ds.avgLayerArea
    avgLayerThickness = ds.avgLayerThickness

    print '  Compute OHC and make plots...'
    for index in range(len(regionIndicesToPlot)):
        regionIndex = regionIndicesToPlot[index]

        # Compute volume of each layer in the region:
        layerArea = sumLayerMaskValue[:, regionIndex, :] * \
            avgLayerArea[:, regionIndex, :]
        layerVolume = layerArea * avgLayerThickness[:, regionIndex, :]

        # Compute OHC:
        ohc = layerVolume * avgLayTemperatureAnomaly[:, regionIndex, :]
        # OHC over 0-bottom depth range:
        ohcTotal = ohc.sum('nVertLevels')
        ohcTotal = factor*ohcTotal

        # OHC over 0-700m depth range:
        ohc700m = factor*ohc[:, 0:k700m].sum('nVertLevels')

        # OHC over 700m-2000m depth range:
        ohc2000m = factor*ohc[:, k700m+1:k2000m].sum('nVertLevels')

        # OHC over 2000m-bottom depth range:
        ohcBottom = ohc[:, k2000m+1:kbtm].sum('nVertLevels')
        ohcBottom = factor*ohcBottom

        title = 'OHC, {}, 0-bottom (thick-), 0-700m (thin-), 700-2000m (--),' \
                ' 2000m-bottom (-.) \n {}'.format(plotTitles[regionIndex],
                                                  mainRunName)

        xLabel = 'Time [years]'
        yLabel = '[x$10^{22}$ J]'

        figureName = '{}/ohc_{}_{}.png'.format(plotsDirectory,
                                               regions[regionIndex],
                                               mainRunName)

        if preprocessedReferenceRunName != 'None':
            ohcPreprocessedTotal = dsPreprocessedTimeSlice.ohc_tot
            ohcPreprocessed700m = dsPreprocessedTimeSlice.ohc_700m
            ohcPreprocessed2000m = dsPreprocessedTimeSlice.ohc_2000m
            ohcPreprocessedBottom = dsPreprocessedTimeSlice.ohc_btm
            title = '{} (r), {} (b)'.format(title,
                                            preprocessedReferenceRunName)
            timeseries_analysis_plot(config, [ohcTotal, ohc700m, ohc2000m,
                                              ohcBottom, ohcPreprocessedTotal,
                                              ohcPreprocessed700m,
                                              ohcPreprocessed2000m,
                                              ohcPreprocessedBottom],
                                     movingAveragePoints, title,
                                     xLabel, yLabel, figureName,
                                     lineStyles=['r-', 'r-', 'r--', 'r-.',
                                                 'b-', 'b-', 'b--', 'b-.'],
                                     lineWidths=[2, 1, 1.5, 1.5, 2, 1, 1.5,
                                                 1.5],
                                     calendar=calendar)

        if (not compareWithObservations and
                preprocessedReferenceRunName == 'None'):
            timeseries_analysis_plot(config, [ohcTotal, ohc700m, ohc2000m,
                                              ohcBottom],
                                     movingAveragePoints, title,
                                     xLabel, yLabel, figureName,
                                     lineStyles=['r-', 'r-', 'r--', 'r-.'],
                                     lineWidths=[2, 1, 1.5, 1.5],
                                     calendar=calendar)
