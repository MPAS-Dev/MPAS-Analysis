import xarray as xr
import pandas as pd
import datetime

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, \
    remove_repeated_time_index

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.io import StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.timekeeping.Date import Date


def sst_timeseries(config, streamMap=None, variableMap=None):
    """
    Performs analysis of the time-series output of sea-surface temperature
    (SST).

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Author: Xylar Asay-Davis, Milena Veneziani
    Last Modified: 02/02/2017
    """

    # Define/read in general variables
    print '  Load SST data...'
    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('timeSeries', 'startDate')
    endDate = config.get('timeSeries', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    inFiles = streams.readpath(streamName, startDate=startDate,
                               endDate=endDate)
    print 'Reading files {} through {}'.format(inFiles[0], inFiles[-1])

    mainRunName = config.get('runs', 'mainRunName')
    preprocessedReferenceRunName = config.get('runs',
                                              'preprocessedReferenceRunName')
    preprocessedInputDirectory = config.get('oceanPreprocessedReference',
                                            'baseDirectory')

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')

    yearOffset = config.getint('time', 'yearOffset')

    movingAveragePoints = config.getint('timeSeriesSST', 'movingAveragePoints')

    regions = config.getExpression('regions', 'regions')
    plotTitles = config.getExpression('regions', 'plotTitles')
    regionIndicesToPlot = config.getExpression('timeSeriesSST',
                                               'regionIndicesToPlot')

    # Load data:
    varList = ['avgSurfaceTemperature']
    ds = xr.open_mfdataset(
        inFiles,
        preprocess=lambda x: preprocess_mpas(x, yearoffset=yearOffset,
                                             timestr='Time',
                                             onlyvars=varList,
                                             varmap=variableMap))
    ds = remove_repeated_time_index(ds)

    # convert the start and end dates to datetime objects using
    # the Date class, which ensures the results are within the
    # supported range
    timeStart = Date(startDate).to_datetime(yearOffset)
    timeEnd = Date(endDate).to_datetime(yearOffset)
    # select only the data in the specified range of years
    ds = ds.sel(Time=slice(timeStart, timeEnd))

    SSTregions = ds.avgSurfaceTemperature

    yearStart = (pd.to_datetime(ds.Time.min().values)).year
    yearEnd = (pd.to_datetime(ds.Time.max().values)).year
    timeStart = datetime.datetime(yearStart, 1, 1)
    timeEnd = datetime.datetime(yearEnd, 12, 31)

    if preprocessedReferenceRunName != 'None':
        print '  Load in SST for a preprocesses reference run...'
        inFilesPreprocessed = '{}/SST.{}.year*.nc'.format(
            preprocessedInputDirectory, preprocessedReferenceRunName)
        dsPreprocessed = xr.open_mfdataset(
            inFilesPreprocessed,
            preprocess=lambda x: preprocess_mpas(x, yearoffset=yearOffset))
        dsPreprocessed = remove_repeated_time_index(dsPreprocessed)
        yearEndPreprocessed = \
            (pd.to_datetime(dsPreprocessed.Time.max().values)).year
        if yearStart <= yearEndPreprocessed:
            dsPreprocessedTimeSlice = \
                dsPreprocessed.sel(Time=slice(timeStart, timeEnd))
        else:
            print '   Warning: Preprocessed time series ends before the ' \
                'timeSeries startYear and will not be plotted.'
            preprocessedReferenceRunName = 'None'

    print '  Make plots...'
    for index in range(len(regionIndicesToPlot)):
        regionIndex = regionIndicesToPlot[index]

        title = plotTitles[regionIndex]
        title = 'SST, %s, %s (r-)' % (title, mainRunName)
        xLabel = 'Time [years]'
        yLabel = '[$^\circ$ C]'

        SST = SSTregions[:, regionIndex]

        if preprocessedReferenceRunName != 'None':
            figureName = '{}/sst_{}_{}_{}.png'.format(
                plotsDirectory, regions[regionIndex], mainRunName,
                preprocessedReferenceRunName)
            SST_v0 = dsPreprocessedTimeSlice.SST

            title = '{}\n {} (b-)'.format(title, preprocessedReferenceRunName)
            timeseries_analysis_plot(config, [SST, SST_v0],
                                     movingAveragePoints,
                                     title, xLabel, yLabel, figureName,
                                     lineStyles=['r-', 'b-'],
                                     lineWidths=[1.2, 1.2])
        else:
            figureName = '{}/sst_{}_{}.png'.format(plotsDirectory,
                                                   regions[regionIndex],
                                                   mainRunName)
            timeseries_analysis_plot(config, [SST], movingAveragePoints, title,
                                     xLabel, yLabel, figureName,
                                     lineStyles=['r-'], lineWidths=[1.2])
