import os

from ..shared.analysis_task import AnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime

from ..shared.time_series import time_series

from ..shared.io import build_config_full_path, make_directories, \
    check_path_exists


class TimeSeriesSST(AnalysisTask):
    """
    Performs analysis of the time-series output of sea-surface temperature
    (SST).

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesSST, self).__init__(
            config=config,
            taskName='timeSeriesSST',
            componentName='ocean',
            tags=['timeSeries', 'sst'])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        OSError
            If files are not present

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(TimeSeriesSST, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config

        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
                check_path_exists(config.get('oceanPreprocessedReference',
                                             'baseDirectory'))
        return  # }}}

    def run(self):  # {{{
        """
        Performs analysis of the time-series output of sea-surface temperature
        (SST).

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting 2-d maps of SST climatologies..."

        print '  Load SST data...'

        simulationStartTime = get_simulation_start_time(self.runStreams)
        config = self.config
        calendar = self.calendar

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = config.get('timeSeries', 'startDate')
        endDate = config.get('timeSeries', 'endDate')
        streamName = \
            self.historyStreams.find_stream(self.streamMap['timeSeriesStats'])
        fileNames = self.historyStreams.readpath(streamName,
                                                 startDate=startDate,
                                                 endDate=endDate,
                                                 calendar=calendar)
        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(fileNames[0]),
                  os.path.basename(fileNames[-1]))

        mainRunName = config.get('runs', 'mainRunName')
        preprocessedReferenceRunName = \
            config.get('runs', 'preprocessedReferenceRunName')
        preprocessedInputDirectory = config.get('oceanPreprocessedReference',
                                                'baseDirectory')

        movingAveragePoints = config.getint('timeSeriesSST',
                                            'movingAveragePoints')

        regions = config.getExpression('regions', 'regions')
        plotTitles = config.getExpression('regions', 'plotTitles')
        regionIndicesToPlot = config.getExpression('timeSeriesSST',
                                                   'regionIndicesToPlot')

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        regionNames = config.getExpression('regions', 'regions')
        regionNames = [regionNames[index] for index in regionIndicesToPlot]

        # Load data:
        varList = ['avgSurfaceTemperature']
        ds = open_multifile_dataset(fileNames=fileNames,
                                    calendar=calendar,
                                    config=config,
                                    simulationStartTime=simulationStartTime,
                                    timeVariableName='Time',
                                    variableList=varList,
                                    variableMap=self.variableMap,
                                    startDate=startDate,
                                    endDate=endDate)

        ds = ds.isel(nOceanRegions=regionIndicesToPlot)

        yearStart = days_to_datetime(ds.Time.min(), calendar=calendar).year
        yearEnd = days_to_datetime(ds.Time.max(), calendar=calendar).year
        timeStart = date_to_days(year=yearStart, month=1, day=1,
                                 calendar=calendar)
        timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                               calendar=calendar)

        if preprocessedReferenceRunName != 'None':
            print '  Load in SST for a preprocesses reference run...'
            inFilesPreprocessed = '{}/SST.{}.year*.nc'.format(
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
                dsPreprocessedTimeSlice = \
                    dsPreprocessed.sel(Time=slice(timeStart, timeEnd))
            else:
                print '   Warning: Preprocessed time series ends before the ' \
                    'timeSeries startYear and will not be plotted.'
                preprocessedReferenceRunName = 'None'

        cacheFileName = '{}/sstTimeSeries.nc'.format(outputDirectory)

        # save ds so it's avaliable in _compute_sst_part
        self.ds = ds
        dsSST = time_series.cache_time_series(ds.Time.values,
                                              self._compute_sst_part,
                                              cacheFileName, calendar,
                                              yearsPerCacheUpdate=10,
                                              printProgress=True)

        print '  Make plots...'
        for index, regionIndex in enumerate(regionIndicesToPlot):

            title = plotTitles[regionIndex]
            title = 'SST, %s, %s (r-)' % (title, mainRunName)
            xLabel = 'Time [years]'
            yLabel = '[$^\circ$ C]'

            SST = dsSST.avgSurfaceTemperature.isel(nOceanRegions=index)

            figureName = '{}/sst_{}_{}.png'.format(self.plotsDirectory,
                                                   regions[regionIndex],
                                                   mainRunName)

            if preprocessedReferenceRunName != 'None':
                SST_v0 = dsPreprocessedTimeSlice.SST

                title = '{}\n {} (b-)'.format(title,
                                              preprocessedReferenceRunName)
                timeseries_analysis_plot(config, [SST, SST_v0],
                                         movingAveragePoints,
                                         title, xLabel, yLabel, figureName,
                                         lineStyles=['r-', 'b-'],
                                         lineWidths=[1.2, 1.2],
                                         calendar=calendar)
            else:
                timeseries_analysis_plot(config, [SST], movingAveragePoints,
                                         title, xLabel, yLabel, figureName,
                                         lineStyles=['r-'], lineWidths=[1.2],
                                         calendar=calendar)
        # }}}

    def _compute_sst_part(self, timeIndices, firstCall):  # {{{
        '''
        Compute part of the SST time series, given time indices to process.
        '''
        dsLocal = self.ds.isel(Time=timeIndices)
        return dsLocal
        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
