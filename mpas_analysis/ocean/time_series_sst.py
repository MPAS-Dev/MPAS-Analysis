from ..shared import AnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime

from ..shared.io.utility import build_config_full_path, make_directories, \
    check_path_exists
from ..shared.html import write_image_xml


class TimeSeriesSST(AnalysisTask):
    """
    Performs analysis of the time-series output of sea-surface temperature
    (SST).

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, config, mpasTimeSeriesTask):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

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

        self.mpasTimeSeriesTask = mpasTimeSeriesTask

        self.run_after(mpasTimeSeriesTask)

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

        config = self.config

        self.startDate = self.config.get('timeSeries', 'startDate')
        self.endDate = self.config.get('timeSeries', 'endDate')

        self.variableList = \
            ['timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']
        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
                check_path_exists(config.get('oceanPreprocessedReference',
                                             'baseDirectory'))

        self.inputFile = self.mpasTimeSeriesTask.outputFile

        mainRunName = config.get('runs', 'mainRunName')
        regions = config.getExpression('regions', 'regions')
        regionIndicesToPlot = config.getExpression('timeSeriesSST',
                                                   'regionIndicesToPlot')

        self.xmlFileNames = []
        self.filePrefixes = {}

        regions = [regions[index] for index in regionIndicesToPlot]

        for region in regions:
            filePrefix = 'sst_{}_{}'.format(region, mainRunName)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[region] = filePrefix

        return  # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of the time-series output of sea-surface temperature
        (SST).

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        self.logger.info("\nPlotting SST time series...")

        self.logger.info('  Load SST data...')

        simulationStartTime = get_simulation_start_time(self.runStreams)
        config = self.config
        calendar = self.calendar

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

        dsSST = open_multifile_dataset(fileNames=self.inputFile,
                                       calendar=calendar,
                                       config=config,
                                       simulationStartTime=simulationStartTime,
                                       timeVariableName=['xtime_startMonthly',
                                                         'xtime_endMonthly'],
                                       variableList=self.variableList,
                                       startDate=self.startDate,
                                       endDate=self.endDate)

        yearStart = days_to_datetime(dsSST.Time.min(), calendar=calendar).year
        yearEnd = days_to_datetime(dsSST.Time.max(), calendar=calendar).year
        timeStart = date_to_days(year=yearStart, month=1, day=1,
                                 calendar=calendar)
        timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                               calendar=calendar)

        if preprocessedReferenceRunName != 'None':
            self.logger.info('  Load in SST for a preprocesses reference '
                             'run...')
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
                self.logger.warning('Preprocessed time series ends before the '
                                    'timeSeries startYear and will not be '
                                    'plotted.')
                preprocessedReferenceRunName = 'None'

        self.logger.info('  Make plots...')
        for regionIndex in regionIndicesToPlot:
            region = regions[regionIndex]

            title = plotTitles[regionIndex]
            title = 'SST, %s \n %s (black)' % (title, mainRunName)
            xLabel = 'Time [years]'
            yLabel = '[$^\circ$C]'

            varName = self.variableList[0]
            SST = dsSST[varName].isel(nOceanRegions=regionIndex)

            filePrefix = self.filePrefixes[region]

            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            if preprocessedReferenceRunName != 'None':
                SST_v0 = dsPreprocessedTimeSlice.SST

                title = '{}\n {} (red)'.format(title,
                                              preprocessedReferenceRunName)
                timeseries_analysis_plot(config, [SST, SST_v0],
                                         movingAveragePoints,
                                         title, xLabel, yLabel, figureName,
                                         lineStyles=['k-', 'r-'],
                                         lineWidths=[3, 1.5],
                                         legendText=[None, None],
                                         calendar=calendar)
            else:
                timeseries_analysis_plot(config, [SST], movingAveragePoints,
                                         title, xLabel, yLabel, figureName,
                                         lineStyles=['k-'], lineWidths=[3],
                                         legendText=[None], calendar=calendar)

            caption = 'Running Mean of {} Sea Surface Temperature'.format(
                    region)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Time Series',
                groupLink='timeseries',
                thumbnailDescription='{} SST'.format(region),
                imageDescription=caption,
                imageCaption=caption)

        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
