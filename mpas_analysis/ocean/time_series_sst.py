# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot.plotting import timeseries_analysis_plot

from mpas_analysis.shared.generalized_reader import open_multifile_dataset
from mpas_analysis.shared.io import open_mpas_dataset

from mpas_analysis.shared.timekeeping.utility import date_to_days, \
    days_to_datetime

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, check_path_exists
from mpas_analysis.shared.html import write_image_xml


class TimeSeriesSST(AnalysisTask):
    """
    Performs analysis of the time-series output of sea-surface temperature
    (SST).

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    refConfig :  ``MpasAnalysisConfigParser``
        Configuration options for a reference run (if any)
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasTimeSeriesTask, refConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesSST, self).__init__(
            config=config,
            taskName='timeSeriesSST',
            componentName='ocean',
            tags=['timeSeries', 'sst', 'publicObs'])

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.refConfig = refConfig

        self.run_after(mpasTimeSeriesTask)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        OSError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

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
        regions = config.getExpression('timeSeriesSST', 'regions')

        self.xmlFileNames = []
        self.filePrefixes = {}

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
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani

        self.logger.info("\nPlotting SST time series...")

        self.logger.info('  Load SST data...')

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
        regionsToPlot = config.getExpression('timeSeriesSST', 'regions')

        regionIndicesToPlot = [regions.index(region) for region in
                               regionsToPlot]

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        dsSST = open_mpas_dataset(fileName=self.inputFile,
                                  calendar=calendar,
                                  variableList=self.variableList,
                                  startDate=self.startDate,
                                  endDate=self.endDate)

        yearStart = days_to_datetime(dsSST.Time.min(), calendar=calendar).year
        yearEnd = days_to_datetime(dsSST.Time.max(), calendar=calendar).year
        timeStart = date_to_days(year=yearStart, month=1, day=1,
                                 calendar=calendar)
        timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                               calendar=calendar)

        if self.refConfig is not None:
            baseDirectory = build_config_full_path(
                self.refConfig, 'output', 'timeSeriesSubdirectory')

            refFileName = '{}/{}.nc'.format(
                    baseDirectory, self.mpasTimeSeriesTask.fullTaskName)

            refStartYear = self.refConfig.getint('timeSeries', 'startYear')
            refEndYear = self.refConfig.getint('timeSeries', 'endYear')
            refStartDate = '{:04d}-01-01_00:00:00'.format(refStartYear)
            refEndDate = '{:04d}-12-31_23:59:59'.format(refEndYear)

            dsRefSST = open_mpas_dataset(
                    fileName=refFileName,
                    calendar=calendar,
                    variableList=self.variableList,
                    startDate=refStartDate,
                    endDate=refEndDate)
        else:
            dsRefSST = None

        if preprocessedReferenceRunName != 'None':
            self.logger.info('  Load in SST for a preprocesses reference '
                             'run...')
            inFilesPreprocessed = '{}/SST.{}.year*.nc'.format(
                preprocessedInputDirectory, preprocessedReferenceRunName)
            dsPreprocessed = open_multifile_dataset(
                fileNames=inFilesPreprocessed,
                calendar=calendar,
                config=config,
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

            title = '{} SST'.format(plotTitles[regionIndex])
            xLabel = 'Time [years]'
            yLabel = '[$\degree$C]'

            varName = self.variableList[0]
            SST = dsSST[varName].isel(nOceanRegions=regionIndex)

            filePrefix = self.filePrefixes[region]

            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            lineColors = ['k']
            lineWidths = [3]

            fields = [SST]
            legendText = [mainRunName]

            if dsRefSST is not None:
                refSST = dsRefSST[varName].isel(nOceanRegions=regionIndex)
                fields.append(refSST)
                lineColors.append('r')
                lineWidths.append(1.5)
                refRunName = self.refConfig.get('runs', 'mainRunName')
                legendText.append(refRunName)

            if preprocessedReferenceRunName != 'None':
                SST_v0 = dsPreprocessedTimeSlice.SST
                fields.append(SST_v0)
                lineColors.append('purple')
                lineWidths.append(1.5)
                legendText.append(preprocessedReferenceRunName)

            if config.has_option(self.taskName, 'firstYearXTicks'):
                firstYearXTicks = config.getint(self.taskName,
                                                'firstYearXTicks')
            else:
                firstYearXTicks = None

            if config.has_option(self.taskName, 'yearStrideXTicks'):
                yearStrideXTicks = config.getint(self.taskName,
                                                 'yearStrideXTicks')
            else:
                yearStrideXTicks = None

            timeseries_analysis_plot(config, fields, movingAveragePoints,
                                     title, xLabel, yLabel, figureName,
                                     calendar=calendar,
                                     lineColors=lineColors,
                                     lineWidths=lineWidths,
                                     legendText=legendText,
                                     firstYearXTicks=firstYearXTicks,
                                     yearStrideXTicks=yearStrideXTicks)

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
