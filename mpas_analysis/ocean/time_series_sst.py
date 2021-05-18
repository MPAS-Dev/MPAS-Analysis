# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig

from mpas_analysis.shared.time_series import combine_time_series_with_ncrcat
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

    controlConfig :  ``MpasAnalysisConfigParser``
        Configuration options for a control run (if any)
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasTimeSeriesTask, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
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
        self.controlConfig = controlConfig

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

        if self.controlConfig is not None:
            baseDirectory = build_config_full_path(
                self.controlConfig, 'output', 'timeSeriesSubdirectory')

            controlFileName = '{}/{}.nc'.format(
                baseDirectory, self.mpasTimeSeriesTask.fullTaskName)

            controlStartYear = self.controlConfig.getint(
                'timeSeries', 'startYear')
            controlEndYear = self.controlConfig.getint('timeSeries', 'endYear')
            controlStartDate = '{:04d}-01-01_00:00:00'.format(controlStartYear)
            controlEndDate = '{:04d}-12-31_23:59:59'.format(controlEndYear)

            dsRefSST = open_mpas_dataset(
                fileName=controlFileName,
                calendar=calendar,
                variableList=self.variableList,
                startDate=controlStartDate,
                endDate=controlEndDate)
        else:
            dsRefSST = None

        if preprocessedReferenceRunName != 'None':
            self.logger.info('  Load in SST for a preprocesses reference '
                             'run...')
            inFilesPreprocessed = '{}/SST.{}.year*.nc'.format(
                preprocessedInputDirectory, preprocessedReferenceRunName)

            outFolder = '{}/preprocessed'.format(outputDirectory)
            make_directories(outFolder)
            outFileName = '{}/sst.nc'.format(outFolder)

            combine_time_series_with_ncrcat(inFilesPreprocessed,
                                            outFileName, logger=self.logger)
            dsPreprocessed = open_mpas_dataset(fileName=outFileName,
                                               calendar=calendar,
                                               timeVariableNames='xtime')
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
            yLabel = r'[$\degree$C]'

            varName = self.variableList[0]
            SST = dsSST[varName].isel(nOceanRegions=regionIndex)

            filePrefix = self.filePrefixes[region]

            outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            lineColors = [config.get('timeSeries', 'mainColor')]
            lineWidths = [3]

            fields = [SST]
            legendText = [mainRunName]

            if dsRefSST is not None:
                refSST = dsRefSST[varName].isel(nOceanRegions=regionIndex)
                fields.append(refSST)
                lineColors.append(config.get('timeSeries', 'controlColor'))
                lineWidths.append(1.5)
                controlRunName = self.controlConfig.get('runs', 'mainRunName')
                legendText.append(controlRunName)

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

            timeseries_analysis_plot(config, fields, calendar=calendar,
                                     title=title, xlabel=xLabel, ylabel=yLabel,
                                     movingAveragePoints=movingAveragePoints,
                                     lineColors=lineColors,
                                     lineWidths=lineWidths,
                                     legendText=legendText,
                                     firstYearXTicks=firstYearXTicks,
                                     yearStrideXTicks=yearStrideXTicks)

            savefig(outFileName)

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
