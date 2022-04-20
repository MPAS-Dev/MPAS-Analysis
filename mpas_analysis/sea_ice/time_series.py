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

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot import timeseries_analysis_plot, \
    timeseries_analysis_plot_polar, savefig

from mpas_analysis.shared.io.utility import build_config_full_path, \
    check_path_exists, make_directories, build_obs_path

from mpas_analysis.shared.timekeeping.utility import date_to_days, \
    days_to_datetime, datetime_to_days, get_simulation_start_time
from mpas_analysis.shared.timekeeping.MpasRelativeDelta import \
    MpasRelativeDelta

from mpas_analysis.shared.time_series import combine_time_series_with_ncrcat
from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.html import write_image_xml


class TimeSeriesSeaIce(AnalysisTask):
    """
    Performs analysis of time series of sea-ice properties.

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    controlconfig : mpas_tools.config.MpasConfigParser
        Configuration options for a control run (if any)

    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasTimeSeriesTask,
                 controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesSeaIce, self).__init__(
            config=config,
            taskName='timeSeriesSeaIceAreaVol',
            componentName='seaIce',
            tags=['timeSeries', 'publicObs', 'arctic', 'antarctic'])

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.controlConfig = controlConfig

        self.run_after(mpasTimeSeriesTask)

    def setup_and_check(self):
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
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(TimeSeriesSeaIce, self).setup_and_check()

        config = self.config

        self.startDate = self.config.get('timeSeries', 'startDate')
        self.endDate = self.config.get('timeSeries', 'endDate')

        self.variableList = ['timeMonthly_avg_iceAreaCell',
                             'timeMonthly_avg_iceVolumeCell']
        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        self.inputFile = self.mpasTimeSeriesTask.outputFile

        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
            check_path_exists(config.get('seaIcePreprocessedReference',
                                         'baseDirectory'))

        # get a list of timeSeriesStatsMonthly output files from the streams
        # file, reading only those that are between the start and end dates
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.startDate = config.get('timeSeries', 'startDate')
        self.endDate = config.get('timeSeries', 'endDate')
        self.inputFiles = \
            self.historyStreams.readpath(streamName,
                                         startDate=self.startDate,
                                         endDate=self.endDate,
                                         calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDate,
                                       self.endDate))

        self.simulationStartTime = get_simulation_start_time(self.runStreams)

        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-SeaIce restart file found: need at least '
                          'one restart file to perform remapping of '
                          'climatologies.')

        # these are redundant for now.  Later cleanup is needed where these
        # file names are reused in run()
        self.xmlFileNames = []

        polarPlot = config.getboolean('timeSeriesSeaIceAreaVol', 'polarPlot')
        mainRunName = config.get('runs', 'mainRunName')
        preprocessedReferenceRunName = \
            config.get('runs', 'preprocessedReferenceRunName')
        compareWithObservations = config.getboolean('timeSeriesSeaIceAreaVol',
                                                    'compareWithObservations')

        polarXMLFileNames = []

        if (not compareWithObservations and
                preprocessedReferenceRunName == 'None'):
            for variableName in ['iceArea', 'iceVolume']:
                filePrefix = '{}.{}'.format(mainRunName,
                                            variableName)

                self.xmlFileNames.append('{}/{}.xml'.format(
                    self.plotsDirectory, filePrefix))
                polarXMLFileNames.append('{}/{}_polar.xml'.format(
                    self.plotsDirectory, filePrefix))
        else:

            for hemisphere in ['NH', 'SH']:
                for variableName in ['iceArea', 'iceVolume']:
                    filePrefix = '{}{}_{}'.format(variableName,
                                                  hemisphere,
                                                  mainRunName)

                    self.xmlFileNames.append('{}/{}.xml'.format(
                        self.plotsDirectory, filePrefix))
                    polarXMLFileNames.append('{}/{}_polar.xml'.format(
                        self.plotsDirectory, filePrefix))

        if polarPlot:
            self.xmlFileNames.extend(polarXMLFileNames)
        return

    def run_task(self):
        """
        Performs analysis of time series of sea-ice properties.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani

        self.logger.info("\nPlotting sea-ice area and volume time series...")

        config = self.config
        calendar = self.calendar

        sectionName = self.taskName

        plotTitles = {'iceArea': 'Sea-ice area',
                      'iceVolume': 'Sea-ice volume',
                      'iceThickness': 'Sea-ice mean thickness'}

        units = {'iceArea': '[km$^2$]',
                 'iceVolume': '[10$^3$ km$^3$]',
                 'iceThickness': '[m]'}

        obsFileNames = {
            'iceArea': {'NH': build_obs_path(
                config, 'seaIce',
                relativePathOption='areaNH',
                relativePathSection=sectionName),
                'SH': build_obs_path(
                config, 'seaIce',
                relativePathOption='areaSH',
                relativePathSection=sectionName)},
            'iceVolume': {'NH': build_obs_path(
                config, 'seaIce',
                relativePathOption='volNH',
                relativePathSection=sectionName),
                'SH': build_obs_path(
                config, 'seaIce',
                relativePathOption='volSH',
                relativePathSection=sectionName)}}

        # Some plotting rules
        titleFontSize = config.get('timeSeriesSeaIceAreaVol', 'titleFontSize')

        mainRunName = config.get('runs', 'mainRunName')
        preprocessedReferenceRunName = \
            config.get('runs', 'preprocessedReferenceRunName')
        preprocessedReferenceDirectory = \
            config.get('seaIcePreprocessedReference', 'baseDirectory')

        compareWithObservations = config.getboolean('timeSeriesSeaIceAreaVol',
                                                    'compareWithObservations')

        movingAveragePoints = config.getint('timeSeriesSeaIceAreaVol',
                                            'movingAveragePoints')

        polarPlot = config.getboolean('timeSeriesSeaIceAreaVol', 'polarPlot')

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        self.logger.info('  Load sea-ice data...')
        # Load mesh

        dsTimeSeries = self._compute_area_vol()

        yearStart = days_to_datetime(dsTimeSeries['NH'].Time.min(),
                                     calendar=calendar).year
        yearEnd = days_to_datetime(dsTimeSeries['NH'].Time.max(),
                                   calendar=calendar).year
        timeStart = date_to_days(year=yearStart, month=1, day=1,
                                 calendar=calendar)
        timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                               calendar=calendar)

        if preprocessedReferenceRunName != 'None':
            # determine if we're beyond the end of the preprocessed data
            # (and go ahead and cache the data set while we're checking)
            outFolder = '{}/preprocessed'.format(outputDirectory)
            make_directories(outFolder)
            inFilesPreprocessed = '{}/icevol.{}.year*.nc'.format(
                preprocessedReferenceDirectory, preprocessedReferenceRunName)
            outFileName = '{}/iceVolume.nc'.format(outFolder)

            combine_time_series_with_ncrcat(inFilesPreprocessed,
                                            outFileName,
                                            logger=self.logger)
            dsPreprocessed = open_mpas_dataset(fileName=outFileName,
                                               calendar=calendar,
                                               timeVariableNames='xtime')
            preprocessedYearEnd = days_to_datetime(dsPreprocessed.Time.max(),
                                                   calendar=calendar).year
            if yearStart <= preprocessedYearEnd:
                dsPreprocessedTimeSlice = \
                    dsPreprocessed.sel(Time=slice(timeStart, timeEnd))
            else:
                self.logger.warning('Preprocessed time series ends before the '
                                    'timeSeries startYear and will not be '
                                    'plotted.')
                preprocessedReferenceRunName = 'None'

        if self.controlConfig is not None:

            dsTimeSeriesRef = {}
            baseDirectory = build_config_full_path(
                self.controlConfig, 'output', 'timeSeriesSubdirectory')

            controlRunName = self.controlConfig.get('runs', 'mainRunName')

            for hemisphere in ['NH', 'SH']:
                inFileName = '{}/seaIceAreaVol{}.nc'.format(baseDirectory,
                                                            hemisphere)

                dsTimeSeriesRef[hemisphere] = xr.open_dataset(inFileName)

        norm = {'iceArea': 1e-6,  # m^2 to km^2
                'iceVolume': 1e-12,  # m^3 to 10^3 km^3
                'iceThickness': 1.}

        xLabel = 'Time [years]'

        galleryGroup = 'Time Series'
        groupLink = 'timeseries'

        obs = {}
        preprocessed = {}
        figureNameStd = {}
        figureNamePolar = {}
        title = {}
        plotVars = {}
        obsLegend = {}
        plotVarsRef = {}

        for hemisphere in ['NH', 'SH']:

            self.logger.info('  Make {} plots...'.format(hemisphere))

            for variableName in ['iceArea', 'iceVolume']:
                key = (hemisphere, variableName)

                # apply the norm to each variable
                plotVars[key] = (norm[variableName] *
                                 dsTimeSeries[hemisphere][variableName])

                if self.controlConfig is not None:
                    plotVarsRef[key] = norm[variableName] * \
                        dsTimeSeriesRef[hemisphere][variableName]

                prefix = '{}/{}{}_{}'.format(self.plotsDirectory,
                                             variableName,
                                             hemisphere,
                                             mainRunName)

                figureNameStd[key] = '{}.png'.format(prefix)
                figureNamePolar[key] = '{}_polar.png'.format(prefix)

                title[key] = '{} ({})'.format(plotTitles[variableName],
                                              hemisphere)

            if compareWithObservations:
                key = (hemisphere, 'iceArea')
                obsLegend[key] = 'SSM/I observations, annual cycle '
                if hemisphere == 'NH':
                    key = (hemisphere, 'iceVolume')
                    obsLegend[key] = 'PIOMAS, annual cycle (blue)'

            if preprocessedReferenceRunName != 'None':
                for variableName in ['iceArea', 'iceVolume']:
                    key = (hemisphere, variableName)

            if compareWithObservations:

                outFolder = '{}/obs'.format(outputDirectory)
                make_directories(outFolder)
                outFileName = '{}/iceArea{}.nc'.format(outFolder, hemisphere)

                combine_time_series_with_ncrcat(
                    obsFileNames['iceArea'][hemisphere],
                    outFileName, logger=self.logger)
                dsObs = open_mpas_dataset(fileName=outFileName,
                                          calendar=calendar,
                                          timeVariableNames='xtime')
                key = (hemisphere, 'iceArea')
                obs[key] = self._replicate_cycle(plotVars[key], dsObs.IceArea,
                                                 calendar)

                key = (hemisphere, 'iceVolume')
                if hemisphere == 'NH':
                    outFileName = '{}/iceVolume{}.nc'.format(outFolder,
                                                             hemisphere)
                    combine_time_series_with_ncrcat(
                        obsFileNames['iceVolume'][hemisphere],
                        outFileName, logger=self.logger)
                    dsObs = open_mpas_dataset(fileName=outFileName,
                                              calendar=calendar,
                                              timeVariableNames='xtime')
                    obs[key] = self._replicate_cycle(plotVars[key],
                                                     dsObs.IceVol,
                                                     calendar)
                else:
                    obs[key] = None

            if preprocessedReferenceRunName != 'None':
                outFolder = '{}/preprocessed'.format(outputDirectory)
                inFilesPreprocessed = '{}/icearea.{}.year*.nc'.format(
                    preprocessedReferenceDirectory,
                    preprocessedReferenceRunName)

                outFileName = '{}/iceArea.nc'.format(outFolder)

                combine_time_series_with_ncrcat(inFilesPreprocessed,
                                                outFileName,
                                                logger=self.logger)
                dsPreprocessed = open_mpas_dataset(fileName=outFileName,
                                                   calendar=calendar,
                                                   timeVariableNames='xtime')
                dsPreprocessedTimeSlice = dsPreprocessed.sel(
                    Time=slice(timeStart, timeEnd))
                key = (hemisphere, 'iceArea')
                preprocessed[key] = dsPreprocessedTimeSlice[
                    'icearea_{}'.format(hemisphere.lower())]

                inFilesPreprocessed = '{}/icevol.{}.year*.nc'.format(
                    preprocessedReferenceDirectory,
                    preprocessedReferenceRunName)
                outFileName = '{}/iceVolume.nc'.format(outFolder)

                combine_time_series_with_ncrcat(inFilesPreprocessed,
                                                outFileName,
                                                logger=self.logger)
                dsPreprocessed = open_mpas_dataset(fileName=outFileName,
                                                   calendar=calendar,
                                                   timeVariableNames='xtime')
                dsPreprocessedTimeSlice = dsPreprocessed.sel(
                    Time=slice(timeStart, timeEnd))
                key = (hemisphere, 'iceVolume')
                preprocessed[key] = dsPreprocessedTimeSlice[
                    'icevolume_{}'.format(hemisphere.lower())]

            for variableName in ['iceArea', 'iceVolume']:
                key = (hemisphere, variableName)
                dsvalues = [plotVars[key]]
                legendText = [mainRunName]
                lineColors = [config.get('timeSeries', 'mainColor')]
                lineWidths = [3]
                if compareWithObservations and key in obsLegend.keys():
                    dsvalues.append(obs[key])
                    legendText.append(obsLegend[key])
                    lineColors.append(config.get('timeSeries', 'obsColor1'))
                    lineWidths.append(1.2)
                if preprocessedReferenceRunName != 'None':
                    dsvalues.append(preprocessed[key])
                    legendText.append(preprocessedReferenceRunName)
                    lineColors.append('purple')
                    lineWidths.append(1.2)

                if self.controlConfig is not None:
                    dsvalues.append(plotVarsRef[key])
                    legendText.append(controlRunName)
                    lineColors.append(config.get('timeSeries',
                                                 'controlColor'))
                    lineWidths.append(1.2)

                if config.has_option(sectionName, 'firstYearXTicks'):
                    firstYearXTicks = config.getint(sectionName,
                                                    'firstYearXTicks')
                else:
                    firstYearXTicks = None

                if config.has_option(sectionName, 'yearStrideXTicks'):
                    yearStrideXTicks = config.getint(sectionName,
                                                     'yearStrideXTicks')
                else:
                    yearStrideXTicks = None

                # separate plots for nothern and southern hemispheres
                timeseries_analysis_plot(
                    config, dsvalues, calendar=calendar, title=title[key],
                    xlabel=xLabel, ylabel=units[variableName],
                    movingAveragePoints=movingAveragePoints,
                    lineColors=lineColors, lineWidths=lineWidths,
                    legendText=legendText, titleFontSize=titleFontSize,
                    firstYearXTicks=firstYearXTicks,
                    yearStrideXTicks=yearStrideXTicks)

                savefig(figureNameStd[key], config)

                filePrefix = '{}{}_{}'.format(variableName,
                                              hemisphere,
                                              mainRunName)
                thumbnailDescription = '{} {}'.format(
                    hemisphere, plotTitles[variableName])
                caption = 'Running mean of {}'.format(
                    thumbnailDescription)
                write_image_xml(
                    config,
                    filePrefix,
                    componentName='Sea Ice',
                    componentSubdirectory='sea_ice',
                    galleryGroup=galleryGroup,
                    groupLink=groupLink,
                    thumbnailDescription=thumbnailDescription,
                    imageDescription=caption,
                    imageCaption=caption)

                if (polarPlot):
                    timeseries_analysis_plot_polar(config, dsvalues, title[key],
                                                   movingAveragePoints,
                                                   lineColors=lineColors,
                                                   lineWidths=lineWidths,
                                                   legendText=legendText,
                                                   titleFontSize=titleFontSize)

                    savefig(figureNamePolar[key], config)

                    filePrefix = '{}{}_{}_polar'.format(variableName,
                                                        hemisphere,
                                                        mainRunName)
                    write_image_xml(
                        config,
                        filePrefix,
                        componentName='Sea Ice',
                        componentSubdirectory='sea_ice',
                        galleryGroup=galleryGroup,
                        groupLink=groupLink,
                        thumbnailDescription=thumbnailDescription,
                        imageDescription=caption,
                        imageCaption=caption)

    def _replicate_cycle(self, ds, dsToReplicate, calendar):
        """
        Replicates a periodic time series `dsToReplicate` to cover the
        timeframe of the dataset `ds`.

        Parameters
        ----------
        ds : dataset used to find the start and end time of the replicated
            cycle

        dsToReplicate : dataset to replicate.  The period of the cycle is the
            length of dsToReplicate plus the time between the first two time
            values (typically one year total).

        calendar : {'gregorian', 'gregorian_noleap'}
            The name of one of the calendars supported by MPAS cores

        Returns:
        --------
        dsShift : a cyclicly repeated version of `dsToReplicte` covering the
            range of time of `ds`.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani

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
        while(dsStartTime > repStartTime + (startIndex + 1) * period):
            startIndex += 1

        endIndex = 0
        while(dsEndTime > repEndTime + endIndex * period):
            endIndex += 1

        dsShift = dsToReplicate.copy()

        times = days_to_datetime(dsShift.Time, calendar=calendar)
        dsShift.coords['Time'] = ('Time',
                                  datetime_to_days(times + startIndex * period,
                                                   calendar=calendar))
        # replicate cycle:
        for cycleIndex in range(startIndex, endIndex):
            dsNew = dsToReplicate.copy()
            dsNew.coords['Time'] = \
                ('Time', datetime_to_days(times + (cycleIndex + 1) * period,
                                          calendar=calendar))
            dsShift = xr.concat([dsShift, dsNew], dim='Time')

        # clip dsShift to the range of ds
        dsStartTime = dsShift.Time.sel(Time=ds.Time.min(),
                                       method=str('nearest')).values
        dsEndTime = dsShift.Time.sel(Time=ds.Time.max(),
                                     method=str('nearest')).values
        dsShift = dsShift.sel(Time=slice(dsStartTime, dsEndTime))

        return dsShift

    def _compute_area_vol(self):
        '''
        Compute part of the time series of sea ice volume and area, given time
        indices to process.
        '''

        outFileNames = {}
        for hemisphere in ['NH', 'SH']:
            baseDirectory = build_config_full_path(
                self.config, 'output', 'timeSeriesSubdirectory')

            make_directories(baseDirectory)

            outFileName = '{}/seaIceAreaVol{}.nc'.format(baseDirectory,
                                                         hemisphere)
            outFileNames[hemisphere] = outFileName

        dsTimeSeries = {}
        dsMesh = xr.open_dataset(self.restartFileName)
        dsMesh = dsMesh[['latCell', 'areaCell']]
        # Load data
        ds = open_mpas_dataset(
            fileName=self.inputFile,
            calendar=self.calendar,
            variableList=self.variableList,
            startDate=self.startDate,
            endDate=self.endDate)

        for hemisphere in ['NH', 'SH']:

            if hemisphere == 'NH':
                mask = dsMesh.latCell > 0
            else:
                mask = dsMesh.latCell < 0

            dsAreaSum = (ds.where(mask) * dsMesh.areaCell).sum('nCells')
            dsAreaSum = dsAreaSum.rename(
                {'timeMonthly_avg_iceAreaCell': 'iceArea',
                 'timeMonthly_avg_iceVolumeCell': 'iceVolume'})
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

            dsTimeSeries[hemisphere] = dsAreaSum

            write_netcdf(dsAreaSum, outFileNames[hemisphere])

        return dsTimeSeries
