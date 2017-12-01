
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr

from .sea_ice_analysis_task import SeaIceAnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot, \
    timeseries_analysis_plot_polar

from ..shared.io.utility import build_config_full_path, check_path_exists, \
    make_directories

from ..shared.timekeeping.utility import date_to_days, days_to_datetime, \
    datetime_to_days
from ..shared.timekeeping.MpasRelativeDelta import MpasRelativeDelta

from ..shared.generalized_reader import open_multifile_dataset
from ..shared.io import open_mpas_dataset
from ..shared.mpas_xarray.mpas_xarray import subset_variables

from ..shared.html import write_image_xml


class TimeSeriesSeaIce(SeaIceAnalysisTask):
    """
    Performs analysis of time series of sea-ice properties.

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
        # first, call the constructor from the base class (SeaIceAnalysisTask)
        super(TimeSeriesSeaIce, self).__init__(
            config=config,
            taskName='timeSeriesSeaIceAreaVol',
            componentName='seaIce',
            tags=['timeSeries'])

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
        # first, call setup_and_check from the base class (SeaIceAnalysisTask),
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
        return  # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of time series of sea-ice properties.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        self.logger.info("\nPlotting sea-ice area and volume time series...")

        config = self.config
        calendar = self.calendar

        plotTitles = {'iceArea': 'Sea-ice area',
                      'iceVolume': 'Sea-ice volume',
                      'iceThickness': 'Sea-ice mean thickness'}

        units = {'iceArea': '[km$^2$]',
                 'iceVolume': '[10$^3$ km$^3$]',
                 'iceThickness': '[m]'}

        obsFileNames = {
            'iceArea': {'NH': build_config_full_path(config,
                                                     'seaIceObservations',
                                                     'areaNH'),
                        'SH': build_config_full_path(config,
                                                     'seaIceObservations',
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
        self.dsMesh = xr.open_dataset(self.restartFileName)
        self.dsMesh = subset_variables(self.dsMesh,
                                       variableList=['lonCell', 'latCell',
                                                     'areaCell'])

        # Load data
        ds = open_mpas_dataset(
            fileName=self.inputFile,
            calendar=calendar,
            variableList=self.variableList,
            startDate=self.startDate,
            endDate=self.endDate)

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
                dsPreprocessedTimeSlice = \
                    dsPreprocessed.sel(Time=slice(timeStart, timeEnd))
            else:
                self.logger.warning('Preprocessed time series ends before the '
                                    'timeSeries startYear and will not be '
                                    'plotted.')
                preprocessedReferenceRunName = 'None'

        norm = {'iceArea': 1e-6,  # m^2 to km^2
                'iceVolume': 1e-12,  # m^3 to 10^3 km^3
                'iceThickness': 1.}

        xLabel = 'Time [years]'

        galleryGroup = 'Time Series'
        groupLink = 'timeseries'

        dsTimeSeries = {}
        obs = {}
        preprocessed = {}
        figureNameStd = {}
        figureNamePolar = {}
        title = {}
        plotVars = {}

        for hemisphere in ['NH', 'SH']:

            dsTimeSeries[hemisphere] = self._compute_area_vol(ds, hemisphere)

            self.logger.info('  Make {} plots...'.format(hemisphere))

            for variableName in ['iceArea', 'iceVolume']:
                key = (hemisphere, variableName)

                # apply the norm to each variable
                plotVars[key] = (norm[variableName] *
                                 dsTimeSeries[hemisphere][variableName])

                prefix = '{}/{}{}_{}'.format(self.plotsDirectory,
                                             variableName,
                                             hemisphere,
                                             mainRunName)

                figureNameStd[key] = '{}.png'.format(prefix)
                figureNamePolar[key] = '{}_polar.png'.format(prefix)

                title[key] = '{} ({}) \n {} (black)'.format(
                    plotTitles[variableName], hemisphere, mainRunName)

            if compareWithObservations:
                key = (hemisphere, 'iceArea')
                title[key] = '{}\nSSM/I observations, annual cycle (blue)'.format(
                    title[key])
                if hemisphere == 'NH':
                    key = (hemisphere, 'iceVolume')
                    title[key] = \
                        '{}\nPIOMAS, annual cycle (blue)'.format(title[key])

            if preprocessedReferenceRunName != 'None':
                for variableName in ['iceArea', 'iceVolume']:
                    key = (hemisphere, variableName)
                    title[key] = '{}\n {} (red)'.format(
                        title[key], preprocessedReferenceRunName)

            if compareWithObservations:
                dsObs = open_multifile_dataset(
                    fileNames=obsFileNames['iceArea'][hemisphere],
                    calendar=calendar,
                    config=config,
                    timeVariableName='xtime')
                key = (hemisphere, 'iceArea')
                obs[key] = self._replicate_cycle(plotVars[key], dsObs.IceArea,
                                                 calendar)

                key = (hemisphere, 'iceVolume')
                if hemisphere == 'NH':
                    dsObs = open_multifile_dataset(
                        fileNames=obsFileNames['iceVolume'][hemisphere],
                        calendar=calendar,
                        config=config,
                        timeVariableName='xtime')
                    obs[key] = self._replicate_cycle(plotVars[key],
                                                     dsObs.IceVol,
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
                        lineStyles = ['k-', 'b-', 'r-']
                        lineWidths = [3, 1.2, 1.2]
                        legendText = [None, None, None]
                    else:
                        # just v1 model and obs
                        plotVars[key] = [plotVars[key], obs[key]]
                        lineStyles = ['k-', 'b-']
                        lineWidths = [3, 1.2]
                        legendText = [None, None]
                elif preprocessedReferenceRunName != 'None':
                    # just v1 and v0 models
                    plotVars[key] = [plotVars[key], preprocessed[key]]
                    lineStyles = ['k-', 'r-']
                    lineWidths = [3, 1.2]
                    legendText = [None, None]

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
                                             legendText=legendText,
                                             titleFontSize=titleFontSize,
                                             calendar=calendar)
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
                        timeseries_analysis_plot_polar(
                            config,
                            plotVars[key],
                            movingAveragePoints,
                            title[key],
                            figureNamePolar[key],
                            lineStyles=lineStyles,
                            lineWidths=lineWidths,
                            legendText=legendText,
                            titleFontSize=titleFontSize,
                            calendar=calendar)

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

        if (not compareWithObservations and
                preprocessedReferenceRunName == 'None'):
            for variableName in ['iceArea', 'iceVolume']:
                # we will combine north and south onto a single graph
                figureNameStd = '{}/{}.{}.png'.format(self.plotsDirectory,
                                                      mainRunName,
                                                      variableName)
                figureNamePolar = \
                    '{}/{}.{}_polar.png'.format(self.plotsDirectory,
                                                mainRunName,
                                                variableName)
                title = \
                    '{}, NH (black), SH (blue)\n{}'.format(plotTitles[variableName],
                                                    mainRunName)
                varList = [plotVars[('NH', variableName)],
                           plotVars[('SH', variableName)]]
                timeseries_analysis_plot(config, varList,
                                         movingAveragePoints,
                                         title, xLabel, units[variableName],
                                         figureNameStd,
                                         lineStyles=['k-', 'b-'],
                                         lineWidths=[2, 2],
                                         legendText=[None, None],
                                         titleFontSize=titleFontSize,
                                         calendar=calendar)
                filePrefix = '{}.{}'.format(mainRunName, variableName)
                thumbnailDescription = plotTitles[variableName]
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
                    timeseries_analysis_plot_polar(config, varList,
                                                   movingAveragePoints,
                                                   title, figureNamePolar,
                                                   lineStyles=['k-', 'b-'],
                                                   lineWidths=[2, 2],
                                                   legendText=[None, None],
                                                   titleFontSize=titleFontSize,
                                                   calendar=calendar)
                    filePrefix = '{}.{}_polar'.format(mainRunName,
                                                      variableName)
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

        # }}}

    def _replicate_cycle(self, ds, dsToReplicate, calendar):  # {{{
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

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
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
            dsNew.coords['Time'] = \
                ('Time', datetime_to_days(times + (cycleIndex+1)*period,
                                          calendar=calendar))
            dsShift = xr.concat([dsShift, dsNew], dim='Time')

        # clip dsShift to the range of ds
        dsStartTime = dsShift.Time.sel(Time=ds.Time.min(),
                                       method=str('nearest')).values
        dsEndTime = dsShift.Time.sel(Time=ds.Time.max(),
                                     method=str('nearest')).values
        dsShift = dsShift.sel(Time=slice(dsStartTime, dsEndTime))

        return dsShift  # }}}

    def _compute_area_vol(self, ds, hemisphere):  # {{{
        '''
        Compute part of the time series of sea ice volume and area, given time
        indices to process.
        '''

        if hemisphere == 'NH':
            mask = self.dsMesh.latCell > 0
        else:
            mask = self.dsMesh.latCell < 0
        ds = ds.where(mask)

        dsAreaSum = (ds*self.dsMesh.areaCell).sum('nCells')
        dsAreaSum = dsAreaSum.rename(
                {'timeMonthly_avg_iceAreaCell': 'iceArea',
                 'timeMonthly_avg_iceVolumeCell': 'iceVolume'})
        dsAreaSum['iceThickness'] = (dsAreaSum.iceVolume /
                                     self.dsMesh.areaCell.sum('nCells'))

        dsAreaSum['iceArea'].attrs['units'] = 'm$^2$'
        dsAreaSum['iceArea'].attrs['description'] = \
            'Total {} sea ice area'.format(hemisphere)
        dsAreaSum['iceVolume'].attrs['units'] = 'm$^3$'
        dsAreaSum['iceVolume'].attrs['description'] = \
            'Total {} sea ice volume'.format(hemisphere)
        dsAreaSum['iceThickness'].attrs['units'] = 'm'
        dsAreaSum['iceThickness'].attrs['description'] = \
            'Mean {} sea ice volume'.format(hemisphere)

        return dsAreaSum  # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
