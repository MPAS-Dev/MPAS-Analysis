import xarray as xr
import os

from .sea_ice_analysis_task import SeaIceAnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot, \
    timeseries_analysis_plot_polar

from ..shared.io.utility import build_config_full_path, check_path_exists, \
    make_directories

from ..shared.timekeeping.utility import date_to_days, days_to_datetime, \
    datetime_to_days
from ..shared.timekeeping.MpasRelativeDelta import MpasRelativeDelta

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset
from ..shared.mpas_xarray.mpas_xarray import subset_variables

from ..shared.time_series import cache_time_series


class TimeSeriesSeaIce(SeaIceAnalysisTask):
    """
    Performs analysis of time series of sea-ice properties.

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
        # first, call the constructor from the base class (SeaIceAnalysisTask)
        super(TimeSeriesSeaIce, self).__init__(
            config=config,
            taskName='timeSeriesSeaIceAreaVol',
            componentName='seaIce',
            tags=['timeSeries'])

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
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(TimeSeriesSeaIce, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config
        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
                check_path_exists(config.get('seaIcePreprocessedReference',
                                             'baseDirectory'))

        # get a list of timeSeriesStatsMonthly output files from the streams
        # file, reading only those that are between the start and end dates
        streamName = self.historyStreams.find_stream(
            self.streamMap['timeSeriesStats'])
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
        return  # }}}

    def run(self):  # {{{
        """
        Performs analysis of time series of sea-ice properties.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting sea-ice area and volume time series..."

        config = self.config
        calendar = self.calendar

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))

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

        print '  Load sea-ice data...'
        # Load mesh
        self.dsMesh = xr.open_dataset(self.restartFileName)
        self.dsMesh = subset_variables(self.dsMesh,
                                       variableList=['lonCell', 'latCell',
                                                     'areaCell'])

        # Load data
        ds = open_multifile_dataset(
            fileNames=self.inputFiles,
            calendar=calendar,
            config=config,
            simulationStartTime=self.simulationStartTime,
            timeVariableName='Time',
            variableList=['iceAreaCell', 'iceVolumeCell'],
            variableMap=self.variableMap,
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

            # store some variables for use in _compute_area_vol_part
            self.hemisphere = hemisphere
            self.ds = ds
            dsTimeSeries[hemisphere] = cache_time_series(
                ds.Time.values, self._compute_area_vol_part, cacheFileName,
                calendar, yearsPerCacheUpdate=10, printProgress=True)

            print '  Make {} plots...'.format(hemisphere)

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

                title[key] = '{} ({}), {} (r)'.format(
                    plotTitles[variableName], hemisphere, mainRunName)

            if compareWithObservations:
                key = (hemisphere, 'iceArea')
                title[key] = '{}\nSSM/I observations, annual cycle (k)'.format(
                    title[key])
                if hemisphere == 'NH':
                    key = (hemisphere, 'iceVolume')
                    title[key] = \
                        '{}\nPIOMAS, annual cycle (k)'.format(title[key])

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
                        timeseries_analysis_plot_polar(
                            config,
                            plotVars[key],
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
                figureNameStd = '{}/{}.{}.png'.format(self.plotsDirectory,
                                                      mainRunName,
                                                      variableName)
                figureNamePolar = \
                    '{}/{}.{}_polar.png'.format(self.plotsDirectory,
                                                mainRunName,
                                                variableName)
                title = \
                    '{}, NH (r), SH (k)\n{}'.format(plotTitles[variableName],
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
                                       method='nearest').values
        dsEndTime = dsShift.Time.sel(Time=ds.Time.max(),
                                     method='nearest').values
        dsShift = dsShift.sel(Time=slice(dsStartTime, dsEndTime))

        return dsShift  # }}}

    def _compute_area_vol_part(self, timeIndices, firstCall):  # {{{
        '''
        Compute part of the time series of sea ice volume and area, given time
        indices to process.
        '''
        dsLocal = self.ds.isel(Time=timeIndices)

        if self.hemisphere == 'NH':
            mask = self.dsMesh.latCell > 0
        else:
            mask = self.dsMesh.latCell < 0
        dsLocal = dsLocal.where(mask)

        dsAreaSum = (dsLocal*self.dsMesh.areaCell).sum('nCells')
        dsAreaSum = dsAreaSum.rename({'iceAreaCell': 'iceArea',
                                      'iceVolumeCell': 'iceVolume'})
        dsAreaSum['iceThickness'] = (dsAreaSum.iceVolume /
                                     self.dsMesh.areaCell.sum('nCells'))

        dsAreaSum['iceArea'].attrs['units'] = 'm$^2$'
        dsAreaSum['iceArea'].attrs['description'] = \
            'Total {} sea ice area'.format(self.hemisphere)
        dsAreaSum['iceVolume'].attrs['units'] = 'm$^3$'
        dsAreaSum['iceVolume'].attrs['description'] = \
            'Total {} sea ice volume'.format(self.hemisphere)
        dsAreaSum['iceThickness'].attrs['units'] = 'm'
        dsAreaSum['iceThickness'].attrs['description'] = \
            'Mean {} sea ice volume'.format(self.hemisphere)

        return dsAreaSum  # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
