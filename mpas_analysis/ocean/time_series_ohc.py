import numpy as np
import netCDF4
import os

from ..shared.analysis_task import AnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime, string_to_datetime

from ..shared.time_series import time_series

from ..shared.io.utility import build_config_full_path, make_directories, \
    check_path_exists


class TimeSeriesOHC(AnalysisTask):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.

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
        super(TimeSeriesOHC, self).__init__(config)

        # name the task, component and category
        self.taskName = 'timeSeriesOHC'
        self.componentName = 'ocean'
        self.tags = ['timeSeries', 'ohc']

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
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(TimeSeriesOHC, self).setup_and_check()

        config = self.config

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
                check_path_exists(config.get('oceanPreprocessedReference',
                                             'baseDirectory'))
        return  # }}}

    def run(self):  # {{{
        """
        Performs analysis of ocean heat content (OHC) from time-series output.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting OHC time series..."

        simulationStartTime = get_simulation_start_time(self.runStreams)
        config = self.config
        calendar = self.calendar

        # read parameters from config file
        mainRunName = config.get('runs', 'mainRunName')
        preprocessedReferenceRunName = \
            config.get('runs', 'preprocessedReferenceRunName')
        preprocessedInputDirectory = config.get('oceanPreprocessedReference',
                                                'baseDirectory')

        compareWithObservations = config.getboolean('timeSeriesOHC',
                                                    'compareWithObservations')

        movingAveragePoints = config.getint('timeSeriesOHC',
                                            'movingAveragePoints')

        regions = config.getExpression('regions', 'regions')
        plotTitles = config.getExpression('regions', 'plotTitles')
        regionIndicesToPlot = config.getExpression('timeSeriesOHC',
                                                   'regionIndicesToPlot')

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        regionNames = config.getExpression('regions', 'regions')
        regionNames = [regionNames[index] for index in regionIndicesToPlot]

        # Note: input file, not a mesh file because we need dycore specific
        # fields such as refBottomDepth and namelist fields such as
        # config_density0, as well as simulationStartTime, that are not
        # guaranteed to be in the mesh file.
        try:
            restartFile = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for OHC calculation')

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = config.get('timeSeries', 'startDate')
        endDate = config.get('timeSeries', 'endDate')
        streamName = self.historyStreams.find_stream(
            self.streamMap['timeSeriesStats'])
        fileNames = self.historyStreams.readpath(streamName,
                                                 startDate=startDate,
                                                 endDate=endDate,
                                                 calendar=calendar)
        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(fileNames[0]),
                  os.path.basename(fileNames[-1]))

        # Define/read in general variables
        print '  Read in depth and compute specific depth indexes...'
        ncFile = netCDF4.Dataset(restartFile, mode='r')
        # reference depth [m]
        depth = ncFile.variables['refBottomDepth'][:]
        ncFile.close()

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
                                    variableMap=self.variableMap,
                                    startDate=startDate,
                                    endDate=endDate)

        ds = ds.isel(nOceanRegionsTmp=regionIndicesToPlot)

        timeStart = string_to_datetime(startDate)
        timeEnd = string_to_datetime(endDate)

        # Select year-1 data and average it (for later computing anomalies)
        timeStartFirstYear = string_to_datetime(simulationStartTime)
        if timeStartFirstYear < timeStart:
            startDateFirstYear = simulationStartTime
            firstYear = int(startDateFirstYear[0:4])
            endDateFirstYear = '{:04d}-12-31_23:59:59'.format(firstYear)
            filesFirstYear = \
                self.historyStreams.readpath(streamName,
                                             startDate=startDateFirstYear,
                                             endDate=endDateFirstYear,
                                             calendar=calendar)
            dsFirstYear = open_multifile_dataset(
                fileNames=filesFirstYear,
                calendar=calendar,
                config=config,
                simulationStartTime=simulationStartTime,
                timeVariableName='Time',
                variableList=['avgLayerTemperature'],
                variableMap=self.variableMap,
                startDate=startDateFirstYear,
                endDate=endDateFirstYear)

            dsFirstYear = \
                dsFirstYear.isel(nOceanRegionsTmp=regionIndicesToPlot)

            firstYearAvgLayerTemperature = dsFirstYear.avgLayerTemperature
        else:
            firstYearAvgLayerTemperature = ds.avgLayerTemperature
            firstYear = timeStart.year

        timeStartFirstYear = date_to_days(year=firstYear, month=1, day=1,
                                          calendar=calendar)
        timeEndFirstYear = date_to_days(year=firstYear, month=12, day=31,
                                        hour=23, minute=59, second=59,
                                        calendar=calendar)

        firstYearAvgLayerTemperature = firstYearAvgLayerTemperature.sel(
            Time=slice(timeStartFirstYear, timeEndFirstYear))

        firstYearAvgLayerTemperature = \
            firstYearAvgLayerTemperature.mean('Time')

        print '  Compute temperature anomalies...'

        ds['avgLayTemperatureAnomaly'] = (ds.avgLayerTemperature -
                                          firstYearAvgLayerTemperature)

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
                dsPreprocessedTimeSlice = \
                    dsPreprocessed.sel(Time=slice(timeStart, timeEnd))
            else:
                print '   Warning: Preprocessed time series ends before the ' \
                    'timeSeries startYear and will not be plotted.'
                preprocessedReferenceRunName = 'None'

        cacheFileName = '{}/ohcTimeSeries.nc'.format(outputDirectory)

        # store fields needed by _compute_ohc_part
        self.ds = ds
        self.regionNames = regionNames
        dsOHC = time_series.cache_time_series(ds.Time.values,
                                              self._compute_ohc_part,
                                              cacheFileName, calendar,
                                              yearsPerCacheUpdate=10,
                                              printProgress=True)

        unitsScalefactor = 1e-22

        print '  Compute OHC and make plots...'
        for index, regionIndex in enumerate(regionIndicesToPlot):

            ohc = dsOHC.ohc.isel(nOceanRegionsTmp=index)

            # OHC over 0-bottom depth range:
            ohcTotal = ohc.sum('nVertLevels')
            ohcTotal = unitsScalefactor*ohcTotal

            # OHC over 0-700m depth range:
            ohc700m = unitsScalefactor*ohc[:, 0:k700m].sum('nVertLevels')

            # OHC over 700m-2000m depth range:
            ohc2000m = \
                unitsScalefactor*ohc[:, k700m+1:k2000m].sum('nVertLevels')

            # OHC over 2000m-bottom depth range:
            ohcBottom = ohc[:, k2000m+1:kbtm].sum('nVertLevels')
            ohcBottom = unitsScalefactor*ohcBottom

            title = 'OHC, {}, 0-bottom (thick-),' \
                    ' 0-700m (thin-), 700-2000m (--),' \
                    ' 2000m-bottom (-.) \n {}'.format(plotTitles[regionIndex],
                                                      mainRunName)

            xLabel = 'Time [years]'
            yLabel = '[x$10^{22}$ J]'

            figureName = '{}/ohc_{}_{}.png'.format(self.plotsDirectory,
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
                                                  ohcBottom,
                                                  ohcPreprocessedTotal,
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
        # }}}

    def _compute_ohc_part(self, timeIndices, firstCall):  # {{{
        '''
        Compute part of the OHC time series, given time indices to process.
        '''

        # specific heat [J/(kg*degC)]
        cp = self.namelist.getfloat('config_specific_heat_sea_water')
        # [kg/m3]
        rho = self.namelist.getfloat('config_density0')

        dsLocal = self.ds.isel(Time=timeIndices)

        dsLocal['ohc'] = rho*cp*dsLocal.sumLayerMaskValue * \
            dsLocal.avgLayerArea * dsLocal.avgLayerThickness * \
            dsLocal.avgLayTemperatureAnomaly
        dsLocal.ohc.attrs['units'] = 'J'
        dsLocal.ohc.attrs['description'] = 'Ocean heat content in each region'
        dsLocal['regionNames'] = ('nOceanRegionsTmp', self.regionNames)

        return dsLocal  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
