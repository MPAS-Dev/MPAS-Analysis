# -*- coding: utf-8 -*-
import numpy as np
import netCDF4
import os

from ..shared.analysis_task import AnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot, \
    plot_vertical_section, setup_colormap

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime, string_to_datetime

from ..shared.time_series import time_series

from ..shared.io.utility import build_config_full_path, make_directories, \
    check_path_exists
from ..shared.html import write_image_xml


class TimeSeriesOHC(AnalysisTask):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani, Greg Streletz
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
        super(TimeSeriesOHC, self).__init__(
            config=config,
            taskName='timeSeriesOHC',
            componentName='ocean',
            tags=['timeSeries', 'ohc'])

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
        Xylar Asay-Davis, Greg Streletz
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(TimeSeriesOHC, self).setup_and_check()

        config = self.config

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
                check_path_exists(config.get('oceanPreprocessedReference',
                                             'baseDirectory'))

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        self.streamName = 'timeSeriesStatsMonthlyOutput'
        self.startDate = self.config.get('timeSeries', 'startDate')
        self.endDate = self.config.get('timeSeries', 'endDate')
        self.inputFiles = self.historyStreams.readpath(
                self.streamName, startDate=self.startDate,
                endDate=self.endDate, calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(self.streamName, self.startDate,
                                       self.endDate))

        mainRunName = config.get('runs', 'mainRunName')
        regions = config.getExpression('regions', 'regions')
        regionIndicesToPlot = config.getExpression('timeSeriesOHC',
                                                   'regionIndicesToPlot')

        self.xmlFileNames = []
        self.filePrefixes = {}

        regions = [regions[index] for index in regionIndicesToPlot]

        for region in regions:
            filePrefix = 'TAnomalyZ_{}_{}'.format(region, mainRunName)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                            filePrefix))
            self.filePrefixes[0, region] = filePrefix

            filePrefix = 'SAnomalyZ_{}_{}'.format(region, mainRunName)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                            filePrefix))
            self.filePrefixes[1, region] = filePrefix

            filePrefix = 'OHCAnomalyZ_{}_{}'.format(region, mainRunName)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                            filePrefix))
            self.filePrefixes[2, region] = filePrefix

            filePrefix = 'OHCAnomaly_{}_{}'.format(region, mainRunName)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[3, region] = filePrefix
        return  # }}}

    def run(self):  # {{{
        """
        Performs analysis of ocean heat content (OHC) from time-series output.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani, Greg Streletz
        """

        print "\nPlotting OHC time series and T, S, and OHC vertical trends..."

        simulationStartTime = get_simulation_start_time(self.runStreams)
        config = self.config
        calendar = self.calendar

        # read parameters from config file
        mainRunName = config.get('runs', 'mainRunName')
        preprocessedReferenceRunName = \
            config.get('runs', 'preprocessedReferenceRunName')
        preprocessedInputDirectory = config.get('oceanPreprocessedReference',
                                                'baseDirectory')

        configSectionName = 'timeSeriesOHC'

        plotOriginalFields = config.getboolean(configSectionName,
                                               'plotOriginalFields')

        compareWithObservations = config.getboolean(configSectionName,
                                                    'compareWithObservations')

        movingAveragePoints = config.getint(configSectionName,
                                            'movingAveragePoints')

        movingAveragePointsHovmoller = config.getint(configSectionName,
                                                     'movingAveragePointsHovmoller')

        regions = config.getExpression('regions', 'regions')
        plotTitles = config.getExpression('regions', 'plotTitles')
        regionIndicesToPlot = config.getExpression(configSectionName,
                                                   'regionIndicesToPlot')

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        regionNames = config.getExpression('regions', 'regions')

        # Note: input file, not a mesh file because we need dycore specific
        # fields such as refBottomDepth and namelist fields such as
        # config_density0, as well as simulationStartTime, that are not
        # guaranteed to be in the mesh file.
        try:
            restartFile = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for OHC calculation')

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))

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
        avgTemperatureVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature'
        avgSalinityVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerSalinity'
        sumMaskVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue'
        avgAreaVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea'
        avgThickVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness'
        variableList = [avgTemperatureVarName, avgSalinityVarName, sumMaskVarName,
                        avgAreaVarName, avgThickVarName]
        ds = open_multifile_dataset(fileNames=self.inputFiles,
                                    calendar=calendar,
                                    config=config,
                                    simulationStartTime=simulationStartTime,
                                    timeVariableName=['xtime_startMonthly',
                                                      'xtime_endMonthly'],
                                    variableList=variableList,
                                    startDate=self.startDate,
                                    endDate=self.endDate)

        timeStart = string_to_datetime(self.startDate)
        timeEnd = string_to_datetime(self.endDate)

        # Select year-1 data and average it (for later computing anomalies)
        timeStartFirstYear = string_to_datetime(simulationStartTime)
        if timeStartFirstYear < timeStart:
            startDateFirstYear = simulationStartTime
            firstYear = int(startDateFirstYear[0:4])
            endDateFirstYear = '{:04d}-12-31_23:59:59'.format(firstYear)
            filesFirstYear = \
                self.historyStreams.readpath(self.streamName,
                                             startDate=startDateFirstYear,
                                             endDate=endDateFirstYear,
                                             calendar=calendar)
            dsFirstYear = open_multifile_dataset(
                fileNames=filesFirstYear,
                calendar=calendar,
                config=config,
                simulationStartTime=simulationStartTime,
                timeVariableName=['xtime_startMonthly', 'xtime_endMonthly'],
                variableList=[avgTemperatureVarName, avgSalinityVarName],
                startDate=startDateFirstYear,
                endDate=endDateFirstYear)

            firstYearAvgLayerTemperature = dsFirstYear[avgTemperatureVarName]
            firstYearAvgLayerSalinity = dsFirstYear[avgSalinityVarName]
        else:
            firstYearAvgLayerTemperature = ds[avgTemperatureVarName]
            firstYearAvgLayerSalinity = ds[avgSalinityVarName]
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

        firstYearAvgLayerSalinity = firstYearAvgLayerSalinity.sel(
            Time=slice(timeStartFirstYear, timeEndFirstYear))
        firstYearAvgLayerSalinity = \
            firstYearAvgLayerSalinity.mean('Time')

        print '  Compute temperature and salinity anomalies...'

        ds['avgLayerTemperatureAnomaly'] = (ds[avgTemperatureVarName] - firstYearAvgLayerTemperature)

        ds['avgLayerSalinityAnomaly'] = (ds[avgSalinityVarName] - firstYearAvgLayerSalinity)

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
        for regionIndex in regionIndicesToPlot:
            region = regions[regionIndex]

            # Plot temperature, salinity and OHC anomalies (and full fields)
            # trends with depth

            #  First T vs depth/time
            x = ds.Time.values
            y = depth
            z = ds.avgLayerTemperatureAnomaly.isel(nOceanRegionsTmp=regionIndex)
            z = z.transpose()

            colorbarLabel = '[$^\circ$ C]'
            xLabel = 'Time [years]'
            yLabel = 'Depth [m]'

            title = 'Temperature Anomaly, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

            figureName = '{}/TAnomalyZ_{}_{}.png'.format(self.plotsDirectory,
                                                         region,
                                                         mainRunName)

            (colormapName, colorbarLevels) = setup_colormap(config,
                                                            configSectionName,
                                                            suffix='TemperatureAnomaly')

            contourLevels = \
                config.getExpression(configSectionName,
                                     'contourLevels{}'.format('TemperatureAnomaly'),
                                     usenumpyfunc=True)
   
            plot_vertical_section(config, x, y, z,
                                  colormapName, colorbarLevels, contourLevels,
                                  colorbarLabel, title, xLabel, yLabel,
                                  figureName, linewidths=1, xArrayIsTime=True,
                                  N=movingAveragePointsHovmoller, calendar=calendar)

            filePrefix = self.filePrefixes[0, region]
            caption = 'Trend of {} Temperature Anomaly vs depth from Year 0001'.format(region)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Trends vs Depth',
                groupLink='trendsvsdepth',
                thumbnailDescription=u'{} ΔT'.format(region),
                imageDescription=caption,
                imageCaption=caption)
             
            if plotOriginalFields:
                z = ds[avgTemperatureVarName].isel(nOceanRegionsTmp=regionIndex)
                z = z.transpose()

                title = 'Temperature, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

                figureName = '{}/TZ_{}_{}.png'.format(self.plotsDirectory,
                                                      region,
                                                      mainRunName)

                (colormap, colorbarLevels) = setup_colormap(config,
                                                            configSectionName,
                                                            suffix='Temperature')

                contourLevels = config.getExpression(configSectionName,
                                                     'contourLevels{}'.format('Temperature'),
                                                     usenumpyfunc=True)

                plot_vertical_section(config, x, y, z,
                                      colormap, colorbarLevels, contourLevels,
                                      colorbarLabel, title, xLabel, yLabel,
                                      figureName, linewidths=1, xArrayIsTime=True,
                                      N=movingAveragePointsHovmoller, calendar=calendar)
            
            #  Second S vs depth/time
            z = ds.avgLayerSalinityAnomaly.isel(nOceanRegionsTmp=regionIndex)
            z = z.transpose()

            title = 'Salinity Anomaly, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

            colorbarLabel = '[PSU]'

            figureName = '{}/SAnomalyZ_{}_{}.png'.format(self.plotsDirectory,
                                                         region,
                                                         mainRunName)

            (colormapName, colorbarLevels) = setup_colormap(config,
                                                            configSectionName,
                                                            suffix='SalinityAnomaly')
   
            contourLevels = \
                config.getExpression(configSectionName,
                                     'contourLevels{}'.format('SalinityAnomaly'),
                                     usenumpyfunc=True)

            plot_vertical_section(config, x, y, z,
                                  colormapName, colorbarLevels, contourLevels,
                                  colorbarLabel, title, xLabel, yLabel,
                                  figureName, linewidths=1, xArrayIsTime=True,
                                  N=movingAveragePointsHovmoller, calendar=calendar)

            filePrefix = self.filePrefixes[1, region]
            caption = 'Trend of {} Salinity Anomaly vs depth from Year 0001'.format(region)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Trends vs Depth',
                groupLink='trendsvsdepth',
                thumbnailDescription=u'{} ΔS'.format(region),
                imageDescription=caption,
                imageCaption=caption)

            if plotOriginalFields:
                z = ds[avgSalinityVarName].isel(nOceanRegionsTmp=regionIndex)
                z = z.transpose()

                title = 'Salinity, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

                figureName = '{}/SZ_{}_{}.png'.format(self.plotsDirectory,
                                                      region,
                                                      mainRunName)

                (colormapName, colorbarLevels) = setup_colormap(config,
                                                                configSectionName,
                                                                suffix='Salinity')
  
                contourLevels = config.getExpression(configSectionName,
                                                     'contourLevels{}'.format('Salinity'),
                                                     usenumpyfunc=True)

                plot_vertical_section(config, x, y, z,
                                      colormapName, colorbarLevels, contourLevels,
                                      colorbarLabel, title, xLabel, yLabel,
                                      figureName, linewidths=1, xArrayIsTime=True,
                                      N=movingAveragePointsHovmoller, calendar=calendar)

            #  Third OHC vs depth/time
            ohcAnomaly = dsOHC.ohcAnomaly.isel(nOceanRegionsTmp=regionIndex)
            ohcAnomalyScaled = unitsScalefactor*ohcAnomaly            
            z = ohcAnomalyScaled.transpose()

            title = 'OHC Anomaly, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

            colorbarLabel = '[x$10^{22}$ J]'

            figureName = '{}/OHCAnomalyZ_{}_{}.png'.format(self.plotsDirectory,
                                                   region,
                                                   mainRunName)
            
            (colormap, colorbarLevels) = setup_colormap(config,
                                                        configSectionName,
                                                        suffix='OHCAnomaly')

            contourLevels = \
                config.getExpression(configSectionName,
                                     'contourLevels{}'.format('OHCAnomaly'),
                                     usenumpyfunc=True)

            plot_vertical_section(config, x, y, z,
                                  colormap, colorbarLevels, contourLevels,
                                  colorbarLabel, title, xLabel, yLabel,
                                  figureName, linewidths=1, xArrayIsTime=True,
                                  N=movingAveragePointsHovmoller, calendar=calendar)

            filePrefix = self.filePrefixes[2, region]
            caption = 'Trend of {} OHC Anomaly vs depth from Year 0001'.format(region)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Trends vs Depth',
                groupLink='trendsvsdepth',
                thumbnailDescription=u'{} ΔOHC'.format(region),
                imageDescription=caption,
                imageCaption=caption)

            if plotOriginalFields:
                ohc = dsOHC.ohc.isel(nOceanRegionsTmp=regionIndex)
                ohcScaled = unitsScalefactor*ohc
                z = ohcScaled.transpose()

                title = 'OHC, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

                figureName = '{}/OHCZ_{}_{}.png'.format(self.plotsDirectory,
                                                        region,
                                                        mainRunName)
            
                (colormap, colorbarLevels) = setup_colormap(config,
                                                            configSectionName,
                                                            suffix='OHC')

                contourLevels = config.getExpression(configSectionName,
                                                     'contourLevels{}'.format('OHC'),
                                                     usenumpyfunc=True)

                plot_vertical_section(config, x, y, z,
                                      colormap, colorbarLevels, contourLevels,
                                      colorbarLabel, title, xLabel, yLabel,
                                      figureName, linewidths=1, xArrayIsTime=True,
                                      N=movingAveragePointsHovmoller, calendar=calendar)

            # Now plot OHC timeseries

            # OHC over 0-bottom depth range:
            ohcAnomalyTotal = unitsScalefactor*ohcAnomaly.sum('nVertLevels')

            # OHC over 0-700m depth range:
            ohcAnomaly700m = unitsScalefactor*ohcAnomaly[:, 0:k700m].sum('nVertLevels')

            # OHC over 700m-2000m depth range:
            ohcAnomaly2000m = \
                unitsScalefactor*ohcAnomaly[:, k700m+1:k2000m].sum('nVertLevels')

            # OHC over 2000m-bottom depth range:
            ohcAnomalyBottom = unitsScalefactor*ohcAnomaly[:, k2000m+1:kbtm].sum('nVertLevels')

            xLabel = 'Time [years]'
            yLabel = '[x$10^{22}$ J]'

            title = 'OHC Anomaly, {}, 0-bottom (thick-),' \
                    ' 0-700m (thin-), 700-2000m (--),' \
                    ' 2000m-bottom (-.) \n {}'.format(plotTitles[regionIndex],
                                                      mainRunName)

            figureName = '{}/OHCAnomaly_{}_{}.png'.format(self.plotsDirectory,
                                                   region,
                                                   mainRunName)

            if preprocessedReferenceRunName != 'None':
                # these preprocessed data are OHC *anomalies*
                ohcPreprocessedTotal = dsPreprocessedTimeSlice.ohc_tot
                ohcPreprocessed700m = dsPreprocessedTimeSlice.ohc_700m
                ohcPreprocessed2000m = dsPreprocessedTimeSlice.ohc_2000m
                ohcPreprocessedBottom = dsPreprocessedTimeSlice.ohc_btm
                title = '{} (r), {} (b)'.format(title,
                                                preprocessedReferenceRunName)
                timeseries_analysis_plot(config, [ohcAnomalyTotal, ohcAnomaly700m, ohcAnomaly2000m,
                                                  ohcAnomalyBottom,
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
                timeseries_analysis_plot(config, [ohcAnomalyTotal, ohcAnomaly700m, ohcAnomaly2000m,
                                                  ohcAnomalyBottom],
                                         movingAveragePoints, title,
                                         xLabel, yLabel, figureName,
                                         lineStyles=['r-', 'r-', 'r--', 'r-.'],
                                         lineWidths=[2, 1, 1.5, 1.5],
                                         calendar=calendar)

            filePrefix = self.filePrefixes[3, region]
            caption = 'Running Mean of the Anomaly in {} Ocean Heat Content ' \
                      'from Year 0001'.format(region)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Time Series',
                groupLink='timeseries',
                thumbnailDescription=u'{} ΔOHC'.format(region),
                imageDescription=caption,
                imageCaption=caption)

            if plotOriginalFields:
                ohcTotal = unitsScalefactor*ohc.sum('nVertLevels')
                ohc700m = unitsScalefactor*ohc[:, 0:k700m].sum('nVertLevels')
                ohc2000m = \
                    unitsScalefactor*ohc[:, k700m+1:k2000m].sum('nVertLevels')
                ohcBottom = unitsScalefactor*ohc[:, k2000m+1:kbtm].sum('nVertLevels')

                title = 'OHC, {}, 0-bottom (thick-),' \
                        ' 0-700m (thin-), 700-2000m (--),' \
                        ' 2000m-bottom (-.) \n {}'.format(plotTitles[regionIndex],
                                                          mainRunName)
                        
                figureName = '{}/OHC_{}_{}.png'.format(self.plotsDirectory,
                                                       region,
                                                       mainRunName)

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


        avgTemperatureVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature'
        sumMaskVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue'
        avgAreaVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea'
        avgThickVarName = \
            'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness'


        dsLocal['ohc'] = rho*cp*dsLocal[sumMaskVarName] * \
            dsLocal[avgAreaVarName] * dsLocal[avgThickVarName] * \
            dsLocal[avgTemperatureVarName]
        dsLocal.ohc.attrs['units'] = 'J'
        dsLocal.ohc.attrs['description'] = 'Ocean heat content in each region'

        dsLocal['ohcAnomaly'] = rho*cp*dsLocal[sumMaskVarName] * \
            dsLocal[avgAreaVarName] * dsLocal[avgThickVarName] * \
            dsLocal.avgLayerTemperatureAnomaly
        dsLocal.ohcAnomaly.attrs['units'] = 'J'
        dsLocal.ohcAnomaly.attrs['description'] = 'Ocean heat content anomaly in each region'

        dsLocal['regionNames'] = ('nOceanRegionsTmp', self.regionNames)

        return dsLocal  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
