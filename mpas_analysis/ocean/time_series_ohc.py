# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import netCDF4

from ..shared import AnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot, \
    plot_vertical_section, setup_colormap

from ..shared.generalized_reader import open_multifile_dataset
from ..shared.io import open_mpas_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime, string_to_datetime

from ..shared.io.utility import build_config_full_path, make_directories, \
    check_path_exists
from ..shared.html import write_image_xml


class TimeSeriesOHC(AnalysisTask):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani, Greg Streletz
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
        super(TimeSeriesOHC, self).__init__(
            config=config,
            taskName='timeSeriesOHC',
            componentName='ocean',
            tags=['timeSeries', 'ohc'])

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
        Xylar Asay-Davis, Greg Streletz
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(TimeSeriesOHC, self).setup_and_check()

        config = self.config

        self.startDate = self.config.get('timeSeries', 'startDate')
        self.endDate = self.config.get('timeSeries', 'endDate')

        self.variables = {}
        for suffix in ['avgLayerTemperature', 'avgLayerSalinity',
                       'sumLayerMaskValue', 'avgLayerArea',
                       'avgLayerThickness']:
            self.variables[suffix] = \
                    'timeMonthly_avg_avgValueWithinOceanLayerRegion_' + suffix

        self.mpasTimeSeriesTask.add_variables(
                variableList=self.variables.values())

        self.inputFile = self.mpasTimeSeriesTask.outputFile

        if config.get('runs', 'preprocessedReferenceRunName') != 'None':
                check_path_exists(config.get('oceanPreprocessedReference',
                                             'baseDirectory'))

        self.startDate = self.config.get('timeSeries', 'startDate')
        self.endDate = self.config.get('timeSeries', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')
        regions = config.getExpression('regions', 'regions')
        regionIndicesToPlot = config.getExpression('timeSeriesOHC',
                                                   'regionIndicesToPlot')

        self.xmlFileNames = []
        self.filePrefixes = {}

        regions = [regions[index] for index in regionIndicesToPlot]

        configSectionName = 'timeSeriesOHC'
        plotOriginalFields = config.getboolean(configSectionName,
                                               'plotOriginalFields')

        if plotOriginalFields:
            plotTypes = ['TZ', 'TAnomalyZ', 'SZ', 'SAnomalyZ', 'OHCZ',
                         'OHCAnomalyZ', 'OHC', 'OHCAnomaly']
        else:
            plotTypes = ['TAnomalyZ', 'SAnomalyZ', 'OHCAnomalyZ', 'OHCAnomaly']

        for region in regions:
            for plotType in plotTypes:
                filePrefix = '{}_{}_{}'.format(plotType, region, mainRunName)
                self.xmlFileNames.append('{}/{}.xml'.format(
                        self.plotsDirectory, filePrefix))
                self.filePrefixes[plotType, region] = filePrefix

        return  # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of ocean heat content (OHC) from time-series output.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani, Greg Streletz
        """

        self.logger.info("\nPlotting OHC time series and T, S, and OHC "
                         "vertical trends...")

        simulationStartTime = get_simulation_start_time(self.runStreams)
        config = self.config
        calendar = self.calendar
        variables = self.variables

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

        movingAveragePointsTimeSeries = config.getint(
                configSectionName, 'movingAveragePointsTimeSeries')

        movingAveragePointsHovmoller = config.getint(
                configSectionName, 'movingAveragePointsHovmoller')

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

        # Define/read in general variables
        self.logger.info('  Read in depth and compute specific depth '
                         'indexes...')
        ncFile = netCDF4.Dataset(restartFile, mode='r')
        # reference depth [m]
        depth = ncFile.variables['refBottomDepth'][:]
        ncFile.close()

        k700m = np.where(depth > 700.)[0][0] - 1
        k2000m = np.where(depth > 2000.)[0][0] - 1

        kbtm = len(depth)-1

        # Load data
        self.logger.info('  Load ocean data...')
        ds = open_mpas_dataset(fileName=self.inputFile,
                               calendar=calendar,
                               variableList=self.variables.values(),
                               startDate=self.startDate,
                               endDate=self.endDate)
        # rename the variables to shorter names for convenience
        renameDict = dict((v, k) for k, v in variables.items())
        ds.rename(renameDict, inplace=True)

        timeStart = string_to_datetime(self.startDate)
        timeEnd = string_to_datetime(self.endDate)

        # Select year-1 data and average it (for later computing anomalies)
        timeStartFirstYear = string_to_datetime(simulationStartTime)
        if timeStartFirstYear < timeStart:
            startDateFirstYear = simulationStartTime
            firstYear = int(startDateFirstYear[0:4])
            endDateFirstYear = '{:04d}-12-31_23:59:59'.format(firstYear)
            dsFirstYear = open_mpas_dataset(
                fileName=self.inputFile,
                calendar=calendar,
                variableList=self.variables.values(),
                startDate=startDateFirstYear,
                endDate=endDateFirstYear)

            dsFirstYear.rename(renameDict, inplace=True)

            firstYearAvgLayerTemperature = dsFirstYear['avgLayerTemperature']
            firstYearavgLayerSalinity = dsFirstYear['avgLayerSalinity']
        else:
            firstYearAvgLayerTemperature = ds['avgLayerTemperature']
            firstYearavgLayerSalinity = ds['avgLayerSalinity']
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

        firstYearavgLayerSalinity = firstYearavgLayerSalinity.sel(
            Time=slice(timeStartFirstYear, timeEndFirstYear))
        firstYearavgLayerSalinity = \
            firstYearavgLayerSalinity.mean('Time')

        self.logger.info('  Compute temperature and salinity anomalies...')

        ds['avgLayerTemperatureAnomaly'] = \
            (ds['avgLayerTemperature'] - firstYearAvgLayerTemperature)

        ds['avgLayerSalinityAnomaly'] = \
            (ds['avgLayerSalinity'] - firstYearavgLayerSalinity)

        yearStart = days_to_datetime(ds.Time.min(), calendar=calendar).year
        yearEnd = days_to_datetime(ds.Time.max(), calendar=calendar).year
        timeStart = date_to_days(year=yearStart, month=1, day=1,
                                 calendar=calendar)
        timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                               calendar=calendar)

        if preprocessedReferenceRunName != 'None':
            self.logger.info('  Load in OHC from preprocessed reference '
                             'run...')
            inFilesPreprocessed = '{}/OHC.{}.year*.nc'.format(
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

        # Add the OHC to the data set
        ds = self._compute_ohc(ds, regionNames)

        unitsScalefactor = 1e-22

        self.logger.info('  Compute OHC and make plots...')
        for regionIndex in regionIndicesToPlot:
            region = regions[regionIndex]

            # Plot temperature, salinity and OHC anomalies (and full fields)
            # trends with depth

            #  First T vs depth/time
            x = ds.Time.values
            y = depth
            z = ds.avgLayerTemperatureAnomaly.isel(
                    nOceanRegionsTmp=regionIndex)
            z = z.transpose()

            colorbarLabel = '[$^\circ$C]'
            xLabel = 'Time [years]'
            yLabel = 'Depth [m]'

            title = 'Temperature Anomaly, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

            figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefixes['TAnomalyZ', region])

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

            caption = 'Trend of {} Temperature Anomaly vs depth from Year ' \
                      '{}'.format(region, simulationStartTime[0:4])
            write_image_xml(
                config=config,
                filePrefix=self.filePrefixes['TAnomalyZ', region],
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Trends vs Depth',
                groupLink='trendsvsdepth',
                thumbnailDescription=u'{} ΔT'.format(region),
                imageDescription=caption,
                imageCaption=caption)

            if plotOriginalFields:
                z = ds['avgLayerTemperature'].isel(nOceanRegionsTmp=regionIndex)
                z = z.transpose()

                title = 'Temperature, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

                figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefixes['TZ', region])

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

                caption = 'Trend of {} Temperature vs depth'.format(region)
                write_image_xml(
                    config=config,
                    filePrefix=self.filePrefixes['TZ', region],
                    componentName='Ocean',
                    componentSubdirectory='ocean',
                    galleryGroup='Trends vs Depth',
                    groupLink='trendsvsdepth',
                    thumbnailDescription=u'{} T'.format(region),
                    imageDescription=caption,
                    imageCaption=caption)


            #  Second S vs depth/time
            z = ds.avgLayerSalinityAnomaly.isel(nOceanRegionsTmp=regionIndex)
            z = z.transpose()

            title = 'Salinity Anomaly, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

            colorbarLabel = '[PSU]'

            figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefixes['SAnomalyZ', region])

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

            caption = 'Trend of {} Salinity Anomaly vs depth from Year ' \
                      '{}'.format(region, simulationStartTime[0:4])
            write_image_xml(
                config=config,
                filePrefix=self.filePrefixes['SAnomalyZ', region],
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Trends vs Depth',
                groupLink='trendsvsdepth',
                thumbnailDescription=u'{} ΔS'.format(region),
                imageDescription=caption,
                imageCaption=caption)

            if plotOriginalFields:
                z = ds['avgLayerSalinity'].isel(nOceanRegionsTmp=regionIndex)
                z = z.transpose()

                title = 'Salinity, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

                figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefixes['SZ', region])

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

                caption = 'Trend of {} Salinity vs depth'.format(region)
                write_image_xml(
                    config=config,
                    filePrefix=self.filePrefixes['SZ', region],
                    componentName='Ocean',
                    componentSubdirectory='ocean',
                    galleryGroup='Trends vs Depth',
                    groupLink='trendsvsdepth',
                    thumbnailDescription=u'{} S'.format(region),
                    imageDescription=caption,
                    imageCaption=caption)


            #  Third OHC vs depth/time
            ohcAnomaly = ds.ohcAnomaly.isel(nOceanRegionsTmp=regionIndex)
            ohcAnomalyScaled = unitsScalefactor*ohcAnomaly
            z = ohcAnomalyScaled.transpose()

            title = 'OHC Anomaly, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

            colorbarLabel = '[x$10^{22}$ J]'

            figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefixes['OHCAnomalyZ', region])

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

            caption = 'Trend of {} OHC Anomaly vs depth from Year ' \
                      '{}'.format(region, simulationStartTime[0:4])
            write_image_xml(
                config=config,
                filePrefix=self.filePrefixes['OHCAnomalyZ', region],
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Trends vs Depth',
                groupLink='trendsvsdepth',
                thumbnailDescription=u'{} ΔOHC'.format(region),
                imageDescription=caption,
                imageCaption=caption)

            if plotOriginalFields:
                ohc = ds.ohc.isel(nOceanRegionsTmp=regionIndex)
                ohcScaled = unitsScalefactor*ohc
                z = ohcScaled.transpose()

                title = 'OHC, {} \n {}'.format(plotTitles[regionIndex], mainRunName)

                figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefixes['OHCZ', region])

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

                caption = 'Trend of {} OHC vs depth'.format(region)
                write_image_xml(
                    config=config,
                    filePrefix=self.filePrefixes['OHCZ', region],
                    componentName='Ocean',
                    componentSubdirectory='ocean',
                    galleryGroup='Trends vs Depth',
                    groupLink='trendsvsdepth',
                    thumbnailDescription=u'{} OHC'.format(region),
                    imageDescription=caption,
                    imageCaption=caption)


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
            yLabel = '$\Delta$OHC  [x$10^{22}$ J]'

            title = 'OHC anomaly, {} \n {}'.format(plotTitles[regionIndex],
                                                   mainRunName)

            figureName = '{}/{}.png'.format(self.plotsDirectory,
                                            self.filePrefixes['OHCAnomaly', region])

            if preprocessedReferenceRunName != 'None':
                # these preprocessed data are OHC *anomalies*
                ohcPreprocessedTotal = dsPreprocessedTimeSlice.ohc_tot
                ohcPreprocessed700m = dsPreprocessedTimeSlice.ohc_700m
                ohcPreprocessed2000m = dsPreprocessedTimeSlice.ohc_2000m
                ohcPreprocessedBottom = dsPreprocessedTimeSlice.ohc_btm
                title = '{} (black lines) \n {} (red lines)'.format(title,
                                                                    preprocessedReferenceRunName)
                timeseries_analysis_plot(config, [ohcAnomalyTotal,
                                                  ohcAnomaly700m,
                                                  ohcAnomaly2000m,
                                                  ohcAnomalyBottom,
                                                  ohcPreprocessedTotal,
                                                  ohcPreprocessed700m,
                                                  ohcPreprocessed2000m,
                                                  ohcPreprocessedBottom],
                                         movingAveragePointsTimeSeries, title,
                                         xLabel, yLabel, figureName,
                                         lineStyles=['k-', 'k-', 'k--', 'k+',
                                                     'r-', 'r-', 'r--', 'r+'],
                                         lineWidths=[5, 3, 3, 3,
                                                     5, 3, 3, 3],
                                         legendText=['0-bottom', '0-700m',
                                                     '700-2000m', '2000m-bottom',
                                                     None, None, None, None],
                                         calendar=calendar)

            if (not compareWithObservations and
                    preprocessedReferenceRunName == 'None'):
                timeseries_analysis_plot(config, [ohcAnomalyTotal,
                                                  ohcAnomaly700m,
                                                  ohcAnomaly2000m,
                                                  ohcAnomalyBottom],
                                         movingAveragePointsTimeSeries, title,
                                         xLabel, yLabel, figureName,
                                         lineStyles=['k-', 'k-', 'k--', 'k+'],
                                         lineWidths=[5, 3, 3, 3],
                                         legendText=['0-bottom', '0-700m',
                                                     '700-2000m', '2000m-bottom'],
                                         calendar=calendar)

            caption = 'Running Mean of the Anomaly in {} Ocean Heat Content ' \
                      'from Year {}'.format(region, simulationStartTime[0:4])
            write_image_xml(
                config=config,
                filePrefix=self.filePrefixes['OHCAnomaly', region],
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

                title = 'OHC, {} \n {}'.format(plotTitles[regionIndex],
                                               mainRunName)

                figureName = '{}/{}.png'.format(self.plotsDirectory,
                                                self.filePrefixes['OHC', region])

                timeseries_analysis_plot(config, [ohcTotal,
                                                  ohc700m,
                                                  ohc2000m,
                                                  ohcBottom],
                                         movingAveragePointsTimeSeries, title,
                                         xLabel, yLabel, figureName,
                                         lineStyles=['k-', 'k-', 'k--', 'k+'],
                                         lineWidths=[5, 3, 3, 3],
                                         legendText=['0-bottom', '0-700m',
                                                     '700-2000m', '2000m-bottom'],
                                         calendar=calendar)

                caption = 'Running Mean of {} Ocean Heat Content'.format(region)
                write_image_xml(
                    config=config,
                    filePrefix=self.filePrefixes['OHC', region],
                    componentName='Ocean',
                    componentSubdirectory='ocean',
                    galleryGroup='Time Series',
                    groupLink='timeseries',
                    thumbnailDescription=u'{} OHC'.format(region),
                    imageDescription=caption,
                    imageCaption=caption)

        # }}}

    def _compute_ohc(self, ds, regionNames):  # {{{
        '''
        Compute the OHC time series.
        '''

        # specific heat [J/(kg*degC)]
        cp = self.namelist.getfloat('config_specific_heat_sea_water')
        # [kg/m3]
        rho = self.namelist.getfloat('config_density0')

        ds['ohc'] = rho*cp*ds['sumLayerMaskValue'] * \
            ds['avgLayerArea'] * ds['avgLayerThickness'] * \
            ds['avgLayerTemperature']
        ds.ohc.attrs['units'] = 'J'
        ds.ohc.attrs['description'] = 'Ocean heat content in each region'

        ds['ohcAnomaly'] = rho*cp*ds['sumLayerMaskValue'] * \
            ds['avgLayerArea'] * ds['avgLayerThickness'] * \
            ds.avgLayerTemperatureAnomaly
        ds.ohcAnomaly.attrs['units'] = 'J'
        ds.ohcAnomaly.attrs['description'] = \
            'Ocean heat content anomaly in each region'

        ds['regionNames'] = ('nOceanRegionsTmp', regionNames)

        return ds  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
