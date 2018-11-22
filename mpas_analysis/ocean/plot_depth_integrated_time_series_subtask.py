# -*- coding: utf-8 -*-
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot.plotting import timeseries_analysis_plot

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.timekeeping.utility import date_to_days, \
    days_to_datetime

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.time_series import compute_moving_avg, \
    combine_time_series_with_ncrcat


class PlotDepthIntegratedTimeSeriesSubtask(AnalysisTask):
    """
    Plots a time series, summed or averaged over various depth ranges

    Attributes
    ----------

    regionName : str
        The name of the region to plot

    inFileName : str
        The file containing the time-depth data set to plot

    outFileLabel : str
        The prefix on each plot and associated XML file

    fieldNameInTitle : str
        The name of the field being plotted, as used in the plot title

    mpasFieldName : str
        The name of the variable in the MPAS timeSeriesStatsMonthly output

    yAxisLabel : str
        the y-axis label of the plotted field (including units)

    sectionName : str
        A section in the config file where the colormap and contour values
        are defined

    thumbnailSuffix : str
        The text to be displayed under the thumbnail image, to which the
        region name will be prepended

    imageCaption : str
        The caption when mousing over the plot or displaying it full
        screen

    galleryGroup : str
        The name of the group of galleries in which this plot belongs

    groupSubtitle : str
        The subtitle of the group in which this plot belongs (or blank
        if none)

    groupLink : str
        A short name (with no spaces) for the link to the gallery group

    galleryName : str
        The name of the gallery in which this plot belongs

    controlConfig : ``MpasAnalysisConfigParser``
        The configuration options for the control run (if any)
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    def __init__(self, parentTask, regionName, inFileName, outFileLabel,
                 fieldNameInTitle, mpasFieldName, yAxisLabel, sectionName,
                 thumbnailSuffix, imageCaption, galleryGroup, groupSubtitle,
                 groupLink, galleryName, subtaskName=None, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``AnalysisTask``
            The parent task of which this is a subtask

        regionName : str
            The name of the region to plot

        inFileName : str
            The file containing the time-depth data set to plot

        outFileLabel : str
            The prefix on each plot and associated XML file

        fieldNameInTitle : str
            The name of the field being plotted, as used in the plot title

        mpasFieldName : str
            The name of the variable in the MPAS timeSeriesStatsMonthly output

        yAxisLabel : str
            the y-axis label of the plotted field

        sectionName : str
            a section in the config file where the colormap and contour values
            are defined

        thumbnailSuffix : str
            The text to be displayed under the thumbnail image, to which the
            region name will be prepended

        imageCaption : str
            the caption when mousing over the plot or displaying it full
            screen

        galleryGroup : str
            the name of the group of galleries in which this plot belongs

        groupSubtitle : str
            the subtitle of the group in which this plot belongs (or blank
            if none)

        groupLink : str
            a short name (with no spaces) for the link to the gallery group

        galleryName : str
            the name of the gallery in which this plot belongs

        subtaskName :  str, optional
            The name of the subtask (``plotTimeSeries<RegionName>`` by default)

        controlConfig : ``MpasAnalysisConfigParser``, optional
            The configuration options for the control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if subtaskName is None:
            suffix = regionName[0].upper() + regionName[1:]
            subtaskName = 'plotDepthIntegratedTimeSeries{}'.format(suffix)

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotDepthIntegratedTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName='ocean',
            tags=parentTask.tags,
            subtaskName=subtaskName)

        self.regionName = regionName
        self.inFileName = inFileName
        self.outFileLabel = outFileLabel
        self.fieldNameInTitle = fieldNameInTitle
        self.mpasFieldName = mpasFieldName
        self.yAxisLabel = yAxisLabel
        self.sectionName = sectionName

        self.controlConfig = controlConfig

        # xml/html related variables
        self.thumbnailSuffix = thumbnailSuffix
        self.imageCaption = imageCaption
        self.galleryGroup = galleryGroup
        self.groupSubtitle = groupSubtitle
        self.groupLink = groupLink
        self.galleryName = galleryName

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Greg Streletz

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotDepthIntegratedTimeSeriesSubtask, self).setup_and_check()

        config = self.config

        if self.controlConfig is not None:
            # we need to know what file to read from the control run so
            # an absolute path won't work
            assert(not os.path.isabs(self.inFileName))

            baseDirectory = build_config_full_path(
                self.controlConfig, 'output', 'timeSeriesSubdirectory')

            self.refFileName = '{}/{}'.format(baseDirectory,
                                              self.inFileName)

        preprocessedReferenceRunName = config.get(
                'runs', 'preprocessedReferenceRunName')
        if preprocessedReferenceRunName != 'None':

            assert(not os.path.isabs(self.inFileName))

            baseDirectory = build_config_full_path(
                config, 'output', 'timeSeriesSubdirectory')

            make_directories('{}/preprocessed'.format(baseDirectory))

            self.preprocessedIntermediateFileName = \
                '{}/preprocessed/intermediate_{}'.format(baseDirectory,
                                                         self.inFileName)
            self.preprocessedFileName = '{}/preprocessed/{}'.format(
                    baseDirectory, self.inFileName)

        if not os.path.isabs(self.inFileName):
            baseDirectory = build_config_full_path(
                config, 'output', 'timeSeriesSubdirectory')

            self.inFileName = '{}/{}'.format(baseDirectory,
                                             self.inFileName)

        mainRunName = self.config.get('runs', 'mainRunName')

        self.filePrefix = '{}_{}_{}'.format(self.outFileLabel,
                                            self.regionName,
                                            mainRunName)
        self.xmlFileNames = ['{}/{}.xml'.format(
                self.plotsDirectory, self.filePrefix)]

        return  # }}}

    def run_task(self):  # {{{
        """
        Compute vertical agregates of the data and plot the time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

        self.logger.info("\nPlotting depth-integrated time series of "
                         "{}...".format(self.fieldNameInTitle))

        config = self.config
        calendar = self.calendar

        mainRunName = config.get('runs', 'mainRunName')

        plotTitles = config.getExpression('regions', 'plotTitles')
        allRegionNames = config.getExpression('regions', 'regions')
        regionIndex = allRegionNames.index(self.regionName)
        regionNameInTitle = plotTitles[regionIndex]

        startDate = config.get('timeSeries', 'startDate')
        endDate = config.get('timeSeries', 'endDate')

        # Load data
        self.logger.info('  Load ocean data...')
        ds = open_mpas_dataset(fileName=self.inFileName,
                               calendar=calendar,
                               variableList=[self.mpasFieldName, 'depth'],
                               timeVariableNames=None,
                               startDate=startDate,
                               endDate=endDate)
        ds = ds.isel(nOceanRegionsTmp=regionIndex)

        depths = ds.depth.values

        divisionDepths = config.getExpression(self.sectionName, 'depths')

        # for each depth interval to plot, determine the top and bottom depth
        topDepths = [0, 0] + divisionDepths
        bottomDepths = [depths[-1]] + divisionDepths + [depths[-1]]

        legends = []
        for top, bottom in zip(topDepths, bottomDepths):
            if bottom == depths[-1]:
                legends.append('{}m-bottom'.format(top))
            else:
                legends.append('{}m-{}m'.format(top, bottom))

        # more possible symbols than we typically use
        lines = ['-', '-', '--', None, None, None, None]
        markers = [None, None, None, '+', 'o', '^', 'v']
        widths = [5, 3, 3, 3, 3, 3, 3]
        points = [None, None, None, 300, 300, 300, 300]

        color = 'k'

        xLabel = 'Time [years]'
        yLabel = self.yAxisLabel

        title = '{}, {} \n {} (black)'.format(self.fieldNameInTitle,
                                              regionNameInTitle,
                                              mainRunName)

        figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefix)

        timeSeries = []
        lineColors = []
        lineStyles = []
        lineMarkers = []
        lineWidths = []
        maxPoints = []
        legendText = []

        for rangeIndex in range(len(topDepths)):
            top = topDepths[rangeIndex]
            bottom = bottomDepths[rangeIndex]
            field = ds[self.mpasFieldName].where(ds.depth > top)
            field = field.where(ds.depth <= bottom)
            timeSeries.append(field.sum('nVertLevels'))

            lineColors.append(color)
            lineStyles.append(lines[rangeIndex])
            lineMarkers.append(markers[rangeIndex])
            lineWidths.append(widths[rangeIndex])
            maxPoints.append(points[rangeIndex])
            legendText.append(legends[rangeIndex])

        preprocessedReferenceRunName = config.get(
                'runs', 'preprocessedReferenceRunName')
        if preprocessedReferenceRunName != 'None':
            preprocessedInputDirectory = config.get(
                    'oceanPreprocessedReference', 'baseDirectory')

            self.logger.info('  Load in preprocessed reference data...')
            preprocessedFilePrefix = config.get(self.sectionName,
                                                'preprocessedFilePrefix')
            inFilesPreprocessed = '{}/{}.{}.year*.nc'.format(
                preprocessedInputDirectory, preprocessedFilePrefix,
                preprocessedReferenceRunName)

            combine_time_series_with_ncrcat(
                    inFilesPreprocessed, self.preprocessedIntermediateFileName,
                    logger=self.logger)
            dsPreprocessed = open_mpas_dataset(
                    fileName=self.preprocessedIntermediateFileName,
                    calendar=calendar,
                    timeVariableNames='xtime')

            yearStart = days_to_datetime(ds.Time.min(), calendar=calendar).year
            yearEnd = days_to_datetime(ds.Time.max(), calendar=calendar).year
            timeStart = date_to_days(year=yearStart, month=1, day=1,
                                     calendar=calendar)
            timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                                   calendar=calendar)

            yearEndPreprocessed = days_to_datetime(dsPreprocessed.Time.max(),
                                                   calendar=calendar).year
            if yearStart <= yearEndPreprocessed:
                dsPreprocessed = dsPreprocessed.sel(Time=slice(timeStart,
                                                               timeEnd))
            else:
                self.logger.warning('Warning: Preprocessed time series ends '
                                    'before the timeSeries startYear and will '
                                    'not be plotted.')
                preprocessedReferenceRunName = 'None'

            # rolling mean seems to have trouble with dask data sets so we
            # write out the data set and read it back as a single-file data set
            # (without dask)
            dsPreprocessed = dsPreprocessed.drop('xtime')
            write_netcdf(dsPreprocessed, self.preprocessedFileName)
            dsPreprocessed = xarray.open_dataset(self.preprocessedFileName)

        if preprocessedReferenceRunName != 'None':
            color = 'purple'
            title = '{} \n {} (purple)'.format(title,
                                               preprocessedReferenceRunName)

            preprocessedFieldPrefix = config.get(self.sectionName,
                                                 'preprocessedFieldPrefix')

            movingAveragePoints = config.getint(self.sectionName,
                                                'movingAveragePoints')

            suffixes = ['tot'] + ['{}m'.format(depth) for depth in
                                  divisionDepths] + ['btm']

            # these preprocessed data are already anomalies
            dsPreprocessed = compute_moving_avg(dsPreprocessed,
                                                movingAveragePoints)
            for rangeIndex in range(len(suffixes)):
                variableName = '{}_{}'.format(preprocessedFieldPrefix,
                                              suffixes[rangeIndex])
                if variableName in list(dsPreprocessed.data_vars.keys()):
                    timeSeries.append(dsPreprocessed[variableName])
                else:
                    self.logger.warning('Warning: Preprocessed variable {} '
                                        'not found. Skipping.'.format(
                                                variableName))
                    timeSeries.extend(None)

                lineColors.append(color)
                lineStyles.append(lines[rangeIndex])
                lineMarkers.append(markers[rangeIndex])
                lineWidths.append(widths[rangeIndex])
                maxPoints.append(points[rangeIndex])
                legendText.append(None)

        if self.controlConfig is not None:

            controlRunName = self.controlConfig.get('runs', 'mainRunName')

            title = '{} \n {} (red)'.format(title, controlRunName)

            self.logger.info('  Load ocean data from control run...')
            controlStartYear = self.controlConfig.getint('timeSeries',
                                                         'startYear')
            controlEndYear = self.controlConfig.getint('timeSeries',
                                                       'endYear')
            controlStartDate = '{:04d}-01-01_00:00:00'.format(controlStartYear)
            controlEndDate = '{:04d}-12-31_23:59:59'.format(controlEndYear)
            dsRef = open_mpas_dataset(fileName=self.refFileName,
                                      calendar=calendar,
                                      variableList=[self.mpasFieldName,
                                                    'depth'],
                                      timeVariableNames=None,
                                      startDate=controlStartDate,
                                      endDate=controlEndDate)
            dsRef = dsRef.isel(nOceanRegionsTmp=regionIndex)

            color = 'r'

            for rangeIndex in range(len(topDepths)):
                top = topDepths[rangeIndex]
                bottom = bottomDepths[rangeIndex]
                field = dsRef[self.mpasFieldName].where(dsRef.depth > top)
                field = field.where(dsRef.depth <= bottom)
                timeSeries.append(field.sum('nVertLevels'))

                lineColors.append(color)
                lineStyles.append(lines[rangeIndex])
                lineMarkers.append(markers[rangeIndex])
                lineWidths.append(widths[rangeIndex])
                maxPoints.append(points[rangeIndex])
                legendText.append(None)

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

        timeseries_analysis_plot(config=config, dsvalues=timeSeries, N=None,
                                 title=title, xlabel=xLabel, ylabel=yLabel,
                                 fileout=figureName, calendar=calendar,
                                 lineColors=lineColors, lineStyles=lineStyles,
                                 markers=lineMarkers, lineWidths=lineWidths,
                                 maxPoints=maxPoints, legendText=legendText,
                                 firstYearXTicks=firstYearXTicks,
                                 yearStrideXTicks=yearStrideXTicks)

        write_image_xml(
            config=config,
            filePrefix=self.filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup=self.galleryGroup,
            groupLink=self.groupLink,
            galleryName=self.galleryName,
            thumbnailDescription='{} {}'.format(self.regionName,
                                                self.thumbnailSuffix),
            imageDescription=self.imageCaption,
            imageCaption=self.imageCaption)

        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
