# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import os

from ..shared import AnalysisTask

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.generalized_reader import open_multifile_dataset
from ..shared.io import open_mpas_dataset

from ..shared.timekeeping.utility import date_to_days, days_to_datetime

from ..shared.io.utility import build_config_full_path

from ..shared.html import write_image_xml

from ..shared.time_series import compute_moving_avg


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

    refConfig : ``MpasAnalysisConfigParser``
        The configuration options for the reference run (if any)

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani, Greg Streletz
    """

    def __init__(self, parentTask, regionName, inFileName, outFileLabel,
                 fieldNameInTitle, mpasFieldName, yAxisLabel, sectionName,
                 thumbnailSuffix, imageCaption, galleryGroup, groupSubtitle,
                 groupLink, galleryName, subtaskName=None, refConfig=None):
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

        refConfig : ``MpasAnalysisConfigParser``, optional
            The configuration options for the reference run (if any)

        Authors
        -------
        Xylar Asay-Davis
        """

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

        self.refConfig = refConfig

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

        Authors
        -------
        Xylar Asay-Davis, Greg Streletz
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotDepthIntegratedTimeSeriesSubtask, self).setup_and_check()

        config = self.config

        if self.refConfig is not None:
            # we need to know what file to read from the reference run so
            # an absolute path won't work
            assert(not os.path.isabs(self.inFileName))

            baseDirectory = build_config_full_path(
                self.refConfig, 'output', 'timeSeriesSubdirectory')

            self.refFileName = '{}/{}'.format(baseDirectory,
                                              self.inFileName)

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

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani, Greg Streletz
        """

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
                               variableList=[self.mpasFieldName],
                               timeVariableNames=None,
                               startDate=startDate,
                               endDate=endDate)
        ds = ds.isel(nOceanRegionsTmp=regionIndex)

        # Note: restart file, not a mesh file because we need refBottomDepth,
        # not in a mesh file
        try:
            restartFile = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for OHC calculation')

        # Define/read in general variables
        self.logger.info('  Read in depth...')
        with xr.open_dataset(restartFile) as dsRestart:
            # reference depth [m]
            depths = dsRestart.refBottomDepth.values

        # add depths as a coordinate to the data set
        ds.coords['depth'] = (('nVertLevels',), depths)

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
        lines = ['-', '-', '--', '+', 'o', '^', 'v']
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
        lineStyles = []
        lineWidths = []
        maxPoints = []
        legendText = []

        for rangeIndex in range(len(topDepths)):
            top = topDepths[rangeIndex]
            bottom = bottomDepths[rangeIndex]
            field = ds[self.mpasFieldName].where(ds.depth > top)
            field = field.where(ds.depth <= bottom)
            timeSeries.append(field.sum('nVertLevels'))

            lineStyles.append('{}{}'.format(color, lines[rangeIndex]))
            lineWidths.append(widths[rangeIndex])
            maxPoints.append(points[rangeIndex])
            legendText.append(legends[rangeIndex])

        preprocessedReferenceRunName = config.get(
                'runs', 'preprocessedReferenceRunName')
        if preprocessedReferenceRunName != 'None':
            preprocessedInputDirectory = config.get(
                    'oceanPreprocessedReference', 'baseDirectory')

            self.logger.info('  Load in OHC from preprocessed reference '
                             'run...')
            preprocessedFilePrefix = config.get(self.sectionName,
                                                'preprocessedFilePrefix')
            inFilesPreprocessed = '{}/{}.{}.year*.nc'.format(
                preprocessedInputDirectory, preprocessedFilePrefix,
                preprocessedReferenceRunName)
            dsPreprocessed = open_multifile_dataset(
                fileNames=inFilesPreprocessed,
                calendar=calendar,
                config=config,
                timeVariableName='xtime')

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

        if preprocessedReferenceRunName != 'None':
            color = 'r'
            title = '{} \n {} (red)'.format(title,
                                            preprocessedReferenceRunName)

            preprocessedFieldPrefix = config.get(self.sectionName,
                                                 'preprocessedFieldPrefix')

            movingAveragePoints = config.getint(self.sectionName,
                                                'movingAveragePoints')

            suffixes = ['tot'] + ['{}m'.format(depth) for depth in
                                  divisionDepths] + ['btm']

            # these preprocessed data are OHC *anomalies*
            for rangeIndex in range(len(suffixes)):
                variableName = '{}_{}'.format(preprocessedFieldPrefix,
                                              suffixes[rangeIndex])
                if variableName in list(dsPreprocessed.data_vars.keys()):
                    field = dsPreprocessed[variableName]
                    field = compute_moving_avg(field, movingAveragePoints)
                    timeSeries.append(field)
                else:
                    self.logger.warning('Warning: Preprocessed variable {} '
                                        'not found. Skipping.'.format(
                                                variableName))
                    timeSeries.extend(None)

                lineStyles.append('{}{}'.format(color, lines[rangeIndex]))
                lineWidths.append(widths[rangeIndex])
                maxPoints.append(points[rangeIndex])
                legendText.append(None)

        if self.refConfig is not None:

            refRunName = self.refConfig.get('runs', 'mainRunName')

            title = '{} \n {} (blue)'.format(title, refRunName)

            self.logger.info('  Load ocean data from reference run...')
            refStartDate = self.refConfig.get('timeSeries', 'startDate')
            refEndDate = self.refConfig.get('timeSeries', 'endDate')
            dsRef = open_mpas_dataset(fileName=self.refFileName,
                                      calendar=calendar,
                                      variableList=[self.mpasFieldName],
                                      timeVariableNames=None,
                                      startDate=refStartDate,
                                      endDate=refEndDate)
            dsRef = dsRef.isel(nOceanRegionsTmp=regionIndex)

            # add depths as a coordinate to the data set
            dsRef.coords['depth'] = (('nVertLevels',), depths)

            color = 'b'

            for rangeIndex in range(len(topDepths)):
                top = topDepths[rangeIndex]
                bottom = bottomDepths[rangeIndex]
                field = dsRef[self.mpasFieldName].where(dsRef.depth > top)
                field = field.where(dsRef.depth <= bottom)
                timeSeries.append(field.sum('nVertLevels'))

                lineStyles.append('{}{}'.format(color, lines[rangeIndex]))
                lineWidths.append(widths[rangeIndex])
                maxPoints.append(points[rangeIndex])
                legendText.append(None)

        timeseries_analysis_plot(config=config, dsvalues=timeSeries, N=None,
                                 title=title, xlabel=xLabel, ylabel=yLabel,
                                 fileout=figureName, lineStyles=lineStyles,
                                 lineWidths=lineWidths, maxPoints=maxPoints,
                                 legendText=legendText, calendar=calendar)

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
