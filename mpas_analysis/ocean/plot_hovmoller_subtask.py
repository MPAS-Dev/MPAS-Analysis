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

import xarray as xr
import os

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot.plotting import plot_vertical_section

from mpas_analysis.shared.io import open_mpas_dataset

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.html import write_image_xml


class PlotHovmollerSubtask(AnalysisTask):
    """
    Plots a time series vs. depth

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

    unitsLabel : str
        The units of the plotted field, to be displayed on color bars

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

    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    def __init__(self, parentTask, regionName, inFileName, outFileLabel,
                 fieldNameInTitle, mpasFieldName, unitsLabel, sectionName,
                 thumbnailSuffix, imageCaption, galleryGroup, groupSubtitle,
                 groupLink, galleryName, subtaskName=None):  # {{{
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

        unitsLabel : str
            the units of the plotted field, to be displayed on color bars

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

        subtaskName :  str
            The name of the subtask (``plotHovmoller<RegionName>`` by default)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if subtaskName is None:
            suffix = regionName[0].upper() + regionName[1:]
            subtaskName = 'plotHovmoller{}'.format(suffix)

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotHovmollerSubtask, self).__init__(
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
        self.unitsLabel = unitsLabel
        self.sectionName = sectionName

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
        super(PlotHovmollerSubtask, self).setup_and_check()

        config = self.config

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
        Make the Hovmoller plot from the time series.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

        self.logger.info("\nPlotting {} trends vs. depth...".format(
                self.fieldNameInTitle))

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')

        plotTitles = config.getExpression('regions', 'plotTitles')
        allRegionNames = config.getExpression('regions', 'regions')
        regionIndex = allRegionNames.index(self.regionName)
        regionNameInTitle = plotTitles[regionIndex]

        startDate = self.config.get('timeSeries', 'startDate')
        endDate = self.config.get('timeSeries', 'endDate')

        # Load data
        self.logger.info('  Load ocean data...')
        ds = open_mpas_dataset(fileName=self.inFileName,
                               calendar=self.calendar,
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
            depth = dsRestart.refBottomDepth.values

        Time = ds.Time.values
        field = ds[self.mpasFieldName].values.transpose()

        xLabel = 'Time [years]'
        yLabel = 'Depth [m]'

        title = '{}, {} \n {}'.format(self.fieldNameInTitle, regionNameInTitle,
                                      mainRunName)

        figureName = '{}/{}.png'.format(self.plotsDirectory, self.filePrefix)

        if config.has_option(self.sectionName, 'firstYearXTicks'):
            firstYearXTicks = config.getint(self.sectionName,
                                            'firstYearXTicks')
        else:
            firstYearXTicks = None

        if config.has_option(self.sectionName, 'yearStrideXTicks'):
            yearStrideXTicks = config.getint(self.sectionName,
                                            'yearStrideXTicks')
        else:
            yearStrideXTicks = None

        plot_vertical_section(config, Time, depth, field, self.sectionName,
                              suffix='', colorbarLabel=self.unitsLabel,
                              title=title, xlabel=xLabel, ylabel=yLabel,
                              fileout=figureName, lineWidth=1,
                              xArrayIsTime=True, calendar=self.calendar,
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
