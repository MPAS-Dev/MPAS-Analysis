"""
An analysis subtasks for plotting comparison of 2D model fields against
observations.

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy as np

from ..shared import AnalysisTask

from ..shared.plot.plotting import plot_global_comparison, \
    setup_colormap, plot_polar_projection_comparison

from ..shared.html import write_image_xml

from ..shared.grid import interp_extrap_corner


def nans_to_numpy_mask(field):
    field = np.ma.masked_array(
        field, np.isnan(field))
    return field


class PlotClimatologyMapSubtask(AnalysisTask):  # {{{
    """
    An analysis task for plotting 2D model fields against observations.

    Attributes
    ----------
    season : str
        A season (key in ``shared.constants.monthDictionary``) to be
        plotted.

    comparisonGridName : {'latlon', 'antarctic'}
        The name of the comparison grid to plot.

    remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
        The subtask for remapping the MPAS climatology that this subtask
        will plot

    remapObsClimatologySubtask : ``RemapObservedClimatologySubtask``
        The subtask for remapping the observational climatology that this
        subtask will plot

    outFileLabel : str
        The prefix on each plot and associated XML file

    fieldNameInTitle : str
        The name of the field being plotted, as used in the plot title

    mpasFieldName : str
        The name of the variable in the MPAS timeSeriesStatsMonthly output

    obsFieldName : str
        The name of the variable to use from the observations file

    observationTitleLabel : str
        the title of the observations subplot

    unitsLabel : str
        the units of the plotted field, to be displayed on color bars

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

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, parentTask, season, comparisonGridName,
                 remapMpasClimatologySubtask, remapObsClimatologySubtask):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        season : str
            A season (key in ``shared.constants.monthDictionary``) to be
            plotted.

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid to plot.

        remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
            The subtask for remapping the MPAS climatology that this subtask
            will plot

        remapObsClimatologySubtask : ``RemapObservedClimatologySubtask``
            The subtask for remapping the observational climatology that this
            subtask will plot

        Authors
        -------
        Xylar Asay-Davis

        '''

        self.season = season
        self.comparisonGridName = comparisonGridName
        self.remapMpasClimatologySubtask = remapMpasClimatologySubtask
        self.remapObsClimatologySubtask = remapObsClimatologySubtask
        subtaskName = 'plot{}_{}'.format(season, comparisonGridName)
        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(PlotClimatologyMapSubtask, self).__init__(
                config=config, taskName=taskName, subtaskName=subtaskName,
                componentName='ocean', tags=tags)

        # this task should not run until the remapping subtasks are done, since
        # it relies on data from those subtasks
        self.run_after(remapMpasClimatologySubtask)
        self.run_after(remapObsClimatologySubtask)
        # }}}

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      obsFieldName, observationTitleLabel, unitsLabel,
                      imageCaption, galleryGroup, groupSubtitle, groupLink,
                      galleryName):  # {{{
        """
        Store attributes related to plots, plot file names and HTML output.

        Parameters
        ----------
        outFileLabel : str
            The prefix on each plot and associated XML file

        fieldNameInTitle : str
            The name of the field being plotted, as used in the plot title

        mpasFieldName : str
            The name of the variable in the MPAS timeSeriesStatsMonthly output

        obsFieldName : str
            The name of the variable to use from the observations file

        observationTitleLabel : str
            the title of the observations subplot

        unitsLabel : str
            the units of the plotted field, to be displayed on color bars

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

        Authors
        -------
        Xylar Asay-Davis
        """
        self.outFileLabel = outFileLabel
        self.fieldNameInTitle = fieldNameInTitle
        self.mpasFieldName = mpasFieldName
        self.obsFieldName = obsFieldName
        self.observationTitleLabel = observationTitleLabel
        self.unitsLabel = unitsLabel

        # xml/html related variables
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
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotClimatologyMapSubtask, self).setup_and_check()

        config = self.config
        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')
        self.startDate = config.get('climatology', 'startDate')
        self.endDate = config.get('climatology', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')

        self.xmlFileNames = []
        self.filePrefixes = {}

        if self.comparisonGridName == 'latlon':
            self.filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
                    self.outFileLabel, mainRunName,
                    self.season, self.startYear, self.endYear)
        else:
            self.filePrefix = '{}_{}_{}_{}_years{:04d}-{:04d}'.format(
                    self.outFileLabel, self.comparisonGridName,
                    mainRunName, self.season, self.startYear, self.endYear)
        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    self.filePrefix))
        # }}}

    def run_task(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations

        Authors
        -------
        Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
        """

        season = self.season
        comparisonGridName = self.comparisonGridName
        self.logger.info("\nPlotting 2-d maps of {} climatologies for {} on "
                         "the {} grid...".format(self.fieldNameInTitle,
                                                 season, comparisonGridName))

        # first read the model climatology
        remappedFileName = self.remapMpasClimatologySubtask.get_file_name(
            season=season, stage='remapped',
            comparisonGridName=comparisonGridName)

        remappedModelClimatology = xr.open_dataset(remappedFileName)

        # now the observations
        remappedFileName = self.remapObsClimatologySubtask.get_file_name(
            stage='remapped', season=season,
            comparisonGridName=comparisonGridName)

        remappedObsClimatology = xr.open_dataset(remappedFileName)

        if self.comparisonGridName == 'latlon':
            self._plot_latlon(remappedModelClimatology, remappedObsClimatology)
        elif self.comparisonGridName == 'antarctic':
            self._plot_antarctic(remappedModelClimatology,
                                 remappedObsClimatology)
        # }}}

    def _plot_latlon(self, remappedModelClimatology, remappedObsClimatology):
        # {{{
        """ plotting a global lat-lon data set """

        season = self.season
        config = self.config
        configSectionName = self.taskName

        mainRunName = config.get('runs', 'mainRunName')

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, configSectionName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, configSectionName, suffix='Difference')

        modelOutput = nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)

        lon = remappedModelClimatology['lon'].values
        lat = remappedModelClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        observations = nans_to_numpy_mask(
            remappedObsClimatology[self.obsFieldName].values)

        bias = modelOutput - observations

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput,
                               observations,
                               bias,
                               colormapResult,
                               colorbarLevelsResult,
                               colormapDifference,
                               colorbarLevelsDifference,
                               fileout=outFileName,
                               title=title,
                               modelTitle='{}'.format(mainRunName),
                               obsTitle=self.observationTitleLabel,
                               diffTitle='Model-Observations',
                               cbarlabel=self.unitsLabel)

        caption = '{} {}'.format(season, self.imageCaption)
        write_image_xml(
            config,
            filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Global {}'.format(self.galleryGroup),
            groupSubtitle=self.groupSubtitle,
            groupLink=self.groupLink,
            gallery=self.galleryName,
            thumbnailDescription=season,
            imageDescription=caption,
            imageCaption=caption)

        # }}}

    def _plot_antarctic(self, remappedModelClimatology,
                        remappedObsClimatology):  # {{{
        """ plotting an Antarctic data set """

        season = self.season
        comparisonGridName = self.comparisonGridName
        config = self.config
        configSectionName = self.taskName

        mainRunName = config.get('runs', 'mainRunName')

        oceanMask = remappedModelClimatology['validMask'].values
        self.landMask = np.ma.masked_array(
            np.ones(oceanMask.shape),
            mask=np.logical_not(np.isnan(oceanMask)))

        modelOutput = nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)

        observations = nans_to_numpy_mask(
            remappedObsClimatology[self.obsFieldName].values)

        bias = modelOutput - observations

        x = interp_extrap_corner(remappedModelClimatology['x'].values)
        y = interp_extrap_corner(remappedModelClimatology['y'].values)

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)

        if config.has_option(configSectionName, 'colormapIndicesResult'):
            colorMapType = 'indexed'
        elif config.has_option(configSectionName, 'normTypeResult'):
            colorMapType = 'norm'
        else:
            raise ValueError('config section {} contains neither the info'
                             'for an indexed color map nor for computing a '
                             'norm'.format(configSectionName))

        plot_polar_projection_comparison(
            config,
            x,
            y,
            self.landMask,
            modelOutput,
            observations,
            bias,
            fileout=outFileName,
            colorMapSectionName=configSectionName,
            colorMapType=colorMapType,
            title=title,
            modelTitle='{}'.format(mainRunName),
            obsTitle=self.observationTitleLabel,
            diffTitle='Model - Observations',
            cbarlabel=self.unitsLabel)

        upperGridName = comparisonGridName[0].upper() + comparisonGridName[1:]
        caption = '{} {}'.format(season, self.imageCaption)
        write_image_xml(
            config,
            filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='{} {}'.format(upperGridName,
                                        self.galleryGroup),
            groupSubtitle=self.groupSubtitle,
            groupLink=self.groupLink,
            gallery=self.galleryName,
            thumbnailDescription=season,
            imageDescription=caption,
            imageCaption=caption)

        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
