
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy.ma as ma
import numpy as np

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot.plotting import plot_polar_comparison

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.climatology import \
    get_remapped_mpas_climatology_file_name


class PlotClimatologyMapSubtask(AnalysisTask):  # {{{
    """
    An analysis task for plotting 2D model fields against observations or a
    reference run.

    Attributes
    ----------
    hemisphere : str
        The hemisphere to plot

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

    refConfig :  ``MpasAnalysisConfigParser``
        Configuration options for a reference run (if any)

    outFileLabel : str
        The prefix on each plot and associated XML file

    fieldNameInTitle : str
        The name of the field being plotted, as used in the plot title

    mpasFieldName : str
        The name of the variable in the MPAS timeSeriesStatsMonthly output

    refFieldName : str
        The name of the variable to use from the observations or reference file

    refTitleLabel : str
        the title of the observations or reference run subplot

    diffTitleLabel : str
        the title of the difference (e.g. bias) subplot

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
    Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, parentTask, hemisphere, season, comparisonGridName,
                 remapMpasClimatologySubtask, remapObsClimatologySubtask=None,
                 refConfig=None, subtaskSuffix=None):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        hemisphere : str
            The hemisphere to plot

        season : str
            A season (key in ``shared.constants.monthDictionary``) to be
            plotted.

        comparisonGridName : {'latlon'}
            The name of the comparison grid to plot.

        remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
            The subtask for remapping the MPAS climatology that this subtask
            will plot

        remapObsClimatologySubtask : ``RemapObservedClimatologySubtask``, optional
            The subtask for remapping the observational climatology that this
            subtask will plot

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)

        subtaskSuffix : str, optional
            A suffix on the subtask to ensure that it is unique (e.g. the
            observations being plotted)

        Authors
        -------
        Xylar Asay-Davis

        '''

        self.hemisphere = hemisphere
        self.season = season
        self.comparisonGridName = comparisonGridName
        self.remapMpasClimatologySubtask = remapMpasClimatologySubtask
        self.remapObsClimatologySubtask = remapObsClimatologySubtask
        self.refConfig = refConfig

        subtaskName = 'plot{}'.format(season)
        if subtaskSuffix is not None:
            subtaskName = '{}_{}'.format(subtaskName, subtaskSuffix)

        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(PlotClimatologyMapSubtask, self).__init__(
                config=config, taskName=taskName, subtaskName=subtaskName,
                componentName='seaIce', tags=tags)

        # this task should not run until the remapping subtasks are done, since
        # it relies on data from those subtasks
        self.run_after(remapMpasClimatologySubtask)
        if remapObsClimatologySubtask is not None:
            self.run_after(remapObsClimatologySubtask)
        # }}}

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      refFieldName, refTitleLabel, diffTitleLabel, unitsLabel,
                      imageDescription, imageCaption, galleryGroup,
                      groupSubtitle, groupLink, galleryName, maskValue=None):
        # {{{
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

        refFieldName : str
            The name of the variable to use from the observations or reference
            run file

        refTitleLabel : str
            the title of the observations or reference run subplot

        diffTitleLabel : str
            the title of the difference (e.g. bias) subplot

        unitsLabel : str
            the units of the plotted field, to be displayed on color bars

        imageDescription : str
            the caption when mousing over the plot

        imageCaption : str
            the caption when displaying the plot full screen

        galleryGroup : str
            the name of the group of galleries in which this plot belongs

        groupSubtitle : str
            the subtitle of the group in which this plot belongs (or blank
            if none)

        groupLink : str
            a short name (with no spaces) for the link to the gallery group

        galleryName : str
            the name of the gallery in which this plot belongs

        maskValue : float, optional
            a value to mask out in plots

        Authors
        -------
        Xylar Asay-Davis
        """
        self.outFileLabel = outFileLabel
        self.fieldNameInTitle = fieldNameInTitle
        self.mpasFieldName = mpasFieldName
        self.refFieldName = refFieldName
        self.refTitleLabel = refTitleLabel
        self.diffTitleLabel = diffTitleLabel
        self.unitsLabel = unitsLabel

        # xml/html related variables
        self.imageDescription = imageDescription
        self.imageCaption = imageCaption
        self.galleryGroup = galleryGroup
        self.groupSubtitle = groupSubtitle
        self.groupLink = groupLink
        self.galleryName = galleryName

        self.maskValue = maskValue
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup
        super(PlotClimatologyMapSubtask, self).setup_and_check()

        config = self.config
        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')
        self.startDate = config.get('climatology', 'startDate')
        self.endDate = config.get('climatology', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')

        self.xmlFileNames = []

        self.filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
                self.outFileLabel, mainRunName,
                self.season, self.startYear, self.endYear)
        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    self.filePrefix))

        # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        config = self.config
        season = self.season
        comparisonGridName = self.comparisonGridName

        self.logger.info("\nPlotting 2-d maps of {} climatologies for "
                         "{} against {}...".format(
                                 self.fieldNameInTitle,
                                 season, self.refTitleLabel))

        mainRunName = config.get('runs', 'mainRunName')
        startYear = self.startYear
        endYear = self.endYear

        hemisphere = self.hemisphere
        sectionName = self.taskName
        vertical = config.getboolean(sectionName, 'vertical')
        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        referenceLongitude = config.getfloat(sectionName,
                                             'referenceLongitude')
        minimumLatitude = config.getfloat(sectionName,
                                          'minimumLatitude')

        remappedFileName = \
            self.remapMpasClimatologySubtask.get_remapped_file_name(
                    season=season, comparisonGridName=comparisonGridName)
        remappedClimatology = xr.open_dataset(remappedFileName)

        modelOutput = remappedClimatology[self.mpasFieldName].values
        if self.maskValue is not None:
            modelOutput = ma.masked_values(modelOutput, self.maskValue)
        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        if self.remapObsClimatologySubtask is not None:
            remappedFileName = self.remapObsClimatologySubtask.get_file_name(
                stage='remapped', season=season,
                comparisonGridName=comparisonGridName)

            remappedRefClimatology = xr.open_dataset(remappedFileName)

            refOutput = remappedRefClimatology[self.refFieldName].values
            if self.maskValue is not None:
                refOutput = ma.masked_values(refOutput, self.maskValue)

            difference = modelOutput - refOutput
        elif self.refConfig is not None:
            climatologyName = self.remapMpasClimatologySubtask.climatologyName
            remappedFileName = \
                get_remapped_mpas_climatology_file_name(
                    self.refConfig, season=season,
                    componentName=self.componentName,
                    climatologyName=climatologyName,
                    comparisonGridName=comparisonGridName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)
            refStartYear = self.refConfig.getint('climatology', 'startYear')
            refEndYear = self.refConfig.getint('climatology', 'endYear')
            if refStartYear != self.startYear or refEndYear != self.endYear:
                self.refTitleLabel = '{}\n(years {:04d}-{:04d})'.format(
                        self.refTitleLabel, refStartYear, refEndYear)
        else:
            remappedRefClimatology = None

        if remappedRefClimatology is None:
            refOutput = None
            difference = None
        else:
            refOutput = remappedRefClimatology[self.refFieldName].values
            if self.maskValue is not None:
                refOutput = ma.masked_values(refOutput, self.maskValue)

            difference = modelOutput - refOutput

        startYear = self.startYear
        endYear = self.endYear
        filePrefix = self.filePrefix
        title = '{} ({}, years {:04d}-{:04d})'.format(
            self.fieldNameInTitle, season, startYear, endYear)
        fileout = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        plot_polar_comparison(
            config,
            lonTarg,
            latTarg,
            modelOutput,
            refOutput,
            difference,
            sectionName,
            title=title,
            fileout=fileout,
            plotProjection=plotProjection,
            latmin=minimumLatitude,
            lon0=referenceLongitude,
            modelTitle=mainRunName,
            refTitle=self.refTitleLabel,
            diffTitle=self.diffTitleLabel,
            cbarlabel=self.unitsLabel,
            vertical=vertical)

        write_image_xml(
            config,
            filePrefix,
            componentName='Sea Ice',
            componentSubdirectory='sea_ice',
            galleryGroup=self.galleryGroup,
            groupSubtitle=self.groupSubtitle,
            groupLink=self.groupLink,
            gallery=self.galleryName,
            thumbnailDescription=season,
            imageDescription=self.imageDescription,
            imageCaption=self.imageCaption)
        # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
