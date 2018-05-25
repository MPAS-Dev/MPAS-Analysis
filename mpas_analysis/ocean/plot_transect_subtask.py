# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
"""
An analysis subtasks for plotting comparison of a model transect against
observations or reference data.

Authors
-------
Xylar Asay-Davis, Greg Streletz
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy

import matplotlib.pyplot as plt
from mpas_analysis.shared.plot.plotting import plot_vertical_section_comparison

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.climatology import compute_climatology, \
    get_remapped_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants

from mpas_analysis.ocean.utility import nans_to_numpy_mask


class PlotTransectSubtask(AnalysisTask):  # {{{
    """
    An analysis task for plotting 2D model fields against observations.

    Attributes
    ----------
    season : str
        A season (key in ``shared.constants.monthDictionary``) to be
        plotted.

    transectName : str
        The name of the transect to plot

    remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
        The subtask for remapping the MPAS climatology that this subtask
        will plot

    plotObs : bool, optional
        Whether to plot against observations.

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

    diffTitleLabel : str, optional
        the title of the difference subplot

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

    configSectionName : str
        the name of the section where the color map and range is defined

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, parentTask, season, transectName, fieldName,
                 remapMpasClimatologySubtask, plotObs=True,
                 refConfig=None):
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

        transectName : str
            The name of the transect to plot

        fieldName : str
            The name of the field to plot (for use in the subtask name only)

        remapMpasClimatologySubtask : ``RemapMpasClimatologySubtask``
            The subtask for remapping the MPAS climatology that this subtask
            will plot

        plotObs : bool, optional
            Whether to plot against observations.

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any), ignored if
            ``plotObs == True``

        Authors
        -------
        Xylar Asay-Davis

        '''

        self.season = season
        self.transectName = transectName
        self.remapMpasClimatologySubtask = remapMpasClimatologySubtask
        self.plotObs = plotObs
        self.refConfig = refConfig
        subtaskName = 'plot{}_{}_{}'.format(season,
                                            transectName.replace(' ', '_'),
                                            fieldName)

        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(PlotTransectSubtask, self).__init__(
                config=config, taskName=taskName, subtaskName=subtaskName,
                componentName='ocean', tags=tags)

        # this task should not run until the remapping subtasks are done, since
        # it relies on data from those subtasks
        self.run_after(remapMpasClimatologySubtask)
        # }}}

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      refFieldName, refTitleLabel, unitsLabel,
                      imageCaption, galleryGroup, groupSubtitle, groupLink,
                      galleryName, configSectionName,
                      diffTitleLabel='Model - Observations'):
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
            file

        refTitleLabel : str
            the title of the observations or reference subplot

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

        configSectionName : str
            the name of the section where the color map and range is defined

        diffTitleLabel : str, optional
            the title of the difference subplot

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
        self.imageCaption = imageCaption
        self.galleryGroup = galleryGroup
        self.groupSubtitle = groupSubtitle
        self.groupLink = groupLink
        self.galleryName = galleryName

        self.fieldNameInTitle = fieldNameInTitle
        self.thumbnailDescription = self.season

        self.configSectionName = configSectionName
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
        super(PlotTransectSubtask, self).setup_and_check()

        config = self.config
        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')
        self.startDate = config.get('climatology', 'startDate')
        self.endDate = config.get('climatology', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')

        self.xmlFileNames = []

        prefixPieces = []
        if self.outFileLabel is not '':
            prefixPieces.append(self.outFileLabel)
        prefixPieces.append(self.transectName.replace(' ', '_'))
        prefixPieces.append(mainRunName)
        years = 'years{:04d}-{:04d}'.format(self.startYear, self.endYear)
        prefixPieces.extend([self.season, years])

        self.filePrefix = '_'.join(prefixPieces)

        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    self.filePrefix))
        # }}}

    def run_task(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations
        or a reference run

        Authors
        -------
        Xylar Asay-Davis, Greg Streletz
        """

        season = self.season
        transectName = self.transectName
        self.logger.info("\nPlotting {} over season {}"
                         "...".format(self.fieldNameInTitle,
                                      season))

        # first read the model climatology
        remappedFileName = \
            self.remapMpasClimatologySubtask.get_remapped_file_name(
                    season=season, comparisonGridName=transectName)

        remappedModelClimatology = xr.open_dataset(remappedFileName)

        # now the observations or reference run
        if self.plotObs:
            verticalComparisonGridName = \
                self.remapMpasClimatologySubtask.verticalComparisonGridName
            remappedFileName = \
                self.remapMpasClimatologySubtask.obsDatasets.get_out_file_name(
                    transectName,
                    verticalComparisonGridName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)

            # if Time is an axis, take the appropriate avarage to get the
            # climatology
            if 'Time' in remappedRefClimatology.dims:
                monthValues = constants.monthDictionary[season]
                remappedRefClimatology = compute_climatology(
                        remappedRefClimatology, monthValues, maskVaries=True)

        elif self.refConfig is not None:
            climatologyName = self.remapMpasClimatologySubtask.climatologyName
            remappedFileName = \
                get_remapped_mpas_climatology_file_name(
                    self.refConfig, season=season,
                    componentName=self.componentName,
                    climatologyName=climatologyName,
                    comparisonGridName=transectName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)
            refStartYear = self.refConfig.getint('climatology', 'startYear')
            refEndYear = self.refConfig.getint('climatology', 'endYear')
            if refStartYear != self.startYear or refEndYear != self.endYear:
                self.refTitleLabel = '{}\n(years {:04d}-{:04d})'.format(
                        self.refTitleLabel, refStartYear, refEndYear)

        else:
            remappedRefClimatology = None

        self._plot_transect(remappedModelClimatology, remappedRefClimatology)

        remappedModelClimatology.close()

        if remappedRefClimatology is not None:
            remappedRefClimatology.close()

        # }}}

    def _plot_transect(self, remappedModelClimatology, remappedRefClimatology):
        # {{{
        """ plotting the transect """

        season = self.season
        config = self.config
        configSectionName = self.configSectionName

        mainRunName = config.get('runs', 'mainRunName')


        # broadcast x and z to have the same dimensions
        x, z = xr.broadcast(remappedModelClimatology.x,
                            remappedModelClimatology.z)

        # convert x and z to numpy arrays, make a copy because they are
        # sometimes read-only (not sure why)
        x = x.values.copy().transpose()
        z = z.values.copy().transpose()

        # z is masked out with NaNs in some locations (where there is land) but
        # this makes pcolormesh unhappy so we'll zero out those locations
        z[numpy.isnan(z)] = 0.

        modelOutput = nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)
        modelOutput = modelOutput.transpose()

        if remappedRefClimatology is None:
            refOutput = None
            bias = None
        else:
            refOutput = remappedRefClimatology[self.refFieldName]
            dims = refOutput.dims
            refOutput = nans_to_numpy_mask(refOutput.values)
            if dims[1] != 'nPoints':
                assert(dims[0] == 'nPoints')
                refOutput = refOutput.transpose()

            bias = modelOutput - refOutput

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{}\n({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)

        # construct a three-panel comparison  plot for the transect

        xLabel = 'Distance [km]'
        yLabel = 'Depth [m]'

        plot_vertical_section_comparison(config,
                                         x,
                                         z,
                                         modelOutput,
                                         refOutput,
                                         bias,
                                         outFileName,
                                         configSectionName,
                                         cbarLabel=self.unitsLabel,
                                         xlabel=xLabel,
                                         ylabel=yLabel,
                                         title=title,
                                         modelTitle='{}'.format(mainRunName),
                                         refTitle=self.refTitleLabel,
                                         diffTitle=self.diffTitleLabel,
                                         invertYAxis=False,
                                         backgroundColor='#c49868')

        caption = '{} {}'.format(season, self.imageCaption)
        write_image_xml(
            config,
            filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup=self.galleryGroup,
            groupSubtitle=self.groupSubtitle,
            groupLink=self.groupLink,
            gallery=self.galleryName,
            thumbnailDescription=self.thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)

        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
