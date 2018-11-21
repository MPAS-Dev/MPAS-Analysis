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
"""
An analysis subtasks for plotting comparison of 2D model fields against
observations.
"""
# Authors
# -------
# Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy as np

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot.plotting import plot_global_comparison, \
    plot_polar_projection_comparison

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.grid import interp_extrap_corner

from mpas_analysis.shared.climatology import \
    get_remapped_mpas_climatology_file_name

from mpas_analysis.ocean.utility import nans_to_numpy_mask


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

    removeMean : bool, optional
        If True, a common mask for the model and reference data sets is
        computed (where both are valid) and the mean over that mask is
        subtracted from both the model and reference results.  This is
        useful for data sets where the desire is to compare the spatial
        pattern but the mean offset is not meaningful (e.g. SSH)

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

    depth : {None, float, 'top', 'bot'}
        Depth at which to perform the comparison, 'top' for the sea surface
        'bot' for the sea floor

    configSectionName : str
        the name of the section where the color map and range is defined
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    def __init__(self, parentTask, season, comparisonGridName,
                 remapMpasClimatologySubtask, remapObsClimatologySubtask=None,
                 refConfig=None, depth=None, removeMean=False,
                 subtaskName=None):  # {{{
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

        remapObsClimatologySubtask : ``RemapObservedClimatologySubtask``, optional
            The subtask for remapping the observational climatology that this
            subtask will plot

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)

        depth : {float, 'top', 'bot'}, optional
            Depth the data is being plotted, 'top' for the sea surface
            'bot' for the sea floor

        removeMean : bool, optional
            If True, a common mask for the model and reference data sets is
            computed (where both are valid) and the mean over that mask is
            subtracted from both the model and reference results.  This is
            useful for data sets where the desire is to compare the spatial
            pattern but the mean offset is not meaningful (e.g. SSH)

        subtaskName : str, optinal
            The name of the subtask.  If not specified, it is
            ``plot<season>_<comparisonGridName>`` with a suffix indicating the
            depth being sliced (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.season = season
        self.depth = depth
        self.comparisonGridName = comparisonGridName
        self.remapMpasClimatologySubtask = remapMpasClimatologySubtask
        self.remapObsClimatologySubtask = remapObsClimatologySubtask
        self.refConfig = refConfig
        self.removeMean = removeMean

        if depth is None:
            self.depthSuffix = ''
        else:
            self.depthSuffix = 'depth_{}'.format(depth)

        if subtaskName is None:
            subtaskName = 'plot{}_{}'.format(season, comparisonGridName)
            if depth is not None:
                subtaskName = '{}_{}'.format(subtaskName, self.depthSuffix)

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
        if remapObsClimatologySubtask is not None:
            self.run_after(remapObsClimatologySubtask)
        # }}}

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      controlFieldName, controlTitleLabel, unitsLabel,
                      imageCaption, galleryGroup, groupSubtitle, groupLink,
                      galleryName, diffTitleLabel='Model - Observations',
                      configSectionName=None):
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

        controlFieldName : str
            The name of the variable to use from the observations or reference
            file

        controlTitleLabel : str
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

        diffTitleLabel : str, optional
            the title of the difference subplot

        configSectionName : str, optional
            the name of the section where the color map and range is defined,
            default is the name of the task
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.outFileLabel = outFileLabel
        self.fieldNameInTitle = fieldNameInTitle
        self.mpasFieldName = mpasFieldName
        self.controlFieldName = controlFieldName
        self.controlTitleLabel = controlTitleLabel
        self.diffTitleLabel = diffTitleLabel
        self.unitsLabel = unitsLabel

        # xml/html related variables
        self.imageCaption = imageCaption
        self.galleryGroup = galleryGroup
        self.groupSubtitle = groupSubtitle
        self.groupLink = groupLink
        self.galleryName = galleryName

        if configSectionName is None:
            self.configSectionName = self.taskName
        else:
            self.configSectionName = configSectionName

        season = self.season
        depth = self.depth
        if depth is None:
            self.fieldNameInTitle = fieldNameInTitle
            self.thumbnailDescription = season
        elif depth == 'top':
            self.fieldNameInTitle = 'Sea Surface {}'.format(fieldNameInTitle)
            self.thumbnailDescription = '{} surface'.format(season)
        elif depth == 'bot':
            self.fieldNameInTitle = 'Sea Floor {}'.format(fieldNameInTitle)
            self.thumbnailDescription = '{} floor'.format(season)
        else:
            self.fieldNameInTitle = '{} at z={} m'.format(fieldNameInTitle,
                                                          depth)
            self.thumbnailDescription = '{} z={} m'.format(season, depth)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

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

        prefixPieces = [self.outFileLabel]
        if self.comparisonGridName != 'latlon':
            prefixPieces.append(self.comparisonGridName)
        prefixPieces.append(mainRunName)
        if self.depth is not None:
            prefixPieces.append(self.depthSuffix)
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
        """
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

        season = self.season
        depth = self.depth
        comparisonGridName = self.comparisonGridName
        self.logger.info("\nPlotting 2-d maps of {} climatologies for {} on "
                         "the {} grid...".format(self.fieldNameInTitle,
                                                 season, comparisonGridName))

        # first read the model climatology
        remappedFileName = \
            self.remapMpasClimatologySubtask.get_remapped_file_name(
                season=season, comparisonGridName=comparisonGridName)

        remappedModelClimatology = xr.open_dataset(remappedFileName)

        if depth is not None:
            if str(depth) not in remappedModelClimatology.depthSlice.values:
                raise KeyError('The climatology you are attempting to perform '
                               'depth slices of was originally created\n'
                               'without depth {}. You will need to delete and '
                               'regenerate the climatology'.format(depth))

            remappedModelClimatology = remappedModelClimatology.sel(
                depthSlice=str(depth), drop=True)

        # now the observations or reference run
        if self.remapObsClimatologySubtask is not None:
            remappedFileName = self.remapObsClimatologySubtask.get_file_name(
                stage='remapped', season=season,
                comparisonGridName=comparisonGridName)

            remappedRefClimatology = xr.open_dataset(remappedFileName)
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
                self.controlTitleLabel = '{}\n(years {:04d}-{:04d})'.format(
                    self.controlTitleLabel, refStartYear, refEndYear)

        else:
            remappedRefClimatology = None

        if remappedRefClimatology is not None and depth is not None:
            if str(depth) not in remappedRefClimatology.depthSlice.values:
                raise KeyError('The climatology you are attempting to perform '
                               'depth slices of was originally created\n'
                               'without depth {}. You will need to delete and '
                               'regenerate the climatology'.format(depth))

            remappedRefClimatology = remappedRefClimatology.sel(
                depthSlice=str(depth), drop=True)

        if self.removeMean:
            if remappedRefClimatology is None:
                remappedModelClimatology[self.mpasFieldName] = \
                    remappedModelClimatology[self.mpasFieldName] - \
                    remappedModelClimatology[self.mpasFieldName].mean()
            else:
                masked = remappedModelClimatology[self.mpasFieldName].where(
                    remappedRefClimatology[self.controlFieldName].notnull())
                remappedModelClimatology[self.mpasFieldName] = \
                    remappedModelClimatology[self.mpasFieldName] - \
                    masked.mean()

                masked = remappedRefClimatology[self.controlFieldName].where(
                    remappedModelClimatology[self.mpasFieldName].notnull())
                remappedRefClimatology[self.controlFieldName] = \
                    remappedRefClimatology[self.controlFieldName] - \
                    masked.mean()

        if self.comparisonGridName == 'latlon':
            self._plot_latlon(remappedModelClimatology, remappedRefClimatology)
        elif self.comparisonGridName == 'antarctic':
            self._plot_antarctic(remappedModelClimatology,
                                 remappedRefClimatology)
        # }}}

    def _plot_latlon(self, remappedModelClimatology, remappedRefClimatology):
        # {{{
        """ plotting a global lat-lon data set """

        season = self.season
        config = self.config
        configSectionName = self.configSectionName

        mainRunName = config.get('runs', 'mainRunName')

        modelOutput = nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)

        lon = remappedModelClimatology['lon'].values
        lat = remappedModelClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        if remappedRefClimatology is None:
            controlOutput = None
            bias = None
        else:
            controlOutput = nans_to_numpy_mask(
                remappedRefClimatology[self.controlFieldName].values)

            bias = modelOutput - controlOutput

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput,
                               controlOutput,
                               bias,
                               configSectionName,
                               fileout=outFileName,
                               title=title,
                               modelTitle='{}'.format(mainRunName),
                               controlTitle=self.controlTitleLabel,
                               diffTitle=self.diffTitleLabel,
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
            thumbnailDescription=self.thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)

        # }}}

    def _plot_antarctic(self, remappedModelClimatology,
                        remappedRefClimatology):  # {{{
        """ plotting an Antarctic data set """

        season = self.season
        comparisonGridName = self.comparisonGridName
        config = self.config
        configSectionName = self.configSectionName

        mainRunName = config.get('runs', 'mainRunName')

        oceanMask = remappedModelClimatology['validMask'].values
        self.landMask = np.ma.masked_array(
            np.ones(oceanMask.shape),
            mask=np.logical_not(np.isnan(oceanMask)))

        modelOutput = nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)

        if remappedRefClimatology is None:
            controlOutput = None
            bias = None
        else:
            controlOutput = nans_to_numpy_mask(
                remappedRefClimatology[self.controlFieldName].values)

            bias = modelOutput - controlOutput

        x = interp_extrap_corner(remappedModelClimatology['x'].values)
        y = interp_extrap_corner(remappedModelClimatology['y'].values)

        filePrefix = self.filePrefix
        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, self.startYear,
                self.endYear)

        plot_polar_projection_comparison(
            config,
            x,
            y,
            self.landMask,
            modelOutput,
            controlOutput,
            bias,
            fileout=outFileName,
            colorMapSectionName=configSectionName,
            title=title,
            modelTitle='{}'.format(mainRunName),
            controlTitle=self.controlTitleLabel,
            diffTitle=self.diffTitleLabel,
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
            thumbnailDescription=self.thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)

        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
