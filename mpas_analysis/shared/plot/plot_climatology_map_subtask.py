# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE

import xarray as xr
import numpy as np

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.plot import plot_global_comparison, \
    plot_projection_comparison

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.climatology import \
    get_remapped_mpas_climatology_file_name
from mpas_analysis.shared.climatology.comparison_descriptors import \
    get_comparison_descriptor

from mpas_analysis.shared.projection import comparison_grid_titles


class PlotClimatologyMapSubtask(AnalysisTask):
    """
    An analysis task for plotting 2D model fields against observations.

    Attributes
    ----------
    season : str
        A season (key in ``shared.constants.monthDictionary``) to be
        plotted.

    comparisonGridName : str
        The name of the comparison grid to plot.

    remapMpasClimatologySubtask : mpas_analysis.shared.climatology.RemapMpasClimatologySubtask
        The subtask for remapping the MPAS climatology that this subtask
        will plot

    remapObsClimatologySubtask : mpas_analysis.shared.climatology.RemapObservedClimatologySubtask
        The subtask for remapping the observational climatology that this
        subtask will plot

    secondRemapMpasClimatologySubtask : mpas_analysis.shared.climatology.RemapMpasClimatologySubtask
        A second subtask for remapping another MPAS climatology to plot
        in the second panel and compare with in the third panel

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

    diffTitleLabel : str, optional
        the title of the difference subplot

    unitsLabel : str
        the units of the plotted field, to be displayed on color bars

    imageCaption : str
        the caption when mousing over the plot or displaying it full
        screen

    galleryGroup : str
        the name of the group of galleries in which this plot belongs

    groupSubtitle : str or None
        the subtitle of the group in which this plot belongs (or blank
        if none)

    groupLink : str
        a short name (with no spaces) for the link to the gallery group

    galleryName : str or None
        the name of the gallery in which this plot belongs

    depth : {None, float, 'top', 'bot'}
        Depth at which to perform the comparison, 'top' for the surface
        'bot' for the base

    configSectionName : str
        the name of the section where the color map and range is defined

    maskMinThreshold : float or None
        a value below which the field is mask out in plots

    maskMaxThreshold : float or None
        a value above which the field is mask out in plots

    extend : {'neither', 'both', 'min', 'max'}
        Determines the ``contourf``-coloring of values that are outside the
        range of the levels provided if using an indexed colormap.
    """

    def __init__(self, parentTask, season, comparisonGridName,
                 remapMpasClimatologySubtask, remapObsClimatologySubtask=None,
                 secondRemapMpasClimatologySubtask=None, controlConfig=None,
                 depth=None, removeMean=False, subtaskName=None):
        """
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  mpas_analysis.shared.AnalysisTask
            The parent (main) task for this subtask

        season : str
            A season (key in ``shared.constants.monthDictionary``) to be
            plotted.

        comparisonGridName : str
            The name of the comparison grid to plot.

        remapMpasClimatologySubtask : mpas_analysis.shared.climatology.RemapMpasClimatologySubtask
            The subtask for remapping the MPAS climatology that this subtask
            will plot

        remapObsClimatologySubtask : mpas_analysis.shared.climatology.RemapObservedClimatologySubtask, optional
            The subtask for remapping the observational climatology that this
            subtask will plot

        secondRemapMpasClimatologySubtask : mpas_analysis.shared.climatology.RemapMpasClimatologySubtask, optional
            A second subtask for remapping another MPAS climatology to plot
            in the second panel and compare with in the third panel

        controlConfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)

        depth : {float, 'top', 'bot'}, optional
            Depth the data is being plotted, 'top' for the sea surface
            'bot' for the sea floor

        removeMean : bool, optional
            If True, a common mask for the model and reference data sets is
            computed (where both are valid) and the mean over that mask is
            subtracted from both the model and reference results.  This is
            useful for data sets where the desire is to compare the spatial
            pattern but the mean offset is not meaningful (e.g. SSH)

        subtaskName : str, optional
            The name of the subtask.  If not specified, it is
            ``plot<season>_<comparisonGridName>`` with a suffix indicating the
            depth being sliced (if any)
        """

        self.season = season
        self.depth = depth
        self.comparisonGridName = comparisonGridName
        self.remapMpasClimatologySubtask = remapMpasClimatologySubtask
        self.remapObsClimatologySubtask = remapObsClimatologySubtask
        self.secondRemapMpasClimatologySubtask = \
            secondRemapMpasClimatologySubtask
        self.controlConfig = controlConfig
        self.removeMean = removeMean

        if depth is None:
            self.depthSuffix = ''
        else:
            self.depthSuffix = f'depth_{depth}'

        if subtaskName is None:
            subtaskName = f'plot{season}_{comparisonGridName}'
            if depth is not None:
                subtaskName = f'{subtaskName}_{self.depthSuffix}'

        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(PlotClimatologyMapSubtask, self).__init__(
            config=config, taskName=taskName, subtaskName=subtaskName,
            componentName=parentTask.componentName, tags=tags)

        # this task should not run until the remapping subtasks are done, since
        # it relies on data from those subtasks
        self.run_after(remapMpasClimatologySubtask)
        if remapObsClimatologySubtask is not None:
            self.run_after(remapObsClimatologySubtask)
        if secondRemapMpasClimatologySubtask is not None:
            self.run_after(secondRemapMpasClimatologySubtask)

        self.outFileLabel = None
        self.fieldNameInTitle = None
        self.mpasFieldName = None
        self.refFieldName = None
        self.refTitleLabel = None
        self.diffTitleLabel = None
        self.unitsLabel = None
        self.imageCaption = None
        self.galleryGroup = None
        self.groupSubtitle = None
        self.groupLink = None
        self.galleryName = None
        self.configSectionName = None
        self.thumbnailDescription = None
        self.startYear = None
        self.endYear = None
        self.startDate = None
        self.endDate = None
        self.filePrefix = None
        self.maskMinThreshold = None
        self.maskMaxThreshold = None
        self.extend = 'both'
        self.prependComparisonGrid = None

    def set_plot_info(self, outFileLabel, fieldNameInTitle, mpasFieldName,
                      refFieldName, refTitleLabel, unitsLabel,
                      imageCaption, galleryGroup, groupSubtitle, groupLink,
                      galleryName, diffTitleLabel='Model - Observations',
                      configSectionName=None, maskMinThreshold=None,
                      maskMaxThreshold=None, extend=None,
                      prependComparisonGrid=True):
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

        groupSubtitle : str or None
            the subtitle of the group in which this plot belongs (or blank
            if none)

        groupLink : str
            a short name (with no spaces) for the link to the gallery group

        galleryName : str or None
            the name of the gallery in which this plot belongs

        diffTitleLabel : str, optional
            the title of the difference subplot

        configSectionName : str or None, optional
            the name of the section where the color map and range is defined,
            default is the name of the task

        maskMinThreshold : float or None, optional
            a value below which the field is mask out in plots

        maskMaxThreshold : float or None, optional
            a value above which the field is mask out in plots

        extend : {'neither', 'both', 'min', 'max'}, optional
            Determines the ``contourf``-coloring of values that are outside the
            range of the levels provided if using an indexed colormap

        prependComparisonGrid : bool, optional
            Whether to prepend the name of the comparison grid to the gallery
            group
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
        self.maskMinThreshold = maskMinThreshold
        self.maskMaxThreshold = maskMaxThreshold

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
            self.fieldNameInTitle = f'Sea Surface {fieldNameInTitle}'
            self.thumbnailDescription = f'{season} surface'
        elif depth == 'bot':
            self.fieldNameInTitle = f'Sea Floor {fieldNameInTitle}'
            self.thumbnailDescription = f'{season} floor'
        else:
            self.fieldNameInTitle = f'{fieldNameInTitle} at z={depth} m'
            self.thumbnailDescription = f'{season} z={depth} m'

        if extend is not None:
            self.extend = extend

        self.prependComparisonGrid = prependComparisonGrid

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
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

        prefixPieces = [self.outFileLabel]
        if self.comparisonGridName != 'latlon':
            prefixPieces.append(self.comparisonGridName)
        prefixPieces.append(mainRunName)
        if self.depth is not None:
            prefixPieces.append(self.depthSuffix)
        years = f'years{self.startYear:04d}-{self.endYear:04d}'
        prefixPieces.extend([self.season, years])

        self.filePrefix = '_'.join(prefixPieces)

        self.xmlFileNames.append(
            f'{self.plotsDirectory}/{self.filePrefix}.xml')

    def run_task(self):
        """
        Plots a comparison of E3SM/MPAS output to SST/TEMP, SSS/SALT or MLD
        observations or a control run
        """

        season = self.season
        depth = self.depth
        comparisonGridName = self.comparisonGridName
        self.logger.info(
            f"\nPlotting 2-d maps of {self.fieldNameInTitle} climatologies "
            f"for {season} on the {comparisonGridName} grid...")

        # first read the model climatology
        remappedFileName = \
            self.remapMpasClimatologySubtask.get_remapped_file_name(
                season=season, comparisonGridName=comparisonGridName)

        remappedModelClimatology = xr.open_dataset(remappedFileName)

        if depth is not None:
            if str(depth) not in remappedModelClimatology.depthSlice.values:
                raise KeyError(f'The climatology you are attempting to '
                               f'perform depth slices of was originally '
                               f'created without depth {depth}. You will need '
                               f'to delete and regenerate the climatology')

            remappedModelClimatology = remappedModelClimatology.sel(
                depthSlice=str(depth), drop=True)

        # now the observations or control run
        if self.remapObsClimatologySubtask is not None:
            remappedFileName = self.remapObsClimatologySubtask.get_file_name(
                stage='remapped', season=season,
                comparisonGridName=comparisonGridName)

            remappedRefClimatology = xr.open_dataset(remappedFileName)
        elif self.secondRemapMpasClimatologySubtask is not None:
            remappedFileName = \
                self.secondRemapMpasClimatologySubtask.get_remapped_file_name(
                        season=season, comparisonGridName=comparisonGridName)
            remappedRefClimatology = xr.open_dataset(remappedFileName)
        elif self.controlConfig is not None:
            climatologyName = self.remapMpasClimatologySubtask.climatologyName
            remappedFileName = \
                get_remapped_mpas_climatology_file_name(
                    self.controlConfig, season=season,
                    componentName=self.componentName,
                    climatologyName=climatologyName,
                    comparisonGridName=comparisonGridName,
                    op=self.remapMpasClimatologySubtask.op)
            remappedRefClimatology = xr.open_dataset(remappedFileName)
            controlStartYear = self.controlConfig.getint('climatology',
                                                         'startYear')
            controlEndYear = self.controlConfig.getint('climatology',
                                                       'endYear')
            if controlStartYear != self.startYear or \
                    controlEndYear != self.endYear:
                self.refTitleLabel = \
                    f'{self.refTitleLabel}\n' \
                    f'(years {controlStartYear:04d}-{controlEndYear:04d})'

        else:
            remappedRefClimatology = None

        if remappedRefClimatology is not None and depth is not None:
            depthIndex = -1
            for index, depthSlice in enumerate(
                    remappedRefClimatology.depthSlice.values):
                try:
                    depthSlice = depthSlice.decode("utf-8")
                except AttributeError:
                    pass
                if depthSlice == str(depth):
                    depthIndex = index
            if depthIndex == -1:
                raise KeyError(f'The climatology you are attempting to '
                               f'perform depth slices of was originally '
                               f'created without depth {depth}. You will need '
                               f'to delete and regenerate the climatology')

            remappedRefClimatology = remappedRefClimatology.isel(
                depthSlice=depthIndex, drop=True)

        if self.removeMean:
            if remappedRefClimatology is None:
                remappedModelClimatology[self.mpasFieldName] = \
                    remappedModelClimatology[self.mpasFieldName] - \
                    remappedModelClimatology[self.mpasFieldName].mean()
            else:
                masked = remappedModelClimatology[self.mpasFieldName].where(
                    remappedRefClimatology[self.refFieldName].notnull())
                remappedModelClimatology[self.mpasFieldName] = \
                    remappedModelClimatology[self.mpasFieldName] - \
                    masked.mean()

                masked = remappedRefClimatology[self.refFieldName].where(
                    remappedModelClimatology[self.mpasFieldName].notnull())
                remappedRefClimatology[self.refFieldName] = \
                    remappedRefClimatology[self.refFieldName] - \
                    masked.mean()

        if self.componentName == 'ocean':
            componentName = 'Ocean'
            componentSubdirectory = 'ocean'
        elif self.componentName == 'seaIce':
            componentName = 'Sea Ice'
            componentSubdirectory = 'sea_ice'
        else:
            raise ValueError(f'Unexpected component: {self.componentName}')

        if self.comparisonGridName == 'latlon':
            self._plot_latlon(remappedModelClimatology, remappedRefClimatology,
                              componentName, componentSubdirectory)
        else:
            self._plot_projection(remappedModelClimatology,
                                  remappedRefClimatology,
                                  componentName, componentSubdirectory)

    def _plot_latlon(self, remappedModelClimatology, remappedRefClimatology,
                     componentName, componentSubdirectory):
        """ plotting a global lat-lon data set """

        season = self.season
        config = self.config
        configSectionName = self.configSectionName

        mainRunName = config.get('runs', 'mainRunName')

        modelOutput = _nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)

        lon = remappedModelClimatology['lon'].values
        lat = remappedModelClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        if remappedRefClimatology is None:
            refOutput = None
            bias = None
        else:
            refOutput = _nans_to_numpy_mask(
                remappedRefClimatology[self.refFieldName].values)

            bias = modelOutput - refOutput

            # mask with thresholds only after taking the diff
            refOutput = self._mask_with_thresholds(refOutput)

        modelOutput = self._mask_with_thresholds(modelOutput)

        if config.has_option(configSectionName, 'titleFontSize'):
            titleFontSize = config.getint(configSectionName, 'titleFontSize')
        else:
            titleFontSize = None

        if config.has_option(configSectionName, 'defaultFontSize'):
            defaultFontSize = config.getint(configSectionName,
                                            'defaultFontSize')
        else:
            defaultFontSize = None

        filePrefix = self.filePrefix
        outFileName = f'{self.plotsDirectory}/{filePrefix}.png'
        title = f'{self.fieldNameInTitle} ({season}, years ' \
                f'{self.startYear:04d}-{self.endYear:04d})'
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput,
                               refOutput,
                               bias,
                               configSectionName,
                               fileout=outFileName,
                               title=title,
                               modelTitle=mainRunName,
                               refTitle=self.refTitleLabel,
                               diffTitle=self.diffTitleLabel,
                               cbarlabel=self.unitsLabel,
                               titleFontSize=titleFontSize,
                               defaultFontSize=defaultFontSize,
                               extend=self.extend)

        caption = f'{season} {self.imageCaption}'
        write_image_xml(
            config,
            filePrefix,
            componentName=componentName,
            componentSubdirectory=componentSubdirectory,
            galleryGroup=f'Global {self.galleryGroup}',
            groupSubtitle=self.groupSubtitle,
            groupLink=f'global_{self.groupLink}',
            gallery=self.galleryName,
            thumbnailDescription=self.thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)

    def _plot_projection(self, remappedModelClimatology,
                         remappedRefClimatology,
                         componentName, componentSubdirectory):
        """ plotting a dataset on a projection grid """

        season = self.season
        comparisonGridName = self.comparisonGridName
        config = self.config
        configSectionName = self.configSectionName

        mainRunName = config.get('runs', 'mainRunName')

        validMask = remappedModelClimatology['validMask'].values
        landMask = np.ma.masked_array(
            np.ones(validMask.shape),
            mask=np.logical_not(np.isnan(validMask)))

        modelOutput = _nans_to_numpy_mask(
            remappedModelClimatology[self.mpasFieldName].values)

        if remappedRefClimatology is None:
            refOutput = None
            bias = None
        else:
            refOutput = _nans_to_numpy_mask(
                remappedRefClimatology[self.refFieldName].values)

            bias = modelOutput - refOutput

            # mask with maskValue only after taking the diff
            refOutput = self._mask_with_thresholds(refOutput)

        modelOutput = self._mask_with_thresholds(modelOutput)

        comparisonDescriptor = get_comparison_descriptor(
            config, comparisonGridName)
        x = comparisonDescriptor.x_corner
        y = comparisonDescriptor.y_corner

        aspectRatio = (x[-1] - x[0])/(y[-1] - y[0])

        # if the plots are even a bit wider than they are tall, make them
        # vertical
        vertical = aspectRatio > 1.2

        filePrefix = self.filePrefix
        outFileName = f'{self.plotsDirectory}/{filePrefix}.png'
        title = f'{self.fieldNameInTitle} ({season}, years ' \
                f'{self.startYear:04d}-{self.endYear:04d})'

        if config.has_option(configSectionName, 'titleFontSize'):
            titleFontSize = config.getint(configSectionName, 'titleFontSize')
        else:
            titleFontSize = None

        if config.has_option(configSectionName, 'defaultFontSize'):
            defaultFontSize = config.getint(configSectionName,
                                            'defaultFontSize')
        else:
            defaultFontSize = None

        if config.has_option(configSectionName, 'cartopyGridFontSize'):
            cartopyGridFontSize = config.getint(configSectionName,
                                                'cartopyGridFontSize')
        else:
            cartopyGridFontSize = None

        plot_projection_comparison(
            config,
            x,
            y,
            landMask,
            modelOutput,
            refOutput,
            bias,
            fileout=outFileName,
            colorMapSectionName=configSectionName,
            projectionName=comparisonGridName,
            title=title,
            modelTitle=mainRunName,
            refTitle=self.refTitleLabel,
            diffTitle=self.diffTitleLabel,
            cbarlabel=self.unitsLabel,
            titleFontSize=titleFontSize,
            cartopyGridFontSize=cartopyGridFontSize,
            defaultFontSize=defaultFontSize,
            vertical=vertical,
            extend=self.extend)

        if self.prependComparisonGrid:
            gridTitle = comparison_grid_titles[comparisonGridName]
            galleryGroup = f'{gridTitle} {self.galleryGroup}'
        else:
            galleryGroup = self.galleryGroup

        caption = f'{season} {self.imageCaption}'
        write_image_xml(
            config,
            filePrefix,
            componentName=componentName,
            componentSubdirectory=componentSubdirectory,
            galleryGroup=galleryGroup,
            groupSubtitle=self.groupSubtitle,
            groupLink=f'{comparisonGridName}_{self.groupLink}',
            gallery=self.galleryName,
            thumbnailDescription=self.thumbnailDescription,
            imageDescription=caption,
            imageCaption=caption)

    def _mask_with_thresholds(self, field):
        if self.maskMinThreshold is not None or \
                self.maskMaxThreshold is not None:
            mask = field.mask
            if self.maskMinThreshold is not None:
                mask = np.logical_or(mask, field <= self.maskMinThreshold)
            if self.maskMaxThreshold is not None:
                mask = np.logical_or(mask, field >= self.maskMaxThreshold)
            field = np.ma.masked_array(field, mask)

        return field


def _nans_to_numpy_mask(field):
    """
    Convert a numpy array with NaNs to a masked numpy array
    """
    field = np.ma.masked_array(
        field, np.isnan(field))
    return field
