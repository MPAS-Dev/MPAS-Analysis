# -*- coding: utf-8 -*-
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
#
import os
import xarray
import numpy
import matplotlib.pyplot as plt

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf_with_fill
from mpas_analysis.shared.io.utility import build_config_full_path, \
    build_obs_path, make_directories, decode_strings
from mpas_analysis.shared.climatology import compute_climatology, \
    get_unmasked_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.plot import histogram_analysis_plot, savefig
from mpas_analysis.shared.html import write_image_xml


class OceanHistogram(AnalysisTask):
    """
    Plots a histogram of a 2-d ocean variable.

    """

    def __init__(self, config, mpasClimatologyTask, regionMasksTask,
                 controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig : mpas_tools.config.MpasConfigParser
            Configuration options for a control run (if any)
        """

        # first, call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config,
            taskName='oceanHistogram',
            componentName='ocean',
            tags=['climatology', 'regions', 'histogram', 'publicObs'])

        self.run_after(mpasClimatologyTask)
        self.mpasClimatologyTask = mpasClimatologyTask

        self.controlConfig = controlConfig
        mainRunName = config.get('runs', 'mainRunName')

        self.regionGroups = config.getexpression(self.taskName, 'regionGroups')
        self.regionNames = config.getexpression(self.taskName, 'regionNames')
        self.seasons = config.getexpression(self.taskName, 'seasons')
        self.variableList = config.getexpression(self.taskName, 'variableList')
        if config.has_option(self.taskName, 'weightList'):
            self.weightList = config.getexpression(self.taskName, 'weightList')
            if not self.weightList:
                self.weightList = None
            elif len(self.weightList) != len(self.variableList):
                raise ValueError('Histogram weightList is not the same '
                                 'length as variableList')
        else:
            self.weightList = None

        baseDirectory = build_config_full_path(
            config, 'output', 'histogramSubdirectory')
        if not os.path.exists(baseDirectory):
            make_directories(baseDirectory)

        self.obsList = config.getexpression(self.taskName, 'obsList')
        obsDicts = {
            'AVISO': {
                'suffix': 'AVISO',
                'gridName': 'Global_1.0x1.0degree',
                'gridFileName': 'SSH/zos_AVISO_L4_199210-201012_20180710.nc',
                'lonVar': 'lon',
                'latVar': 'lat',
                'sshVar': 'zos',
                'pressureAdjustedSSHVar': 'zos'}}

        for regionGroup in self.regionGroups:
            groupObsDicts = {}
            mpasMasksSubtask = regionMasksTask.add_mask_subtask(
                regionGroup=regionGroup)
            regionNames = mpasMasksSubtask.expand_region_names(
                self.regionNames)

            regionGroupSuffix = regionGroup.replace(' ', '_')
            filePrefix = f'histogram_{regionGroupSuffix}'

            # Add mask subtasks for observations and prep groupObsDicts
            # groupObsDicts is a subsetted version of localObsDicts with an
            # additional attribute for the maskTask
            for obsName in self.obsList:
                localObsDict = dict(obsDicts[obsName])
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=localObsDict['gridFileName'])
                obsMasksSubtask = regionMasksTask.add_mask_subtask(
                    regionGroup, obsFileName=obsFileName,
                    lonVar=localObsDict['lonVar'],
                    latVar=localObsDict['latVar'],
                    meshName=localObsDict['gridName'])
                localObsDict['maskTask'] = obsMasksSubtask
                groupObsDicts[obsName] = localObsDict

            for regionName in regionNames:
                sectionName = None

                # Compute weights for histogram
                if self.weightList is not None:
                    computeWeightsSubtask = ComputeHistogramWeightsSubtask(
                        self, regionName, mpasMasksSubtask, filePrefix,
                        self.variableList, self.weightList)
                    self.add_subtask(computeWeightsSubtask)

                for season in self.seasons:

                    # Generate histogram plots
                    plotRegionSubtask = PlotRegionHistogramSubtask(
                        self, regionGroup, regionName, controlConfig,
                        sectionName, filePrefix, mpasClimatologyTask,
                        mpasMasksSubtask, obsMasksSubtask, groupObsDicts,
                        self.variableList, self.weightList, season)
                    plotRegionSubtask.run_after(mpasMasksSubtask)
                    plotRegionSubtask.run_after(obsMasksSubtask)
                    if self.weightList is not None:
                        plotRegionSubtask.run_after(computeWeightsSubtask)
                    self.add_subtask(plotRegionSubtask)

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        OSError
            If files are not present
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super().setup_and_check()

        # Add variables and seasons to climatology task
        variableList = []
        for var in self.variableList:
            variableList.append(f'timeMonthly_avg_{var}')

        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=self.seasons)

        if len(self.obsList) > 1:
            raise ValueError('Histogram analysis does not currently support'
                             'more than one observational product')


class ComputeHistogramWeightsSubtask(AnalysisTask):
    """
    Fetches weight variables from MPAS output files for each variable in
    variableList.

    """
    def __init__(self, parentTask, regionName, mpasMasksSubtask, fullSuffix,
                 variableList, weightList):
        """
        Initialize weights task

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionName : str
            Name of the region to plot

        mpasMasksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask MPAS files for each region to plot, used
            to get the mask file name

        obsDicts : dict of dicts
            Information on the observations to compare agains

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix

        variableList: list of str
            List of variables which will be weighted

        weightList: list of str
            List of variables by which to weight the variables in
            variableList, of the same length as variableList

        """

        super(ComputeHistogramWeightsSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=f'weights_{fullSuffix}_{regionName}')

        self.mpasMasksSubtask = mpasMasksSubtask
        self.regionName = regionName
        self.filePrefix = fullSuffix
        self.variableList = variableList
        self.weightList = weightList

    def run_task(self):
        """
        Apply the region mask to each weight variable and save in a file common
        to that region.
        """

        config = self.config
        base_directory = build_config_full_path(
            config, 'output', 'histogramSubdirectory')

        # Get cell mask for the region
        region_mask_filename = self.mpasMasksSubtask.maskFileName
        ds_region_mask = xarray.open_dataset(region_mask_filename)
        mask_region_names = decode_strings(ds_region_mask.regionNames)
        region_index = mask_region_names.index(self.regionName)
        ds_mask = ds_region_mask.isel(nRegions=region_index)
        cell_mask = ds_mask.regionCellMasks == 1

        # Open the mesh file, which contains unmasked weight variables
        mesh_filename = self.get_mesh_filename()
        ds_mesh = xarray.open_dataset(mesh_filename)
        ds_mesh = ds_mesh.isel(Time=0)

        # Save the cell mask only for the region in its own file, which may be
        # referenced by future analysis (i.e., as a control run)
        new_region_mask_filename = \
            f'{base_directory}/{self.filePrefix}_{self.regionName}_mask.nc'
        write_netcdf_with_fill(ds_mask, new_region_mask_filename)

        if self.weightList is not None:
            ds_weights = xarray.Dataset()
            # Fetch the weight variables and mask them for each region
            for index, var in enumerate(self.variableList):
                weight_var_name = self.weightList[index]
                if weight_var_name in ds_mesh.keys():
                    var_name = f'timeMonthly_avg_{var}'
                    ds_weights[f'{var_name}_weight'] = \
                        ds_mesh[weight_var_name].where(cell_mask, drop=True)
                else:
                    self.logger.warning(
                        f'Weight variable {weight_var_name} is '
                        f'not in the mesh file, skipping'
                    )

        weights_filename = \
            f'{base_directory}/{self.filePrefix}_{self.regionName}_weights.nc'
        write_netcdf_with_fill(ds_weights, weights_filename)


class PlotRegionHistogramSubtask(AnalysisTask):
    """
    Plots a histogram diagram for a given ocean region

    Attributes
    ----------
    regionGroup : str
        Name of the collection of region to plot

    regionName : str
        Name of the region to plot

    sectionName : str
        The section of the config file to get options from

    controlConfig : mpas_tools.config.MpasConfigParser
        The configuration options for the control run (if any)

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    mpasMasksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask MPAS files for each region to plot, used
        to get the mask file name

    obsDicts : dict of dicts
        Information on the observations to compare against

    variableList: list of str
        list of variables to plot

    season : str
        The season to compute the climatology for
    """

    def __init__(self, parentTask, regionGroup, regionName, controlConfig,
                 sectionName, fullSuffix, mpasClimatologyTask,
                 mpasMasksSubtask, obsMasksSubtask, obsDicts, variableList,
                 weightList, season):

        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionGroup : str
            Name of the collection of region to plot

        regionName : str
            Name of the region to plot

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)

        sectionName : str
            The config section with options for this regionGroup

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        mpasMasksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask MPAS files for each region to plot, used
            to get the mask file name

        obsDicts : dict of dicts
            Information on the observations to compare agains

        season : str
            The season to comput the climatogy for
        """

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotRegionHistogramSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=f'plot_{fullSuffix}_{regionName}_{season}')

        self.run_after(mpasClimatologyTask)
        self.regionGroup = regionGroup
        self.regionName = regionName
        self.sectionName = sectionName
        self.controlConfig = controlConfig
        self.mpasClimatologyTask = mpasClimatologyTask
        self.mpasMasksSubtask = mpasMasksSubtask
        self.obsMasksSubtask = obsMasksSubtask
        self.obsDicts = obsDicts
        self.variableList = variableList
        self.weightList = weightList
        self.season = season
        self.filePrefix = fullSuffix

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present
        """

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(PlotRegionHistogramSubtask, self).setup_and_check()

        self.xmlFileNames = []
        for var in self.variableList:
            self.xmlFileNames.append(
                f'{self.plotsDirectory}/{self.filePrefix}_{var}_'
                f'{self.regionName}_{self.season}.xml')

    def run_task(self):
        """
        Plots histograms of properties in an ocean region.
        """

        self.logger.info(f"\nPlotting {self.season} histograms for "
                         f"{self.regionName}")

        config = self.config
        sectionName = self.sectionName

        calendar = self.calendar

        main_run_name = config.get('runs', 'mainRunName')

        region_mask_filename = self.mpasMasksSubtask.maskFileName

        ds_region_mask = xarray.open_dataset(region_mask_filename)

        mask_region_names = decode_strings(ds_region_mask.regionNames)
        region_index = mask_region_names.index(self.regionName)

        ds_mask = ds_region_mask.isel(nRegions=region_index)
        cell_mask = ds_mask.regionCellMasks == 1

        if len(self.obsDicts) > 0:
            obs_region_mask_filename = self.obsMasksSubtask.maskFileName
            ds_obs_region_mask = xarray.open_dataset(obs_region_mask_filename)
            mask_region_names = decode_strings(ds_region_mask.regionNames)
            region_index = mask_region_names.index(self.regionName)

            ds_obs_mask = ds_obs_region_mask.isel(nRegions=region_index)
            obs_cell_mask = ds_obs_mask.regionMasks == 1

        in_filename = get_unmasked_mpas_climatology_file_name(
            config, self.season, self.componentName, op='avg')
        ds = xarray.open_dataset(in_filename)

        base_directory = build_config_full_path(
            config, 'output', 'histogramSubdirectory')

        if self.weightList is not None:
            weights_filename = \
                f'{base_directory}/{self.filePrefix}_{self.regionName}_' \
                'weights.nc'
            ds_weights = xarray.open_dataset(weights_filename)

        if self.controlConfig is not None:
            control_run_name = self.controlConfig.get('runs', 'mainRunName')
            control_filename = get_unmasked_mpas_climatology_file_name(
                self.controlConfig, self.season, self.componentName, op='avg')
            ds_control = xarray.open_dataset(control_filename)
            base_directory = build_config_full_path(
                self.controlConfig, 'output', 'histogramSubdirectory')
            control_region_mask_filename = \
                f'{base_directory}/{self.filePrefix}_{self.regionName}_mask.nc'
            ds_control_region_masks = xarray.open_dataset(
                control_region_mask_filename)
            control_cell_mask = ds_control_region_masks.regionCellMasks == 1
            if self.weightList is not None:
                control_weights_filename = f'{base_directory}/' \
                    f'{self.filePrefix}_{self.regionName}_weights.nc'
                ds_control_weights = xarray.open_dataset(
                    control_weights_filename)

        if config.has_option(self.taskName, 'mainColor'):
            mainColor = config.get(self.taskName, 'mainColor')
        else:
            mainColor = 'C0'
        if config.has_option(self.taskName, 'obsColor'):
            obsColor = config.get(self.taskName, 'obsColor')
        else:
            obsColor = 'C1'
        if config.has_option(self.taskName, 'controlColor'):
            controlColor = config.get(self.taskName, 'controlColor')
        else:
            controlColor = 'C2'

        if config.has_option(self.taskName, 'lineWidth'):
            lineWidth = config.getfloat(self.taskName, 'lineWidth')
        else:
            lineWidth = None

        if config.has_option(self.taskName, 'titleFontSize'):
            titleFontSize = config.getint(self.taskName,
                                          'titleFontSize')
        else:
            titleFontSize = None
        if config.has_option(self.taskName, 'axisFontSize'):
            axisFontSize = config.getint(self.taskName,
                                         'axisFontSize')
        else:
            axisFontSize = None

        if config.has_option(self.taskName, 'defaultFontSize'):
            defaultFontSize = config.getint(self.taskName,
                                            'defaultFontSize')
        else:
            defaultFontSize = None
        if config.has_option(self.taskName, 'bins'):
            bins = config.getint(self.taskName, 'bins')
        else:
            bins = None

        yLabel = 'normalized Probability Density Function'

        for index, var in enumerate(self.variableList):

            fields = []
            weights = []
            legendText = []
            lineColors = []

            var_name = f'timeMonthly_avg_{var}'

            title = f'{self.regionName.replace("_", " ")}, {self.season}'

            caption = f'Normalized probability density function for ' \
                      f'{self.season} {var} climatologies in ' \
                      f'{self.regionName.replace("_", " ")}'

            # Note: consider modifying this for more professional headings
            varTitle = var

            fields.append(ds[var_name].where(cell_mask, drop=True))
            if self.weightList is not None:
                if f'{var_name}_weight' in ds_weights.keys():
                    weights.append(ds_weights[f'{var_name}_weight'].values)
                    caption = f'{caption} weighted by {self.weightList[index]}'
                else:
                    weights.append(None)
            else:
                weights.append(None)

            legendText.append(main_run_name)
            lineColors.append(mainColor)

            xLabel = f"{ds[var_name].attrs['long_name']} " \
                     f"({ds[var_name].attrs['units']})"

            for obs_name in self.obsDicts:
                localObsDict = dict(self.obsDicts[obs_name])
                obs_filename = build_obs_path(
                    config, component=self.componentName,
                    relativePath=localObsDict['gridFileName'])
                if f'{var}Var' not in localObsDict.keys():
                    self.logger.warn(
                        f'{var}Var is not present in {obs_name}, skipping '
                        f'{obs_name}')
                    continue
                obs_var_name = localObsDict[f'{var}Var']
                ds_obs = xarray.open_dataset(obs_filename)
                ds_obs = ds_obs.where(obs_cell_mask, drop=True)
                fields.append(ds_obs[obs_var_name])
                legendText.append(obs_name)
                lineColors.append(obsColor)
                weights.append(None)
            if self.controlConfig is not None:
                fields.append(ds_control[var_name].where(control_cell_mask,
                                                         drop=True))
                control_run_name = self.controlConfig.get('runs',
                                                          'mainRunName')
                legendText.append(control_run_name)
                lineColors.append(controlColor)
                weights.append(ds_control_weights[f'{var_name}_weight'].values)

            if lineWidth is not None:
                lineWidths = [lineWidth for i in fields]
            else:
                lineWidths = None

            histogram_analysis_plot(config, fields, calendar=calendar,
                                    title=title, xLabel=xLabel, yLabel=yLabel,
                                    bins=bins, weights=weights,
                                    lineColors=lineColors,
                                    lineWidths=lineWidths,
                                    legendText=legendText,
                                    titleFontSize=titleFontSize,
                                    defaultFontSize=defaultFontSize)

            out_filename = f'{self.plotsDirectory}/{self.filePrefix}_{var}_' \
                           f'{self.regionName}_{self.season}.png'
            savefig(out_filename, config)

            write_image_xml(
                config=config,
                filePrefix=f'{self.filePrefix}_{var}_{self.regionName}_'
                           f'{self.season}',
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup=f'{self.regionGroup} Histograms',
                groupLink=f'histogram{var}',
                gallery=varTitle,
                thumbnailDescription=f'{self.regionName.replace("_", " ")} '
                                     f'{self.season}',
                imageDescription=caption,
                imageCaption=caption)
