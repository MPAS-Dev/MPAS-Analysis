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
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#
import os
import xarray
import numpy
import matplotlib.pyplot as plt

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io import open_mpas_dataset
from mpas_analysis.shared.io.utility import build_config_full_path, build_obs_path, decode_strings
from mpas_analysis.shared.climatology import compute_climatology, \
    get_unmasked_mpas_climatology_file_name

from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.plot import histogram_analysis_plot, savefig
from mpas_analysis.shared.html import write_image_xml

class OceanHistogram(AnalysisTask):
    """
    Plots a histogram of a 2-d ocean variable.

    Attributes
    ----------
    variableDict : dict
        A dictionary of variables from the time series stats monthly output
        (keys), together with shorter, more convenient names (values)

    histogramFileName : str
        The name of the file where the histogram is stored

    controlConfig : mpas_tools.config.MpasConfigParser
        Configuration options for a control run (if one is provided)

    filePrefix : str
        The basename (without extension) of the PNG and XML files to write out
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, regionMasksTask, controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasHistogram: ``MpasHistogramTask``
            The task that extracts the time series from MPAS monthly output

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig : mpas_tools.config.MpasConfigParser
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

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
        self.filePrefix = f'histogram_{mainRunName}'


        obsList = config.getexpression(self.taskName, 'obsList')
        #TODO add gridName
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
            print(f'Add mask subtask for {regionGroup}, MPAS')
            mpasMasksSubtask = regionMasksTask.add_mask_subtask(
                regionGroup=regionGroup)
            regionNames = mpasMasksSubtask.expand_region_names(self.regionNames)
            if controlConfig is None:
                for obsName in obsList:
                    localObsDict = dict(obsDicts[obsName])
                    obsFileName = build_obs_path(
                        config, component=self.componentName,
                        relativePath=localObsDict['gridFileName'])
                    print(f'Add mask subtask for {regionGroup}, {obsName}, {localObsDict["gridName"]}')
                    obsMasksSubtask = regionMasksTask.add_mask_subtask(
                        regionGroup, obsFileName=obsFileName,
                        lonVar=localObsDict['lonVar'],
                        latVar=localObsDict['latVar'],
                        meshName=localObsDict['gridName'])
                    regionNames = obsMasksSubtask.expand_region_names(self.regionNames)
                    obsDicts[obsName]['maskTask'] = obsMasksSubtask

                    localObsDict['maskTask'] = obsMasksSubtask
                    groupObsDicts[obsName] = localObsDict
            for regionName in regionNames:
                sectionName = None
                fullSuffix = self.filePrefix
                for season in self.seasons:
                    plotRegionSubtask = PlotRegionHistogramSubtask(
                        self, regionGroup, regionName, controlConfig,
                        sectionName, fullSuffix, mpasClimatologyTask,
                        mpasMasksSubtask, obsMasksSubtask, groupObsDicts, self.variableList, season)
                    plotRegionSubtask.run_after(mpasMasksSubtask)
                    plotRegionSubtask.run_after(obsMasksSubtask)
                    self.add_subtask(plotRegionSubtask)
                    print(f'Add regional histogram subtask for {regionName}, {season}')

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        OSError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super().setup_and_check()

        config = self.config

        regionGroups = config.getexpression(self.taskName, 'regionGroups')
        variableList = []
        for var in self.variableList:
            variableList.append(f'timeMonthly_avg_{var}')

        # Specify variables and seasons to compute climology over
        print(f'add climatology variables')
        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=self.seasons)


    def run_task(self):
        """
        Performs histogram analysis of the output of variables in variableList.
        """
        # Authors
        # -------
        # Carolyn Begeman, Adrian Turner, Xylar Asay-Davis


class PlotRegionHistogramSubtask(AnalysisTask):
    """
    Plots a T-S diagram for a given ocean region

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

    varList: list of str
        list of variables to plot

    season : str
        The season to compute the climatology for
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, regionName, controlConfig,
                 sectionName, fullSuffix, mpasClimatologyTask,
                 mpasMasksSubtask, obsMasksSubtask, obsDicts, varList, season):

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
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        print(f'Initialize histogram task')
        super(PlotRegionHistogramSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=f'plot{fullSuffix}_{regionName}_{season}')

        self.run_after(mpasClimatologyTask)
        self.regionGroup = regionGroup
        self.regionName = regionName
        self.sectionName = sectionName
        self.controlConfig = controlConfig
        self.mpasClimatologyTask = mpasClimatologyTask
        self.mpasMasksSubtask = mpasMasksSubtask
        self.obsMasksSubtask = obsMasksSubtask
        self.obsDicts = obsDicts
        self.varList = varList
        self.season = season
        self.filePrefix = fullSuffix

        #TODO
        #parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        #self.subprocessCount = min(parallelTaskCount,
        #                           self.config.getint(self.taskName,
        #                                              'subprocessCount'))
        #self.daskThreads = min(
        #    multiprocessing.cpu_count(),
        #    self.config.getint(self.taskName, 'daskThreads'))

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(PlotRegionHistogramSubtask, self).setup_and_check()

        self.xmlFileNames = []
        for var in self.varList:
            print(f'add xml from subtask: {self.plotsDirectory}/{self.filePrefix}_{var}_{self.regionName}_{self.season}.xml')
            # Add xml file names for each season
            self.xmlFileNames.append(f'{self.plotsDirectory}/{self.filePrefix}_{var}_{self.regionName}_{self.season}.xml')

        print(f'end subtask setup and check')

    def run_task(self):
        """
        Plots time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nPlotting histogram for {}"
                         "...".format(self.regionName))

        config = self.config
        sectionName = self.sectionName

        self.logger.info('  Make plots...')

        calendar = self.calendar

        mainRunName = config.get('runs', 'mainRunName')

        baseDirectory = build_config_full_path(
            config, 'output', 'histogramSubdirectory')
        print(f'baseDirectory={baseDirectory}')
        print(f'plotsDirectory={self.plotsDirectory}')

        regionMaskFileName = self.mpasMasksSubtask.maskFileName
        print(f'Open {regionMaskFileName}')

        dsRegionMask = xarray.open_dataset(regionMaskFileName)

        maskRegionNames = decode_strings(dsRegionMask.regionNames)
        regionIndex = maskRegionNames.index(self.regionName)

        dsMask = dsRegionMask.isel(nRegions=regionIndex)
        cellMask = dsMask.regionCellMasks == 1

        inFileName = get_unmasked_mpas_climatology_file_name(
            config, self.season, self.componentName, op='avg')

        #TODO support control run

        #TODO: currently does not support len(obsList) > 1
        if len(self.obsDicts) > 0:
            obsRegionMaskFileName = self.obsMasksSubtask.maskFileName
            dsObsRegionMask = xarray.open_dataset(obsRegionMaskFileName)
            maskRegionNames = decode_strings(dsRegionMask.regionNames)
            regionIndex = maskRegionNames.index(self.regionName)

            dsObsMask = dsObsRegionMask.isel(nRegions=regionIndex)
            obsCellMask = dsObsMask.regionMasks == 1
        ds = xarray.open_dataset(inFileName)
        if config.has_option(self.taskName, 'weightByVariable'):
            weightVarName = config.get(self.taskName, 'weightByVariable')
        else:
            weightVarName = 'areaCell'
        weights = []
        restartFileName = self.runStreams.readpath('restart')[0]
        dsRestart = xarray.open_dataset(restartFileName)
        dsRestart = dsRestart.isel(Time=0)
        ds = ds.where(cellMask, drop=True)
        dsRestart = dsRestart.where(cellMask, drop=True)

        if config.has_option(self.taskName, 'lineColors'):
            lineColors = [config.get(self.taskName, 'mainColor')]
        else:
            lineColors = None
        if config.has_option(self.taskName, 'obsColor'):
            obsColor = config.get_expression(self.taskName, 'obsColor')
            if lineColors is None:
                lineColors = ['b']
        else:
            if lineColors is not None:
                obsColor = 'k'

        if config.has_option(self.taskName, 'lineWidths'):
            lineWidths = [config.get(self.taskName, 'lineWidths')]
        else:
            lineWidths = None
        legendText = [mainRunName]

        title = mainRunName
        if config.has_option(self.taskName, 'titleFontSize'):
            titleFontSize = config.getint(self.taskName,
                                          'titleFontSize')
        else:
            titleFontSize = None
        if config.has_option(self.taskName, 'titleFontSize'):
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

        for var in self.varList:

            varname = f'timeMonthly_avg_{var}'

            #TODO title as attribute or dict of var
            varTitle = var

            fields = [ds[varname]]
            weights.append(dsRestart[weightVarName].values)
            xLabel = f"{ds[varname].attrs['long_name']} ({ds[varname].attrs['units']})"
            for obsName in self.obsDicts:
                localObsDict = dict(self.obsDicts[obsName])
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=localObsDict['gridFileName'])
                varnameObs = localObsDict[f'{var}Var']
                print(f'{var} in obs is {varnameObs}')
                dsObs = xarray.open_dataset(obsFileName)
                dsObs = dsObs.where(obsCellMask, drop=True)
                fields.append(dsObs[varnameObs])
                legendText.append(obsName)
                if lineColors is not None:
                    lineColors.append(obsColor)
                if lineWidths is not None:
                    lineWidths.append([lineWidths[0]])
                weights.append(None)
            histogram_analysis_plot(config, fields, calendar=calendar,
                                    title=title, xlabel=xLabel, ylabel=yLabel, bins=bins, weights=weights,
                                    lineColors=lineColors, lineWidths=lineWidths,
                                    legendText=legendText,
                                    titleFontSize=titleFontSize, defaultFontSize=defaultFontSize)

            outFileName = f'{self.plotsDirectory}/{self.filePrefix}_{var}_{self.regionName}_{self.season}.png'
            savefig(outFileName, config)

            caption = f'Normalized probability density function for SSH climatologies in {self.regionName}'
            write_image_xml(
                config=config,
                filePrefix=f'{self.filePrefix}_{var}_{self.regionName}_{self.season}',
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup=f'{self.regionGroup} Histograms',
                groupLink=f'histogram{var}',
                gallery=varTitle,
                thumbnailDescription=self.regionName.replace('_', ' '),
                imageDescription=caption,
                imageCaption=caption)
